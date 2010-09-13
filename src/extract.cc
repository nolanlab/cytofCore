#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

	// Replacement for assert
	#define assert(X) { if (!(X)) { error("Assertion Failure %s [%s:%d]: %s\n",__FUNCTION__,__FILE__,__LINE__,#X); } }

    /* Get the list element named by str, or return NULL */
	SEXP getListElement(SEXP list, const char *str) {
		SEXP elmt = NULL_USER_OBJECT, names = getAttrib(list, R_NamesSymbol);
		for (int i = 0; i < length(list); i++) {
			if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
				elmt = VECTOR_ELT(list, i);
				break;
			}
		}
		return elmt;
	}

	struct Options {
		size_t analytes;
		
		double *conf;
		double  pulse_threshold;
		
		double *smooth;
		size_t  smooth_len;
		double  cell_threshold;

		size_t  min_cell_length;
		size_t  max_cell_length;
		bool    noise_subtraction;
		size_t  min_noise_length;

		bool    slope_filter;

		size_t  num_pushes;
		ssize_t blk_size;
	};

#define SLOPE(CONF, A, ANA) (CONF[(A)+ANA])
#define INTERCEPT(CONF, A, ANA) (CONF[(A)])

#define INTENSITY(PUSH, A) (PUSH[2*(A)])
#define PULSE(PUSH, A) (PUSH[2*(A)+1])

	class Dueler {
	public:
		const double* Conf;
		size_t        Analytes;

		double        PulseThreshold;

		double*       Buf;
		size_t        BufMask;
		size_t        Push;

	public:
		Dueler(const double* conf, size_t analytes, double pulse_threshold, size_t max_pushes=256) 
			: Conf(conf), Analytes(analytes), PulseThreshold(pulse_threshold)
		{
			// Create a cicular buffer for dual counts
			size_t p = std::ceil(std::log(max_pushes)/std::log(2.));  // Next largest power of two
			
			Buf = Calloc((1<<p)*Analytes, double);
			BufMask = ~(std::numeric_limits<size_t>::max() << p);
			Push = 0;
		}
		~Dueler() { Free(Buf); }

		size_t getPush() const { return Push; }
		void advPush() { ++Push; }
		
		void fromIMD(const uint16_t* imd) {
			double* dual = Buf + (Push & BufMask) * Analytes;
			for (size_t a=0; a<Analytes; a++) {	
				double d_d = fround((SLOPE(Conf,a,Analytes) * INTENSITY(imd,a) + INTERCEPT(Conf,a,Analytes)),0);
				dual[a] = (d_d < PulseThreshold) ? std::max((double)PULSE(imd,a), d_d) : d_d;
			}
		}

		double pushTotal() const {
			double* dual = Buf + (Push & BufMask) * Analytes;
			double acc = 0.;
			for (size_t a=0; a<Analytes; a++)
				acc += dual[a];
			return acc;
		}

		void integratePushes(double* integral, size_t begin, size_t end) const {
			assert((Push - begin) < BufMask && end < Push);
			// Integrate by analyte across specified pushes
			for (size_t p=begin; p<end; p++) {
				double* dual = Buf + (p & BufMask) * Analytes;
				for (size_t a=0; a<Analytes; a++) {
					integral[a] += dual[a]; 
				}
			}
		}

	};

	class Smoother {
	public:
		const double* Coeffs;
		size_t		  Len;

		double		  *RawBuf, *OutBuf;
		size_t        BufPtr;
		size_t        BufMask;

		int64_t       Push;

	public:
		Smoother(const double* coeffs, size_t len, size_t max_samples=256) : Coeffs(coeffs), Len(len) {
			// Create a circular buffer for filter
			size_t p = std::ceil(std::log(max_samples)/std::log(2.));  // Next largest power of two
			
			RawBuf = Calloc(1<<p, double);
			OutBuf = Calloc(1<<p, double);
			BufPtr = 0;
			BufMask = ~(std::numeric_limits<size_t>::max() << p);

			Push = -lag();
		}
		~Smoother() { Free(RawBuf); Free(OutBuf); }

		size_t lag() const { return Len/2; }

		double fromSample(double new_sample) {
			RawBuf[BufPtr & BufMask] = new_sample;			
			double acc = 0.0;
			for (size_t s=0, p=(BufPtr+1)-Len; s<Len; s++)
				acc += RawBuf[(p+s) & BufMask] * Coeffs[s];
			OutBuf[BufPtr & BufMask] = acc;
			return acc;
		}

		int64_t getPush() const { return Push; }
		void advPush() { ++BufPtr; ++Push; }

		class Iterator {
			double *Buf;
			size_t Mask, P;
		public:
			Iterator(double *buf, size_t mask, size_t p) : Buf(buf), Mask(mask), P(p) {}
			double   operator* () const { return *(Buf + (P & Mask)); }
			bool     operator!=(const Iterator& rhs) { return P != rhs.P; }
			Iterator operator++() { ++P; return *this; } 
			Iterator operator+(size_t i) const { return Iterator(Buf, Mask, P+i); }
		};

		Iterator iterator(size_t idx) const { return Iterator(OutBuf, BufMask, (size_t)((int64_t)idx-Push+(int64_t)BufPtr)); }

	};

	enum ExtractState {
		InNoise,
		ACell
	};

	class Observation {
	public:
		static size_t Analytes;
		static void Observation_global(size_t analytes) { Analytes = analytes; }

	public:
		size_t Start, End; // [Start, End)
		double *Obs;
		double Quality;

		Observation(size_t start) : Start(start), End(start), Obs(NULL), Quality(1.) {}
		~Observation() { if (Obs) Free(Obs); }

		void reset(size_t start) { Start = start; std::fill(Obs,Obs+Analytes,0.); }

		void initObs() { if (!Obs) Obs = Calloc(Analytes, double); }
		double* getObs() { this->initObs(); return Obs; }

		void setStart(size_t start) { Start = start; }
		size_t getStart() const { return Start; }

		void setEnd(size_t end) { End = end; }
		size_t getEnd() const { return End; }

		double getQuality() const { return Quality; }
		void setQuality(double q) { Quality = q; }

		size_t length() const { return End - Start; }
	};
	size_t Observation::Analytes;

	bool length_filter(Observation* cell, const Smoother& smoother, const Dueler& dueler, const Options& opt) {
		return cell->length() >= opt.min_cell_length && cell->length() <= opt.max_cell_length;
	}

	bool slope_filter(Observation* cell, const Smoother& smoother, const Dueler& dueler, const Options& opt) {
		assert(cell->getEnd() > 2);
		size_t cross = 0;			
		for (Smoother::Iterator b=smoother.iterator(cell->getStart()), e=smoother.iterator(cell->getEnd()-2); b!=e; ++b) {
			double sm0 = *b, sm1 = *(b+1), sm2 = *(b+2);
			double slope1 = sm1 - sm0, slope2 = sm2 - sm1;
			if (sign(slope2) != 0 && sign(slope1) != sign(slope2))
				cross++;
		}
		if (cross > 1)
			cell->setQuality(0.);
		return true;  // We don't block cells with multiple crossing, but do mark them as low quality
	}

	void
	subtract_noise(Observation* cell, Observation* noise) {
		double *cell_d = cell->getObs(), *noise_d = noise->getObs();
		for (size_t a=0; a<Observation::Analytes; a++)
			cell_d[a] -= (noise_d[a] / noise->length() * cell->length()); 
	}
	
	void
	extract(int fh, const Options& opt, std::vector<Observation*>& cells) {

		// Initialize state of extraction
		Dueler   dueler(opt.conf, opt.analytes, opt.pulse_threshold);
		Smoother smoother(opt.smooth, opt.smooth_len); 
		ExtractState state = InNoise;
	
		Observation::Observation_global(opt.analytes);
		Observation *cell = NULL;
		Observation *noise_act = new Observation(0); noise_act->initObs();  // Initialize "initial" baseline noise
		Observation *noise_tmp = new Observation(0); noise_tmp->initObs();  // Initialize "in progress" noise

		// Initialize buffers to read from the disk in "opt.blk_size" chunks 
		size_t   read_buf_len = opt.blk_size + 4*opt.analytes;
		uint8_t *read_buf = Calloc(read_buf_len, uint8_t), 
				*read_end = read_buf + read_buf_len, 
				*read_ptr = read_buf;

		while (dueler.getPush() < opt.num_pushes) {
	
			assert(read_ptr + opt.blk_size < read_buf + read_buf_len);
			ssize_t bytes_read = read(fh, read_ptr, opt.blk_size); 
			if (bytes_read < 0)
				error("File read failed\n");
			else if (bytes_read == 0)
				break; // Break out of extraction loop
			read_end = read_ptr + bytes_read;

			uint16_t *imd_push; 
			for (imd_push = (uint16_t*)read_buf; (imd_push+2*opt.analytes)<=(uint16_t*)read_end; imd_push+=2*opt.analytes) {

				dueler.fromIMD(imd_push);			
				double sm = smoother.fromSample( dueler.pushTotal() );  
			
				if (smoother.getPush() >= 0) {
					switch (state) {
						case InNoise: {
							if (sm >= opt.cell_threshold) { // We have found a potential cell ...
								assert(cell == NULL);
								
								state = ACell;  // Transition to "cell" state
								cell = new Observation(smoother.getPush());

								// Noise filters
								noise_tmp->setEnd(smoother.getPush());
								if (noise_tmp->length() >= opt.min_noise_length) {
									std::swap(noise_tmp, noise_act);
								}	
																
							} else
								dueler.integratePushes(noise_tmp->getObs(), smoother.getPush(), smoother.getPush()+1);
							break;
						}
						case ACell: {
							// We have found the end of a cell ...
							if (sm < opt.cell_threshold) {
								assert(cell);

								state = InNoise;  // Transition to "noise" state
								noise_tmp->reset(smoother.getPush());
								dueler.integratePushes(noise_tmp->getObs(), smoother.getPush(), smoother.getPush()+1);
								
								// Cell filters
								cell->setEnd(smoother.getPush());
								if (!length_filter(cell, smoother, dueler, opt) ||
									(opt.slope_filter && !slope_filter (cell, smoother, dueler, opt)))
									delete cell;
								else {
									dueler.integratePushes(cell->getObs(), cell->getStart(), cell->getEnd());
									// Subtract out most recently valid noise
									if (opt.noise_subtraction && noise_act->length() > 0)
										subtract_noise(cell, noise_act);
									cells.push_back(cell);
								}
								cell = NULL;
								break;						
							}
						}
					}  // switch(state) ...
				}  // if (smoother.getPush() ...

				smoother.advPush();
				dueler.advPush();
			}

			{  // Shift remaining data to beginning of buffer
				size_t remainder = read_end - (uint8_t*)imd_push;
				read_ptr = (uint8_t*)memcpy(read_buf, imd_push, remainder) + remainder;
			}
		}

		Free(read_buf);
	}


}

extern "C" {

    SEXP
	CC_extract(SEXP file, SEXP conf, SEXP options) {

		Options opt;
		SEXP arg;

		/*
		 * Required arguments
		 */

		const char* file_name = CHAR(STRING_ELT(file,0));

		if (conf) {
			opt.conf = REAL(conf);
			opt.analytes = static_cast<size_t>(INTEGER(GET_DIM(conf))[0]);
		} else
			error("Must supply conf matrix\n");
	
		if (NULL_USER_OBJECT != (arg = getListElement(options, "smooth"))) {
			opt.smooth = REAL(arg);
			opt.smooth_len = static_cast<size_t>(GET_LENGTH(arg));
		} else 
			error("Must supply smoothing filter\n");

		/*
		 * Optional parameters
		 */
	
#define ARG( NAME, OPT, SET, DEFAULT) \
		if (NULL_USER_OBJECT != (arg = getListElement(options, NAME))) { \
			opt.OPT = SET; \
		} else { \
			opt.OPT = DEFAULT; \
		}
		
		ARG("pulse.threshold", pulse_threshold, *(REAL(arg)), 3.0);	
		ARG("cell.threshold", cell_threshold, *(REAL(arg)), 10.);
		ARG("cell.min.length", min_cell_length, (size_t)*(INTEGER(arg)), 10);
		ARG("cell.max.length", max_cell_length, (size_t)*(INTEGER(arg)), 75);
		ARG("noise.subtraction", noise_subtraction, *(LOGICAL(arg)), true);
		ARG("noise.min.length", min_noise_length, (size_t)*(INTEGER(arg)), 30);
		ARG("slope.filter", slope_filter, *(LOGICAL(arg)), false);
		ARG("num.pushes", num_pushes, (size_t)*(INTEGER(arg)), std::numeric_limits<size_t>::max());

	
#undef ARG

		/*
		 * Open file for binary reading
		 */

		int fh;
		{
			fh = open(file_name, O_RDONLY);
			if (fh == -1)
				error("Error %d: open failed on file %s\n",errno,file_name);
		}

		{  // Extract the optimal block size for reading from this storage device
			struct stat stat_buf;
			if (fstat(fh, &stat_buf) != 0)
				error("Error %d: fstat failed on file %s\n",errno,file_name);
			opt.blk_size = stat_buf.st_blksize;
		}


		/*
		 * Extract cells from IMD file
		 */ 
		
		std::vector<Observation*> cells;

		extract(fh, opt, cells);

		close(fh);


		/*
		 * Build return R data structure
		 */

		SEXP ans, names, cells_extent, cells_obs, cells_qual; int n_protected = 0;

		PROTECT(ans = allocVector(VECSXP,3)); n_protected++;
		PROTECT(names = allocVector(STRSXP,3)); n_protected++;

		size_t num_cells = cells.size();
		
		PROTECT(cells_extent = allocMatrix(INTSXP, num_cells, 2)); n_protected++;
		PROTECT(cells_obs = allocMatrix(REALSXP, num_cells, opt.analytes)); n_protected++;
		PROTECT(cells_qual = allocVector(REALSXP, num_cells)); n_protected++;
		for (size_t i=0; i<num_cells; i++) {
			Observation *c = cells[i];
			// Recall R data structures are column major
			INTEGER(cells_extent)[i] = (int)(c->getStart());
			INTEGER(cells_extent)[i+num_cells] = (int)(c->length());
			for (size_t a=0; a<opt.analytes; a++)
				REAL(cells_obs)[a*num_cells+i] = c->getObs()[a];
			REAL(cells_qual)[i] = c->getQuality();
			// Destroy cell now that we're done
			delete c;
		}

		SET_VECTOR_ELT(ans, 0, cells_extent);
		SET_STRING_ELT(names, 0, mkChar("cells.extent"));
	
		SET_VECTOR_ELT(ans, 1, cells_obs);
		SET_STRING_ELT(names, 1, mkChar("cells.obs"));

		SET_VECTOR_ELT(ans, 2, cells_qual);
		SET_STRING_ELT(names, 2, mkChar("cells.quality"));

		setAttrib(ans, R_NamesSymbol, names);
		UNPROTECT(n_protected);

		return ans;
	}

}
