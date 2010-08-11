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

namespace {

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

		ssize_t blk_size;
		double  pulse_threshold;
	};

	class Filter {
	public:
		const double* Coeffs;
		size_t Len;


	public:
		Filter(const double* coeffs, size_t len) : Coeffs(coeffs), Len(len) {
			
		}


	};


#define SLOPE(CONF, A, ANA) (CONF[(A)+ANA])
#define INTERCEPT(CONF, A, ANA) (CONF[(A)])

#define INTENSITY(PUSH, A) (PUSH[2*(A)])
#define PULSE(PUSH, A) (PUSH[2*(A)+1])
#define DUAL(PUSH, A) (PUSH[(A)])

	void
	imdToDual(const uint16_t* imd, uint16_t* dual, const double* conf, size_t analytes, double pulse_threshold) {
		for (size_t a=0; a<analytes; a++) {	
			double   d_d = fround((SLOPE(conf, a, analytes) * INTENSITY(imd,a) + INTERCEPT(conf, a, analytes)),0);
			DUAL(dual,a) = (d_d < pulse_threshold) ? std::max(PULSE(imd,a), (uint16_t)d_d) : (uint16_t)d_d; 
		}
	}

	uint32_t
	totalIntensity(const uint16_t* dual, size_t analytes) {
		uint32_t res = 0;
		for (size_t a=0; a<analytes; a++)
			res += dual[a];
		return res;
	}


	void
	extract(int fh, const Options& opt) {

		// Setup disk read lengths
		size_t  read_buf_len = opt.blk_size + (2*opt.analytes-1) * 2;
		uint8_t *read_buf = Calloc(read_buf_len, uint8_t), 
				*read_end = read_buf + read_buf_len, 
				*read_ptr = read_buf;
	
		uint16_t *dual_buf = Calloc(read_buf_len / (2 * 2), uint16_t);

		ssize_t bytes_read = read(fh, read_ptr, opt.blk_size); 
		if (bytes_read < 0) {
			error("File read failed\n");
		}
		read_end = read_ptr + bytes_read;

		uint16_t *imd_push = (uint16_t*)read_ptr, *dual_push = dual_buf;	
		while ((imd_push+2*opt.analytes) <= (uint16_t*)read_end) {

			imdToDual(imd_push, dual_push, opt.conf, opt.analytes, opt.pulse_threshold);
	
			uint32_t ti = totalIntensity(dual_push, opt.analytes);

			imd_push+=2*opt.analytes, dual_push+=opt.analytes;
		}
	
		{  // Shift remaining data to beginning of buffer
			size_t remainder = read_end - (uint8_t*)imd_push;
			read_ptr = (uint8_t*)memcpy(read_buf, imd_push, remainder) + remainder;
		}

		Free(read_buf);
		Free(dual_buf);
	}

}

extern "C" {

    SEXP
	CC_extract(SEXP file, SEXP conf, SEXP options) {

		Options opt;

		const char* file_name = CHAR(STRING_ELT(file,0));

		opt.conf = REAL(conf);
		
		opt.analytes = static_cast<size_t>(INTEGER(GET_DIM(conf))[0]);
	
		Rprintf("Analytes: %zu\n",opt.analytes);

		opt.pulse_threshold = 3.0;

		double *filter = NULL;
		size_t  filter_len = 0;
		{
			SEXP filter_s = getListElement(options, "filter");
			if (filter_s != NULL_USER_OBJECT) {
				filter = REAL(filter_s);
				filter_len = static_cast<size_t>(GET_LENGTH(filter_s));
				Rprintf("Filter Length: %zu\n",filter_len);
			}
		} 

		int fh;
		{
			fh = open(file_name, O_RDONLY);
			if (fh == -1)
				error("Error %d: open failed on file %s\n",errno,file_name);
		}

		{
			struct stat stat_buf;
			if (fstat(fh, &stat_buf) != 0)
				error("Error %d: fstat failed on file %s\n",errno,file_name);
			opt.blk_size = stat_buf.st_blksize;
		}
		Rprintf("Block size: %zu\n",opt.blk_size);

		extract(fh, opt);
	
		close(fh);
		return (NULL_USER_OBJECT);

	}

}
