#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

namespace {

    /* get the list element named str, or return NULL */
     
    SEXP getListElement(SEXP list, const char *str) {
	SEXP elmt = NULL_USER_OBJECT, names = getAttrib(list, R_NamesSymbol);
	for (int i = 0; i < length(list); i++)
	    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
		elmt = VECTOR_ELT(list, i);
	    break;
	}
	return elmt;
    }


}


extern "C" {

    SEXP
    CC_extract(SEXP file, SEXP conf, SEXP options) {

	const char* file_name = CHAR(STRING_ELT(file,0));
	
	size_t analytes = static_cast<size_t>(INTEGER(GET_DIM(conf))[0]);
	Rprintf("Analytes: %zu\n",analytes);

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
	    fh = open(file_name, O_RDONLY|O_LARGEFILE);
	    if (fh == -1)
		error("Error %d: open failed on file %s\n",errno,file_name);
	}

	blksize_t in_blk_size;
	{
	    struct stat stat_buf;
	    if (fstat(fh, &stat_buf) != 0)
		error("Error %d: fstat failed on file %s\n",errno,file_name);
	    in_blk_size = stat_buf.st_blksize;
	}
	Rprintf("Block size: %zu\n",in_blk_size);

	close(fh);
	return (NULL_USER_OBJECT);

    }

}
