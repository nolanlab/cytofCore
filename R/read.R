cytofCore.read.conf <- function(file) {
    read.table(file,header=TRUE)
}

cytofCore.read.imd <- function(file, analytes=NULL, conf=NULL, pulse_thresh=3.0, start_push=0, num_pushes=NULL) {
    if (!is.null(analytes)) {
	num_analytes <- length(analytes)
    } else if (!is.null(conf)) {
	num_analytes <- nrow(conf)
	analytes <- paste(round(conf$Mass),conf$Symbol,sep="")
    } else {
	stop("One of analytes or conf must be specified to know how many analytes are present")
    }

    if (is.character(file)) {
	file <- file(file, "rb")
	on.exit(close(file))
    }
    if (!inherits(file, "connection")) {
	stop("'file' must be a character string or connection")
    }
    if (!isOpen(file, "rb")) {
	open(file, "rb")
	on.exit(close(file))
    }

    # Read the IMD file in chunks
    if (start_push > 0) {
	if (!isSeekable(file)) {
	    stop("Cannot seek in specified file or connection")
	}
	# Each analyte consumes 4 bytes, 2 each for intensity and pulse values
	seek(file, where=start_push*num_analytes*4, origin="current")
    }
   
    I <- matrix(nrow=0,ncol=num_analytes)
    P <- matrix(nrow=0,ncol=num_analytes)
 
    colnames(I) <- analytes
    colnames(P) <- analytes

    I_cols <- seq(from=1, by=2, length.out=num_analytes)
    P_cols <- seq(from=2, by=2, length.out=num_analytes)
    
    current_push = 0
    while (is.null(num_pushes) || current_push < num_pushes) { 
	desired_n <- min(1024, num_pushes - current_push) * num_analytes * 2 
	
	n <- readBin(file, integer(), size=2, n=desired_n)		     # Records are 16-bit unsigned integers
	
	IP <- matrix(n,ncol=2*num_analytes,byrow=TRUE)
	I  <- rbind(I, IP[,I_cols])
	P  <- rbind(P, IP[,P_cols])
	
	if (length(n) < desired_n) {
	    break;
	}
	current_push <- current_push + nrow(IP)
    }

    if (is.null(conf)) {
	return (list(intensity=I,pulse=P))
    } else {	
	# If conf file is provided, compute the dual counts
	D <- c()
	for (i in 1:num_analytes) {
	    d <- round(I[,i]*conf$Slope[i]+conf$Intercept[i])
	    use_pulse = d < pulse_thresh
	    d[use_pulse] <- P[use_pulse,i]
	    D <- cbind(D, d)
	}
	colnames(D) <- analytes
     
	return (list(intensity=I,pulse=P,dual=D))
    }
}
