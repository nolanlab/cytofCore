cytofCore.read.conf <- function(file) {
    read.table(file,header=TRUE)
}

cytofCore.read.imd <- function(file, analytes=NULL, conf=NULL, pulse_thresh=3.0, start_push=0, num_pushes=NULL) {
    if (!is.null(analytes)) {
		num_analytes <- length(analytes)
    } else if (!is.null(conf)) {
		if (is.character(conf) && file.exists(conf)) {
			conf <- cytofCore.read.conf(conf)
		}
		num_analytes <- nrow(conf)
		analytes <- paste(round(conf$Mass),conf$Symbol,sep="")
    } else {
		stop("One of analytes or conf must be specified to know how many analytes are present")
    }

    if (is.character(file)) {
		if (is.null(num_pushes) && !is.na(file.info(file)$size)) {
			num_pushes <- as.integer(file.info(file)$size / num_analytes / 4)
		}
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
    current_push <- 0
	if (start_push > 0) {
		if (!isSeekable(file)) {
			stop("Cannot seek in specified file or connection")
		}
		# Each analyte consumes 4 bytes, 2 each for intensity and pulse values
		seek(file, where=start_push*num_analytes*4, origin="current")
    }
      
	while (is.null(num_pushes) || current_push < num_pushes) { 
		desired_n <- ifelse(is.null(num_pushes), 2^16, num_pushes-current_push) * num_analytes * 2
		if (exists("N")) {
			before_n <- length(N)
			N <- c(N,readBin(file, integer(), size=2, n=desired_n, signed=FALSE))  # Records are 16-bit unsigned integers
			actual_n <- length(N)-before_n
		} else {
			N <- readBin(file, integer(), size=2, n=desired_n, signed=FALSE)  # Records are 16-bit unsigned integers
			actual_n <- length(N)
		}
                
		if (actual_n < desired_n) {
            break;
		}
		current_push <- as.integer(length(N) / num_analytes / 2)
    }
  
    I_cols <- seq(from=1, by=2, length.out=num_analytes)
    P_cols <- seq(from=2, by=2, length.out=num_analytes)  
    IP <- matrix(N,ncol=2*num_analytes,byrow=TRUE) 
    
    if (is.null(conf)) {
		l <- list(intensity=IP[,I_cols],pulse=IP[,P_cols])
		colnames(l$intensity) <- analytes
		colnames(l$pulse)     <- analytes
		return(invisible(l))
    } else {
		D <- c()
		for (i in 1:num_analytes) {
			d <- round(IP[,I_cols[i]]*conf$Slope[i]+conf$Intercept[i])
			use_pulse = d < pulse_thresh
			d[use_pulse] <- pmax(d[use_pulse],IP[use_pulse,P_cols[i]])
			D <- cbind(D, d)
		}
		colnames(D) <- analytes
		return(invisible(list(dual=D)))
    }
}
