cytofCore.read.conf <- function(file) {
    read.table(file,header=TRUE)
}

cytofCore.read.imd.xml <- function(file) {
  # extract xml from tail of imd
  imd <- file(file, "rb")
  on.exit(close(imd))
  if (!isSeekable(imd)) {
    stop("Cannot seek in specified file or connection")
  }
  
  endTag="</ExperimentSchema>"
  seek(imd,where=-2*nchar(endTag),origin="end")
  fileEndString=paste(readBin(imd, what="character", size=2, n=2*nchar(endTag), signed=FALSE),collapse="")
  
  if (fileEndString != endTag) {
    stop("The xml tail is either missing or irregular.")
  }
  
  #read in chunks until the beginning of the xml is reached
  chunkLength=1024
  seek(imd,where=-2*chunkLength,origin="current")
  xmlChunk=c()
  while (!(0 %in% xmlChunk)) {
    xmlChunk=c(readBin(imd, what="integer", size=2, n=chunkLength, signed=FALSE),xmlChunk)
    seek(imd,where=-4*chunkLength,origin="current")
  }
  
  # convert to char and trim extra before the xml
  imdString=sub(".*<ExperimentSchema","<ExperimentSchema",rawToChar(as.raw(xmlChunk[which(xmlChunk!=0)])))
  
  # parse xml
  xmlList=xmlToList(xmlInternalTreeParse(imdString))
  
  # analyte info
  analyteList=xmlList[names(xmlList)=="AcquisitionMarkers"]
  analytes=c()
  for (metal in analyteList) {
    analytes=rbind(analytes,c(metal$Mass,metal$Description,metal$MassSymbol ))
  }
  analytes[,3]= paste(analytes[,3],as.character(round(as.numeric(analytes[,1]))),sep="")
  
  colnames(analytes)=c("mass","description","channel")
  
  # dual calibration info
  dualList=xmlList[names(xmlList)=="DualAnalytesSnapshot"]
  dualCalibration=c()
  for (metal in dualList) {
    dualCalibration=rbind(dualCalibration,as.numeric(c(metal$Mass,metal$DualSlope, metal$DualIntercept)))
  }
  colnames(dualCalibration)=c("mass","slope","intercept")
  
  imd.xml=list()
  imd.xml$analytes=analytes
  imd.xml$dualCalibration=dualCalibration
  
  return(imd.xml)
  
}

cytofCore.read.imd <- function(file, analytes=NULL, conf=NULL, pulse_thresh=3.0, start_push=0, num_pushes=NULL) {
    
    
    
    # change to allow input analytes for when imd is missing tail
    if (!is.null(analytes)) {
		num_analytes <- length(analytes)
    } else if (!is.null(conf)) {
		if (is.character(conf) && file.exists(conf)) {
			conf <- cytofCore.read.conf(conf)
		}
		num_analytes <- nrow(conf)
		analytes <- paste(round(conf$Mass),conf$Symbol,sep="")
    } else {
      # Change this to read xml from tail of IMD file.
		stop("One of analytes or conf must be specified to know how many analytes are present")
    }

    if (is.character(file)) {
		if (is.null(num_pushes) && !is.na(file.info(file)$size)) {
			# File size should be multiple of 2 * # analytes * 2 bytes
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
		# If you supplied conf file, we assume you only want the dual counts
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
