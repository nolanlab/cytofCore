cytofCore.read.conf <- function(file) {
  read.table(file,header=TRUE)
}

cytofCore.read.imd.xml <- function(file) {
  # file is the name of an IMD file
  # this function extracts the XML tail from the IMD file
  # returns a list with elements "analytes", "dualCalibration" and "rawText"
  # analytes is a data frame with columns of the mass, description and channel
  # dualCalibration is a data frame with columns of the dual mass, slope and intercept
  # rawText is a character string of the entire XML tail
  
  library("XML")
  
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
    if (length(xmlChunk)>2^17) {stop("Didn't find zeros preceding xml.")}
  }

  # convert to char and trim extra before the xml
  startTag="<ExperimentSchema"
  matchInd=match(utf8ToInt("<"),xmlChunk)
  if (is.na(matchInd)) {stop("XML start tag not found")}
  
  imdString=sub(".*<ExperimentSchema","<ExperimentSchema",intToUtf8(xmlChunk[matchInd:length(xmlChunk)]))
  
  if (substr(imdString,1,nchar(startTag)) != startTag) {stop("XML start not aligned.")}
  
  # parse xml
  xmlList=xmlToList(xmlInternalTreeParse(imdString))
  
  # analyte info
  analyteList=xmlList[names(xmlList)=="AcquisitionMarkers"]
  analytes=c()
  for (metal in analyteList) {
    if is.null(metal$Description) {
      metal$Description = paste0(metal$MassSymbol, as.character(round(metal$Mass)))
    }
    analytes=rbind(analytes,c(metal$Mass,metal$Description,metal$MassSymbol))
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
  imd.xml$analytes=data.frame(analytes)
  imd.xml$dualCalibration=data.frame(dualCalibration)
  imd.xml$rawText=imdString
  
  return(imd.xml)
  
}

cytofCore.read.imd <- function(file, conf=NULL, pulse_thresh=1.0, start_push=0, num_pushes=2^16) {
  # returns a list with elements "intensity", "pulse" and "dual" which each have rows corresponding to pushes and columns corresponding to analytes
  # the dual counts are computed using the Di slopes and intercepts from the XML tail of the IMD, unless a conf file is provided, in which case the Dd slopes and intercepts are used
  # pulse_thresh is the dual count start value
  
  # change to allow input analytes from conf file
  if (!is.null(conf)) {
    if (is.character(conf) && file.exists(conf)) {
      conf <- cytofCore.read.conf(conf)
    }
    
    analytes <- paste(conf$Symbol,round(conf$Mass),sep="")
    slopes=conf$Slope
    intercepts=conf$Intercept
  } else {
    xml=cytofCore.read.imd.xml(file)
    analytes=xml$analytes$description
    masses=as.numeric(as.character(xml$analytes$mass))
    slopes=approx(xml$dualCalibration$mass,xml$dualCalibration$slope,masses,method="linear",rule=2)$y
    intercepts=approx(xml$dualCalibration$mass,xml$dualCalibration$intercept,masses,method="linear",rule=2)$y
  } 
  
  num_analytes <- length(analytes)
  
  if (is.character(file)) {
    # 		if (is.null(num_pushes) && !is.na(file.info(file)$size)) {
    # 			# File size should be multiple of 2 * # analytes * 2 bytes
    # 			num_pushes <- as.integer(file.info(file)$size / num_analytes / 4)
    # 		}
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
  
#   while (is.null(num_pushes) || current_push < num_pushes) { 
#     desired_n <- ifelse(is.null(num_pushes), 2^16, num_pushes-current_push) * num_analytes * 2
#     if (exists("N")) {
#       before_n <- length(N)
#       N <- c(N,readBin(file, integer(), size=2, n=desired_n, signed=FALSE))  # Records are 16-bit unsigned integers
#       actual_n <- length(N)-before_n
#     } else {
#       N <- readBin(file, integer(), size=2, n=desired_n, signed=FALSE)  # Records are 16-bit unsigned integers
#       actual_n <- length(N)
#     }
#     
#     if (actual_n < desired_n) {
#       break;
#     }
#     current_push <- as.integer(length(N) / num_analytes / 2)
#   }

  #simplifying: just read in num_pushes number of pushes
  N=readBin(file,integer(),size=2,n=num_pushes*2*num_analytes,signed=F)
  
  I_cols <- seq(from=1, by=2, length.out=num_analytes)
  P_cols <- seq(from=2, by=2, length.out=num_analytes)  
  IP <- matrix(N,ncol=2*num_analytes,byrow=TRUE) 
  
  # 20140709 updated dual conditions to those used by CyTOF software
  D <- c()
  for (i in 1:num_analytes) {
    d <- round(IP[,I_cols[i]]*slopes[i]+intercepts[i])
    use_pulse = d < IP[,P_cols[i]] & IP[,P_cols[i]]<= pulse_thresh 
    d[use_pulse] <- IP[use_pulse,P_cols[i]]
    D <- cbind(D, d)
  }
  colnames(D) <- analytes
  
  # Changing from earlier version and calculating I and P even if conf supplied     
  l <- list(intensity=IP[,I_cols],pulse=IP[,P_cols],dual=D)
  colnames(l$intensity) <- analytes
  colnames(l$pulse)     <- analytes
  
  return(l)
}

cytofCore.read.rd8.cytof1 <- function(file, segmentSize=NULL,  massCalibration=NULL, start_push=0, num_pushes=2^8) {
  #returns a list with raw intensities where each row is a push and each column is the tof trace of that push.
  #optionally returns which column corresponds to which mass
  #segmentSize defaults to 3200 if not specified but if this is incorrect the results won't make sense
  #including the optional input argument of the massCalibration adds the massPeak locations to the list
  #the massCalibration parameter is a list with names "time", "mass", and "triggerDelay". 
  #"time" and "mass" are 2-element vectors whose values are visible in the CyTOF software mass calibration window
  # Example:  massCalibration=list(mass=c(132.905,192.963),time=c(9597,11537),triggerDelay=8416)

  if (is.character(file)) {
    file <- file(file, "rb")
    on.exit(close(file))
  }
  
  results=list()
  
  # cytof1 has 1 uint8 value per intensity sample
  
    if(is.null(segmentSize)) {
      segmentSize=3200 }
    seek(file, where=start_push*segmentSize, origin="start")
    N=readBin(file,integer(),size=1,n=num_pushes*segmentSize,signed=F)
    intensity <- matrix(N,ncol=segmentSize,byrow=TRUE,dimnames=list(as.character(start_push:(start_push+num_pushes-1)))) 
    
  results$intensity=intensity  
  
  #calibrate the mass windows to the raw data
  if (!is.null(massCalibration)) { 
    A=(massCalibration$time[1]-massCalibration$time[2])/(sqrt(massCalibration$mass[1]) - sqrt(massCalibration$mass[2]))
    T0=massCalibration$time[1]-A*sqrt(massCalibration$mass[1])
    masses=102:195
    massPeaks=T0+A*sqrt(masses) - massCalibration$triggerDelay;
    names(massPeaks)=as.character(masses)
    
    # the mass peaks correspond to the columns of the intensities
    results$massPeaks=massPeaks 
  }

  return(results)
    
}

cytofCore.read.rd8.cytof2 <- function(file, referenceFCS, start_push=0, num_pushes=2^8) {
  #returns a list with two tables of raw intensities, where each row is a push and each column is the tof trace of that push.
  #also returns which column corresponds to which mass
  #the mass calibration values are read in from the xml after the data segment in a cytof2 fcs file acquired with the same settings as the rd8
  
  if (is.character(file)) {
    file <- file(file, "rb")
    on.exit(close(file))
  }
  
  results=list()
  
  # read segment size, trigger delay, A and T0 from an fcs file aquired with the same settings
  tofCalibration=cytofCore.read.tof.calibration(referenceFCS)
  segmentSize=tofCalibration$segmentSize
  # cytof2 has 2 uint8 values for each sample
 

    seek(file, where=start_push*segmentSize*2, origin="current")
    N=readBin(file,integer(),size=1,n=num_pushes*segmentSize*2,signed=F)
    board1_cols <- seq(from=1, by=2, length.out=segmentSize)
    board2_cols <- seq(from=2, by=2, length.out=segmentSize)  
    bothBoards <- matrix(N,ncol=2*segmentSize,byrow=TRUE,dimnames=list(as.character(start_push:(start_push+num_pushes-1))))  
    intensity=list(board1=bothBoards[,board1_cols], board2=bothBoards[,board2_cols] )
  
  
  results$intensity=intensity  
  
  #calibrate the mass windows to the raw data
  
    masses=89:195
    massPeaks=tofCalibration$T0+tofCalibration$A*sqrt(masses) - tofCalibration$triggerDelay;
    names(massPeaks)=as.character(masses)
    
    # the mass peaks correspond to the columns of the intensities
    results$massPeaks=massPeaks 
  
  return(results)
  
}

cytofCore.read.tof.calibration <- function(filename) {
  # finds the segment size, trigger delay, A and T0 from the xml written after the data segment of a Cytof2 fcs file
  
  if (is.character(file)) {
    file <- file(filename, "rb")
    on.exit(close(file))
  }
  
  seek(file, where=10, origin="start")
  headerStartPos=as.integer(intToUtf8(readBin(file,"integer",n=8,size=1,signed=FALSE)))
  headerStopPos=as.integer(intToUtf8(readBin(file,"integer",n=8,size=1,signed=FALSE)))
  
  seek(file, where=headerStartPos,origin="start")
  fcsHeader=intToUtf8(readBin(file,"integer",n=headerStopPos-headerStartPos+1,size=1,signed=FALSE))
  delimiter=substr(fcsHeader,1,1)
  
  # parsing of header to get enddata position
  pattern=paste0(".*",delimiter,"\\$ENDDATA",delimiter,"(\\d+)",delimiter,".*")
  endPos=sub(pattern,"\\1",fcsHeader)
  
  # move to position just after end of data
  seek(file,where=as.integer(endPos)+1,origin="start")
  
  #read in chunks until the end of the xml is reached    
  chunkLength=1024
  xmlChunk=intToUtf8(readBin(file,what="integer",size=1,n=2*chunkLength,signed=F))
  while (length(grep("</FCSHeaderSchema>",substr(xmlChunk,nchar(xmlChunk)-2*chunkLength+1,nchar(xmlChunk))))==0) {
    xmlChunk=paste0(xmlChunk,intToUtf8(readBin(file, what="integer", size=1, n=chunkLength, signed=FALSE)))
    if (nchar(xmlChunk)>2^17) {stop("Didn't find end of xml.")}
  }
  
  # trim to just the xml
  xml=sub("(.*</FCSHeaderSchema>).*","\\1",xmlChunk)
  
  result=list()
  
  # extract the tof calibration values from the xml
  xmlList=xmlToList(xml)
  acquisitionBoardParams=xmlList[["Acquisition"]][["AcquisitionBoardParams"]]
  boardList=xmlToList(acquisitionBoardParams)
  result$triggerDelay=as.integer(boardList[["Params"]][["TriggerDelay"]])
  result$segmentSize=as.integer(boardList[["Params"]][["SegmentSize"]])
  result$A=as.numeric(xmlList[["CalibrationResults"]][["A"]])
  result$T0=as.numeric(xmlList[["CalibrationResults"]][["T0"]])
  
  return(result)
}