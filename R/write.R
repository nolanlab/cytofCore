cytofCore.write.FCS <- function(x, filename, what="numeric",channelDescriptions=NULL,referenceDescription=NULL, oldDescription=NULL) {
  if (is.matrix(x)) {
    # Don't write "Leading_Push" to FCS files		
    x <- x[,which(colnames(x) != "Leading_Push")]
    
    # Build metadata for FCS file
    pd <- c()  # 'params' phenoData
    dl <- list()  # 'description' list
    
    dl[["$DATATYPE"]] <- "F"
    
    if (!is.null(referenceDescription)) {
      if (!is.null(referenceDescription[["$DATE"]])){
        dl[["$DATE"]] <- referenceDescription[["$DATE"]]
      }
      if (!is.null(referenceDescription[["$BTIM"]])){
        dl[["$BTIM"]] <- referenceDescription[["$BTIM"]]
      }
      if (!is.null(referenceDescription[["$ETIM"]])){
        dl[["$ETIM"]] <- referenceDescription[["$ETIM"]]
      }
      if (!is.null(referenceDescription[["$CYT"]])){
        dl[["$CYT"]] <- referenceDescription[["$CYT"]]
      }
      if (!is.null(referenceDescription[["$CYTSN"]])){
        dl[["$CYTSN"]] <- referenceDescription[["$CYTSN"]]
      }

    }
    
    for (c in 1:ncol(x)) {
      c_name <- colnames(x)[c]
      c_desc <- colnames(x)[c]
      if (!is.null(channelDescriptions)){
        if (!is.na(channelDescriptions[c])){
          c_desc <- channelDescriptions[c]  
        }
      }
      
      c_min <- floor(min(0,min(x[,c])))  # Hack to prevent flowCore from shifting range
      c_max <- ceiling(max(x[,c]))
      c_rng <- c_max - c_min + 1
      
      pl <- matrix(c(c_name, c_desc, c_rng, c_min, c_max),nrow=1)
      colnames(pl) <- c("name", "desc", "range", "minRange", "maxRange")
      rownames(pl) <- paste("$P",c,sep="") 
      pd <- rbind(pd, pl)
      
      if (!is.null(referenceDescription[[paste0("P",c,"DISPLAY")]])){
        dl[[paste0("P",c,"DISPLAY")]] <- referenceDescription[[paste0("P",c,"DISPLAY")]]
      } 
      
      if (!is.null(referenceDescription[[paste0("$P",c,"G")]])){
        dl[[paste0("$P",c,"G")]] <- referenceDescription[[paste0("$P",c,"G")]]
      } 
      
      if (!is.null(referenceDescription[[paste0("$P",c,"R")]])){
        dl[[paste0("$P",c,"R")]] <- referenceDescription[[paste0("$P",c,"R")]]
      } else {
        dl[[paste("$P",c,"R",sep="")]] <- toString(c_rng); # Range
      }
      
      if (!is.null(referenceDescription[[paste0("$P",c,"B")]])){
        dl[[paste0("$P",c,"B")]] <- referenceDescription[[paste0("$P",c,"B")]]
      } else {
        dl[[paste("$P",c,"B",sep="")]] <- "32";      # Number of bits
      }
      
      if (!is.null(referenceDescription[[paste0("$P",c,"E")]])){
        dl[[paste0("$P",c,"E")]] <- referenceDescription[[paste0("$P",c,"E")]]
      } else { 
        dl[[paste("$P",c,"E",sep="")]] <- "0,0";      # Exponent
      }
      
      dl[[paste("$P",c,"N",sep="")]] <- c_name;	    # Name
      dl[[paste("$P",c,"S",sep="")]] <- c_desc;	    # Desc	
    }	
    
    
    if (!is.null(oldDescription)) {
      if (!is.null(oldDescription[["$CYT"]])){
        dl[["$CYT"]] <- oldDescription[["$CYT"]]
      }		  
      if (!is.null(oldDescription[["$DATE"]])){
        dl[["$DATE"]] <- oldDescription[["$DATE"]]
      }
      if (!is.null(oldDescription[["$BTIM"]])){
        dl[["$BTIM"]] <- oldDescription[["$BTIM"]]
      }
      if (!is.null(oldDescription[["$ETIM"]])){
        dl[["$ETIM"]] <- oldDescription[["$ETIM"]]
      }
    }
    
    x <- flowFrame(x, as(data.frame(pd), "AnnotatedDataFrame"), description=dl) 
  }
  write.FCS(x, filename, what)
}

cytofCore.write.IMD.stream = function(imdFile,confFile,windowSize=10000,startTime=0,endTime=NULL){
  
  timeConstant = 1/76800
  continue=T
  counter=0;
  first=T
  
  #FILL THIS OUT
  num_pushes <- as.integer(file.info(imdFile)$size / nrow(cytofCore.read.conf(confFile)) / 4)
  max_time = (num_pushes-1)*timeConstant*1000
  if (startTime>max_time){
    stop(paste("Specified startTime is greater than IMD max time of",max_time))
  }
  if (!is.null(endTime)&&(endTime>max_time)){
    stop(paste("Specified endTime is greater than IMD max time of",max_time,". Leave endTime as NULL (default) to read to end of file."))
  }
  
  while (continue){
    windowStartTime = (counter)*(windowSize)*timeConstant*1000
    windowEndTime = (counter*windowSize+windowSize-1)*timeConstant*1000
    windowTime = seq(from=windowStartTime,to=windowEndTime,by=timeConstant*1000)
    #print(paste(windowStartTime,windowEndTime))
    windowKeep=F
    if (any(windowTime>=startTime)){
      windowKeep = windowTime>=startTime
    }
    
    if (!is.null(endTime)){
      windowKeep = windowKeep & (windowTime<=endTime)
      if (windowStartTime>endTime){
        continue=F	
      }
    }
    
    
    if (sum(windowKeep)>0){
      startPush=counter*windowSize
      if ((startPush+windowSize)>num_pushes){
        pushes=cytofCore.read.imd(file=imdFile,conf=confFile,start_push=startPush)$dual
        windowTime = windowTime[1:nrow(pushes)]
        windowKeep = windowKeep[1:nrow(pushes)]
        continue=F
      } else {
        pushes = cytofCore.read.imd(file=imdFile,conf=confFile,start_push=startPush,num_pushes=windowSize)$dual
      }
      pushes = cbind(time=sprintf("%.3f",windowTime[windowKeep]),pushes[windowKeep,])
      print(paste("Writing",nrow(pushes),"pushes"))
      if (first){
        write.table(pushes,file=paste(imdFile,".txt",sep=""),quote=F,row.names=F,col.names=T,append=F)
        first=F
      } else {
        write.table(pushes,file=paste(imdFile,".txt",sep=""),quote=F,row.names=F,col.names=F,append=T)	
      }
    }
    
    counter=counter+1;
    
  }
}
