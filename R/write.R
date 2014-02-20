cytofCore.write.FCS <- function(x, filename, what="numeric",channelDescriptions=NULL) {
    if (is.matrix(x)) {
		# Don't write "Leading_Push" to FCS files		
		x <- x[,which(colnames(x) != "Leading_Push")]

		# Build metadata for FCS file
		pd <- c()  # 'params' phenoData
		dl <- list()  # 'description' list
		
		dl[["$DATATYPE"]] <- "F"
		for (c in 1:ncol(x)) {
		    c_name <- colnames(x)[c]
		    if (!is.null(channelDescriptions)){
          if (!is.na(channelDescriptions[c])){
            c_name <- channelDescriptions[c]  
          }
		    }
		    
		    c_min <- min(0,min(x[,c]))  # Hack to prevent flowCore from shifting range
		    c_max <- max(x[,c])
		    c_rng <- c_max - c_min + 1

		    pl <- matrix(c(c_name, c_name, c_rng, c_min, c_max),nrow=1)
		    colnames(pl) <- c("name", "desc", "range", "minRange", "maxRange")
		    rownames(pl) <- paste("$P",c,sep="") 
		    pd <- rbind(pd, pl)
		
		    dl[[paste("$P",c,"B",sep="")]] <- "32";	    # Number of bits
		    dl[[paste("$P",c,"R",sep="")]] <- toString(c_rng); # Range
		    dl[[paste("$P",c,"E",sep="")]] <- "0,0";	    # Exponent
		    dl[[paste("$P",c,"N",sep="")]] <- c_name;	    # Name
		    dl[[paste("$P",c,"S",sep="")]] <- c_name;	    # Desc	
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
