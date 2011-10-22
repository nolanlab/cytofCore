cytofCore.updateFlowFrameKeywords = function(flowFrame){
  
	params = parameters(flowFrame)
 	pdata = pData(params)
	for (i in 1:ncol(flowFrame)){
		s = paste("$P",i,"S",sep="");
		n = paste("$P",i,"N",sep=""); 
		r = paste("$P",i,"R",sep=""); 
		b = paste("$P",i,"B",sep=""); 
		e = paste("$P",i,"E",sep=""); 
		keyval=list();
		keyval[[s]] = colnames(flowFrame)[i];
		keyval[[n]] = colnames(flowFrame)[i];
		keyval[[r]] = ceiling(max(exprs(flowFrame)[,i])-min(exprs(flowFrame)[,i]))
		keyval[[b]] = 32;
		keyval[[e]] = "0,0";
		keyword(flowFrame) = keyval;
	
		
	  pdata[i,"minRange"]=min(exprs(flowFrame)[,i])
		pdata[i,"maxRange"]=max(exprs(flowFrame)[,i])
	  
	}
 	pData(params)=pdata
	 parameters(flowFrame)=params
  
  # keyval[["$DATATYPE"]] <- "F"
	return(flowFrame)
}


cytofCore.combineChannels = function(flowFrame,channelList,newName=NULL){
	if (!all(channelList %in% colnames(flowFrame))){
		offending = which(!(channelList %in% colnames(flowFrame)))
		stop( paste("Channel(s)",paste(channelList[offending],collapse=","),"not found in frame.") )
	}
 	
  # Add data in combined channels
	combinedValues = apply(exprs(flowFrame)[,channelList],1,sum)
 	
 	channelNumber = ncol(flowFrame)+1
 	channelId = paste("$P",channelNumber,sep="");
  if (is.null(newName)){
  	newName = paste("channel",paste(channelList,collapse="-"),"combined",sep="_");
  }
 	
  
  params = parameters(flowFrame)
 	pdata = pData(params)
 	plist = matrix(c(newName,newName,ceiling(max(combinedValues)-min(combinedValues)),floor(min(combinedValues)),ceiling(max(combinedValues))));
 	rownames(plist) = c("name","desc","range","minRange","maxRange");
	colnames(plist) = c(channelId);
	pdata = rbind(pdata,t(plist))
 	pData(params) = pdata 
 	
  return(cytofCore.updateFlowFrameKeywords(flowFrame(exprs=cbind(exprs(flowFrame), `colnames<-`(cbind(combinedValues), newName)),parameters=params,description=description(flowFrame))))
}

cytofCore.subtract = function(flowFrame,value=100,exclude=c("Time","Cell Length")){
	
	if (!all(exclude %in% colnames(flowFrame))){
		offending = which(!(exclude %in% colnames(flowFrame)))
		stop( paste("Channel(s)",paste(exclude[offending],collapse=","),"not found in frame.") )
	}
	subtractCols = setdiff(colnames(flowFrame),exclude)
	exprs(flowFrame)[,subtractCols] = exprs(flowFrame)[,subtractCols]-value
 
  return(cytofCore.updateFlowFrameKeywords(flowFrame))
}

cytofCore.concatenateFiles = function(inputDir,outputDir=NULL,pattern=NULL,overwrite=F,timeCol="Time"){
	currentwd = getwd();
 
	if (is.null(outputDir)){
		outputDir = inputDir;
	}
	if (!file.exists(inputDir)){
		stop(paste("Input Directory",inputDir,"not found."))
	}
 	if (!file.exists(outputDir)){
		stop(paste("Output Directory",outputDir,"not found."))
	}
 
 	setwd(inputDir);
 	if (is.null(pattern)){
  	outputFileName="combined.fcs";
  } else {
  	outputFileName = paste(paste(pattern,collapse="-"),"_combined.fcs",sep="")	
  }
  
  if (file.exists(outputFileName) && !overwrite){
  	stop(paste("File ",outputFileName," already exists in dir ",outputDir,". Set overwrite argument to TRUE or remove file first.",sep=""))
  }
 
	fcsFiles = list.files(inputDir,pattern=".fcs",ignore.case=F)
  
  matched = c(1:length(fcsFiles));
  for (match in pattern){
  	matched=intersect(matched,grep(match,fcsFiles))
  }
	combineFiles = fcsFiles[matched]
 
  begin=T
 	for (file in combineFiles){
 	  cat(paste("Adding",file,"\n"));
 		if (begin){
 			combinedData = exprs(read.FCS(file));
			if (!(timeCol %in% colnames(combinedData))){
				stop(paste("Time Column",timeCol,"not found for file",file));		
			}
		
			begin=F
		
 		} else {
			tmpData = exprs(read.FCS(file));
 			if (!(timeCol %in% colnames(tmpData))){
				stop(paste("Time Column",timeCol,"not found for file",file));		
			}
	 		tmpData[,timeCol] = tmpData[,timeCol]+(max(combinedData[,timeCol])+round((max(combinedData[,timeCol])-min(combinedData[,timeCol]))/nrow(combinedData)))
 			combinedData = rbind(combinedData,tmpData)
 		}
 	}
  
 	write.FCS(cytofCore.updateFlowFrameKeywords(flowFrame(exprs=combinedData)),outputFileName);
  rm(combinedData,tmpData);
  gc()
  setwd(currentwd);
  return(paste(outputDir,"->",outputFileName))
}