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

cytofCore.subtract = function(flowFrame,value=100,exclude=c(1,2)){
	
	if (length(exclude)>0){
		if (!is.numeric(exclude)){
			if (!all(exclude %in% colnames(flowFrame))){
				offending = which(!(exclude %in% colnames(flowFrame)))
				stop( paste("Channel(s)",paste(exclude[offending],collapse=","),"not found in frame.") )
			}
		} else {
			exclude = colnames(flowFrame)[exclude];
		}	
	}
	
	subtractCols = setdiff(colnames(flowFrame),exclude)
	exprs(flowFrame)[,subtractCols] = exprs(flowFrame)[,subtractCols]-value
 
  return(cytofCore.updateFlowFrameKeywords(flowFrame))
}

cytofCore.concatenateDirectoryFiles = function(inputDir,outputDir=NULL,pattern=NULL,overwrite=F,timeCol="Time"){
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

cytofCore.concatenateFiles = function(fcsFiles,timeCol="Time"){
	begin=T
	for (file in fcsFiles){
		if (!file.exists(file)){
			stop(paste(file,"Not Found! Stopping."))
		}
		cat(paste("Adding",file,"\n"));
		if (begin){
			combinedData = exprs(read.FCS(file));
			if (!(timeCol %in% colnames(combinedData))){
				stop(paste("Time Column",timeCol,"not found for file",file));		
			}
			begin=F
		} else {
			tmpData = exprs(read.FCS(file));
			if (!all(colnames(tmpData)==colnames(combinedData))){
				stop(paste("File names in files do not match.",colnames(tmpData),colnames(combinedData),sep="\n"))
			}
			if (!(timeCol %in% colnames(tmpData))){
				stop(paste("Time Column",timeCol,"not found for file",file));		
			}
			tmpData[,timeCol] = tmpData[,timeCol]+(max(combinedData[,timeCol])+round((max(combinedData[,timeCol])-min(combinedData[,timeCol]))/nrow(combinedData)))
			combinedData = rbind(combinedData,tmpData)
		}
	}
	outputFileName = paste(dirname(fcsFiles[1]),"/concatenated_",gsub("\\s+","_",date()),".fcs",sep="")
	write.FCS(cytofCore.updateFlowFrameKeywords(flowFrame(exprs=combinedData)),outputFileName);
	return(outputFileName)
}

cytofCore.updatePanel = function(){
  library("tcltk")
  #select template file, which has list of metals you want to include and their desired marker names
  Filters <- matrix(c("template", ".txt","template", ".TXT", "template",".fcs","template",".FCS","template",".csv","template","CSV"), 6, 2, byrow = TRUE)
  templateFile<- tk_choose.files(caption="Select template file", multi=FALSE, filters=Filters)
  templateExt=substr(templateFile,nchar(templateFile)-2,nchar(templateFile))
  
  #read in metal and marker list from template file
  if (templateExt=="csv" || templateExt=="CSV") {
    data =read.csv(templateFile, header=FALSE)
    channels=as.character(data[,1])
    markers=as.character(data[,2])
  } else if (templateExt=="txt" || templateExt=="TXT") {
    data=read.table(templateFile, header=FALSE)
    channels=as.character(data[,1])
    markers=as.character(data[,2])
  } else if (templateExt=="fcs" || templateExt=="FCS") {
    data=read.FCS(templateFile)
    channels=as.character(parameters(data)$name)
    markers=as.character(parameters(data)$desc)
  }
  
  
  #for simplicity, copy over the channel values when the markers are empty
  for (i in 1:length(markers)) {
    if ((nchar(markers[i]) == 0) || is.na(markers[i])) {
      markers[i]=channels[i]
    }
  }
  
  #select folder of FCS files to relabel
  fcsFolder<-tk_choose.dir()
  newFolder<-file.path(fcsFolder,"relabeled")
  dir.create(newFolder)
  fcsFiles<-list.files(path=fcsFolder, pattern=".fcs$")
  
  addZeroCol<-matrix(0,ncol=length(channels))  #keep track of whether to add missing column for all files
  for (file in fcsFiles) {
    oldFile=read.FCS(file.path(fcsFolder,file))
    oldChannels=as.character(parameters(oldFile)$name)
    oldData=exprs(oldFile)
    newData=matrix(data=0,nrow=nrow(oldData),ncol=length(channels))
    
    keepCols=matrix(TRUE,ncol=length(channels)) #keep track of columns to include in final matrix
    for (i in 1:length(channels)) {
      ind=grep(paste("\\Q",channels[i],"\\E",sep=""),oldChannels)
      if (length(ind)==1) { #found this template metal in this fcs file
        newData[,i]=oldData[,ind]
      } else if (addZeroCol[i]==0) { #haven't decided yet whether to create dummy column
        print(paste("Channel", channels[i], "was not found in",file,"!")) 
        yn=readline("Do you want to create a column of zeros in its place? y/n: ") 
        
        allone=readline("Do you want to apply this to all FCS files in this folder? y/n: ")
        
        if (yn=="y" && allone=="y") { #create dummy column in this and all files
          addZeroCol[i]=1
        } else if (yn=="n" && allone=="y") { #don't create dummy column in any file
          addZeroCol[i]=-1
          keepCols[i]=FALSE
        } else if (yn=="n") {#don't create dummy column but only in this file
          keepCols[i]=FALSE
        }       
      } else if (addZeroCol[i]==-1) { #previously decided not to create this dummy column
        keepCols[i]=FALSE
      }
      
    }
    
    newData=newData[,keepCols]
    newChannels=channels[keepCols]
    newMarkers=markers[keepCols]
    
    colnames(newData) = newChannels
    
    print(paste("Writing",file))
    cytofCore.write.FCS(newData,file.path(newFolder,file),channelDescriptions=newMarkers,referenceDescription=description(oldFile))
  }
  
}