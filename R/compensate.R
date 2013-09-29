applyCompensation = function(flowFrame,compMatrix){
  return(compensate(flowFrame,compMatrix))
}

cytofCore.applyFileCompensation = function(inputFCSFile,inputCompFile,outputFile,row.names=1,skip=0){
  flowFrame = read.FCS(inputFCSFile)
  compMatrix = read.csv(inputCompFile,check.names=F,header=T,row.names=row.names,skip=skip)
  compedFrame = applyCompensation(flowFrame,compMatrix)
  write.FCS(compedFrame,filename=outputFile)
}

cytofCore.batchFileCompensate = function(inputFCSDirectory,inputCompDirectory,outputDirectory,fcsFiles,compFiles){
  if ((length(compFiles)=!1)||(length(compFiles)!=length(fcsFiles))){
    stop("Must provide a single comp file or a comp file for each fcs file.")
  }
  if (length(compFiles)==1){
    compFiles = rep(compFile,length(fcsFiles))
  }
  for (i in 1:length(fcsFiles)){
    fcsFile = fcsFiles[i]
    compFile = compFiles[i]
    print("Comping",fcsFile,"with",compFile)
    cytofCore.applyFileCompensation(inputFCSFile=file.path(inputFCSDirectory,fcsFile),
                                    inputCompFile=file.path(inputCompDirectory,compFile),
                                    outputFile=file.path(outputDirectory,fcsFile))
  }
}

