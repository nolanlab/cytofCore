library("cytofCore")
args = commandArgs();
for (filename in args){
	if (sum(grep(".fcs",filename,ignore.case=T))==0){
		cat(paste("Ignoring argument",filename,"\n"))
	} else {
		cat(paste("Subtracting from ",filename,".\n",sep=""))
		outputFile = write.FCS(cytofCore.subtract(flowFrame=read.FCS(filename),value=100),filename=paste(filename,".subtracted.fcs",sep=""))
		cat(paste("Output File: ",outputFile,".\n\n",sep=""))
	}
}