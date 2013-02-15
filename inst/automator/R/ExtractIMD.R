library("cytofCore")
args = commandArgs()
print(args)
if (sum(grep("imd$",args))<1 || sum(grep("conf$",args))<1){
	stop("Can't find IMD or Conf missing")
}
imdFile = args[grep("imd$",args)]
confFile = args[grep("conf$",args)]
startTime = as.numeric(args[5])
endTime=as.numeric(args[6])
if (endTime=="EndOfFile"){
	endTime=NULL
} else {
	endTime=as.numeric(endTime)
}
print(paste("Extracting from",startTime,"to",endTime));
cytofCore.write.IMD.stream(imdFile=imdFile,confFile=confFile,startTime=startTime,endTime=endTime)

