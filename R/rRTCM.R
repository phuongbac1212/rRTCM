require(parallel)
require(dplyr)
file = "data/data.rtcm"

source("R/RcppExports.R")
rtcmDecodeMsm7 <-function(message) {
	return(rtcmDecodeSigMsm7(message, NULL))
}

rtcmRead <- function(file) {
	data = readBin(file, "raw", n= file.info(file)$size)
	part =c(seq(1, length(data), by =length(data)/8), length(data))
	start = part[1:8]
	end = part[2:9]


	cut.list = pbmcapply::pbmcmapply(
		FUN = rtcmCut,
		from = start,
		to = end,
		MoreArgs = list(msg = data),
		mc.cores = detectCores()
	) %>% unlist() %>% unique()
	closeAllConnections()

	msg.list =DescTools::SplitAt(x = data, cut.list)
	all= unlist(lapply(msg.list, rtcmGetMsgNumber))

	x = msg.list[[426]]
	rtcmDecodeMsm7(x)
}


