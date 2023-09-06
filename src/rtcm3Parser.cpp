#include <Rcpp.h>
#include <vector>
#include "crc24q.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::export]]
List rtcmCut(IntegerVector msg, int from =0, int to=NA_INTEGER, bool display_progress=true) {
	int index = from;
	if (to == NA_INTEGER)
		to = msg.length();
	List res;
	Progress p(msg.length(), display_progress);
	while (msg[index] != 0xd3) index++;
	while (index < to) {
		if (Progress::check_abort())
			return -1.0;
		res.push_back(index+1);
		int len = ((msg[index+1] & 0x03) << 8) | msg[index+2];
		index += len + 6;
		p.increment(len+6);
	}
	return(res);
}

// [[Rcpp::export]]
int rtcmGetMsgNumber(IntegerVector msg) {
	if (checkMessage(msg))
		return(((msg[3] << 8) | msg[4]) >> 4);
	else
		return NA_INTEGER;
}

unsigned long getbitu(IntegerVector buff, int pos, int len)
{
	unsigned long bits=0;
	int i;
	for (i=pos;i<pos+len;i++) bits=(bits<<1)+((buff[i/8]>>(7-i%8))&1u);
	return bits;
}
Rcpp::CharacterVector intToChar(Rcpp::IntegerVector x){
	int n = x.size();
	Rcpp::CharacterVector BASEL_SEG(n);
	for(int i = 0; i < n; i++){
		BASEL_SEG[i] = std::to_string(x[i]);
	}
	return BASEL_SEG;
}

//[[Rcpp::export]]
List rtcmMsmDecodeHeader(IntegerVector msg){
	int index = 24;
	unsigned int msgNum =getbitu(msg, index, 12); index+=12;
	unsigned int stationId = getbitu(msg, index, 12); index +=12;
	unsigned int epoch = getbitu(msg, index, 30); index+= 30;
	bool MMB =getbitu(msg, index, 1); index+= 1;
	unsigned int IODS = getbitu(msg, index, 3); index+= 3;
	index+=7;// reversed
	unsigned int SCI = getbitu(msg, index, 2); index+= 2;
	unsigned int ECI = getbitu(msg, index, 2); index+= 2;
	bool GDfSI =getbitu(msg, index, 1); index+= 1;
	unsigned int GSI = getbitu(msg, index, 3); index+= 3;
	unsigned long GSM = getbitu(msg, index, 64); index+= 64;
	unsigned long GSgM =getbitu(msg, index, 32); index+= 32;
	unsigned int satCount = 0;
	IntegerVector avaiableSat = IntegerVector::create();
	unsigned int sigCount = 0;
	IntegerVector avaiableSig = IntegerVector::create();
	for (int i=63; i>=0; i--) {
		if ((GSM >>i) & 1)
		{
			satCount += 1;
			avaiableSat.push_back(64 -i);
		}
		if ((i < 32) && ((GSgM >> i) & 1)) {
			sigCount += 1;
			avaiableSig.push_back(32-i);
		}
	}
	unsigned long cellMask =getbitu(msg, index, satCount*sigCount);
	IntegerVector cellVal = IntegerVector::create();
	for (int i=0; i<satCount*sigCount; i++) {
		cellVal.push_back(cellMask & 1);
		cellMask >>= 1;
	}
	NumericMatrix cellValues(sigCount, satCount, cellVal.begin());
	rownames(cellValues) = intToChar(avaiableSig);
	colnames(cellValues) = intToChar(avaiableSat);
	List L = List::create(
			Named("msgNum") = msgNum ,
                       _["stationId"] = stationId,
                       _["epoch"] = epoch,
                       _["MMB"] = MMB,
                       _["IODS"] = IODS,
                       _["SCI"] = SCI,
                       _["ECI"] = ECI,
                       _["GDfSI"] = GDfSI,
                       _["GSI"] = GSI,
                       _["GSM"] = GSM,
                       _["GSgM"] = GSgM,
                       _["avaiableSat"] = avaiableSat,
                       _["avaiableSig"] = avaiableSig,
                       _["cellValues"] = cellValues,
                       _["index"] = index
                       );
	return L;
}

// [[Rcpp::export]]
List rtcmDecodeSatMsm7(IntegerVector msg, List header = R_NilValue) {
	if (header.length() == 0)
		header = rtcmMsmDecodeHeader(msg);
	int index = header["index"];
	IntegerVector roughRange =IntegerVector();
	IntegerVector avaiableSat =header["avaiableSat"];
	for (int i =0; i<avaiableSat.length(); i++){
		roughRange.push_back(getbitu(msg, index, 8)); index+=8;
	}

	IntegerVector ESI =IntegerVector();
	for (int i =0; i<avaiableSat.length(); i++){
		ESI.push_back(getbitu(msg, index, 4)); index+=4;
	}

	IntegerVector roughRangeM1m =IntegerVector();
	for (int i =0; i<avaiableSat.length(); i++){
		roughRangeM1m.push_back(getbitu(msg, index, 10)); index+=10;
	}

	IntegerVector roughPhaseRangeRate =IntegerVector();
	for (int i =0; i<avaiableSat.length(); i++){
		roughPhaseRangeRate.push_back(getbitu(msg, index, 14)); index+=14;
	}
	header.push_back(roughRange,"roughRange");
	header.push_back(ESI,"ESI");
	header.push_back(roughRangeM1m,"roughRangeM1m");
	header.push_back(roughPhaseRangeRate,"roughPhaseRangeRate");
	header["index"] = index;
	return(header);
}

// [[Rcpp::export]]
List rtcmDecodeSigMsm7(IntegerVector msg, List header = R_NilValue) {
	if (header.length() == 0) {
		header = rtcmDecodeSatMsm7(msg, NULL);
	}
	int index = header["index"];
	IntegerVector avaiableSat =header["avaiableSat"];
	IntegerVector avaiableSig =header["avaiableSig"];
	int ncell = avaiableSat.length() * avaiableSig.length();

	IntegerVector Pseudoranges =IntegerVector();
	for (int i =0; i<ncell; i++){
		Pseudoranges.push_back(getbitu(msg, index, 20)); index+=20;
	}

	IntegerVector PhaseRange =IntegerVector();
	for (int i =0; i<ncell; i++){
		PhaseRange.push_back(getbitu(msg, index, 24)); index+=24;
	}

	IntegerVector PhaseRangeLti =IntegerVector();
	for (int i =0; i<ncell; i++){
		PhaseRangeLti.push_back(getbitu(msg, index, 10)); index+=10;
	}

	IntegerVector HAI =IntegerVector();
	for (int i =0; i<ncell; i++){
		HAI.push_back(getbitu(msg, index, 1)); index+=1;
	}

	IntegerVector CNR =IntegerVector();
	for (int i =0; i<ncell; i++){
		CNR.push_back(getbitu(msg, index, 10)); index+=10;
	}

	IntegerVector finePhaseRangeRates =IntegerVector();
	for (int i =0; i<ncell; i++){
		finePhaseRangeRates.push_back(getbitu(msg, index, 15)); index+=15;
	}
	Pseudoranges.attr("dim") = Dimension(avaiableSig.length(), avaiableSat.length());
	PhaseRange.attr("dim") = Dimension(avaiableSig.length(), avaiableSat.length());
	PhaseRangeLti.attr("dim") = Dimension(avaiableSig.length(), avaiableSat.length());
	HAI.attr("dim") = Dimension(avaiableSig.length(), avaiableSat.length());
	CNR.attr("dim") = Dimension(avaiableSig.length(), avaiableSat.length());

	finePhaseRangeRates.attr("dim") = Dimension(avaiableSig.length(), avaiableSat.length());
	header.push_back(Pseudoranges,"Pseudoranges");
	header.push_back(PhaseRange,"PhaseRange");
	header.push_back(PhaseRangeLti,"PhaseRangeLti");
	header.push_back(HAI,"HAI");
	header.push_back(CNR,"CNR");
	header.push_back(finePhaseRangeRates,"finePhaseRangeRates");
	header["index"] = index;
	return header;
}

