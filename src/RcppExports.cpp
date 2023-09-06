// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// calculate_crc24q
unsigned int calculate_crc24q(IntegerVector data, int len);
RcppExport SEXP _rRTCM3_calculate_crc24q(SEXP dataSEXP, SEXP lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_crc24q(data, len));
    return rcpp_result_gen;
END_RCPP
}
// checkMessage
bool checkMessage(IntegerVector msg);
RcppExport SEXP _rRTCM3_checkMessage(SEXP msgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type msg(msgSEXP);
    rcpp_result_gen = Rcpp::wrap(checkMessage(msg));
    return rcpp_result_gen;
END_RCPP
}
// rtcmCut
List rtcmCut(IntegerVector msg, int from, int to, bool display_progress);
RcppExport SEXP _rRTCM3_rtcmCut(SEXP msgSEXP, SEXP fromSEXP, SEXP toSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type msg(msgSEXP);
    Rcpp::traits::input_parameter< int >::type from(fromSEXP);
    Rcpp::traits::input_parameter< int >::type to(toSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(rtcmCut(msg, from, to, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// rtcmGetMsgNumber
int rtcmGetMsgNumber(IntegerVector msg);
RcppExport SEXP _rRTCM3_rtcmGetMsgNumber(SEXP msgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type msg(msgSEXP);
    rcpp_result_gen = Rcpp::wrap(rtcmGetMsgNumber(msg));
    return rcpp_result_gen;
END_RCPP
}
// rtcmMsmDecodeHeader
List rtcmMsmDecodeHeader(IntegerVector msg);
RcppExport SEXP _rRTCM3_rtcmMsmDecodeHeader(SEXP msgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type msg(msgSEXP);
    rcpp_result_gen = Rcpp::wrap(rtcmMsmDecodeHeader(msg));
    return rcpp_result_gen;
END_RCPP
}
// rtcmDecodeSatMsm7
List rtcmDecodeSatMsm7(IntegerVector msg, List header);
RcppExport SEXP _rRTCM3_rtcmDecodeSatMsm7(SEXP msgSEXP, SEXP headerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type msg(msgSEXP);
    Rcpp::traits::input_parameter< List >::type header(headerSEXP);
    rcpp_result_gen = Rcpp::wrap(rtcmDecodeSatMsm7(msg, header));
    return rcpp_result_gen;
END_RCPP
}
// rtcmDecodeSigMsm7
List rtcmDecodeSigMsm7(IntegerVector msg, List header);
RcppExport SEXP _rRTCM3_rtcmDecodeSigMsm7(SEXP msgSEXP, SEXP headerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type msg(msgSEXP);
    Rcpp::traits::input_parameter< List >::type header(headerSEXP);
    rcpp_result_gen = Rcpp::wrap(rtcmDecodeSigMsm7(msg, header));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rRTCM3_calculate_crc24q", (DL_FUNC) &_rRTCM3_calculate_crc24q, 2},
    {"_rRTCM3_checkMessage", (DL_FUNC) &_rRTCM3_checkMessage, 1},
    {"_rRTCM3_rtcmCut", (DL_FUNC) &_rRTCM3_rtcmCut, 4},
    {"_rRTCM3_rtcmGetMsgNumber", (DL_FUNC) &_rRTCM3_rtcmGetMsgNumber, 1},
    {"_rRTCM3_rtcmMsmDecodeHeader", (DL_FUNC) &_rRTCM3_rtcmMsmDecodeHeader, 1},
    {"_rRTCM3_rtcmDecodeSatMsm7", (DL_FUNC) &_rRTCM3_rtcmDecodeSatMsm7, 2},
    {"_rRTCM3_rtcmDecodeSigMsm7", (DL_FUNC) &_rRTCM3_rtcmDecodeSigMsm7, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_rRTCM3(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}