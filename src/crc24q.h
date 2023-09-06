#ifndef __CRC24Q_H__
#define __CRC24Q_H__

#include <Rcpp.h>
using namespace Rcpp;

unsigned int calculate_crc24q(IntegerVector data, int len);
bool checkMessage(IntegerVector msg);
#endif //__CRC24Q_H__
