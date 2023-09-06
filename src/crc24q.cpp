#include <Rcpp.h>
using namespace Rcpp;
// CRC-24Q polynomial: X^24 + X^23 + X^18 + X^17 + X^14 + X^11 + X^10 + X^7 + X^6 + X^5 + X^4 + X^3 + X + 1
#define CRC_POLYNOMIAL 0x1864CFB

// Function to calculate CRC-24Q checksum
// [[Rcpp::export]]
unsigned int calculate_crc24q(IntegerVector data, int len) {
	uint32_t crc = 0x000000;

	for (size_t i = 0; i < len; i++) {
		crc ^= (data[i] << 16);

		for (int j = 0; j < 8; j++) {
			crc <<= 1;
			if (crc & 0x1000000) {
				crc ^= CRC_POLYNOMIAL;
			}
		}
	}

	return crc & 0xFFFFFF;
}
//[[Rcpp::export]]
bool checkMessage(IntegerVector msg) {
	int len = msg.length();
	return (calculate_crc24q(msg, len-3) == (msg[len-3]<<16| msg[len-2]<<8 | msg[len-1]));
}

