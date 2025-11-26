// strobemer_chop.cpp
#include <Rcpp.h>
#include "strobemer_extract.h"  // Include the strobemer header

// [[Rcpp::export]]
std::vector<std::string> chopStrobemer(std::string seq, int n, int k, int w_min, int w_max, int strobemerType) {
  std::vector<std::string> result;
  
  // Initialize the strobemer library with parameters
  strobemer::init(n, k, w_min, w_max, static_cast<strobemer_type>(strobemerType));
  
  int number = seq.size() - strobemer::strobmer_span() + 1;
  strobemer* buff = new strobemer[number];
  
  // Chop the strobemer
  strobemer::chop_strobemer(seq.c_str(), seq.size(), buff);
  
  // Collect the valid strobemers into the result vector
  for (int i = 0; i < number; i++) {
    if (buff[i].valid) {
      result.push_back(buff[i].to_string());
    }
  }
  
  delete[] buff; // Clean up
  return result;
}
