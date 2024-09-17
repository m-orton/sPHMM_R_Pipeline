#include <Rcpp.h>
#include <vector>
#include <string>

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame kmer_posFindCpp(List sequences, CharacterVector seq_ids, int klen) {
  int num_sequences = sequences.size();
  
  // Determine the maximum number of k-mers any sequence can have
  int max_kmers = 0;
  for (int i = 0; i < num_sequences; ++i) {
    std::string seq = as<std::string>(sequences[i]);
    int seq_length = seq.size();
    int num_kmers = seq_length - klen + 1;
    if (num_kmers > max_kmers) {
      max_kmers = num_kmers;
    }
  }
  
  // Create a matrix to store the k-mers
  CharacterMatrix kmers_matrix(num_sequences, max_kmers);
  
  for (int i = 0; i < num_sequences; ++i) {
    std::string seq = as<std::string>(sequences[i]);
    int seq_length = seq.size();
    int num_kmers = seq_length - klen + 1;
    
    for (int j = 0; j < max_kmers; ++j) {
      if (j < num_kmers) {
        kmers_matrix(i, j) = seq.substr(j, klen);
      } else {
        kmers_matrix(i, j) = NA_STRING; // Fill with NA if there are not enough k-mers
      }
    }
  }
  
  // Create column names for the kmers
  CharacterVector col_names(max_kmers + 1);
  col_names[0] = "SeqID";
  for (int i = 1; i <= max_kmers; ++i) {
    col_names[i] = std::to_string(i);
  }
  
  // Convert the matrix to a DataFrame with sequence IDs and set column names
  DataFrame df = DataFrame::create(_["SeqID"] = seq_ids, _["Kmers"] = kmers_matrix);
  df.attr("names") = col_names;
  return df;
}
