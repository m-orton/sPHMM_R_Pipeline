#include <Rcpp.h>
#include <fstream>
using namespace Rcpp;

inline const short int V(const char x) {
  switch (x) {
    case 'A': case 'a':
      return 0;
    break;
    case 'C': case 'c':
      return 1;
    break;
    case 'G': case 'g':
      return 2;
    break;
    default:
      return 3;
    break;
  }
}

inline unsigned int X0(const std::string A, const int k, const int n) {
  unsigned int result = 0;
  int j = k;
  for (int i = n - 1; i > n - k - 1; i--) {
    result += pow(4, k - j) * V(A[i]);
    j--;
  }
  return result;
}

// Function to write a single sequence's result to a TSV file
void writeSequenceResultToTSV(std::ofstream& outfile,
                              const int sequenceId,
                              const IntegerVector& counts,
                              const List& positions,
                              const CharacterVector& kmers,
                              const NumericVector& frequencies,
                              const CharacterVector& taxonomy) {
  int size = counts.size();
  for (int i = 0; i < size; ++i) {
    outfile << sequenceId << "\t";  // Include sequence ID
    outfile << counts[i] << "\t";
    outfile << as<IntegerVector>(positions[i])[0] << "\t";
    outfile << kmers[i] << "\t";
    outfile << frequencies[i] << "\t";
    outfile << taxonomy[i] << "\n"; // Include taxonomy
  }
}

// Modified function to include sequence IDs in the output
// [[Rcpp::export]]
inline int kmer_countToFile(const CharacterVector sequences, const CharacterVector taxonomy, const int k, const std::string outfilePath) {
  try {
    int num_sequences = sequences.size();
    std::ofstream outfile(outfilePath);

    // Write column headers
    outfile << "SequenceID\tCount\tPosition\tKmer\tFrequency\tTaxonomy\n";
    
    for (int s = 0; s < num_sequences; ++s) {
      std::string sequence = as<std::string>(sequences[s]);
      std::string taxon = as<std::string>(taxonomy[s]); 
      
      int n = sequence.size();
      
      IntegerVector counts(pow(4, k), 0);
      List positions(pow(4, k));
      CharacterVector kmers(pow(4, k));
      NumericVector frequencies(pow(4, k), 0);
      CharacterVector taxonomyVector(pow(4, k), taxon);
      
      int x = X0(sequence, k, n);
      counts[x]++;
      positions[x] = IntegerVector::create(n - k + 1);
      kmers[x] = sequence.substr(n - k, k);
      frequencies[x] = 1;
      
      const int N = pow(4, k - 1);
      for (int i = n - k - 1; i > -1; i--) {
        x = N * V(sequence[i]) + x / 4 - x % 4 / 4;
        counts[x]++;
        frequencies[x]++;
        positions[x] = (positions[x], i + 1);
        kmers[x] = (kmers[x], sequence.substr(i, k));
      }
      
      // Calculate frequencies
      for (int j = 0; j < frequencies.size(); ++j) {
        frequencies[j] /= n - k + 1;
      }
      
      // Filter out entries with count 0
      IntegerVector non_zero_counts = counts[counts > 0];
      List non_zero_positions = positions[counts > 0];
      CharacterVector non_zero_kmers = kmers[counts > 0];
      NumericVector non_zero_frequencies = frequencies[counts > 0];
      
      // Write the results to the file
      writeSequenceResultToTSV(outfile, s + 1, non_zero_counts, non_zero_positions, non_zero_kmers, non_zero_frequencies, taxonomyVector);
    }
    
    outfile.close();
    return 0;  // Return 0 for success
  } catch (std::exception& e) {
    Rcerr << "Error: " << e.what() << std::endl;
    return 1;  // Return 1 for failure
  }
}
