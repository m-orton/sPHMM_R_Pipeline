#include <Rcpp.h>
#include <unordered_set>
#include <unordered_map>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame strobemer_filterUniqueSp(DataFrame strobeDT) {
  CharacterVector spID = strobeDT["spID"];
  CharacterVector Strobemer = strobeDT["Strobemer"];
  CharacterVector recordID = strobeDT["recordID"];
  int n = strobeDT.nrows();
  
  // Prepare to store results
  std::unordered_set<std::string> seenSequences;
  std::vector<std::string> uniqueSeqs;
  std::vector<std::string> uniqueSpecies;
  std::vector<std::string> uniqueRecordIDs;
  
  // Loop through the sequences and species
  for (int i = 0; i < n; ++i) {
    std::string seq = Rcpp::as<std::string>(Strobemer[i]);
    std::string sp = Rcpp::as<std::string>(spID[i]);
    std::string rec = Rcpp::as<std::string>(recordID[i]);
    
    // Check if the sequence has been seen for the current species
    if (seenSequences.find(seq) == seenSequences.end()) {
      seenSequences.insert(seq);
      uniqueSeqs.push_back(seq);
      uniqueSpecies.push_back(sp);
      uniqueRecordIDs.push_back(rec);
    }
  }
  
  return DataFrame::create(Named("spID") = uniqueSpecies,
                           Named("Strobemer") = uniqueSeqs,
                           Named("recordID") = uniqueRecordIDs);
}