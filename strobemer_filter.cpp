#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame strobemer_FilterSp(DataFrame strobeDT, String foldName, String spIDName, int maxS) {
  // Extract columns from the data frame
  CharacterVector recordID = strobeDT["recordID"];
  CharacterVector strobemer = strobeDT["Strobemer"];
  
  // Calculate the number of unique record IDs
  int numRecs = unique(recordID).size();
  
  // Map to store unique recordID counts per Strobemer
  std::map<std::string, int> strobemerMap;
  for (int i = 0; i < recordID.size(); ++i) {
    strobemerMap[as<std::string>(strobemer[i])]++; 
  }
  
  // Prepare unique Strobemers and sCoverage_Group
  std::vector<std::string> uniqueStrobemer;
  std::vector<double> sCoverage_Group;
  for (auto& entry : strobemerMap) {
    uniqueStrobemer.push_back(entry.first);
    sCoverage_Group.push_back(static_cast<double>(entry.second) / numRecs);
  }
  
  // Sort the data by sCoverage_Group in descending order
  std::vector<size_t> indices(uniqueStrobemer.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&sCoverage_Group](size_t i, size_t j) {
    return sCoverage_Group[i] > sCoverage_Group[j];
  });
  
  // Reorder the vectors based on sorted indices
  std::vector<std::string> sortedStrobemer;
  std::vector<double> sortedSCoverageGroup;
  for (size_t i : indices) {
    sortedStrobemer.push_back(uniqueStrobemer[i]);
    sortedSCoverageGroup.push_back(sCoverage_Group[i]);
  }
  
  // Limit to maxS rows
  size_t limit = std::min(static_cast<size_t>(maxS), sortedStrobemer.size());
  sortedStrobemer.resize(limit);
  sortedSCoverageGroup.resize(limit);
  
  // Add spID and folds columns
  std::vector<std::string> spID(limit, spIDName);
  std::vector<std::string> folds(limit, foldName);
  
  // Create the resulting data frame
  return DataFrame::create(
    Named("Strobemer") = sortedStrobemer,
    Named("sCoverage_Group") = sortedSCoverageGroup,
    Named("spID") = spID,
    Named("folds") = folds
  );
}
