#include <Rcpp.h>
using namespace Rcpp;

//' Finds the median
//'
//' Returns the median of a list of values with corresponding frequencies. This
//'  is not intended to be called directly by the user.
//'
//' @param values A vector of the unique values that occur
//' @param freqs A vector of the number of occurrences of each value
//' @return Returns the median value of the data comprising each entry in values
//'  repeated the corresponding entry in freqs number of times, as a numeric.
//' @author Tom Mayo \email{t.mayo@@ed.ac.uk}
// [[Rcpp::export]]
double median_freq(NumericVector values, IntegerVector freqs) {
    const int len = freqs.size();
    std::vector<std::pair<double, int> > allDat;
    int freqSum = 0;
    for (int i=0; i < len; ++i) {
        allDat.push_back(std::pair<double, int>(values[i], freqs[i]));
        freqSum += freqs[i];
    }
    std::sort(allDat.begin(), allDat.end());
    int accum = 0;
    for (int i=0; i < len; ++i) {
        accum += allDat[i].second;
        if (freqSum % 2 == 0) {
            if (accum > freqSum / 2) {
                return allDat[i].first;
            } else if (accum == freqSum / 2) {
                return (allDat[i].first + allDat[i+1].first) / 2;
            }
        } else {
            if (accum >= (freqSum+1)/2) {
                return allDat[i].first;
            }
        }
    }
    return NA_REAL;  // Should not be reached
}
