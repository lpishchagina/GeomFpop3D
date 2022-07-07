#include "FPOP.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//converting parameters("approximation", "intersection", "exclusion") to a numeric value.
int NmbOfapproxFPOP(std::string approximation, std::string intersection,  std::string exclusion ) {
  int type_approx = 0;
  if (approximation == "rectangle") {
    if ((intersection == "empty") && (exclusion == "empty")) { type_approx = 1; }
    if (intersection == "all" ) {
      if (exclusion == "all") {  type_approx = 11; }
      if (exclusion == "empty") {  type_approx = 22; }
      if (exclusion == "random") {  type_approx = 33; }
    }
    if ((intersection == "last") && (exclusion == "all")) { type_approx = 44; }
    if ((intersection == "random") && (exclusion == "all")) { type_approx = 55; }
    //if ((intersection == "sph_all") && (exclusion == "all")) { type_approx = 66; }
  }
  if ((approximation == "sphere") && (intersection == "last") && (exclusion == "all"))  { type_approx = 77; }
  return type_approx;
}

//Comparison of two FPOP-methods
bool TestOfComparisonTwoFPOP(Rcpp::NumericMatrix data, double penalty, unsigned int type_approx2, double UnpenalizedCost1, unsigned int* LastChpt1) {
  bool res = false;
  if (type_approx2 == 1) {
    FPOP<rec_empty_empty> X = FPOP<rec_empty_empty>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false );
    res = X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 11) {
    FPOP<rec_all_all> X = FPOP<rec_all_all>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false );
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 22) {
    FPOP<rec_all_empty> X = FPOP<rec_all_empty>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false );
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 33) {
    FPOP<rec_all_rand> X = FPOP<rec_all_rand>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false );
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 44) {
    FPOP<rec_rand_all> X = FPOP<rec_rand_all>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false );
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  } else if (type_approx2 == 55) {
    FPOP<rec_rand_all> X = FPOP<rec_rand_all>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false );
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  }/* else if (type_approx2 == 66) {
    FPOP<rec_sph_all_all> X = FPOP<rec_sph_all_all>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false );
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  }*/ else if (type_approx2 == 77) {
    FPOP<sph_last_all> X = FPOP<sph_last_all>(data, penalty);
    X.algoFPOP(data, type_approx2, false, false );
    res =  X.TestFPOP(UnpenalizedCost1, LastChpt1);
  }
  
  return res;
}

//' @title approxFpop
//'
//' @description FPOP method (using the rectangle approximation of the sets) for  the multiple changepoint detection.
//' @param data is a matrix of data (p-rows x n-columns).
//' @param penalty is a value of penalty (a non-negative real number).
//' @param approximation is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle').
//' @param intersection is the type of intersection : 'empty', 'all', 'last' or 'random' (by default, 'last').
//' @param exclusion is the type of intersection : 'empty', 'all'or 'random'(by default, 'all').
//' @param NbOfCands is the logical parameter (if NbOfCands = TRUE, than the file "NbOfCands.txt" contains the number of change candidates for each iteration.
//'
//' @return a list of  elements  = (changes, means, UnpenalizedCost, NumberOfCandidates).
//'
//' \describe{
//' \item{\code{changes}}{is the changepoint vector that gives the last index of each segment for the p-variate time series.}
//' \item{\code{means}}{is the list of successive means for the p-variate time series.}
//' \item{\code{UnpenalizedCost}}{is a number equal to the global cost.}
//' \item{\code{NumberOfCandidates}}{is a number of candidates at each iteration (vector).}
//' }
//'
//' @examples
//' N <- 100
//' Chpt <-50
//' Means <-  matrix(c(0,1,1,10), nrow = 2)
//' Noise <- 1
//' Dim <- 2
//' Penalty <- 6*log(100)
//' time_series <- rnormChanges(p = Dim, n = N, changes = Chpt, means = Means, noise = Noise)
//' time_series <- rnormChanges(p=3, n = 100, changes = NULL, means = matrix(0, ncol = 1, nrow = 3), noise = 1)
//' approxFpop(data = time_series, penalty = Penalty, approximation = 'rectangle', intersection = 'all', exclusion = 'all',  NbOfCands = TRUE)
//' approxFpop(data = time_series, penalty = Penalty, approximation = 'rectangle', intersection = 'all', exclusion = 'alllines',  NbOfCands = TRUE)
//' approxFpop(data = time_series, penalty = Penalty, approximation = 'rectangle', intersection = 'all', exclusion = 'allh',  NbOfCands = TRUE)

// [[Rcpp::export]]
List approxFpop(Rcpp::NumericMatrix data, double penalty, std::string approximation = "rectangle", std::string intersection = "all",  std::string exclusion = "all", bool NbOfCands = false, bool NbOfEmpirCands = false) {
  List res;
  int type_approx = NmbOfapproxFPOP(approximation, intersection, exclusion);
  //----------stop--------------------------------------------------------------
  if (penalty < 0) {throw std::range_error("Penalty should be a non-negative number!");}
  if(type_approx == 0){throw std::range_error("This combination of parameters 'intersection' and 'exclusion' is not available. ");}
  //----------------------------------------------------------------------------
  if (type_approx == 1) {
    FPOP<rec_empty_empty> X = FPOP<rec_empty_empty>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfEmpirCands);
    res = X.ResAlgoFPOP(NbOfCands, NbOfEmpirCands);
  }  else if (type_approx == 11) {
    FPOP<rec_all_all> X = FPOP<rec_all_all>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfEmpirCands);
    res = X.ResAlgoFPOP(NbOfCands, NbOfEmpirCands);
  } else if (type_approx == 22) {
    FPOP<rec_all_empty> X = FPOP<rec_all_empty>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfEmpirCands);
    res = X.ResAlgoFPOP(NbOfCands, NbOfEmpirCands);
  } else if (type_approx == 33) {
    FPOP<rec_all_rand> X = FPOP<rec_all_rand>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfEmpirCands);
    res = X.ResAlgoFPOP(NbOfCands, NbOfEmpirCands);
  } else if (type_approx == 44) {
    FPOP<rec_rand_all> X = FPOP<rec_rand_all>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfEmpirCands);
    res = X.ResAlgoFPOP(NbOfCands, NbOfEmpirCands);
  } else if (type_approx == 55) {
    FPOP<rec_rand_all> X = FPOP<rec_rand_all>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfEmpirCands);
    res = X.ResAlgoFPOP(NbOfCands, NbOfEmpirCands);
  }/* else if (type_approx == 66) {
    FPOP<rec_sph_all_all> X = FPOP<rec_sph_all_all>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfEmpirCands);
    res = X.ResAlgoFPOP(NbOfCands, NbOfEmpirCands);
  }*/ else if (type_approx == 77) {
    FPOP<sph_last_all> X = FPOP<sph_last_all>(data, penalty);
    X.algoFPOP(data, type_approx, NbOfCands, NbOfEmpirCands);
    res = X.ResAlgoFPOP(NbOfCands, NbOfEmpirCands);
  }
  return res;
}

//'@title TestTwoApproxFpop
//'
//' @description Сomparing the parameters ("UnpenalizedCost", "LastChpt") for two different FPOP methods (using the rectangle approximation of the sets) .
//' @param data is a matrix of data (p-rows x n-columns).
//' @param penalty is a value of penalty (a non-negative real number).
//' @param approximation1 is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle').
//' @param intersection1 is the type of intersection : 'empty', 'all', 'last', 'random' or 'sphere'(by default, 'last').
//' @param exclusion1 is the type of intersection : 'empty', 'all', 'random' or 'sphere'(by default, 'all').
//' @param str1 is the structure of a repository for exclusion disks : 'l' (all) or 'v'(reduction) (by default, 'v').
//' @param approximation2 is the type of approximation : 'rectangle' or 'sphere' (by default, 'rectangle').
//' @param intersection2 is the type of intersection : 'empty', 'all', 'last', 'random' or 'sphere'(by default, 'empty').
//' @param exclusion2 is the type of intersection : 'empty', 'all', 'random' or 'sphere'(by default, 'empty').
//' @param str2 is the structure of a repository for exclusion disks : 'l' (all) or 'v'(reduction) (by default, 'v').
//'
//' @return TRUE or FALSE
//'
//' \describe{
//' \item{\code{TRUE}}{ 'TRUE' if parameters are the same.}
//' \item{\code{FALSE}}{'TRUE' if parameters are different.}
//' }
//'
//' @examples
//' time_series <- rnormChanges(p = 2, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = 2), noise = 1)
  //' TestTwoApproxFpop(data = time_series, penalty = Penality, approximation1 = 'rectangle', intersection1 = 'all', exclusion1 = 'all')

bool TestTwoApproxFpop(Rcpp::NumericMatrix data, double penalty, std::string approximation1 = "rectangle", std::string intersection1 = "all",  std::string exclusion1 = "all",  std::string approximation2 = "rectangle", std::string intersection2 = "empty",  std::string exclusion2 = "empty") {
  int type_approx1 = NmbOfapproxFPOP(approximation1, intersection1, exclusion1 );
  int type_approx2 = NmbOfapproxFPOP(approximation2, intersection2, exclusion2);
  
  //----------stop--------------------------------------------------------------
  if (penalty < 0) {
    throw std::range_error("Penalty should be a non-negative number!");
    return false;
  }
  if((type_approx1 == 0)||(type_approx2 == 0)) {
    throw std::range_error("These combinations of parameters 'intersection','exclusion' and 'str' are  not available. ");
    return false;
  }
  if(type_approx1 == type_approx2) {
    throw std::range_error("These combinations have same parameters 'intersection', 'exclusion' and 'str'. ");
    return true;
  }
  //----------------------------------------------------------------------------
  bool res = false;
  if (type_approx1 == 1) {
    FPOP<rec_empty_empty> X1 = FPOP<rec_empty_empty>(data, penalty);
    X1.algoFPOP(data, type_approx1, false ,false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 11) {
    FPOP<rec_all_all> X1 = FPOP<rec_all_all>(data, penalty);
    X1.algoFPOP(data, type_approx1, false ,false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 22) {
    FPOP<rec_all_empty> X1 = FPOP<rec_all_empty>(data, penalty);
    X1.algoFPOP(data, type_approx1, false ,false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 33) {
    FPOP<rec_all_rand> X1 = FPOP<rec_all_rand>(data, penalty);
    X1.algoFPOP(data, type_approx1, false ,false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 44) {
    FPOP<rec_rand_all> X1 = FPOP<rec_rand_all>(data, penalty);
    X1.algoFPOP(data, type_approx1, false ,false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  } else if (type_approx1 == 55) {
    FPOP<rec_rand_all> X1 = FPOP<rec_rand_all>(data, penalty);
    X1.algoFPOP(data, type_approx1, false ,false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  }/* else if (type_approx1 == 66) {
    FPOP<rec_sph_all_all> X1 = FPOP<rec_sph_all_all>(data, penalty);
    X1.algoFPOP(data, type_approx1, false ,false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  }*/ else if (type_approx1 == 77) {
    FPOP<sph_last_all> X1 = FPOP<sph_last_all>(data, penalty);
    X1.algoFPOP(data, type_approx1, false ,false);
    res = TestOfComparisonTwoFPOP(data,penalty, type_approx2, X1.GetUnpenalizedCost(), X1.GetLastChpt());
  }
  return res;
}
