#include "d3cost.h"
#include <iostream>
#include "math.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
//+
d3cost::d3cost(unsigned int i, unsigned int t, std::array<double, 3>* &csdY, std::array<double, 3>* &csdY2, double* &lCosts){
  k = t - i + 1;
  mi1beta = lCosts[i];
  double sumEYit2 = 0;
  double sumYit2 = 0;
  for (unsigned int p = 0; p < 3; p++) {
    EYit[p] = (csdY[t+1][p] - csdY[i][p])/k;
    
    sumEYit2 = sumEYit2 + EYit[p] * EYit[p];
    sumYit2 = sumYit2 + (csdY2[t+1][p] - csdY2[i][p]);
  }
  kVYit = sumYit2 - k * sumEYit2;
}

d3cost::d3cost(const d3cost &cost){
  k = cost.k;
  mi1beta = cost.mi1beta;
  kVYit =cost.kVYit;
  EYit = cost.EYit;
}

void d3cost::idCost(unsigned int i, unsigned int t, std::array<double, 3>* &csdY, std::array<double, 3>* &csdY2, double* &lCosts){
  k = t - i + 1;
  mi1beta = lCosts[i];
  double sumEYit2 = 0;
  double sumYit2 = 0;
  for (unsigned int p = 0; p < 3; p++){
    EYit[p] = (csdY[t+1][p] - csdY[i][p])/k;
    
    sumEYit2 = sumEYit2 + EYit[p] * EYit[p];
    sumYit2 = sumYit2 + (csdY2[t+1][p] - csdY2[i][p]);
  }
  kVYit = sumYit2 - k * sumEYit2;
}

unsigned int d3cost::get_k() const {return k;}
double d3cost::get_kVYit() const {return kVYit;}
double d3cost::get_mi1beta() const {return mi1beta;}
std::array<double,3> d3cost::get_EYit() {return EYit;}
double d3cost::get_min() { return (kVYit + mi1beta);}

double d3cost::CostOfMu(std::array<double, 3> muFix, unsigned int i, unsigned int t, std::array<double, 3>* &csdY, std::array<double, 3>* &csdY2, double* &lCosts) {
  double mi_1_p = lCosts[i];
  double sum_mu2 = 0;
  double sum_x2 = 0;
  double sum_scal = 0;
  for (unsigned int k = 0; k < 3; k++) {
    sum_mu2 = sum_mu2 + muFix[k] * muFix[k];
    sum_x2 = sum_x2 + (csdY2[t+1][k] - csdY2[i][k]);
    sum_scal = sum_scal - 2 * (csdY[t+1][k] - csdY[i][k]) * muFix[k];
  }
  return(sum_x2 + sum_scal+(t-i+1) * sum_mu2 + mi_1_p);
      
}


