#include "rec_all_empty.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

rec_all_empty::rec_all_empty(const rec_all_empty & candidate) {
  tau = candidate.tau;
  rectangle = d3rectangle();
  csY = candidate.csY;
  csY2 = candidate.csY2;
  locCosts = candidate.locCosts;
}

rec_all_empty::~rec_all_empty() { csY = NULL; csY2 = NULL; locCosts = NULL; }

unsigned int rec_all_empty::get_tau()const { return tau; }

void rec_all_empty::CleanOfCandidate() { csY = NULL; csY2 = NULL; locCosts = NULL;}

bool rec_all_empty::EmptyOfCandidate() { return (rectangle.IsEmptyRect()); }

void rec_all_empty::idCandidate(unsigned int t, std::array<double,3>* &csy, std::array<double,3>* &csy2, double* &loccosts) {
  tau = t;
  csY = csy;
  csY2 = csy2;
  locCosts = loccosts;
}

double rec_all_empty::get_dist(std::array<double,3>pnt1, std::array<double,3>pnt2) {
  double res = 0;
  for (unsigned int k = 0; k < 3; k++) {
    res = res + (pnt2[k] - pnt1[k]) * (pnt2[k] - pnt1[k]);
  }
  return sqrt(res);    
}

d3rectangle rec_all_empty::getRectangle() const {return rectangle;}

void rec_all_empty::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<rec_all_empty>::iterator> &vectlinktocands) {
  d3sphere sphere = d3sphere();
  //intersection approximation    //intersection set:
  for (unsigned int i = IndexToLinkOfUpdCand; i < vectlinktocands.size(); i++) {
    sphere.createSphere(tau, vectlinktocands[i] -> get_tau(), csY, csY2, locCosts);
    if (sphere.get_r()  == 0) { 
      rectangle.DoEmptyRect();
      return; 
    }     //pelt
    rectangle.IntersectionSphere(sphere);
    if (rectangle.IsEmptyRect()) {
      return;
    }
  }
}