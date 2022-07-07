#include "rec_all_rand.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//constructor copy, destructor--------------------------------------------------
rec_all_rand::rec_all_rand(const rec_all_rand & candidate) {
  tau = candidate.tau;
  rectangle = d3rectangle();
  csY = candidate.csY;
  csY2 = candidate.csY2;
  locCosts = candidate.locCosts;
  indexSpheresBefore.clear();
  indexSpheresBefore = candidate.indexSpheresBefore;
  flCreate = candidate.flCreate;
}

rec_all_rand::~rec_all_rand() {csY = NULL; csY2 = NULL; locCosts = NULL; }

//accessory---------------------------------------------------------------------
unsigned int rec_all_rand::get_tau() const { return tau; }

//tools-------------------------------------------------------------------------
double rec_all_rand::get_dist(std::array<double,3> pnt1, std::array<double,3> pnt2) {
  double res = 0;
  for (unsigned int k = 0; k < 3; k++) {
    res = res + (pnt2[k] - pnt1[k]) * (pnt2[k] - pnt1[k]);
  }
  return sqrt(res);
}

int rec_all_rand::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

bool rec_all_rand::EmptyOfCandidate() { return rectangle.IsEmptyRect(); }

void rec_all_rand::idCandidate(unsigned int t, std::array<double,3>* &csy, std::array<double,3>* &csy2, double* &loccosts) {
  tau = t;
  csY = csy;
  csY2 = csy2;
  locCosts = loccosts;
  indexSpheresBefore.clear();
  flCreate = true;
}

void rec_all_rand::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<rec_all_rand>::iterator> &vectlinktocands) {
  d3sphere sphere = d3sphere();
  //labels of the elements from exclusion set
  if (flCreate) {//flCreate = true =>1 iteration : create labels of the elements from exclusion set
    flCreate = false;
    if (IndexToLinkOfUpdCand > 0) {
      for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
        indexSpheresBefore.push_back(vectlinktocands[i] -> get_tau());
      }
    }
  }
  //intersection approximation    //intersection set:
  for (int i = IndexToLinkOfUpdCand; i < vectlinktocands.size(); i++) {
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
  //exclusion approximation with one random sphere from exclusion set :
  if ((!rectangle.IsEmptyRect()) && (indexSpheresBefore.size() > 0)) {
    unsigned int IndexRandBeforeTau = get_Number(indexSpheresBefore.size()) - 1;
    sphere.createSphere(indexSpheresBefore[IndexRandBeforeTau],  tau-1, csY, csY2, locCosts);
    if (!rectangle.EmptyIntersection(sphere)) {
      rectangle.ExclusionSphere(sphere);
    } else {
      if (IndexRandBeforeTau < (indexSpheresBefore.size() - 1 )) {
        indexSpheresBefore[IndexRandBeforeTau] = indexSpheresBefore.back();
      }
      indexSpheresBefore.pop_back();
    }
  }
}