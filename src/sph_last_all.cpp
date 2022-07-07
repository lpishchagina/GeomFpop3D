#include "sph_last_all.h"
#include <Rcpp.h>
//+
using namespace Rcpp;
using namespace std;

sph_last_all::sph_last_all(const sph_last_all & candidate) {
  tau = candidate.tau;
  csY = candidate.csY;
  csY2 = candidate.csY2;
  locCosts = candidate.locCosts;
  spheresBefore.clear();
  spheresBefore = candidate.spheresBefore;
  flCreate = candidate.flCreate;
  fl_empty = candidate.fl_empty;
  rectangle = candidate.rectangle;
}

sph_last_all::~sph_last_all() { csY = NULL; csY2 = NULL; locCosts = NULL; }

unsigned int sph_last_all::get_tau()const { return tau; }

bool sph_last_all::EmptyOfCandidate() { return (fl_empty); }

void sph_last_all::idCandidate(unsigned int t, std::array<double,3>* &csy, std::array<double,3>* &csy2, double* &loccosts) {
  tau = t;
  csY = csy;
  csY2 = csy2;
  locCosts = loccosts;
  spheresBefore.clear();
  flCreate = true;
  fl_empty = false;
  rectangle.DoEmptyRect();
}

d3rectangle sph_last_all::getRectangle() const {return rectangle;}

double sph_last_all::get_dist(std::array<double,3>pnt1, std::array<double,3>pnt2) {
  double res = 0;
  for (unsigned int k = 0; k < 3; k++) {
    res = res + (pnt2[k] - pnt1[k]) * (pnt2[k] - pnt1[k]);
  }
  return sqrt(res);    
}

void sph_last_all::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<sph_last_all>::iterator> &vectlinktocands) {
  std::list<unsigned int> spheresAfter;
  std::list<d3sphere>::iterator iter;
  typename std::list<unsigned int>::reverse_iterator riter;
  
  d3sphere sphere = d3sphere();
  d3sphere sphereTest = d3sphere();
  //flCreate = true =>1 iteration : Creation of sphereListBefore
  //exclusion spheres
  if (flCreate) {
    flCreate = false;
    if (IndexToLinkOfUpdCand > 0) {
      for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
        sphere.createSphere(vectlinktocands[i] -> get_tau(), tau-1, csY, csY2, locCosts);
        if (sphere.get_r() == 0) { fl_empty = true; return;} 
        spheresBefore.push_back(sphere);
      }
    }
  }
  //intersection spheres
  for (unsigned int i = IndexToLinkOfUpdCand; i < vectlinktocands.size(); i++) {
    spheresAfter.push_back(vectlinktocands[i] -> get_tau());
  }
  sphere.createSphere(tau, vectlinktocands[vectlinktocands.size() - 1] -> get_tau(), csY, csY2, locCosts);
  if (sphere.get_r() == 0) { fl_empty = true; return; }
  //spheres (check exclusion)
  if (spheresBefore.size() > 0) {
    iter = spheresBefore.begin();
    while (iter != spheresBefore.end()) {
      if (sphere.isInclusion(*iter)) {
        fl_empty = true;
        return;
      } else if (sphere.isnotIntersection(*iter)) {
        iter = spheresBefore.erase(iter);
      } else {
        ++iter;
      }
    }
  }
  //spheres (check intersection)
  if (spheresAfter.size() > 1) {
    riter = spheresAfter.rbegin();
    while (riter != (spheresAfter.rend())) {
      sphereTest.createSphere(tau, (*riter), csY, csY2, locCosts);
      if (sphereTest.get_r() == 0) { fl_empty = true; return; }
      if (sphere.isnotIntersection(sphereTest)) {
        fl_empty = true; return;
      } else {
        ++riter;
      }
    }
  }
}
  
