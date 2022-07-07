#include "rec_all_all.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

rec_all_all::rec_all_all(const rec_all_all & candidate) {
  tau = candidate.tau;
  rectangle = d3rectangle();
  csY = candidate.csY;
  csY2 = candidate.csY2;
  locCosts = candidate.locCosts;
  spheresBefore.clear();
  spheresBefore = candidate.spheresBefore;
  flCreate = candidate.flCreate;
}

rec_all_all::~rec_all_all() { csY = NULL; csY2 = NULL; locCosts = NULL; }

std::list<d3sphere> rec_all_all::get_spheresBefore()const { return spheresBefore; }
unsigned int rec_all_all::get_tau()const { return tau; }

void rec_all_all::CleanOfCandidate() { csY = NULL; csY2 = NULL; locCosts = NULL;  spheresBefore.clear(); }

bool rec_all_all::EmptyOfCandidate() { return (rectangle.IsEmptyRect()); }

void rec_all_all::idCandidate(unsigned int t, std::array<double,3>* &csy, std::array<double,3>* &csy2, double* &loccosts) {
  tau = t;
  csY = csy;
  csY2 = csy2;
  locCosts = loccosts;
  spheresBefore.clear();
  flCreate = true;
}

double rec_all_all::get_dist(std::array<double,3>pnt1, std::array<double,3>pnt2) {
  double res = 0;
  for (unsigned int k = 0; k < 3; k++) {
    res = res + (pnt2[k] - pnt1[k]) * (pnt2[k] - pnt1[k]);
  }
  return sqrt(res);    
}

d3rectangle rec_all_all::getRectangle() const {return rectangle;}

void rec_all_all::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<rec_all_all>::iterator> &vectlinktocands) {
    std::list<d3sphere>::iterator iter;
    d3sphere sphere = d3sphere();

  //flCreate = true =>1 iteration : Creation of sphereListBefore
  if (flCreate) {
    flCreate = false;
    if (IndexToLinkOfUpdCand > 0) {
      for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
        sphere.createSphere( vectlinktocands[i] -> get_tau(), tau-1, csY, csY2, locCosts);
        spheresBefore.push_back(sphere);
      }
    }
  }
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
  //exclusion approximation:
  if ((spheresBefore.size() > 0) && (!rectangle.IsEmptyRect())) {
    iter = spheresBefore.begin();
    while ( (iter != spheresBefore.end()) && (!rectangle.IsEmptyRect())) {
      if (rectangle.EmptyIntersection(*iter)) {
        iter = spheresBefore.erase(iter);
      } else {
        rectangle.ExclusionSphere(*iter);
        ++iter;
      }
    }
  }
}