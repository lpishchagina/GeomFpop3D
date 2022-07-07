#include "rec_empty_empty.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

rec_empty_empty::rec_empty_empty(const rec_empty_empty & candidate) {
  tau = candidate.tau;
  csY = candidate.csY;
  csY2 = candidate.csY2;
  locCosts = candidate.locCosts;
}

rec_empty_empty::~rec_empty_empty() { csY = NULL; csY2 = NULL;  locCosts = NULL; }

unsigned int rec_empty_empty::get_tau()const { return tau; }

d3rectangle rec_empty_empty::getRectangle() const {return rectangle;}

bool rec_empty_empty::EmptyOfCandidate() { return fl_empty; }

void rec_empty_empty::idCandidate(unsigned int t, std::array<double,3>* &csy, std::array<double,3>* &csy2, double* &loccosts) {
  tau = t;
  csY = csy;
  csY2 = csy2;
  locCosts = loccosts;
  fl_empty = false;
}

void rec_empty_empty::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<rec_empty_empty>::iterator> &vectlinktocands) {
  fl_empty = false;
  d3cost cost = d3cost();
  //last t : vectlinktocands[vectlinktocands.size()-1] -> get_tau();
  cost.idCost(tau, vectlinktocands[vectlinktocands.size() - 1] -> get_tau(), csY, csY2, locCosts);
  double r2 = (locCosts[(vectlinktocands[vectlinktocands.size() - 1] -> get_tau()) + 1] - locCosts[tau] - cost.get_kVYit())/cost.get_k();
  if (r2 < 0) {
    fl_empty = true;
  }
}
  
