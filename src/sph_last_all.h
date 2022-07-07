#ifndef SPH_LAST_ALL_H
#define SPH_LAST_ALL_H
//+
#include "d3cost.h"
#include "d3rectangle.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class sph_last_all {
private:
  unsigned int tau;
  std::array<double,3>* csY;
  std::array<double,3>* csY2;
  double* locCosts;
  std::list<d3sphere> spheresBefore;
  bool flCreate;
  bool fl_empty;
  
  d3rectangle rectangle;

public:
  sph_last_all():tau(0), csY(NULL),csY2(NULL), locCosts(NULL), flCreate(true), fl_empty(false),   rectangle(d3rectangle()){ }
  sph_last_all(unsigned int t):  tau(t), csY(NULL),csY2(NULL), locCosts(NULL), flCreate(true), fl_empty(false) ,   rectangle(d3rectangle()){ }
  sph_last_all(const sph_last_all & candidate);
  ~sph_last_all();
  
  d3rectangle getRectangle() const;
  
  
  double get_dist(std::array<double,3> pnt1, std::array<double,3> pnt2);
  unsigned int get_tau()const;

  bool EmptyOfCandidate();
  void idCandidate(unsigned int t, std::array<double,3>* &csy, std::array<double,3>* &csy2, double* &loccosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<sph_last_all>::iterator> &vectlinktocands);
};
#endif //SPH_LAST_ALL_H