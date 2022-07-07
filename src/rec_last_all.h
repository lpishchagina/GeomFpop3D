#ifndef REC_LAST_ALL_H
#define REC_LAST_ALL_H

#include "d3rectangle.h"
#include "d3cost.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class rec_last_all {
private:
  unsigned int tau;
  d3rectangle rectangle;
  std::array<double,3>* csY;
  std::array<double,3>* csY2;
  double* locCosts;
  std::list<d3sphere> spheresBefore;
  bool flCreate;
  
public:
  rec_last_all():tau(0), rectangle(d3rectangle()),  csY(NULL),csY2(NULL), locCosts(NULL), flCreate(true) { }
  rec_last_all(unsigned int t):  tau(t), rectangle(d3rectangle()), csY(NULL),csY2(NULL), locCosts(NULL), flCreate(true) { }
  rec_last_all(const rec_last_all & candidate);
  ~rec_last_all();
  
  double get_dist(std::array<double,3> pnt1, std::array<double,3> pnt2);
  std::list<d3sphere> get_spheresBefore() const;
  unsigned int get_tau()const;
  d3rectangle getRectangle() const;
  
  bool EmptyOfCandidate();
  void idCandidate(unsigned int t, std::array<double,3>* &csy, std::array<double,3>* &csy2, double* &loccosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<rec_last_all>::iterator> &vectlinktocands);
};
#endif //REC_LAST_ALL_H