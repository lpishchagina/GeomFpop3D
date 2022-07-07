#ifndef REC_ALL_EMPTY_H
#define REC_ALL_EMPTY_H

#include "d3rectangle.h"
#include "d3cost.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class rec_all_empty {
private:
  unsigned int tau;
  d3rectangle rectangle;
  std::array<double,3>* csY;
  std::array<double,3>* csY2;
  double* locCosts;

public:
  rec_all_empty():tau(0), rectangle(d3rectangle()),  csY(NULL),csY2(NULL), locCosts(NULL)  { }
  rec_all_empty(unsigned int t):  tau(t), rectangle(d3rectangle()), csY(NULL),csY2(NULL), locCosts(NULL){ }
  rec_all_empty(const rec_all_empty & candidate);
  ~rec_all_empty();
  
  double get_dist(std::array<double,3> pnt1, std::array<double,3> pnt2);
  unsigned int get_tau()const;
  d3rectangle getRectangle() const;
  
  void CleanOfCandidate();
  bool EmptyOfCandidate();
  void idCandidate(unsigned int t, std::array<double,3>* &csy, std::array<double,3>* &csy2, double* &loccosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<rec_all_empty>::iterator> &vectlinktocands);
};
#endif //REC_ALL_EMPTY_H