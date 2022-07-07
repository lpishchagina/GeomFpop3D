#ifndef REC_EMPTY_EMPTY_H
#define REC_EMPTY_EMPTY_H

#include "d3rectangle.h"
#include "d3cost.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

class rec_empty_empty{
private:
  unsigned int tau;
  d3rectangle rectangle;// empty 
  std::array<double,3>* csY;
  std::array<double,3>* csY2;
  double* locCosts;
  bool fl_empty;
public:
  rec_empty_empty(): tau(0), csY(NULL), csY2(NULL), locCosts(NULL),    fl_empty(false) { }
  rec_empty_empty(unsigned int t):
  tau(t), csY(NULL), csY2(NULL), locCosts(NULL),  fl_empty(false) { }
  rec_empty_empty(const rec_empty_empty & candidate);
  ~rec_empty_empty();
  
  unsigned int get_tau()const;
  d3rectangle getRectangle() const;
  
  bool EmptyOfCandidate();
  void idCandidate(unsigned int t, std::array<double,3>* &csy, std::array<double,3>* &csy2, double* &loccosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<rec_empty_empty>::iterator> &vectlinktocands);
};
#endif //REC_EMPTY_EMPTY_H