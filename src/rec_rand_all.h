#ifndef REC_RAND_ALL_H
#define REC_RAND_ALL_H

#include "d3rectangle.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

/*+++
 Class rec_rand_all
 -------------------------------------------------------------------------------
 Description:
 Rectangular approximation (Zone = AI(ramdom sphere from intersection set) + AE(all))
 Parameters:
 "tau" - change-point candidate;
 "csY" - cumsum of data;
 "csY2" - cumsum data^2;
 "locCosts" - min value of costs
 "rectangle" - approximation zone
 "spheresBefore" - exclusion set;
 "flCreate" - indicator, if "true" => new change-point candidate => create labels of the elements from exclusion set
 -------------------------------------------------------------------------------
 */
class rec_rand_all {
private:
  unsigned int tau;
  d3rectangle rectangle;
  std::array<double,3>* csY;
  std::array<double,3>* csY2;
  double* locCosts;
  std::list<d3sphere> spheresBefore;
  bool flCreate;
  
public:
  rec_rand_all(): tau(0), rectangle(d3rectangle()),  csY(NULL),csY2(NULL), locCosts(NULL), flCreate(true) {}
  rec_rand_all(unsigned int t):tau(t), rectangle(d3rectangle()), csY(NULL), csY2(NULL), locCosts(NULL), flCreate(true) {}
  rec_rand_all(const rec_rand_all & candidate);
  ~rec_rand_all();
  
  std::list<d3sphere> get_spheresBefore() const;
  unsigned int get_tau()const;
  
  int get_Number(int N);
  double get_dist(std::array<double,3> pnt1, std::array<double,3> pnt2);
  bool EmptyOfCandidate();
  void idCandidate(unsigned int t, std::array<double,3>* &csy, std::array<double,3>* &csy2, double* &loccosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<rec_rand_all>::iterator> &vectlinktocands);
};
#endif // REC_RAND_ALL_H