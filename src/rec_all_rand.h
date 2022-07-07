#ifndef REC_ALL_RAND_H
#define REC_ALL_RAND_H

#include "d3rectangle.h"
#include <vector>
#include <list>
#include <iterator>
#include <stdio.h>

/*+++
 Class rec_all_rand
 -------------------------------------------------------------------------------
 Description:
 Rectangular approximation (Zone = AI(ramdom sphere from intersection set) + AE(all))
 Parameters:
 "tau" - change-point candidate;
 "csY" - cumsum of data;
 "csY2" - cumsum data^2;
 "locCosts" - min value of costs
 "rectangle" - approximation zone
 "indexSpheresBefore" - labels of the elements from exclusion set
 "flCreate" - indicator, if "true" => new change-point candidate => create labels of the elements from exclusion set
 -------------------------------------------------------------------------------
 */
class rec_all_rand {
private:
  unsigned int tau;
  d3rectangle rectangle;
  std::array<double,3>* csY;
  std::array<double,3>* csY2;
  double* locCosts;
  std::vector<unsigned int> indexSpheresBefore;
  bool flCreate;
  
public:
  rec_all_rand(): tau(0), rectangle(d3rectangle()),  csY(NULL),csY2(NULL), locCosts(NULL), flCreate(true) {}
  rec_all_rand(unsigned int t):tau(t), rectangle(d3rectangle()), csY(NULL), csY2(NULL), locCosts(NULL), flCreate(true) {}
  rec_all_rand(const rec_all_rand & candidate);
  ~rec_all_rand();
  
  unsigned int get_tau()const;
  
  int get_Number(int N);
  double get_dist(std::array<double,3> pnt1, std::array<double,3> pnt2);
  
  bool EmptyOfCandidate();
  void idCandidate(unsigned int t, std::array<double,3>* &csy, std::array<double,3>* &csy2, double* &loccosts);
  void UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<rec_all_rand>::iterator> &vectlinktocands);
};
#endif // REC_ALL_RAND_H
