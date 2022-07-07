#ifndef D3RECTANGLE_H
#define D3RECTANGLE_H

#include "math.h"
#include "d3sphere.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

class d3rectangle{
private:
  std::array<std::array<double,2>, 3> borders;//matrix(2x2) of constraints for x ,each xi =(xi0,xi1)  i = 0, p-1
  std::array<std::array<bool,2>, 3> idUpdateBorders;
  
public:
  d3rectangle();
  d3rectangle(const d3rectangle &rect);

  std::array<std::array<double,2>, 3> get_borders()const;
  std::array<std::array<bool,2>, 3> get_idUpdateBorders()const;
  
  double get_dist(std::array<double,3> pnt1, std::array<double,3> pnt2);
  double min_ab(double a, double b, bool & flminb);
  double max_ab(double a, double b, bool & flmaxb);

  bool EmptyIntersection(const d3sphere &disk);
  bool IsEmptyRect();
  void DoEmptyRect();
  void ExclusionSphere(const d3sphere &sphere);
  void IntersectionSphere(const d3sphere &sphere);
  void SphereApproximation(const d3sphere &sphere);
 };

#endif //D3RECTANGLE_H