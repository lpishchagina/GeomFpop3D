#ifndef D3SPHERE_H
#define D3SPHERE_H

#include <vector>
#include <array>
#include "math.h"
#include "d3cost.h"
//+
class d3sphere
{
private:
  double r; 
   std::array<double,3> c;                           
                                      
  
public:
  d3sphere(): r(0) {}
  d3sphere(double radius, std::array<double,3> center);
  d3sphere(const d3sphere & sphere);
  void idSphere(double radius, std::array<double,3> center);
  double get_r() const;
  std::array<double,3> get_c() const;
  double get_dist(std::array<double,3> pnt1, std::array<double,3> pnt2);
  
  void createSphere(unsigned int i, unsigned int t, std::array<double,3>* &csY, std::array<double,3>* &csY2,  double* &locCosts);
  
  bool isIntersection(const d3sphere & sphere);
  bool isnotIntersection(const d3sphere & sphere);
  bool isInclusion(const d3sphere & sphere);
}; 

#endif //D3SPHERE_H
//------------------------------------------------------------------------------