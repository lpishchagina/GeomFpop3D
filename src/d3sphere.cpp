#include "d3sphere.h"
//+

d3sphere::d3sphere(double radius,  std::array<double,3> center) {
  r = radius;
  c = center;
}

d3sphere::d3sphere(const d3sphere &sphere){
  r = sphere.r;
  c = sphere.c;
}

void d3sphere::idSphere(double radius, std::array<double,3> center) {
  r = radius ;
  c = center;
}

double d3sphere::get_r() const{return r;}
std::array<double,3> d3sphere::get_c()const{return c;}

double d3sphere::get_dist(std::array<double,3>pnt1, std::array<double,3>pnt2) {
  double res = 0;
  for (unsigned int p = 0; p < 3; p++) {
    res = res + (pnt2[p] - pnt1[p]) * (pnt2[p] - pnt1[p]);
  }
  return sqrt(res);    
}

bool d3sphere::isIntersection(const d3sphere &sphere) {
  bool res = true; 
  double d = get_dist(c, sphere.get_c());
  if (d >= (r + sphere.get_r())) {
    res = false;
  }
  return res;
}


bool d3sphere::isnotIntersection(const d3sphere &sphere) {
  bool res = true; 
  double d = get_dist(c, sphere.get_c());
  if (d < (r + sphere.get_r())) {
    res = false;
  }
  return res;
}

bool d3sphere::isInclusion(const d3sphere &sphere) {
  bool res = false; 
  double d = get_dist(c, sphere.get_c());
  if (d <= (sphere.get_r() - r)) {
    res = true;
  }
  return res;
}

void d3sphere::createSphere(unsigned int i, unsigned int t, std::array<double,3>* &csY, std::array<double,3>* &csY2,  double* &locCosts) {
  d3cost cost = d3cost();
  cost.idCost( i, t, csY, csY2, locCosts);
  double r2 = (locCosts[t+1] - locCosts[i] - cost.get_kVYit()) / cost.get_k();
  if (r2 > 0) {
    idSphere( sqrt(r2), cost.get_EYit());
  } else {
    idSphere(0, cost.get_EYit());
  } 
}
