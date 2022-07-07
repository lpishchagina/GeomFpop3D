#include "d3rectangle.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
//+
d3rectangle::d3rectangle(){
  for (unsigned int i = 0; i < 3; i++) {
    borders[i] = {-INFINITY, INFINITY};
    idUpdateBorders [i]= {false, false};
  }
}

d3rectangle::d3rectangle(const d3rectangle &rect){
  borders = rect.get_borders();
}

std::array<std::array<double,2>, 3> d3rectangle::get_borders()const { return borders;}

std::array<std::array<bool,2>, 3> d3rectangle::get_idUpdateBorders()const {return idUpdateBorders;}


void d3rectangle::DoEmptyRect() {
  borders[0] = {0,0};
  borders[1] = borders[0];
  borders[2] = borders[0];
}

bool d3rectangle::IsEmptyRect() {
  for (unsigned int k = 0; k < 3; k++) {
    if (borders[k][0] >= borders[k][1]) {
      return true;
    }
  }
  return false;
}

double d3rectangle::get_dist(std::array<double,3>pnt1, std::array<double,3>pnt2) {
  double res = 0;
  for (unsigned int k = 0; k < 3; k++) {
    res = res + (pnt2[k] - pnt1[k]) * (pnt2[k] - pnt1[k]);
  }
  return sqrt(res);    
}

double d3rectangle::min_ab(double a, double b, bool & flminb) {
  if (a <= b) {
    flminb = false;
    return a;
  } else {
    flminb = true;
    return b;
  }
}

double d3rectangle::max_ab(double a, double b, bool & flmaxb) {
  if (a >= b) {
    flmaxb = false;
    return a;
  } else {
    flmaxb = true;
    return b;
  }
}

bool d3rectangle::EmptyIntersection(const d3sphere &sphere){
  std::array<double,3> c_s = sphere.get_c();
  //point_min-------------------------------------------------------------------
  std::array<double,3> pnt_min = c_s;
  for (unsigned int k = 0; k < 3; k++){
    if (pnt_min[k] <= borders[k][0]){ pnt_min[k] = borders[k][0];}
    if (pnt_min[k] >= borders[k][1]){ pnt_min[k] = borders[k][1];}
  }
  //check-----------------------------------------------------------------------
  if (get_dist(pnt_min, c_s) >= sphere.get_r()) {
    for (unsigned int i = 0; i < 3; i++) {
      idUpdateBorders [i]= {false, false};
    }
    return true;
  } else {
    return false;
  }
}

//Intersection_sphere
void d3rectangle::IntersectionSphere(const d3sphere &sphere){
  double r_s = sphere.get_r();
  std::array<double,3> c_s = sphere.get_c();
  //point_min-------------------------------------------------------------------
  std::array<double,3> pnt_min = c_s;
  for (unsigned int k = 0; k < 3; k++){
    if (pnt_min[k] <= borders[k][0]){ pnt_min[k] = borders[k][0];}
    if (pnt_min[k] >= borders[k][1]){ pnt_min[k] = borders[k][1];}
  }
  //discriminant----------------------------------------------------------------
  std::array<double,3> dx2;
  double dx2i = 1;
  unsigned int i = 0;
  while ((dx2i > 0)&&(i < 3)){
    dx2i = 0;
    for (unsigned int j = 0; j < 3; j++){if (j != i){dx2i = dx2i + (pnt_min[j] - c_s[j]) * (pnt_min[j] - c_s[j]);}}
    dx2i = r_s * r_s - dx2i;
    dx2[i] = dx2i;
    ++i;
  }
  //----------------------------------------------------------------------------
  if (i != 3) {
    borders[0][0] =  borders[0][1];
    idUpdateBorders [0]= {false, false};
    idUpdateBorders [1]= {false, false};
    idUpdateBorders [2]= {false, false};
  } else {
    for (unsigned int k = 0; k < 3; k++){
      borders[k][0] = max_ab(borders[k][0], c_s[k] - sqrt(dx2[k]), idUpdateBorders[k][0]);
      borders[k][1] = min_ab(borders[k][1], c_s[k] + sqrt(dx2[k]), idUpdateBorders[k][1]);
    }
  }
}

//ExclusionSphere
void d3rectangle::ExclusionSphere(const d3sphere &sphere) {
  double r_s = sphere.get_r();
  std::array<double,3>  c_s = sphere.get_c();
  double dx2;
  //-point_max------------------------------------------------------------------
  std::array<double,3>  pnt_max;
  for (unsigned int k = 0; k < 3; k++) {
    if (abs(c_s[k] - borders[k][1]) >= abs(c_s[k] - borders[k][0])) {
      pnt_max[k] = borders[k][1];
    }
    else{
      pnt_max[k] = borders[k][0];
    }
  }
  //discriminant----------------------------------------------------------------
  idUpdateBorders [0]= {false, false};
  idUpdateBorders [1]= {false, false};
  idUpdateBorders [2]= {false, false};
  for (unsigned int k = 0; k < 3; k++){
    dx2 = 0;
    for (unsigned int j = 0; j < 3; j++){if (j != k){dx2 = dx2 + (pnt_max[j] - c_s[j]) * (pnt_max[j] - c_s[j]);}}
    dx2 = r_s * r_s - dx2;
    if (dx2 > 0){
      if ((pnt_max[k] == borders[k][0]) && (borders[k][1] <=  c_s[k] + sqrt(dx2))) {borders[k][1] = min_ab(borders[k][1], c_s[k] - sqrt(dx2), idUpdateBorders[k][1]);}
      if ((pnt_max[k] == borders[k][1]) && (borders[k][0] >=  c_s[k] - sqrt(dx2))) {borders[k][0] = max_ab(borders[k][0], c_s[k] + sqrt(dx2), idUpdateBorders[k][0]);}
    }
  }
}

void d3rectangle::SphereApproximation(const d3sphere &sphere){
  double radius = sphere.get_r();
  std::array<double,3> c_s = sphere.get_c();
  for (unsigned int k = 0; k < 3; k++){
    borders[k][0] = c_s[k] - radius;
    borders[k][1] = c_s[k] + radius;
  }
}