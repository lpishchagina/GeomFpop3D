#ifndef D3COST_H
#define D3COST_H

#include <array>
//+
class d3cost{
private:
  unsigned int k;
  double kVYit;
  double mi1beta;
  std::array<double, 3>  EYit;
public:
  d3cost() : k(0), kVYit(0), mi1beta(0) {}
  d3cost(unsigned int i, unsigned int t, std::array<double, 3>* &csdY, std::array<double, 3>* &csdY2, double* &lCosts);
  d3cost(const d3cost &d3cost);

  void idCost(unsigned int i, unsigned int t, std::array<double, 3>* &csdY, std::array<double, 3>* &csdY2, double* &lCosts);
  unsigned int get_k() const;
  double get_kVYit() const;
  double get_mi1beta() const;
  std::array<double, 3> get_EYit();
  double get_min();
  
  double CostOfMu(std::array<double, 3> muFix, unsigned int i, unsigned int t, std::array<double, 3>* &csdY, std::array<double, 3>* &csdY2, double* &lCosts);
    
};

#endif // D3COST_H

