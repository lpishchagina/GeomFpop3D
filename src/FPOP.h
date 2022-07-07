#ifndef FPOP_H
#define FPOP_H
//PELT
#include "rec_empty_empty.h"//

#include "sph_last_all.h"//

#include "rec_all_all.h"//

#include "rec_last_all.h"//
#include "rec_rand_all.h"//

#include "rec_all_empty.h"//
#include "rec_all_rand.h"//

#include <Rcpp.h>
#include <array>
#include "math.h"

using namespace Rcpp;
using namespace std;

template <class CandidateOfChange>
class FPOP {
private:
  unsigned int N;
  double Penalty;
  std::array<double,3>* CumSumData;
  std::array<double,3>* CumSumData2;
  std::vector <unsigned int> Changes;
  std::vector <std::vector <double>> SegmentMeans;
  double UnpenalizedCost;
  
  double* VectOfCosts;                    //UnpenalizedCost = VectOfCosts[n] - Changes.size()*Penality
  unsigned int* LastChpt;       //vector of the best last changepoints
  std::vector <unsigned int> NbOfCandidats;
  std::vector <unsigned int> NbOfEmpirCandidats;
  std::vector <unsigned int> times;
public:
  FPOP<CandidateOfChange>() { }
  
  FPOP<CandidateOfChange> (Rcpp::NumericMatrix data, double penalty) {
    N = (unsigned int)data.ncol();
    Penalty = penalty;
    
    VectOfCosts = new double[N + 1];
    LastChpt = new unsigned int[N];
    CumSumData = new std::array<double,3> [N + 1];
    CumSumData2 = new std::array<double,3> [N + 1];
  }
  
  FPOP<CandidateOfChange> (const FPOP<CandidateOfChange> &candidate) {
    N = candidate.N;
    Penalty = candidate.Penalty;
    Changes = candidate.Changes;
    SegmentMeans = candidate.SegmentMeans;
    UnpenalizedCost = candidate.UnpenalizedCost;
    NbOfCandidats = candidate.NbOfCandidats;
    NbOfEmpirCandidats = candidate.NbOfEmpirCandidats;
    times = candidate.times;
    
    VectOfCosts = new double[N + 1];
    LastChpt = new unsigned int[N];
    CumSumData2 = new std::array<double,3> [N + 1];
    CumSumData2 = new std::array<double,3> [N + 1];
    for (unsigned int i = 0; i < N + 1; i++) {
      VectOfCosts[i] = candidate.VectOfCosts[i];
      for (unsigned int k = 0; k < 3; k++) {
        CumSumData[i][k] = candidate.CumSumData[i][k];
        CumSumData2[i][k] = candidate.CumSumData2[i][k];
      }
    }
    for (unsigned int i = 0; i < N; i++) {
      LastChpt[i] = candidate.LastChpt[i];
    }
  }
  
  ~FPOP<CandidateOfChange>() {
    delete [] CumSumData;
    delete [] CumSumData2;
    delete [] VectOfCosts;
    delete [] LastChpt;
    CumSumData = NULL;
    CumSumData2 = NULL;
    VectOfCosts = NULL;
    LastChpt = NULL;
  }
  
  std::vector <unsigned int> GetChanges() const { return Changes; }
  std::vector <std::vector <double>> GetSegmentMeans() const { return SegmentMeans; }
  double GetUnpenalizedCost() const { return UnpenalizedCost; }
  unsigned int GetN() const { return N; }
  double GetPenalty() const { return Penalty; }
  double* GetVectOfCosts() const { return VectOfCosts; }
  unsigned int* GetLastChpt()  { return LastChpt; }
  std::vector <unsigned int> GetNbOfCandidats()  { return NbOfCandidats; }
  std::vector <unsigned int> GetNbOfEmpirCandidats()  { return NbOfEmpirCandidats; }
  std::vector <unsigned int> GetTimes()  { return times; }
  
  void CalcCumSumData(Rcpp::NumericMatrix data) {
    for (unsigned int k = 0; k < 3; k++) {
      CumSumData[0][k] = 0;
    }
    for (unsigned int j = 1; j < (N + 1); j++) {
      for (unsigned int k = 0; k < 3; k++) {
        CumSumData[j][k] = CumSumData[j - 1][k] + data(k, j-1);
      }
    };
  }
  
  void CalcCumSumData2(Rcpp::NumericMatrix data) {
    for (unsigned int k = 0; k < 3; k++) {
      CumSumData2[0][k] = 0;
    }
    for (unsigned int j = 1; j < (N + 1); j++) {
      for (unsigned int k = 0; k < 3; k++) {
        CumSumData2[j][k] = CumSumData2[j - 1][k] + data(k, j-1) * data(k, j - 1);
      }
    }
  }
  
  void algoFPOP(Rcpp::NumericMatrix data, int type_approx, bool NbOfCands, bool NbOfEmpirCands) {
    if (!NbOfCands && (!NbOfEmpirCands || !NbOfEmpirCands) ) {
      NbOfCandidats.push_back(0);
     // NbOfEmpirCandidats.push_back(0);
    }
    unsigned int deltaN;
    unsigned int RealNbExclus;//del
    unsigned int nbOfObservs = 10;
    if (NbOfCands && NbOfEmpirCands) {deltaN = N/nbOfObservs ;}
    VectOfCosts[0] = 0;
    CalcCumSumData(data);
    CalcCumSumData2(data);
    CandidateOfChange candidate = CandidateOfChange();
    d3cost cost = d3cost();
    std::list<CandidateOfChange> ListOfCandidates;                    
    std::vector< typename std::list<CandidateOfChange>::iterator> VectLinkToCandidates ;
    double min_val;
    unsigned int label;
    unsigned int u;
    //Algorithm-----------------------------------------------------------------
    for (unsigned int t = 0; t < N; t++) {
      cost.idCost(t, t, CumSumData, CumSumData2, VectOfCosts);
      min_val = cost.get_min();                       //min value of cost
      label = t;                                 //best last position
      //First run: searching min
      typename std::list<CandidateOfChange>::reverse_iterator rit_candidate = ListOfCandidates.rbegin();
      while (rit_candidate != ListOfCandidates.rend()) {
        u = rit_candidate -> get_tau();
        cost.idCost(u, t, CumSumData, CumSumData2, VectOfCosts);
        if (min_val >= cost.get_min()) {
          min_val = cost.get_min();
          label = u;
        }
        ++rit_candidate;
      }
      //new min, best last changepoint and SegmentMeans--------------------------------
      VectOfCosts[t + 1] = min_val + Penalty;
      LastChpt[t] = label;
      //Candidate of Change.Initialisation.
      candidate.idCandidate(t, CumSumData, CumSumData2, VectOfCosts);
      ListOfCandidates.push_back(candidate);
      
      //Generate vector of link
      VectLinkToCandidates.clear();
      typename std::list<CandidateOfChange>::iterator VecIt = ListOfCandidates.begin();
      while (VecIt != ListOfCandidates.end()) {
        VectLinkToCandidates.push_back(VecIt);
        VecIt++;
      }
      //Second run:
      //Update ListOfCandidates
      unsigned int SizeVectLink = VectLinkToCandidates.size();
      for (unsigned int IndexOfCandVectLink = 0; IndexOfCandVectLink < SizeVectLink; IndexOfCandVectLink++) {
        VectLinkToCandidates[IndexOfCandVectLink] -> UpdateOfCandidate(IndexOfCandVectLink,VectLinkToCandidates);
      }
      //Remove empty candidates
      typename std::list<CandidateOfChange>::iterator it_candidate = ListOfCandidates.begin();
      while (it_candidate != ListOfCandidates.end()) {
        if (it_candidate -> EmptyOfCandidate()) {
          it_candidate = ListOfCandidates.erase(it_candidate);
          --it_candidate;
        }
        ++it_candidate;
      }
      if (NbOfCands) { NbOfCandidats.push_back(VectLinkToCandidates.size());}
    /*  if (NbOfCands && NbOfEmpirCands && ((t % deltaN) == (deltaN - 1) )) { 
        NbOfEmpirCandidats.push_back(calculNbOfEmpirCands(ListOfCandidates, t));
        times.push_back(t+1);
    }*/
    }
    //Result
    //vector of Changes
    unsigned int chp = N;
    while (chp > 0) {
      Changes.push_back(chp);
      chp = LastChpt[chp-1];
    }
    Changes.push_back(0);
    unsigned int j = 1;
    std::vector<double> MeanOneSegment;
    chp = N - 1;
    while (chp > 0) {
      MeanOneSegment.clear();
      for (unsigned int k = 0; k < 3; k++) {
        MeanOneSegment.push_back((CumSumData[chp + 1][k] - CumSumData[Changes[j]][k])/(chp - Changes[j] + 1));
      }
      SegmentMeans.push_back(MeanOneSegment);
      chp = Changes[j];
      j++;
    }
    reverse(SegmentMeans.begin(), SegmentMeans.end());
    Changes.pop_back();//remove 0
    reverse(Changes.begin(), Changes.end());
    Changes.pop_back();//remove N
    
    UnpenalizedCost = VectOfCosts[N] - Penalty * (Changes.size());
  }
  
  List ResAlgoFPOP(bool NbOfCands,bool NbOfEmpirCands){
    List res;
    res["changes"] = GetChanges();
    res["means"] = GetSegmentMeans();
    res["UnpenalizedCost"] = GetUnpenalizedCost();
    if (NbOfCands) {
      res["NumberOfCandidats"] = GetNbOfCandidats();
      if (NbOfEmpirCands) { 
        res["NumberOfEmpirCandidats"] = GetNbOfEmpirCandidats();  
        res["TimesForEmpirCandidats"] = GetTimes();
      }
    }
    return res;
  }
  
  //Test of  the  UnpenalizedCost  and LastChpt
  bool TestLastChpt(unsigned int* LastChptOfTestX) {
    bool res  = true;
    unsigned int  IndexOfPoint = 0;
    while (res && (IndexOfPoint < N)) {
      if (LastChpt[IndexOfPoint] != LastChptOfTestX[IndexOfPoint] ) {
        res = false;
      }
      IndexOfPoint++;
    }
    return res;
  }
  
  bool TestUnpenalizedCost(double UnpenalizedCostOfTestX) {
    bool res  = (UnpenalizedCost == UnpenalizedCostOfTestX);
    return res;
  }
  
  bool TestFPOP(double UnpenalizedCostOfTestX, unsigned int* LastChptOfTestX) {
    bool res  = (UnpenalizedCost == UnpenalizedCostOfTestX);
    if (res) {
      res = TestLastChpt(LastChptOfTestX);
    }
    return res;
  }
  /*
  ////////////////////////////////////////////////////////////////////////////////
  unsigned int calculNbOfEmpirCands(std::list<CandidateOfChange> &listOfCands, unsigned int t) {
    //Rcpp::Rcout<<"t = "<<t<<endl;
    std::vector<unsigned int> nbCandsEmpir;
    unsigned int count = 0;
    typename std::list<CandidateOfChange>::reverse_iterator rit_candidate = listOfCands.rbegin();
    typename std::list<CandidateOfChange>::iterator it_candidate;
    
    std::array<double,2> x0y0;// [x, y]
    std::array<double,2> dxdy; //[dx, dy]
    std::array<double,2> muFix;//[mu_x, mu_y]
    double delta = 512;
    //parameters for min
    d2cost cost = d2cost();
    double min_val;
    unsigned int tauOpt;
    unsigned int u;
    it_candidate = listOfCands.begin();
    while (it_candidate != listOfCands.end()) {//Begin : point [x0,y0] ,calcul deltaOfGrid
      for (unsigned int p = 0; p < 2; p++) {
        x0y0[p] = (*it_candidate).getRectangle().get_borders()[p][0];
        dxdy[p] = ((*it_candidate).getRectangle().get_borders()[p][1] - (*it_candidate).getRectangle().get_borders()[p][0])/delta;
      }
      for (unsigned int i = 0; i < delta; i++) {
        for (unsigned int j = 0; j < delta; j++) {
          muFix[0] = x0y0[0] + i * dxdy[0]; // point into the grid
          muFix[1] = x0y0[1] + j * dxdy[1];
          rit_candidate = listOfCands.rbegin();
          min_val = INFINITY; 
          double cost_value;
          while (rit_candidate != listOfCands.rend()) {
            //u = rit_candidate->get_tau();
            u = (*rit_candidate).get_tau();
            cost_value = cost.CostOfMu(muFix, u, t, CumSumData, CumSumData2, VectOfCosts);
           //!!! if (min_val >= cost_value) {
              if (min_val > cost_value) {
              min_val = cost_value;
              tauOpt = u;
            }
            rit_candidate++;
          }//new min in tauOpt
        //  if (tauOpt == it_candidate->get_tau()) { 
          if (tauOpt == (*it_candidate).get_tau()) { 
           // Rcpp::Rcout<<"muFix = ("<<muFix[0]<<","<<muFix[1]<<")"<<endl;
            count++;
            i = delta;
            j = delta;
          }  
        }
      }
      ++it_candidate;
    }
    return count;
  }
*/
};
#endif //FPOP_H

