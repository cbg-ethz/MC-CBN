//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <assert.h>     /* assert */
#include <vector>
#include <algorithm>
#include <string>

using namespace std;
using namespace Rcpp;

/********************** Topological sorts - not necessary for the MCEM **********************/

bool zeroInDegree(const NumericMatrix &poset, int i) {
  for(int j = 0; j < poset.nrow(); j++) {
    if(poset(j, i) == 1) {
      return false;
    }
  }
  return true;
}

void allTopos(NumericMatrix &poset,  vector<int> &topo_path, long &nrOfSorts, vector<bool> &visited) {
  
  if(topo_path.size() == poset.nrow()) {
    nrOfSorts++;
    cout<<" New path: " << nrOfSorts<<"      :";
    
    for(unsigned int i = 0; i < topo_path.size(); i++) {
      cout<<  topo_path[i] << " ";
    }
    cout<<endl;
    
    return;
  }
  
  for(int i = 0; i < poset.nrow(); i++) {
    if(visited[i] == false && zeroInDegree(poset, i)) { 

      topo_path.push_back(i);
      visited[i] = true;
      vector<int> toNodes;
      for(int k = 0; k < poset.nrow(); k++) {
        if(poset(i, k) == 1) {
          toNodes.push_back(k);
          poset(i, k ) = 0;
        }
      }

      allTopos(poset, topo_path, nrOfSorts, visited);
      topo_path.pop_back();
      
      for(unsigned int k = 0; k < toNodes.size(); k++) {
        poset(i, toNodes[k]) = 1;  
      }
      
      visited[i] = false;
    }
  }
}

RcppExport SEXP allTopoSorts(SEXP poset_){
  NumericMatrix poset = clone(as<NumericMatrix>(poset_) );
  long nrOfSorts = 0;
  vector<int> topo_path;
  vector<bool> visited(poset.nrow(), false);
  
  allTopos(poset, topo_path,  nrOfSorts, visited);
  
  return wrap(nrOfSorts);
}


/////////////////////////////////////////
//
//double logsumexp_positive(const NumericVector &nums) {
//  double max_exp = nums[0], sum = 0.0;
//  int ct = nums.length();
//  for (int i = 1 ; i < ct ; i++)
//    if (nums[i] > max_exp)
//      max_exp = nums[i];
//
//  for (int i = 0; i < ct ; i++)
//    sum += exp(nums[i] - max_exp);
//
//  return log(sum) + max_exp;
//}
