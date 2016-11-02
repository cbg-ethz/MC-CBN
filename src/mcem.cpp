//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <assert.h>     /* assert */
#include <vector>
#include <algorithm>
#include <string>



using namespace std;
using namespace Rcpp;

void print_vec(const IntegerVector &vec) {
  for(int i = 0; i <vec.length(); i++) {
    cout<<vec[i]<< " ";
  }
  cout<<endl;
}

void print_matrix(const NumericMatrix &mat) {
  for(int i = 0; i <mat.nrow(); i++) {
    for(int j = 0; j <mat.ncol(); j++) {
      cout<<mat(i, j) << " ";
    }
    cout<<endl;
  }
}
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


//////////////////////////////////// Internal functions   ///////////////////////////////////

double rtexp( double m, double t) {
   double u = as<double>(runif(1) );
   return(-log(1-u*(1-exp(-t*m)))/m );
} 

double my_dexp(double x, double rate, bool log_scale=false) {
  
  
  NumericVector ret = dexp( NumericVector::create(x), rate, log_scale);
  
  return(as<double>(ret));
}


double my_pexp(double x, double  rate, bool  log_scale=false) {
  NumericVector ret = pexp( NumericVector::create(x), rate, true, log_scale);
  
  return(as<double>(ret));
}


double my_rexp(double n, double rate ) {
  NumericVector ret = rexp(n, rate);
  return(as<double>(ret));
}

NumericVector drawHiddenVarsSample(const NumericVector &O_vec,  double time, const NumericVector &lambda, 
    const NumericVector &topo_path, const vector< vector<int> > &parents, double &dens) {
  int p = lambda.length();
  NumericVector T(p), T_sum(p);
  dens = 0.0;
  for(int k = 0; k < p; k++) {
      int e = topo_path[k];  
      
      double max_parents_t = 0.0;
      
      for(unsigned int j = 0; j < parents[e].size(); j++) {
        if(T_sum[ parents[e][j] ] > max_parents_t) {
          max_parents_t = T_sum[parents[e][j]];
        }
      }
  
      if( O_vec[e] == 1 ) {
        double tmp = rtexp(lambda[e],  time-max_parents_t);
        
        dens = dens +   my_dexp( tmp, lambda[e], true) - my_pexp(time-max_parents_t, lambda[e], true);
        
        T_sum[e] =  max_parents_t + (tmp);
      } else {
        
        double  tmp = my_rexp(1, lambda[e]) ;
        T_sum[e] =  std::max(time, max_parents_t) + tmp;    
        
        double tmp2 = my_dexp( tmp, lambda[e], true) ;
        dens +=  tmp2;
      }
      
      T[e] = T_sum[e] - max_parents_t ;      
  }
  return(T);
//  return Rcpp::List::create(Rcpp::Named("T") = T, Rcpp::Named("density")=dens);
}

double log_cbn_density_( const NumericVector &t,  const NumericVector &rates)
{
  double ret = 0.0;
  
  for(int i = 0; i < rates.length(); i++) {
    ret += my_dexp(t[i], rates[i], true);
  }
  
  return(ret);
}


NumericVector _expectedTimeDifference(const NumericVector &O_vec,  double time, const NumericVector &lambda, 
    const NumericVector &topo_path, const vector< vector<int> > &parents, int nrOfSamples) {
  int p = lambda.length();
  
  double proposal_density, org_density;
  NumericMatrix timeSamples(nrOfSamples, p); 
  NumericVector importanceWeights(nrOfSamples);
  
  for(int i = 0; i < nrOfSamples; i++) {
    timeSamples(i, _) = drawHiddenVarsSample(O_vec , time, lambda, topo_path, parents, proposal_density);  
    org_density = log_cbn_density_(timeSamples(i, _), lambda);
    importanceWeights[i] = exp(org_density - proposal_density);
  }

  if(sum(importanceWeights) == 0) { // just for rare cases in the first iterations of the MC-EM
    importanceWeights = rep(1.0/nrOfSamples, nrOfSamples);
  } else {
    importanceWeights = importanceWeights / sum(importanceWeights);  
  }
  
 /*  cout<<"weights:";
  for(int i = 0 ; i < weights.length(); i++) {
    cout<< weights[i] << "   ";
  }
  cout<<endl;
    for(int i = 0 ; i < O_vec.length(); i++) {
    cout<< O_vec[i] << "   ";
  }
  cout<<endl; */


  NumericVector timeDifference(p);
  for(int i = 0; i < p; i++) {
    timeDifference[i] = sum(timeSamples(_, i) * importanceWeights);
  }
  
  return(timeDifference);
}


//////////////////////////////////// Rcpp Export functions   ///////////////////////////////////

RcppExport SEXP drawHiddenVarsSamples(SEXP N_, SEXP O_vec_,  SEXP time_, SEXP lambda_, 
     SEXP poset_, SEXP topo_path_) {
  int N = as<int>(N_);
  NumericVector O_vec = as<NumericVector>(O_vec_), lambda = as<NumericVector>(lambda_), 
      topo_path = as<NumericVector>(topo_path_);
  double time = as<double>(time_);
  NumericMatrix poset = as<NumericMatrix>(poset_);
  
  int p = lambda.length();
  NumericMatrix T(N, p), suffStat(N, p);
  NumericVector densities(N);
  
  
   vector< vector<int> > parents(p);
  
  for(int k = 0; k < p; k++) {
    for(int i = 0; i < p; i++ ) {
      if(poset(i, k) == 1) {
        parents[k].push_back(i);
      }
    }
  }
  
  for(int i = 0; i < N; i++) {
    T(i, _) = drawHiddenVarsSample(O_vec , time, lambda, topo_path, parents, densities[i]);
  }

//  return(T);
    return Rcpp::List::create(Rcpp::Named("T") = T, Rcpp::Named("densities")=densities);
}


RcppExport SEXP expectedTimeDifference(SEXP nrOfSamples_, SEXP O_vec_,  SEXP time_, SEXP lambda_, 
     SEXP poset_, SEXP topo_path_) {
  int nrOfSamples = as<int>(nrOfSamples_);
  NumericVector O_vec = as<NumericVector>(O_vec_), lambda = as<NumericVector>(lambda_), 
      topo_path = as<NumericVector>(topo_path_);
  double time = as<double>(time_);
  NumericMatrix poset = as<NumericMatrix>(poset_);
  
  int p = lambda.length();

  vector< vector<int> > parents(p);
  
  for(int k = 0; k < p; k++) {
    for(int i = 0; i < p; i++ ) {
      if(poset(i, k) == 1) {
        parents[k].push_back(i);
      }
    }
  }
  return(_expectedTimeDifference(O_vec,  time, lambda, topo_path, parents, nrOfSamples));
 
}


NumericMatrix gernerateMutationTimes(const NumericVector &lambda, const NumericMatrix &poset,  NumericMatrix &O, 
      const vector<double> &times, const NumericVector &topo_path, int nrOfSamples) {
  
  int p = poset.nrow(), n = times.size();
  
  NumericMatrix T = NumericMatrix(n, p); 

  vector< vector<int> > parents(p);
  
//   NumericVector densities(n);
  for(int k = 0; k < p; k++) {
    for(int i = 0; i < p; i++ ) {
      if(poset(i, k) == 1) {
        parents[k].push_back(i);
      }
    }
  }
 
//  for(int i = 0; i < O.nrow(); i++ ) {
//    T(i, _) = drawHiddenVarsSample(O(i, _) , times[i], lambda, topo_path, parents, densities[i]);
//  }


  for(int i = 0; i < O.nrow(); i++ ) {
    T(i, _) =  _expectedTimeDifference(O(i, _) , times[i], lambda, topo_path, parents, nrOfSamples);
  }

  return(T); 
}


RcppExport SEXP MCEM( SEXP ilambda_, SEXP poset_, SEXP O_, SEXP times_, SEXP max_iter_, SEXP alpha_, SEXP topo_path_, SEXP weights_,
  SEXP nrOfSamples_, SEXP verbose_, SEXP maxLambda_) {
  int nrOfSamples = as<int>(nrOfSamples_);
  NumericMatrix poset = as<NumericMatrix>(poset_);
  NumericMatrix O = as<NumericMatrix>(O_);
  vector<double> times = as< vector<double> >(times_);
  NumericVector ilambda = as<NumericVector >( ilambda_);
  NumericVector lambda = clone(ilambda);
  NumericVector weights = as<NumericVector >(weights_);
  int max_iter = as<int>(max_iter_),  p = poset.nrow();
  double maxLambda = as<double>(maxLambda_);
  
  NumericVector average_lambda = clone(lambda), T_colsum(p);
  NumericVector topo_path = as<NumericVector >(topo_path_);
  int verbose = as<int>(verbose_);
  double average_ll = 0, ll = 0;
  
  double N = sum(weights);
  
  double alpha = as<double>(alpha_);
  
  int iter = 0;
  
  int start_iter = max( int((1-alpha) * max_iter), 1)  ;
  
  average_lambda = rep(0, p);
  
    // tempr
  NumericMatrix newMutationTimes = O;
  while(iter < max_iter) {
    
      // E step
    newMutationTimes = gernerateMutationTimes(lambda, poset, O, times, topo_path, nrOfSamples);
    
      // M-step
      // TODO: To update the  rate parameter of a mutation, only use those genotypes that have some
      // evidence for estimation of this mutation, i.e., all of their parent mutations have happened!
    for(int i = 0; i < p; i++) {
      T_colsum[i] = sum(weights*newMutationTimes(_, i));

//      T_colsum[i] = exp( logsumexp_positive(log(weights) + log(newMutationTimes(_, i)) ) ) ;
      lambda[i] = N/T_colsum[i];
      if(lambda[i] > maxLambda) {
        lambda[i] = maxLambda;
      }
    }
    
    ll = N * sum(log(lambda)) - sum(lambda*T_colsum);
    
    
    if(iter >= start_iter) {
      average_lambda =  average_lambda  +  lambda;
      average_ll = average_ll + ll;
    }
    
    iter =  iter + 1;
  }
  average_lambda = average_lambda/ (max_iter-start_iter);
  average_ll = average_ll/ (max_iter-start_iter);

    return Rcpp::List::create(Rcpp::Named("par") = average_lambda, Rcpp::Named("ll")=average_ll);
}

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
    
//      if(nrOfSorts > pow(10, 20)) {
//        
//      }
    
    return;
  }
  
  for(int i = 0; i < poset.nrow(); i++) {
//    cout<<"********:"<< i<<endl;
    if(visited[i] == false && zeroInDegree(poset, i)) { 
//      cout<<"********: zero"<< i<<endl;
//      print_matrix(poset);

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
//      cout<<"********: zero---after"<< i<<endl;
//      print_matrix(poset);

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




///////////////// lattice size


double logsumexp(double a, double b) {
  double max_exp = a;
  
  if(b > a) {
    max_exp = b;
  }
  double sum = exp(a-max_exp) + exp(b-max_exp);
  
  return log(sum) + max_exp;
}

int get_next_node(const IntegerVector &visited, const int &p) {
  for(int i = 0; i < p; i++) {
    if(visited[i] == 0) {
      return i;
    }
  }
  cout<<"error, in get_next_node"<<endl;
  return -1;
}

double genoLatticeSizeRecursion(const NumericMatrix &poset, const IntegerVector &visited_, bool verbose) {
  IntegerVector visited = clone(visited_);
  int p = poset.ncol();

  if(sum(visited) == (p-1) ) {
    return(log(2));
  }
  
  
  int next_node = get_next_node(visited, p);

  
  visited[next_node] = 1;

  double log_size1 = genoLatticeSizeRecursion(poset, visited, verbose);

  bool one_node = true;
  for(int i = 0; i < p; i++) {
    
    if(visited[i] == 0) {
      if(poset(next_node,i) == 1 || poset(i, next_node) == 1) {
        visited[i] = 1;
        one_node = false;
      }
    }
  }
  if( one_node ) {
      // 2 * size1
    return(log(2) + log_size1); 
  } else if(sum(visited) == poset.ncol() ){
    return(logsumexp(log_size1, log(1)));
  }
 
  double log_size2 = genoLatticeSizeRecursion( poset, visited, verbose);
  
  return logsumexp(log_size1, log_size2);
}


RcppExport SEXP genoLatticeSize(SEXP poset_, bool verbose=false){
  NumericMatrix poset = clone(as<NumericMatrix>(poset_) );
  
  IntegerVector visited = rep(0, poset.ncol());
  return wrap(genoLatticeSizeRecursion(poset, visited, verbose));
}
