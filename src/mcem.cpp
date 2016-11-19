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

NumericVector drawHiddenVarsSample(const NumericVector &O_vec,  double time, const NumericVector &lambda,  const NumericVector &topo_path, 
                                   const vector< vector<int> > &parents, double &dens, double lambda_s, bool sampling_time_available) {
  int p = lambda.length();
  NumericVector T(p), T_sum(p);
  dens = 0.0;
  
  if(sampling_time_available == false) {
    time = my_rexp(1, lambda_s);
  }
  
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


NumericVector _expectedTimeDifference(const NumericVector &O_vec,  double time, const NumericVector &lambda,  const NumericVector &topo_path, 
                                  const vector< vector<int> > &parents, int nrOfSamples, double lambda_s, bool sampling_time_available) {
  int p = lambda.length();
  

  double proposal_density, org_density;
  NumericMatrix timeSamples(nrOfSamples, p); 
  NumericVector importanceWeights(nrOfSamples);
  
  for(int i = 0; i < nrOfSamples; i++) {
    timeSamples(i, _) = drawHiddenVarsSample(O_vec , time, lambda, topo_path, parents, proposal_density, lambda_s, sampling_time_available);  
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
                                      SEXP poset_, SEXP topo_path_, SEXP lambda_s_, SEXP sampling_time_available_) {
  int N = as<int>(N_);
  NumericVector O_vec = as<NumericVector>(O_vec_), lambda = as<NumericVector>(lambda_), 
    topo_path = as<NumericVector>(topo_path_);
  double time = as<double>(time_);
  bool sampling_time_available = as<bool>(sampling_time_available_);
  double lambda_s = as<double>(lambda_s_);
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
    T(i, _) = drawHiddenVarsSample(O_vec , time, lambda, topo_path, parents, densities[i], lambda_s, sampling_time_available);
  }
  
  //  return(T);
  return Rcpp::List::create(Rcpp::Named("T") = T, Rcpp::Named("densities")=densities);
}


RcppExport SEXP expectedTimeDifference(SEXP nrOfSamples_, SEXP O_vec_,  SEXP time_, SEXP lambda_, 
                                       SEXP poset_, SEXP topo_path_, SEXP lambda_s_, SEXP sampling_time_available_) {
  int nrOfSamples = as<int>(nrOfSamples_);
  NumericVector O_vec = as<NumericVector>(O_vec_), lambda = as<NumericVector>(lambda_), 
    topo_path = as<NumericVector>(topo_path_);
  double time = as<double>(time_);
  bool sampling_time_available = as<bool>(sampling_time_available_);
  double lambda_s = as<double>(lambda_s_);
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
  return(_expectedTimeDifference(O_vec,  time, lambda, topo_path, parents, nrOfSamples, lambda_s, sampling_time_available));
  
}


NumericMatrix gernerateMutationTimes(const NumericVector &lambda, const NumericMatrix &poset,  NumericMatrix &O, const vector<double> &times, 
              const NumericVector &topo_path, int nrOfSamples, double lambda_s, bool sampling_time_available) {
  int p = poset.nrow(), n = times.size();
  
  NumericMatrix T = NumericMatrix(n, p); 
  
  vector< vector<int> > parents(p);
  
  for(int k = 0; k < p; k++) {
    for(int i = 0; i < p; i++ ) {
      if(poset(i, k) == 1) {
        parents[k].push_back(i);
      }
    }
  }
  
  for(int i = 0; i < O.nrow(); i++ ) {
    T(i, _) =  _expectedTimeDifference(O(i, _) , times[i], lambda, topo_path, parents, nrOfSamples, lambda_s, sampling_time_available);
  }
  
  return(T); 
}


RcppExport SEXP MCEM( SEXP ilambda_, SEXP poset_, SEXP O_, SEXP times_, SEXP max_iter_, SEXP alpha_, SEXP topo_path_, SEXP weights_,
                      SEXP nrOfSamples_, SEXP verbose_, SEXP maxLambda_, SEXP lambda_s_, SEXP sampling_time_available_) {
  int nrOfSamples = as<int>(nrOfSamples_);
  NumericMatrix poset = as<NumericMatrix>(poset_);
  NumericMatrix O = as<NumericMatrix>(O_);
  vector<double> times = as< vector<double> >(times_);
  NumericVector ilambda = as<NumericVector >( ilambda_);
  NumericVector lambda = clone(ilambda);
  NumericVector weights = as<NumericVector >(weights_);
  int max_iter = as<int>(max_iter_),  p = poset.nrow();
  double maxLambda = as<double>(maxLambda_);
  
  bool sampling_time_available = as<bool>(sampling_time_available_);
  double lambda_s = as<double>(lambda_s_);
  
  
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
    newMutationTimes = gernerateMutationTimes(lambda, poset, O, times, topo_path, nrOfSamples, lambda_s, sampling_time_available);
    
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
