//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <assert.h>     /* assert */
#include <vector>
#include <algorithm>
#include <string>


using namespace std;
using namespace Rcpp;


RcppExport SEXP build_transition_matrix(SEXP lambda_, SEXP G_) {
  vector<double> lambda = as< vector<double> > (lambda_); 
	NumericMatrix G = as<NumericMatrix>(G_);
  
	int N = G.nrow();
  int p = lambda.size();
  
  vector<long int> indexes_in_lattice;
  for (long int i = 0; i < N; i++) {
    long int row_sum = 0;
    
    for(long int j = 0; j < p; j++) {
      if(G(i, j) == 1) {
        row_sum = row_sum + pow(2, p-j-1);
      }
    }
    indexes_in_lattice.push_back(row_sum);
  }
  
		// very large number of genotypes
	if( N > 30000)  {
		cout<<"Error, the number of genotypes, "<< N<<", is more than 30000"<< endl;
			// return a matrix
		return(NumericMatrix(1, 1));
	};
	 
	NumericMatrix tr_matrix(N, N);
	
	 
	for (long int i = 0; i < N; i++) {
	//	long int geno_idx = possible_geno_indexes[i];
		//long int tr_geno_idx = tr_indexes[geno_idx];
		
		double row_sum = 0;
		for(long int j = 0; j < p; j++) {
			if(G(i, j) == 0) {
			
				long int next_geno_idx = indexes_in_lattice[i] + pow(2, p-j-1);
				
         std::vector<long int>::iterator next_geno_idx_itr = find(indexes_in_lattice.begin(), indexes_in_lattice.end(), next_geno_idx);
         
					// if the next genotype is in the lattice
				if(next_geno_idx_itr != indexes_in_lattice.end()) {
					tr_matrix(i, next_geno_idx_itr-indexes_in_lattice.begin()) = lambda[j];
					row_sum += lambda[j];
				}
			}
		}
		tr_matrix(i, i) = -row_sum;
	} 

	return wrap(tr_matrix);
}
