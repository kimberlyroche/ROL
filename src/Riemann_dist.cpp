#include <Rcpp.h>
#include <RcppEigen.h>
#include <omp.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

using Eigen::MatrixXd;
using Eigen::VectorXd;

//' Calculate distances between a matrix of positive-definite covariance matrices concatenated
//' horizontally
//'
//' @param samples posterior samples matrix, arranged by host and concatenated across columns
//' @param n_hosts number of hosts
//' @param n_samples_per number of samples per host
//' @details Dimensions of the samples matrix are D rows x DN columns where D = no. of features
//' and N = no. of samples
//' @return distance
// [[Rcpp::export]]
Eigen::MatrixXd Riemann_dist_samples(Eigen::MatrixXd samples, int n_hosts, int n_samples_per) {
  // Eigen parallelization doesn't appear to fight with this
  int P = samples.rows();
  MatrixXd distances((n_hosts*n_samples_per), (n_hosts*n_samples_per));
  #pragma omp parallel shared(samples, P, distances) num_threads(omp_get_max_threads())
  {
  #pragma omp for collapse(2)
  for(int i=0; i < n_hosts*n_samples_per; i++) {
    for (int j=0; j < n_hosts*n_samples_per; j++) {
      if(j >= i) {
        //if(i % 100 == 0 && j % 100 == 0) {
        //  Rcout << "Calculating distance between samples " << i << ", " << j << std::endl;
        //}
        // get samples A and B
        int a_idx = i*P;
        MatrixXd A = samples.middleCols(a_idx, P);
        int b_idx = j*P;
        MatrixXd B = samples.middleCols(b_idx, P);
        // calculate distance
        Eigen::LLT<MatrixXd> LLT_A;
        LLT_A.compute(A);
        MatrixXd L = LLT_A.matrixL();
        MatrixXd invL = L.inverse();
        MatrixXd temp = invL*B*(invL.transpose());
        Eigen::EigenSolver<MatrixXd> es(temp);
        double d = sqrt(es.eigenvalues().array().log().abs2().sum());
        distances(i,j) = d;
      }
    }
  }
  }
  // copy the lower triangular from the upper
  // these distances are demonstrably symmetric
  for(int i=0; i < n_hosts*n_samples_per; i++) {
    for (int j=0; j < n_hosts*n_samples_per; j++) {
      if(j < i) {
        distances(i,j) = distances(j,i);
      }
    }
  }
  return(distances);
}

//' Calculate distances between all pairs of samples from two sample sets
//'
//' @param A posterior sample matrix 1, ordered by sample number, the host concatenated across columns
//' @param B posterior sample matrix 2, ordered by sample number, the host concatenated across columns
//' @param n_hosts number of hosts in set A and B (must be equal)
//' @param host_samples_A number of samples per host in set A
//' @param host_samples_B number of samples per host in set B
//' @details Dimensions of the samples matrix are D rows x DN columns where D = no. of features
//' and N = no. of samples. For two hosts "DUD" and "OMO" with 3 samples each, these are arranged
//' column-wise in the matrix as (DUD 1, OMO 1, DUD 2, OMO 2, DUD 3, OMO 3).
//' @return distance
// [[Rcpp::export]]
Eigen::MatrixXd Riemann_dist_sets(Eigen::MatrixXd A, Eigen::MatrixXd B, int n_hosts,
                                  int host_samples_A, int host_samples_B) {
  // Eigen parallelization doesn't appear to fight with this
  int P = A.rows();
  if(B.rows() != P) {
    Rcpp::stop("Number of features in sets does not match!");
  }
  MatrixXd distances((n_hosts*host_samples_A), (n_hosts*host_samples_B));
  #pragma omp parallel shared(A, B, P, distances) num_threads(omp_get_max_threads())
  {
  #pragma omp for collapse(2)
    for(int i=0; i < n_hosts*host_samples_A; i++) {
      for (int j=0; j < n_hosts*host_samples_B; j++) {
        // get individual samples a and b
        int a_idx = i*P;
        MatrixXd a = A.middleCols(a_idx, P);
        int b_idx = j*P;
        MatrixXd b = B.middleCols(b_idx, P);
        // calculate distance
        Eigen::LLT<MatrixXd> LLT_A;
        LLT_A.compute(a);
        MatrixXd L = LLT_A.matrixL();
        MatrixXd invL = L.inverse();
        MatrixXd temp = invL*b*(invL.transpose());
        Eigen::EigenSolver<MatrixXd> es(temp);
        double d = sqrt(es.eigenvalues().array().log().abs2().sum());
        distances(i,j) = d;
      }
    }
  }
  return(distances);
}

//' Calculate distances between positive-definite covariance matrices
//'
//' @param A matrix 1
//' @param B matrix 2
//' @return distance
// [[Rcpp::export]]
double Riemann_dist_pair(Eigen::MatrixXd A, Eigen::MatrixXd B) {
  Eigen::LLT<MatrixXd> LLT_A;
  LLT_A.compute(A);
  MatrixXd L = LLT_A.matrixL();
  MatrixXd invL = L.inverse();
  MatrixXd temp = invL*B*(invL.transpose());
  Eigen::EigenSolver<MatrixXd> es(temp);
  double d = sqrt(es.eigenvalues().array().log().abs2().sum());
  // Frobenius difference alternative
  // MatrixXd C = A - B;
  // double d = sqrt(C.array().abs2().sum());
  return(d);
}

