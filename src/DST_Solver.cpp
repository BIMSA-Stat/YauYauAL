#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
// [[Rcpp::depends(RcppEigen)]]

Eigen::VectorXd mydstn(const Eigen::VectorXd &input, int Ns) {
 Eigen::VectorXd output = Eigen::VectorXd::Zero(Ns);
 for (int k = 0; k < Ns; ++k) {
   double sum = 0.0;
   for (int n = 0; n < Ns; ++n) {
     sum += input[n] * std::sin(M_PI * (n + 1) * (k + 1) / (Ns + 1));
   }
   output[k] = sum * std::sqrt(2.0 / (Ns + 1));
 }
 return output;
}

//' @name DST_Solver
//' @title Solve a System Using the Discrete Sine Transform (DST-II)
//' @description
//' This function solves a system using the Discrete Sine Transform (DST-II).
//' It computes transformations, applies scaling with eigenvalues, and normalizes the result.
//'
//' @usage
//' DST_Solver(Dim, Lambda, B, U, Ns)
//'
//' @param Dim Integer. The dimension of the problem. For example, \code{Dim = 1} solves a one-dimensional problem.
//' @param Lambda A numeric vector of type \code{Eigen::VectorXd} representing scaling factors (eigenvalues).
//' @param B A sparse matrix of type \code{Eigen::SparseMatrix<double>} representing a linear operator in the system.
//' @param U A numeric vector of type \code{Eigen::VectorXd} representing the initial state vector.
//' @param Ns Integer. The size of each segment to be processed by DST-II.
//'
//' @details
//' The function performs the following steps:
//' \enumerate{
//'   \item Computes the product \eqn{B \cdot U}, where \eqn{B} is a sparse matrix and \eqn{U} is the input vector.
//'   \item Applies the Discrete Sine Transform (DST-II) to the product.
//'   \item Scales the transformed values by dividing by the vector \eqn{\Lambda}.
//'   \item Applies DST-II again to the scaled values.
//'   \item Ensures all values in the result are non-negative and normalizes the vector so that its sum is 1.
//' }
//' It handles one-dimensional problems as a special case and generalizes to higher dimensions by processing segments independently.
//'
//' @return
//' A numeric vector of type \code{Eigen::VectorXd} representing the final normalized result after applying DST-II twice and scaling.
//'
//' @examples
//' library(Matrix)
//' Dim <- 1
//' Lambda <- c(1, 2, 3, 4)
//' B <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(0.5, 0.75), dims = c(4, 4))
//' U <- c(1, 2, 3, 4)
//' Ns <- 4
//' result <- DST_Solver(Dim, Lambda, B, U, Ns)
//' print(result)
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd DST_Solver(int Dim, const Eigen::VectorXd &Lambda, const Eigen::SparseMatrix<double> &B, const Eigen::VectorXd &U, int Ns) {
 Eigen::VectorXd result;

 switch (Dim) {
 case 1: {
   // One-dimensional case
   Eigen::VectorXd BU = B * U;
   result = Eigen::VectorXd::Zero(Ns);

   // Apply DST-II to B * U
   result = mydstn(BU, Ns);
   result = result.array() / Lambda.array();

   // Apply DST-II again
   Eigen::VectorXd temp = result;
   result = mydstn(temp, Ns);
   break;
 }
 default: {
   // Higher dimensions
   int total_size = std::pow(Ns, Dim);
   Eigen::VectorXd BU = B * U;

   result = Eigen::VectorXd::Zero(total_size);
   for (int i = 0; i < total_size / Ns; ++i) {
     Eigen::VectorXd temp = BU.segment(i * Ns, Ns);
     Eigen::VectorXd dst_result(Ns);

     dst_result = mydstn(temp, Ns);
     result.segment(i * Ns, Ns) = dst_result;
   }
   result = result.array() / Lambda.array();

   // Apply DST-II again
   for (int i = 0; i < total_size / Ns; ++i) {
     Eigen::VectorXd temp = result.segment(i * Ns, Ns);
     Eigen::VectorXd dst_result(Ns);

     dst_result = mydstn(temp, Ns);
     result.segment(i * Ns, Ns) = dst_result;
   }
   break;
 }
 }

 // Ensure non-negative values and normalize
 result = result.cwiseMax(0.0);
 result = result / result.sum();

 return result;
}
