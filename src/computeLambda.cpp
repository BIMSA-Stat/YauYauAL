#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <vector>
#include <cmath>
#include <Eigen/Sparse>

//' @name computeLambda
//' @title Compute Lambda Vector
//' @description
//' Computes the Lambda vector using Kronecker products to construct a sparse matrix.
//' The Lambda vector represents scaled eigenvalues derived from the diagonal of the Kronecker sum
//' of a 1D Laplacian operator in multiple dimensions.
//'
//' @usage
//' computeLambda(Dim, Ns, Dt, Ds)
//'
//' @param Dim Integer. The dimension of the Kronecker product.
//' @param Ns Integer. The size of the grid in each dimension.
//' @param Dt Numeric. The time step size.
//' @param Ds Numeric. The space step size.
//'
//' @details
//' The function first computes a vector of eigenvalues \eqn{d} for a 1D Laplacian operator.
//' It then constructs a sparse matrix \eqn{D_\text{kron}} by summing the Kronecker products
//' of these eigenvalues across all dimensions. Finally, it computes the Lambda vector by scaling
//' the diagonal elements of \eqn{D_\text{kron}} using the parameters \code{Dt} and \code{Ds}.
//'
//' @return
//' A numeric vector of type \code{Eigen::VectorXd} containing the scaled eigenvalues (\code{Lambda}).
//' The length of the vector is \eqn{Ns^\text{Dim}}.
//'
//' @examples
//' Dim <- 3
//' Ns <- 4
//' Dt <- 0.01
//' Ds <- 0.5
//' Lambda <- computeLambda(Dim, Ns, Dt, Ds)
//' print(Lambda)
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd computeLambda(int Dim, int Ns, double Dt, double Ds) {
 // Step 1: Compute d vector
 Eigen::VectorXd d(Ns);
 for (int i = 0; i < Ns; ++i) {
   d[i] = -4 * std::pow(std::sin((i + 1) * M_PI / (2 * Ns + 2)), 2);
 }

 // Step 2: Initialize D_kron as a sparse matrix
 int NsDim = std::pow(Ns, Dim);
 Eigen::SparseMatrix<double> D_kron(NsDim, NsDim);
 D_kron.setZero();

 // Step 3: Compute D_kron using Kronecker products and summation
 for (int i = 1; i <= Dim; ++i) {
   int power1 = std::pow(Ns, Dim - i);
   int power2 = std::pow(Ns, i - 1);

   Eigen::SparseMatrix<double> I1(power1, power1);
   I1.setIdentity();

   Eigen::SparseMatrix<double> I2(power2, power2);
   I2.setIdentity();

   Eigen::SparseMatrix<double> D_sparse(Ns, Ns);
   std::vector<Eigen::Triplet<double>> triplets;
   for (int j = 0; j < Ns; ++j) {
     triplets.push_back(Eigen::Triplet<double>(j, j, d[j]));
   }
   D_sparse.setFromTriplets(triplets.begin(), triplets.end());

   Eigen::SparseMatrix<double> kron_result = kroneckerProduct(I1, kroneckerProduct(D_sparse, I2));
   D_kron += kron_result;
 }

 // Step 4: Compute Lambda
 Eigen::VectorXd diag_D_kron = D_kron.diagonal();
 Eigen::VectorXd Lambda = Eigen::VectorXd::Ones(NsDim) - 0.5 * Dt / (Ds * Ds) * diag_D_kron;

 return Lambda;
}

Eigen::SparseMatrix<double> kroneckerProduct(const Eigen::SparseMatrix<double> &A, const Eigen::SparseMatrix<double> &B) {
 int rowsA = A.rows(), colsA = A.cols();
 int rowsB = B.rows(), colsB = B.cols();

 Eigen::SparseMatrix<double> result(rowsA * rowsB, colsA * colsB);
 std::vector<Eigen::Triplet<double>> triplets;

 for (int k = 0; k < A.outerSize(); ++k) {
   for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
     for (int l = 0; l < B.outerSize(); ++l) {
       for (Eigen::SparseMatrix<double>::InnerIterator jt(B, l); jt; ++jt) {
         int row = it.row() * rowsB + jt.row();
         int col = it.col() * colsB + jt.col();
         double value = it.value() * jt.value();
         if (value != 0.0) {
           triplets.push_back(Eigen::Triplet<double>(row, col, value));
         }
       }
     }
   }
 }

 result.setFromTriplets(triplets.begin(), triplets.end());
 return result;
}
