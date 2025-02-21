#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <numeric>  // For std::accumulate and std::inner_product

//' @name computeB
//' @title Compute Sparse Matrix B
//' @description
//' This function computes a sparse matrix \eqn{B} using input matrices, user-defined functions,
//' and sparse Kronecker-product-related operations. The matrix \eqn{B} is constructed by applying
//' transformations to input grid points and updating a diagonal matrix \eqn{Q} based on user-defined
//' functions and sparse matrices.
//'
//' @usage computeB(s, D_list, Dt, Ds, f, df, h)
//'
//' @param s A numeric matrix (\code{Eigen::MatrixXd}) representing grid points.
//'          Each row corresponds to a grid point, and each column corresponds to a dimension.
//' @param D_list A list of sparse matrices (\code{Eigen::SparseMatrix}) representing
//'               differential operators for each dimension.
//' @param Dt A numeric value specifying the time step parameter.
//' @param Ds A numeric value specifying the grid spacing parameter.
//' @param f A user-defined function applied to each row of \eqn{s}. It should return
//'          a numeric vector.
//' @param df A user-defined derivative function applied to each row of \eqn{s}.
//'           It should return a numeric vector.
//' @param h A user-defined function applied to each row of \eqn{s}. It contributes
//'          to the diagonal updates in matrix \eqn{Q}. It should return a numeric vector.
//'
//' @details
//' The function calculates the matrix \eqn{B} as follows:
//' \enumerate{
//'   \item Apply the user-defined function \eqn{f} to each row of \eqn{s} to compute
//'         the negated values stored in \eqn{p}.
//'   \item Construct a diagonal matrix \eqn{Q}, where each diagonal element is computed
//'         as the sum of the user-defined derivative function \eqn{df} and the inner product
//'         of the function \eqn{h} applied to the corresponding row of \eqn{s}.
//'   \item Update \eqn{Q} by adding contributions from \eqn{p} and the input list of
//'         sparse matrices \eqn{D_list}.
//'   \item Compute \eqn{B} as \eqn{I + Dt \cdot Q}, where \eqn{I} is an identity matrix.
//' }
//'
//' @return
//' A sparse matrix (\code{Eigen::SparseMatrix}) representing the computed \eqn{B} matrix.
//'
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> computeB(const Eigen::MatrixXd &s,
                                    const std::vector<Eigen::SparseMatrix<double>> &D_list,
                                    double Dt, double Ds,
                                    Rcpp::Function f, Rcpp::Function df, Rcpp::Function h) {

 int NsDim = s.rows();    // Total number of elements (dim * Ns)
 int Dim = s.cols();      // Number of dimensions

 // Step 1: Compute p = -t(apply(s, 1, f))
 Eigen::MatrixXd p(NsDim, Dim);
 for (int i = 0; i < NsDim; ++i) {
   Rcpp::NumericVector fx = f(s.row(i));
   for (int j = 0; j < Dim; ++j) {
     p(i, j) = -fx[j];
   }
 }

 // Step 2: Initialize Q as diagonal matrix with values calculated from df and h
 Eigen::VectorXd Q_diag(NsDim);
 for (int i = 0; i < NsDim; ++i) {
   Rcpp::NumericVector dfx = df(s.row(i));
   Rcpp::NumericVector hx = h(s.row(i));

   double df_sum = std::accumulate(dfx.begin(), dfx.end(), 0.0);
   double h_sum = 0.5 * std::inner_product(hx.begin(), hx.end(), hx.begin(), 0.0);
   Q_diag[i] = -(df_sum + h_sum);
 }

 Eigen::SparseMatrix<double> Q(NsDim, NsDim);
 Q.setIdentity();
 Q = Q * Q_diag.asDiagonal();

 // Step 3: Update Q with the contributions of D_list and p
 for (int i = 0; i < Dim; ++i) {
   Eigen::SparseMatrix<double> P_sparse(NsDim, NsDim);
   std::vector<Eigen::Triplet<double>> triplets;

   for (int j = 0; j < NsDim; ++j) {
     if (p(j, i) != 0) {
       triplets.push_back(Eigen::Triplet<double>(j, j, p(j, i)));
     }
   }

   P_sparse.setFromTriplets(triplets.begin(), triplets.end());
   Q += P_sparse * D_list[i];
 }

 // Step 4: Compute B = Diagonal(NsDim) + Dt * Q
 Eigen::SparseMatrix<double> B(NsDim, NsDim);
 B.setIdentity();
 B += Dt * Q;

 return B;
}
