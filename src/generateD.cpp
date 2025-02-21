#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <vector>
#include <cmath>
#include <Eigen/Sparse>

//' @name generateD
//' @title Generate Sparse Matrices Using Kronecker Products
//' @description
//' This function generates a list of sparse matrices based on the specified dimension (\code{Dim}),
//' grid size (\code{Ns}), and spacing parameter (\code{Ds}). These matrices are constructed using 
//' Kronecker products of identity matrices and sparse matrices representing differential operators. 
//' They are useful for discretizing high-dimensional differential operators on regular grids.
//' @usage
//' generateD(Dim, Ns, Ds)
//' @param Dim Integer. The number of dimensions of the problem. This determines the number of matrices generated.
//' @param Ns Integer. The size of the grid in each dimension.
//' @param Ds Numeric. The grid spacing parameter, which scales the off-diagonal values in the matrices.
//' @details
//' The function constructs one sparse matrix for each dimension. For each dimension, the following steps are performed:
//' \enumerate{
//'   \item Two identity matrices (\code{I1} and \code{I2}) are created based on the grid size and dimension index.
//'   \item A sparse matrix is generated to represent the discretized 1D operator, with off-diagonal values scaled by \code{0.5 / Ds}.
//'   \item The final sparse matrix for the dimension is constructed as a Kronecker product of these components.
//' }
//' The output is a list of sparse matrices, one for each dimension, which can be used to discretize high-dimensional problems.
//' @return
//' A list of sparse matrices of type \code{Eigen::SparseMatrix<double>}. Each matrix corresponds to one dimension
//' and represents a discretized operator for that dimension.
//' @examples
//' # Example: Generate sparse matrices for a 3D grid
//' Dim <- 3
//' Ns <- 4
//' Ds <- 0.5
//' D_list <- generateD(Dim, Ns, Ds)
//' # Check the structure of the first sparse matrix
//' print(D_list[[1]])
//' @export
// [[Rcpp::export]]
std::vector<Eigen::SparseMatrix<double>> generateD(int Dim, int Ns, double Ds) {
 std::vector<Eigen::SparseMatrix<double>> D_list(Dim);
 
 for (int i = 0; i < Dim; ++i) {
   int power1 = std::pow(Ns, Dim - i - 1); // Size of I1
   int power2 = std::pow(Ns, i);           // Size of I2
   
   // Create identity matrices I1 and I2
   Eigen::SparseMatrix<double> I1(power1, power1);
   I1.setIdentity();
   
   Eigen::SparseMatrix<double> I2(power2, power2);
   I2.setIdentity();
   
   // Create sparse matrix for 1D operator
   Eigen::SparseMatrix<double> D_sparse(Ns, Ns);
   std::vector<Eigen::Triplet<double>> tripletList;
   for (int j = 0; j < Ns - 1; ++j) {
     tripletList.emplace_back(j, j + 1, 0.5 / Ds);  // Upper diagonal
     tripletList.emplace_back(j + 1, j, -0.5 / Ds); // Lower diagonal
   }
   D_sparse.setFromTriplets(tripletList.begin(), tripletList.end());
   
   // Compute the Kronecker product for the current dimension
   Eigen::SparseMatrix<double> kronResult(power1 * Ns, power1 * Ns);
   std::vector<Eigen::Triplet<double>> tripletsKron;
   
   for (int k = 0; k < I1.outerSize(); ++k) {
     for (Eigen::SparseMatrix<double>::InnerIterator it(I1, k); it; ++it) {
       for (int l = 0; l < D_sparse.outerSize(); ++l) {
         for (Eigen::SparseMatrix<double>::InnerIterator jt(D_sparse, l); jt; ++jt) {
           int row = it.row() * Ns + jt.row();
           int col = it.col() * Ns + jt.col();
           double value = it.value() * jt.value();
           if (value != 0.0) {
             tripletsKron.emplace_back(row, col, value);
           }
         }
       }
     }
   }
   kronResult.setFromTriplets(tripletsKron.begin(), tripletsKron.end());
   
   // Add the Kronecker product result to the list
   D_list[i] = kronResult;
 }
 
 return D_list;
}
