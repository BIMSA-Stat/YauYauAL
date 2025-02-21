#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @name ExpandGrid
//' @title Generate a Grid of All Combinations of Input Sequence
//' @description
//' This function generates a grid by expanding the input sequence \code{s} across the specified 
//' number of dimensions (\code{Dim}). The resulting grid contains all possible combinations of the 
//' elements in \code{s}, repeated across \code{Dim} dimensions. The grid is filled column by column, 
//' with the last column cycling through the values in \code{s} the fastest.
//' @usage
//' ExpandGrid(Dim, s)
//' @param Dim Integer. The number of dimensions for the grid. For example, \code{Dim = 2} creates a 2D grid.
//' @param s A numeric vector (\code{NumericVector}) containing the values to expand into the grid.
//' @details
//' This function creates a numeric matrix where each row represents a unique combination of the elements 
//' in \code{s} distributed across \code{Dim} dimensions. The total number of rows in the resulting matrix 
//' is \eqn{\text{length}(s)^\text{Dim}}, and the number of columns is equal to \code{Dim}. 
//' The combinations are generated such that the first column cycles through all combinations the slowest, 
//' while the last column cycles the fastest.
//' @return
//' A numeric matrix (\code{NumericMatrix}) where each row represents a combination of values from 
//' \code{s} across \code{Dim} dimensions. 
//' rows and \code{Dim} columns.
//' @examples
//' # Generate a 3D grid using values from 1 to 4
//' s <- seq(1, 4, 1)
//' Dim <- 3
//' grid <- ExpandGrid(Dim, s)
//' print(grid)
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix ExpandGrid(int Dim, Rcpp::NumericVector s) {
 int Ns = s.size();                      // Length of the input sequence
 int totalRows = std::pow(Ns, Dim);      // Total number of rows in the output grid
 
 arma::mat result(totalRows, Dim);       // Initialize result matrix using Armadillo
 
 // Fill the grid column by column
 for (int col = Dim - 1; col >= 0; col--) {
   int blockSize = std::pow(Ns, col);    // Determine block size for the current column
   for (int i = 0; i < totalRows; i++) {
     result(i, col) = s[(i / blockSize) % Ns];
   }
 }
 
 return Rcpp::wrap(result);              // Convert Armadillo matrix to Rcpp::NumericMatrix
}
