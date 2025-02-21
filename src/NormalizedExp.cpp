#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @name NormalizedExp
//' @title Normalized Exponential Function
//' @description
//' Computes the normalized exponential of a numeric vector. To prevent numerical overflow,
//' the input values are shifted by subtracting the maximum value of the vector. The function
//' then calculates the exponential of each shifted value and normalizes the result so that
//' the sum of the output vector equals 1. This is commonly used for probability distributions
//' or softmax operations in machine learning and numerical analysis.
//'
//' @usage
//' NormalizedExp(x)
//'
//' @param x A numeric vector of type \code{arma::vec}. Contains the values to be processed
//'          by the normalized exponential function.
//'
//' @details
//' The function performs the following steps:
//' \enumerate{
//'   \item Finds the maximum value in the input vector \code{x}, denoted as \code{max_x}.
//'   \item Shifts the input vector by subtracting \code{max_x}, creating a new vector \code{x_shifted}.
//'   \item Computes the element-wise exponential of \code{x_shifted}.
//'   \item Normalizes the exponential values by dividing by their sum, ensuring the output
//'         lies in the range \eqn{[0, 1]} and the sum of the vector equals 1.
//' }
//' This normalization is particularly useful in scenarios where probabilities or relative
//' importance need to be derived from raw scores or unnormalized values.
//'
//' @return
//' A numeric vector of type \code{arma::vec}, where each element is the normalized exponential
//' value of the corresponding input element. The sum of the output vector is guaranteed to be 1.
//'
//' @examples
//' # Example: Normalizing a vector using the normalized exponential function
//' x <- c(1, 2, 3, 4)
//' result <- NormalizedExp(x)
//' print(result)
//'
//' # Example: Handling a vector with large values
//' y <- c(1000, 1001, 1002)
//' result <- NormalizedExp(y)
//' print(result)
//'
//' @export
// [[Rcpp::export]]
arma::vec NormalizedExp(const arma::vec &x) {
 // Step 1: Find the maximum value in the vector
 double max_x = x.max();
 
 // Step 2: Subtract the maximum value from each element to prevent overflow
 arma::vec x_shifted = x - max_x;
 
 // Step 3: Compute the exponential of each element
 arma::vec exp_x = arma::exp(x_shifted);
 
 // Step 4: Normalize by dividing by the sum of the exponentials
 arma::vec result = exp_x / arma::sum(exp_x);
 
 return result;  // Return the normalized vector in the range [0, 1]
}
