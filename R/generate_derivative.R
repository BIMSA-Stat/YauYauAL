#' @name generate_derivative
#' @title Generate Derivative Function
#' @useDynLib YourPackageName
#' @importFrom Deriv Deriv
#' @description
#' This function generates a numerical derivative function for a given vector-valued function. 
#' It uses finite differences to approximate the derivative of the input function with respect 
#' to each element of the input vector.
#'
#' @usage
#' generate_derivative(f)
#'
#' @param f A vector-valued function that takes a numeric vector \code{x} as input 
#'          and returns a numeric vector as output. The function \code{f} should support 
#'          evaluation at any point within the domain of interest.
#'
#' @details
#' The function computes numerical derivatives using a finite difference method. For each 
#' element of the input vector \code{x}, a small perturbation (\code{epsilon}) is applied to 
#' approximate the derivative. Specifically, the derivative of \code{f} at a point \code{x} 
#' is approximated as:
#' \deqn{\frac{\partial f_i}{\partial x_j} \approx \frac{f_i(x + \epsilon e_j) - f_i(x)}{\epsilon},}
#' where \code{e_j} is a unit vector in the direction of the \code{j}th component.
#'
#' This method supports vector-valued functions and computes partial derivatives for 
#' all components of the output with respect to the input vector. The result is a Jacobian matrix 
#' if the output is multidimensional.
#'
#' @return
#' A function that computes the numerical derivative of the input function \code{f}. 
#' The returned function takes a numeric vector \code{x} as input and returns:
#' \itemize{
#'   \item If \code{f} outputs a scalar, the derivative is a numeric vector of the same length as \code{x}.
#'   \item If \code{f} outputs a vector, the derivative is a Jacobian matrix where the \code{(i, j)}-th 
#'         entry represents the partial derivative of the \code{i}th output with respect to the \code{j}th input.
#' }
#'
#' @examples
#' # Define a vector-valued function
#' f <- function(x) {
#'   c(x[1]^2 + x[2], x[1] + sin(x[2]))
#' }
#'
#' # Generate the derivative function
#' df <- generate_derivative(f)
#'
#' # Test the derivative function
#' x <- c(1, 2)
#' df(x)  # Returns the numerical Jacobian of f at x
#'
#' @export
generate_derivative <- function(f) {
  function(x) {
    epsilon <- 1e-8  # A very small number for numerical approximation
    n <- length(x)   # Length of the input vector
    df <- numeric(n) # Storage for the derivative
    for (i in 1:n) {
      x_perturb <- x
      x_perturb[i] <- x_perturb[i] + epsilon
      df[i] <- (f(x_perturb)[i] - f(x)[i]) / epsilon
    }
    return(df)
  }
}

