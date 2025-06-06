% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{DST_Solver}
\alias{DST_Solver}
\title{Solve a System Using the Discrete Sine Transform (DST-II)}
\usage{
DST_Solver(Dim, Lambda, B, U, Ns)
}
\arguments{
\item{Dim}{Integer. The dimension of the problem. For example, \code{Dim = 1} solves a one-dimensional problem.}

\item{Lambda}{A numeric vector of type \code{Eigen::VectorXd} representing scaling factors (eigenvalues).}

\item{B}{A sparse matrix of type \code{Eigen::SparseMatrix<double>} representing a linear operator in the system.}

\item{U}{A numeric vector of type \code{Eigen::VectorXd} representing the initial state vector.}

\item{Ns}{Integer. The size of each segment to be processed by DST-II.}
}
\value{
A numeric vector of type \code{Eigen::VectorXd} representing the final normalized result after applying DST-II twice and scaling.
}
\description{
This function solves a system using the Discrete Sine Transform (DST-II).
It computes transformations, applies scaling with eigenvalues, and normalizes the result.
}
\details{
The function performs the following steps:
\enumerate{
\item Computes the product \eqn{B \cdot U}, where \eqn{B} is a sparse matrix and \eqn{U} is the input vector.
\item Applies the Discrete Sine Transform (DST-II) to the product.
\item Scales the transformed values by dividing by the vector \eqn{\Lambda}.
\item Applies DST-II again to the scaled values.
\item Ensures all values in the result are non-negative and normalizes the vector so that its sum is 1.
}
It handles one-dimensional problems as a special case and generalizes to higher dimensions by processing segments independently.
}
\examples{
library(Matrix)
Dim <- 1
Lambda <- c(1, 2, 3, 4)
B <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(0.5, 0.75), dims = c(4, 4))
U <- c(1, 2, 3, 4)
Ns <- 4
result <- DST_Solver(Dim, Lambda, B, U, Ns)
print(result)

}
