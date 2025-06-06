% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{computeLambda}
\alias{computeLambda}
\title{Compute Lambda Vector}
\usage{
computeLambda(Dim, Ns, Dt, Ds)
}
\arguments{
\item{Dim}{Integer. The dimension of the Kronecker product.}

\item{Ns}{Integer. The size of the grid in each dimension.}

\item{Dt}{Numeric. The time step size.}

\item{Ds}{Numeric. The space step size.}
}
\value{
A numeric vector of type \code{Eigen::VectorXd} containing the scaled eigenvalues (\code{Lambda}).
The length of the vector is \eqn{Ns^\text{Dim}}.
}
\description{
Computes the Lambda vector using Kronecker products to construct a sparse matrix.
The Lambda vector represents scaled eigenvalues derived from the diagonal of the Kronecker sum
of a 1D Laplacian operator in multiple dimensions.
}
\details{
The function first computes a vector of eigenvalues \eqn{d} for a 1D Laplacian operator.
It then constructs a sparse matrix \eqn{D_\text{kron}} by summing the Kronecker products
of these eigenvalues across all dimensions. Finally, it computes the Lambda vector by scaling
the diagonal elements of \eqn{D_\text{kron}} using the parameters \code{Dt} and \code{Ds}.
}
\examples{
Dim <- 3
Ns <- 4
Dt <- 0.01
Ds <- 0.5
Lambda <- computeLambda(Dim, Ns, Dt, Ds)
print(Lambda)

}
