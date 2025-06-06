% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{generateD}
\alias{generateD}
\title{Generate Sparse Matrices Using Kronecker Products}
\usage{
generateD(Dim, Ns, Ds)
}
\arguments{
\item{Dim}{Integer. The number of dimensions of the problem. This determines the number of matrices generated.}

\item{Ns}{Integer. The size of the grid in each dimension.}

\item{Ds}{Numeric. The grid spacing parameter, which scales the off-diagonal values in the matrices.}
}
\value{
A list of sparse matrices of type \code{Eigen::SparseMatrix<double>}. Each matrix corresponds to one dimension
and represents a discretized operator for that dimension.
}
\description{
This function generates a list of sparse matrices based on the specified dimension (\code{Dim}),
grid size (\code{Ns}), and spacing parameter (\code{Ds}). These matrices are constructed using
Kronecker products of identity matrices and sparse matrices representing differential operators.
They are useful for discretizing high-dimensional differential operators on regular grids.
}
\details{
The function constructs one sparse matrix for each dimension. For each dimension, the following steps are performed:
\enumerate{
\item Two identity matrices (\code{I1} and \code{I2}) are created based on the grid size and dimension index.
\item A sparse matrix is generated to represent the discretized 1D operator, with off-diagonal values scaled by \code{0.5 / Ds}.
\item The final sparse matrix for the dimension is constructed as a Kronecker product of these components.
}
The output is a list of sparse matrices, one for each dimension, which can be used to discretize high-dimensional problems.
}
\examples{
# Example: Generate sparse matrices for a 3D grid
Dim <- 3
Ns <- 4
Ds <- 0.5
D_list <- generateD(Dim, Ns, Ds)
# Check the structure of the first sparse matrix
print(D_list[[1]])
}
