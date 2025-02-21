#' @name wrap_outiu_function
#' @title Function to Wrap and Update Output `Iu` with Progress Bar
#' @description
#' This function computes and updates a matrix `Iu` based on a sequence of input vectors, iterative updates, 
#' and transformations. It uses a progress bar to display the computation progress and performs differential 
#' equation solving through the `DST_Solver` function while normalizing values using the `NormalizedExp` function.
#'
#' @usage
#' wrap_outiu_function(s, NtNtau, Ntau, Nt, Dim, y, h, Lambda, B, Ns, NormalizedExp, DST_Solver)
#'
#' @param s A numeric matrix of size \code{(number of sequences x Dim)} representing input state sequences.
#' @param NtNtau Integer. The total number of time steps, including the \code{Ntau} steps.
#' @param Ntau Integer. The number of tau steps in the iterative process.
#' @param Nt Integer. The number of regular time steps within each tau interval.
#' @param Dim Integer. The dimension of the state space.
#' @param y A numeric matrix representing observed values at each time step.
#' @param h A user-defined function that computes the relationship between `s` and the observations `y`. 
#'          It should accept a numeric vector and return a numeric vector.
#' @param Lambda A numeric vector of parameters used within the \code{DST_Solver} function.
#' @param B A numeric matrix used as a coefficient matrix in the \code{DST_Solver} function.
#' @param Ns Integer. The sample size or scaling parameter used in the \code{DST_Solver}.
#' @param NormalizedExp A function that normalizes exponential values for the input sequences.
#' @param DST_Solver A function that solves a differential equation given \code{Dim}, \code{Lambda}, \code{B}, 
#'        and other inputs. It should return a numeric vector of the same length as the input.
#'
#' @details
#' The function iteratively updates the `Iu` matrix as follows:
#' \enumerate{
#'   \item Initializes the state vector \code{U} using an exponential function and normalizes it.
#'   \item Updates \code{U} at each time step by applying the `DST_Solver` and the observation function `h`.
#'   \item Computes weighted sums of `U` and `s` to populate the \code{Iu} matrix.
#'   \item Utilizes a progress bar to display the computation progress for user feedback.
#' }
#'
#' The function ensures numerical stability by normalizing \code{U} at each step and handles multi-dimensional 
#' state spaces by processing each dimension independently.
#'
#' @return
#' A numeric matrix \code{Iu} of size \code{NtNtau x Dim}, where each row corresponds to a time step and each 
#' column corresponds to a dimension of the state space. The values represent the updated state trajectory.
#'
#' @examples
#' # Example setup
#' s <- matrix(runif(100), nrow = 10, ncol = 10)
#' NtNtau <- 50
#' Ntau <- 10
#' Nt <- 5
#' Dim <- 10
#' y <- matrix(runif(100), nrow = Ntau, ncol = Dim)
#'
#' # Define example functions
#' h <- function(x) x^2
#' NormalizedExp <- function(x) exp(x - max(x)) / sum(exp(x - max(x)))
#' DST_Solver <- function(Dim, Lambda, B, U, Ns) {
#'   U * Lambda / sum(Lambda)
#' }
#'
#' # Define parameters
#' Lambda <- runif(Dim)
#' B <- diag(Dim)
#' Ns <- 10
#'
#' # Compute Iu
#' Iu <- wrap_outiu_function(s, NtNtau, Ntau, Nt, Dim, y, h, Lambda, B, Ns, NormalizedExp, DST_Solver)
#' print(Iu)
#'
#' @export
wrap_outiu_function <- function(s, NtNtau, Ntau, Nt, Dim, y, h, Lambda, B, Ns, NormalizedExp, DST_Solver) {
  outiu <- function(order) {
    ss <- s
    ss[, 1] <- s[, order]
    ss[, order] <- s[, 1]
    s <- ss
    
    # Initialize sigma0 and iterate to update Iu
    sigma0 <- exp(-10 * rowSums(s^2))
    
    Iu <- matrix(0, nrow = NtNtau, ncol = Dim)
    Idx <- 1
    U <- sigma0
    U <- U / sum(U)
    
    Iu[Idx, ] <- colSums(U * s)
    
    # Initialize progress bar
    pb <- utils::txtProgressBar(min = 0, max = NtNtau, style = 3)
    
    for (jj in 1:Ntau) {
      if (Idx >= NtNtau) break
      Idx <- Idx + 1
      tmp <- if (jj == 1) y[jj, ] else y[jj, ] - y[jj - 1, ]
      if (Dim == 1) {
        U <- NormalizedExp(rowSums(matrix(t(apply(s, 1, h))) * matrix(tmp, nrow = nrow(s), ncol = length(tmp), byrow = TRUE))) * U
      } else {
        U <- NormalizedExp(rowSums(t(apply(s, 1, h)) * matrix(tmp, nrow = nrow(s), ncol = length(tmp), byrow = TRUE))) * U
      }
      U <- DST_Solver(Dim, Lambda, B, U, Ns)
      U <- U / sum(U)
      U_matrix <- matrix(U, nrow = nrow(s), ncol = Dim, byrow = FALSE)
      Iu[Idx, ] <- colSums(U_matrix * s)
      utils::setTxtProgressBar(pb, Idx)
      
      for (ii in 2:Nt) {
        if (Idx >= NtNtau) break
        Idx <- Idx + 1
        U <- DST_Solver(Dim, Lambda, B, U, Ns)
        U_matrix <- matrix(U, nrow = nrow(s), ncol = Dim, byrow = FALSE)
        Iu[Idx, ] <- colSums(U_matrix * s)
        utils::setTxtProgressBar(pb, Idx)
      }
    }
    
    close(pb)  # Close progress bar
    return(Iu[, order])
  }
  
  Iu <- sapply(1:Dim, outiu)
  Iu <- matrix(Iu, nrow = NtNtau, ncol = Dim)
  return(Iu)
}
