#' Run the state simulation from C++ side
#'
#' @param Dt numeric, time step
#' @param Ntau integer, number of observation points
#' @param NtNtau integer, total steps
#' @param f function, dynamics
#' @param h function, observation
#' @param Dim integer, dimension
#' @param seed integer, random seed
#'
#' @return list with elements \code{x}, \code{y}, and \code{seed}
#' @export
simulate_wrapper <- function(Dt, Ntau, NtNtau, f, h, Dim, seed = NULL) {
  # 直接调用 Rcpp函数
  Simulate_State_Obser(Dt, Ntau, NtNtau, f, h, Dim, seed)
}