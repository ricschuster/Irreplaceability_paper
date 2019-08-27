#' Time an expression or code chunk
#'
#' Time an expression or chunk of R code with `system.time()` and return both
#' the results of evaluating the expression and the execution time in a list.
#'
#' @param code An expression or code chunk (wrapped in curly brackets) to time.
#' @param time_as character; return the execution time and results as list
#'   elements ("list") or return results with time as an attribute ("attr").
#'
#' @return A list with elements `result` and `time`, containing the results of
#'   evaluating `code` and the execution time, respectively. If `time_as =
#'   "attr"`, the result will be returned directly with the time stored as an
#'   attribute.
#' @export
#' @examples
#' # time a simple expression
#' time_this(mean(rnorm(1e6)))
#'
#' # time a code chunk
#' time_this({
#'   x <- (1:1e6) / 1e6
#'   y = 2 * x + 1 + rnorm(x)
#'   lm(y ~ x)
#' })
#'
#' # time can also be stored as an attribute
#' time_this(mean(rnorm(1e6)), time_as = "attr")
time_this <- function(code, time_as = c("list", "attr")) {
  time_as <- match.arg(time_as)
  systime <- system.time({result = code})
  if (time_as == "list") {
    return(list(result = result, time = systime))
  } else {
    attr(result, "time") <- systime
    return(result)
  }
}

#' Wrap a function to add a timer
#'
#' Modify a function so its execution is timed with `system.time()` and both
#' the results of the funciton call and the execution time are returned as a
#' list.
#'
#' @param .f A function to modify to add a timer.
#' @param time_as character; return the execution time and results as list
#'   elements ("list") or return results with time as an attribute ("attr").
#'
#' @return Wrapped function that returns a list with components `result` and
#'   `time`. If `time_as = "attr"`, the result will be returned directly with
#'   the time stored as an attribute.
#' @export
#' @examples
#' mean_timed <- add_timer(mean)
#' mean_timed(rnorm(1e6))
#'
#' lm_timed <- add_timer(lm)
#' x <- (1:1e6) / 1e6
#' y = 2 * x + 1 + rnorm(x)
#' lm_timed(y ~ x)
#'
#' # time can also be stored as an attribute
#' mean_timed_a <- add_timer(mean, time_as = "attr")
#' mean_timed_a(rnorm(1e6))
add_timer <- function(.f, time_as = c("list", "attr")) {
  time_as <- match.arg(time_as)
  function(...) {
    time_this(.f(...), time_as = time_as)
  }
}
