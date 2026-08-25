#' Demonstration dataset with known careless respondents
#'
#' A simulated multi-construct questionnaire (15 factors, 6 items each, 90 items)
#' with 400 respondents, 20\% of whom were made careless with
#' \code{\link{inject_careless}}. Because it is simulated it carries an exact
#' ground-truth label, which makes it convenient for examples and for
#' \code{\link{benchmark_indices}}.
#'
#' @format A list with four elements:
#' \describe{
#'   \item{responses}{A 400 by 90 integer matrix of 1--5 Likert responses.}
#'   \item{careless}{A 0/1 vector: 1 for the injected careless respondents.}
#'   \item{severity}{The corrupted fraction of each careless respondent (0 for clean).}
#'   \item{pattern}{The careless pattern applied to each respondent (\code{NA} for clean).}
#' }
#' @examples
#' fit <- reread(demo_careless$responses)
#' table(flagged = fit$flagged, truth = demo_careless$careless)
"demo_careless"
