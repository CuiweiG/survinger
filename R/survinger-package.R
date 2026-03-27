#' @keywords internal
"_PACKAGE"

#' @importFrom rlang .data .env := %||%
#' @importFrom stats optim qnorm pnorm dnbinom pnbinom rpois rbinom
#'   rnbinom var median quantile sd weighted.mean dpois ppois plnorm
#'   dlnorm na.omit setNames
#' @importFrom methods is
#' @importFrom utils head tail
#' @importFrom generics tidy glance
NULL

#' @export
generics::tidy

#' @export
generics::glance
