#' Color palettes used in outbreaker
#'
#' These functions are different color palettes (color-generating functions) used in outbreaker.
#'
#' @rdname palettes
#'
#' @aliases outbreaker.palettes chains.pal
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @param n a number of colors to be created
#'
#' @export
#'
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' plot(1:8, col=chains.pal(8), cex=10, pch=20)
#'
chains.pal <- function(n) {
    colorRampPalette(c("#660033", "#339966", "#cccc00", "#333399"))(n)
}
