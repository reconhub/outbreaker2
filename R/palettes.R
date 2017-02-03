#' Color palettes used in outbreaker
#'
#' These functions are different color palettes (color-generating functions)
#' used in outbreaker.
#'
#' @rdname palettes
#'
#' @aliases outbreaker_palettes chains_pal
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail_com})
#'
#' @param n a number of colors to be created
#'
#' @export
#'
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' plot(1:8, col = chains_pal(8), cex = 10, pch = 20)
#'
chains_pal <- function(n) {
    colorRampPalette(c("#660033", "#339966", "#cccc00", "#333399"))(n)
}

#' @rdname palettes
#' @export
#' @aliases cases_pal
cases_pal <- function(n) {
    ## This was the viridis palette
    ## cols <- c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF",
    ##           "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF",
    ##           "#B4DE2CFF", "#FDE725FF")


    ## This one is taken from epicontacts
    cols <- c("#ccddff", "#79d2a6", "#ffb3b3", "#a4a4c1","#ffcc00", "#ff9f80",
              "#ccff99", "#df9fbf","#ffcc99", "#cdcdcd")
    
    colorRampPalette(cols)(n)
                
}
