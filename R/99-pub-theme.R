# ============================================================
# Publication-quality ggplot2 theme for survinger
# ============================================================

#' Publication-quality ggplot2 theme
#'
#' A clean, high-contrast theme suitable for journal publication.
#' Follows Nature/Lancet figure guidelines.
#'
#' @param base_size Base font size. Default 11.
#' @param base_family Font family. Default `""`.
#'
#' @return A ggplot2 theme object.
#'
#' @seealso [surv_compare_estimates()]
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_survinger()
#'
#' @export
theme_survinger <- function(base_size = 11, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "#E8E8E8", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "#333333", fill = NA,
                                           linewidth = 0.5),
      panel.background = ggplot2::element_rect(fill = "#FFFFFF"),
      axis.title = ggplot2::element_text(size = base_size, color = "#333333"),
      axis.text = ggplot2::element_text(size = base_size - 1, color = "#555555"),
      axis.ticks = ggplot2::element_line(color = "#999999", linewidth = 0.3),
      axis.ticks.length = ggplot2::unit(2, "pt"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = base_size - 1, face = "bold"),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.key.size = ggplot2::unit(14, "pt"),
      legend.background = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(t = 4, b = 0),
      plot.title = ggplot2::element_text(size = base_size + 2, face = "bold",
                                         color = "#1a1a1a", hjust = 0,
                                         margin = ggplot2::margin(b = 8)),
      plot.subtitle = ggplot2::element_text(size = base_size, color = "#666666",
                                            hjust = 0,
                                            margin = ggplot2::margin(b = 10)),
      plot.caption = ggplot2::element_text(size = base_size - 2, color = "#999999",
                                           hjust = 1, face = "italic"),
      strip.text = ggplot2::element_text(size = base_size, face = "bold",
                                         color = "#333333"),
      strip.background = ggplot2::element_rect(fill = "#F5F5F5", color = NA),
      plot.margin = ggplot2::margin(t = 10, r = 15, b = 10, l = 10)
    )
}

#' survinger color palette
#' @return A named character vector of hex colors.
#' @keywords internal
.surv_colors <- function() {
  c(
    primary    = "#2166AC",
    secondary  = "#D6604D",
    tertiary   = "#4DAF4A",
    quaternary = "#FF7F00",
    quinary    = "#984EA3",
    neutral    = "#878787",
    light_blue = "#92C5DE",
    light_red  = "#F4A582",
    light_green = "#A6D96A",
    bg_highlight = "#FFF7BC"
  )
}
