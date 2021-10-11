# From "clauswilke/colorblindr" --- copied because the package does not work anymore with R 4.1 (Aug 2021)
palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

palette_OkabeIto_black <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

scale_colour_OkabeIto <- function(aesthetics = "colour", ...) {
  scale_OkabeIto(aesthetics, ...)
}

scale_color_OkabeIto <- scale_colour_OkabeIto

scale_fill_OkabeIto <- function(aesthetics = "fill", ...) {
  scale_OkabeIto(aesthetics, ...)
}

scale_OkabeIto <- function(aesthetics, use_black = FALSE, order = 1:8, darken = 0, alpha = NA, ...) {
  if (use_black) {
    values <- palette_OkabeIto_black[order]
  }
  else {
    values <- palette_OkabeIto[order]
  }

  n <- length(values)
  darken <- rep_len(darken, n)
  alpha <- rep_len(alpha, n)

  di <- darken > 0
  if (sum(di) > 0) { # at least one color needs darkening
    values[di] <- colorspace::darken(values[di], amount = darken[di])
  }

  li <- darken < 0
  if (sum(li) > 0) { # at least one color needs lightening
    values[li] <- colorspace::lighten(values[li], amount = -1*darken[li])
  }

  ai <- !is.na(alpha)
  if (sum(ai) > 0) { # at least one color needs alpha
    values[ai] <- scales::alpha(values[ai], alpha[ai])
  }

  pal <- function(n) {
    if (n > length(values)) {
      warning("Insufficient values in manual scale. ", n, " needed but only ",
           length(values), " provided.", call. = FALSE)
    }
    values
  }
  ggplot2::discrete_scale(aesthetics, "manual", pal, ...)
}
