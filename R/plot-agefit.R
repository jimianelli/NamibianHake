#' Plot Age Fit
#'
#' This function creates a plot comparing observed and predicted age compositions for a specified type.
#'
#' @param x A list containing the observed and predicted age compositions. Default is `bc`.
#' @param title A character string for the plot title. Default is NULL.
#' @param type A character string specifying the type of composition to plot. Default is "fishery".
#' @param fage An integer specifying the first age to include in the plot. Default is 2.
#' @param lage An integer specifying the last age to include in the plot. Default is 7.
#'
#' @return A ggplot object comparing observed and predicted age compositions.
#'
#' @examples
#' \dontrun{
#' PlotAgeFit(x = bc, title = "Age Composition", type = "survey", fage = 1, lage = 8)
#' }
#'
#' @export
PlotAgeFit <- function(x = bc, title = NULL, type = "fishery", fage = 2, lage = 7) {
  obs <- x[[paste0(type, "_Pobs")]]
  pred <- x[[paste0(type, "_Phat")]]
  incl <- rowSums(data.frame(obs, src = "Obs")[, 2:7]) > 0
  dftmp <- rbind(
    data.frame(obs, src = "Obs")[incl, ],
    data.frame(pred, src = "Pred")[incl, ]
  )
  names(dftmp) <- c("Year", fage:lage, "type")
  x <- pivot_longer(dftmp,
    cols = 2:(lage + 2 - fage),
    names_to = "Age", values_to = "proportion"
  )
  ggplot(x |> filter(type == "Obs"), aes(x = Age, y = proportion)) +
    geom_bar(stat = "Identity", fill = "salmon") +
    geom_point(
      data = x |> filter(type == "Pred"),
      aes(x = Age, y = proportion), size = 2, shape = 3
    ) +
    facet_wrap(Year ~ .) +
    ggtitle(paste0(title, ", ", type))
}
