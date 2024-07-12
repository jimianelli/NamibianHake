#' Plot Selectivity
#'
#' This function creates a plot of selectivity across years and ages.
#'
#' @param Year A vector of years. Default is `M$Year`.
#' @param sel A matrix of selectivity values. Default is `M[,2:10]`.
#' @param styr An integer specifying the start year for the plot. Default is 1964.
#' @param fage An integer specifying the first age to include in the plot. Default is NULL.
#' @param lage An integer specifying the last age to include in the plot. Default is NULL.
#' @param alpha A numeric value specifying the transparency level of the fill. Default is 0.2.
#' @param scale A numeric value specifying the scale for the density ridges. Default is 3.8.
#' @param fill A character string specifying the fill color for the density ridges. Default is "purple".
#'
#' @return A ggplot object showing selectivity across years and ages.
#'
#' @examples
#' \dontrun{
#' plot_sel(Year = M$Year, sel = M[,2:10], styr = 1970, fage = 1, lage = 8)
#' }
#'
#' @export
plot_sel <- function(Year = M$Year, sel = M[,2:10], styr = 1964, fage = NULL, lage = NULL, alpha = 0.2, scale = 3.8, fill = "purple") {
  df <- data.frame(Year = Year, sel = sel)
  if (is.null(fage)) fage <- 0
  if (is.null(lage)) lage <- length(sel[1, ]) - 1
  df <- df |> select(1:(lage - fage + 2))
  names(df) <- c("Year", fage:lage)
  nages <- length(fage:lage)
  sdf <- pivot_longer(df, names_to = "age", values_to = "sel", cols = 2:(nages + 1)) %>% 
    filter(Year >= styr) %>% 
    mutate(age = as.numeric(age))
  p1 <- ggplot(sdf, aes(x = age, y = as.factor(Year), height = sel)) + 
    geom_density_ridges(stat = "identity", scale = scale, alpha = alpha, fill = fill, color = "black") + 
    ggthemes::theme_few() +
    ylab("Year") + xlab("Age (years)") +
    scale_x_continuous(limits = c(fage, lage), breaks = fage:lage) +
    scale_y_discrete(limits = rev(levels(as.factor(sdf$Year))))
  return(p1)
}
