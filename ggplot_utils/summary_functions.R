#' Summary function that returns a mean, and a min/max defined by user specified quanitles
#' Defaults to return values along the y-axis, but can be changed to x if y_axis=FALSE.

dist_mean_and_quantiles <- function(x, lower_quantile=.025, upper_quantile=.975, y_axis=TRUE){
  m <- mean(x)

  if(y_axis){
    ymin <- quantile(x, lower_quantile)
    ymax <- quantile(x, upper_quantile)
    return(c(y=m, ymin=ymin, ymax=ymax))
  }

  else{
    xmin <- quantile(x, lower_quantile)
    xmax <- quantile(x, upper_quantile)
    return(c(x=m, xmin=xmin, xmax=xmax))
  }
}


