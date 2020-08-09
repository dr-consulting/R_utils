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

#' Generates a simple set of ouptuts from a bivariate regression model that can be added to a plot. Useful for 
#' generating r-squared values, or information about the fit line. 
regression <- function(df, x, y) {
  # setting the regression function. 
  form<-as.formula(
    paste(y, "~", x)
  )
  reg_fun<-lm(data = df, 
              formula = form)  #regression function
  # getting the slope, intercept, R square and adjusted R squared of 
  # the regression function (with 3 decimals).
  slope<-round(coef(reg_fun)[2],3)  
  intercept<-round(coef(reg_fun)[1],3) 
  R2<-round(as.numeric(summary(reg_fun)[8]),3)
  r<-round(sqrt(R2)*slope/abs(slope), digits = 3)
  R2.Adj<-round(as.numeric(summary(reg_fun)[9]),3)
  p_val<-round(as.numeric(summary(reg_fun)$coefficients[2,4]), digits = 3)
  tmp.DF<-data.frame(slope, 
                     intercept, 
                     R2, 
                     r, 
                     R2.Adj, 
                     p_val)
  return(tmp.DF)
}

#' Generates a label_position based on the min and max of each axis and the relative location along that axis
#' 0 = min, 1 = max values of each axis. 
label_pos <- function(min_val, max_val, rel_location){
  range <- max_val - min_val
  delta <- range*rel_location
  return(min_val + delta)
}
