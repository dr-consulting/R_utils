#' Function designed to extract model effects from glm models. Takes a model 
#' object, a target variable/effect name, and an boolean argument to return 
#' odds ratio or not. Default for odds ratio is FALSE, and function will print
#' out preformatted *b* = <coef>, *p* = <p-value> construction that renders in 
#' RMarkdown format. When odds_ratio=TRUE is enable will print out 
#' OR = <odds ratio> 

extract_glm_logistic_effects <- function(model, effect_name, odds_ratio=FALSE){
  b <- round(coef(summary(model))[effect_name, 1], digits = 2)
  p <- round(coef(summary(model))[effect_name, 4], digits = 3)
  OR <- round(exp(b), digits = 2)

  if(!odds_ratio){
    if(p >= .001){
      glue("*b* = {b}, *p* = p")
    }

    else{
      glue("*b* = {b}, *p* < {p}")
    }
  }

  else{
    glue("OR = {OR}")
  }
}
