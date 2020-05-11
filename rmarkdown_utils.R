#' Function designed to extract model effects from glm models. Takes a model 
#' object, a target variable/effect name, and an boolean argument to return 
#' odds ratio or not. Default for odds ratio is FALSE, and function will print
#' out preformatted *b* = <coef>, *p* = <p-value> construction that renders in 
#' RMarkdown format. When odds_ratio=TRUE is enable will print out 
#' OR = <odds ratio> 

extract_glm_logistic_effects <- function(model, effect_name, odds_ratio=FALSE){
  require(glue)
  b <- round(coef(summary(model))[effect_name, 1], digits = 2)
  p <- round(coef(summary(model))[effect_name, 4], digits = 3)
  OR <- round(exp(b), digits = 2)
  
  if(p >= .001){
    out <- glue("*b* = {b}, *p* = {p}")
  }
  
  else{
    out <- glue("*b* = {b}, *p* < .001")
  }
  
  if(odds_ratio){
    out <- paste(out, glue("(OR = {OR})"))
  }
  
  out
}


#' Function for returning list of glm model objects for each x_var value. The 
#' x_var argument takes a vector of variable names. The y_var argument takes 
#' the name of a single dependent variable. Defaults to running a combined model
#' including all x_vars as well. Also defaults to printing out latex version of 
#' a stargazer model summary table. 
run_glm_w_multiple_vars <- function(data, y_var, x_vars, family, combine=TRUE, 
				    y_name=NULL, x_names=NULL, table_title=NULL,
				    print_stargazer=TRUE){
  model_list <- list()
  for(v in 1:length(x_vars)){
    mod_form <- as.formula(
      paste(y_var, x_vars[v], sep = "~")
    ) 
  
    model_list[[x_names[v]]] <- glm(mod_form, data=data, family=family)
  }

  if(combine){
    mod_form <- as.formula(
      paste(y_var, paste(x_vars, collapse = "+"), sep="~")
    )

    model_list[["combined"]] <- glm(mod_form, data=data, family=family)
  }

  if(print_stargazer){
    stargazer::stargazer(model_list, title = table_title,
                         star.char = c("*", "**", "***"), 
  	   	         star.cutoffs = c(.05, .01, .001),
                         dep.var.labels = y_name,
                         covariate.labels = x_names,
                         header=FALSE)
  }
  return(model_list)
}
