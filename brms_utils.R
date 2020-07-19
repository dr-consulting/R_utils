#'  cumulativemodelfit identifies optimal link functions properties for cumulative prior
#'  
#'  Adapted from: 
#'  https://kevinstadler.github.io/blog/bayesian-ordinal-regression-with-random-effects-using-brms/
cumulativemodelfit <- function(y_var, data, links=c("logit", "probit", "cloglog", "cauchit"),
                               thresholds=c("flexible", "equidistant"), verbose=TRUE) {
    require(glue)
    require(ordinal)
    
    # prep data and formula
    data[[y_var]] <- ordered(data[[y_var]])
    form <- formula(glue("{y_var} ~ 1"))
    
    # generate table structure and LL
    names(links) <- links
    names(thresholds) <- thresholds
    llks <- outer(links, thresholds,
                  Vectorize(function(link, threshold)
                      # catch error for responses with 2 levels
                      tryCatch(clm(form, data=data, link=link, threshold=threshold)$logLik,
                               error = function(e) NA)))
    print(llks)
    if (verbose) {
        bestfit <- which.max(llks)
        cat("\nThe best link function is ", links[bestfit %% length(links)], " with a ",
            thresholds[1 + bestfit %/% length(thresholds)], " threshold (logLik ", llks[bestfit],
            ")\n", sep="")
    }
    invisible(llks)
    # save optimal settings for use in defining model priors
    outputs <- list(link = links[bestfit %% length(links)], 
                    threshold = thresholds[1 + bestfit %/% length(thresholds)])
    return(outputs)
}


#' Lightweight wrapper around brms binomial model
brms_binomial_wrapper <- function(y_var, rhs, data, iter, warmup, chains, control_list, priors, out_dir, out_base_name){
    require(glue)
    require(brms)
    form <- glue("{y_var} ~ {rhs}")
    
    model <- bf(form)+bernoulli()
    
    fit <- brm(model, data=data, prior=priors, iter=iter, warmup=warmup, chains=chains, cores=chains, 
               control=control_list)
    print(fit)
    
    sink(glue("{out_dir}/{out_base_name}_{y_var}_binomial_model_summary_{Sys.Date()}.txt"))
    print(fit, digits = 5)
    sink()
    
    return(fit)
}


#' Lighweight wrapper around brms poisson model
brms_poisson_wrapper <- function(y_var, rhs, data, iter, warmup, chains, control_list, priors, out_dir, out_base_name){
    require(glue)
    require(brms)
    form <- glue("{y_var} ~ {rhs}")
    
    model <- bf(form)+poisson()
    
    fit <- brm(model, data=data, prior=priors, iter=iter, warmup=warmup, chains=chains, cores=chains, 
               control=control_list)
    print(fit)
    
    sink(glue("{out_dir}/{out_base_name}_{y_var}_poisson_model_summary_{Sys.Date()}.txt"))
    print(fit, digits = 5)
    sink()
    
    return(fit)
}


#' Lightweight wrapper around brms ordinal model
brms_ordinal_wrapper <- function(y_var, rhs, data, iter, warmup, chains, control_list, priors, out_dir, out_base_name){
    require(glue)
    require(brms)
    
    # find the appropriate cumulative link function settings: 
    cum_link_list <- cumulativemodelfit(y_var, data=data)
    
    form <- glue("{y_var} ~ {rhs}")
    
    model <- bf(form)+cumulative(link = cum_link_list$link[[1]], threshold = cum_link_list$threshold[[1]])
    
    # repetitive for now but easiest fix: 
    data[[y_var]] <- ordered(data[[y_var]])
    
    fit <- brm(model, data=data, prior=priors, iter=iter, warmup=warmup, chains=chains, cores=chains, 
               control=control_list)
    print(fit)
    
    sink(glue("{out_dir}/{out_base_name}_{y_var}_ordinal_model_summary_{Sys.Date()}.txt"))
    print(fit, digits = 5)
    sink()
    
    return(fit)
}


#' generates a grand-intercept to "ground" the time 1 model and speed up estimate times
estimate_intercept_in_logits <- function(y_var, y_hit, data, cond_var=NULL, cond_value=NULL){
    if(is.null(cond_var) | is.null(cond_value)){
       print("Did not pass both a conditional value and variable. Estimating intercept in logits for grand mean.")
       y_obs <- data[[y_var]] 
    }
    else if(!is.null(cond_var) & !is.null(cond_value)){ 
        y_obs <- data[[y_var]][data[[cond_var]] == cond_value]
    }
    else { 
        warning("Incorrect settings passed to function")
    }
    
    y_prop <- sum(y_obs == y_hit, na.rm = TRUE)/sum(!is.na(y_obs))
    y_logit <- log(y_prop/(1-y_prop))
    return(y_logit)
}