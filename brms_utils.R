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
brms_binomial_wrapper <- function(y_var, rhs, data, iter, warmup, chains, control_list, priors, out_dir, 
                                  out_base_name, refresh) {
    
    require(glue)
    require(brms)
    form <- glue("{y_var} ~ {rhs}")
    
    model <- bf(form)+bernoulli()
    
    fit <- brm(model, data=data, iter=iter, warmup=warmup, chains=chains, cores=chains, 
               control=control_list, refresh=refresh, backend="cmdstanr",
               file = glue("{out_dir}/{out_base_name}_{y_var}_binomial_model_summary_{Sys.Date()}"), 
               save_model = glue("{out_dir}/{out_base_name}_{y_var}_binomial_model_summary_{Sys.Date()}"))
    print(fit)
    
    sink(glue("{out_dir}/{out_base_name}_{y_var}_binomial_model_summary_{Sys.Date()}.txt"))
    print(fit, digits = 5)
    sink()
    
    return(fit)
}


#' Lighweight wrapper around brms poisson model
brms_poisson_wrapper <- function(y_var, rhs, data, iter, warmup, chains, control_list, priors, out_dir, out_base_name, 
                                 refresh) {
    require(glue)
    require(brms)
    form <- glue("{y_var} ~ {rhs}")
    
    model <- bf(form)+poisson()
    
    fit <- brm(model, data=data, prior=priors, iter=iter, warmup=warmup, chains=chains, cores=chains, 
               control=control_list, refresh=refresh, 
               file = glue("{out_dir}/{out_base_name}_{y_var}_poisson_model_summary_{Sys.Date()}"), 
               save_model = glue("{out_dir}/{out_base_name}_{y_var}_poisson_model_summary_{Sys.Date()}"))
    print(fit)
    
    sink(glue("{out_dir}/{out_base_name}_{y_var}_poisson_model_summary_{Sys.Date()}.txt"))
    print(fit, digits = 5)
    sink()
    
    return(fit)
}


#' Lightweight wrapper around brms ordinal model
brms_ordinal_wrapper <- function(y_var, rhs, data, iter, warmup, chains, control_list, priors, out_dir, out_base_name, 
                                 refresh) {
    require(glue)
    require(brms)
    
    # find the appropriate cumulative link function settings: 
    cum_link_list <- cumulativemodelfit(y_var, data=data)
    
    form <- glue("{y_var} ~ {rhs}")
    
    model <- bf(form)+cumulative(link = cum_link_list$link[[1]], threshold = cum_link_list$threshold[[1]])
    
    # repetitive for now but easiest fix: 
    data[[y_var]] <- ordered(data[[y_var]])
    
    fit <- brm(model, data=data, prior=priors, iter=iter, warmup=warmup, chains=chains, cores=chains, 
               control=control_list, refresh=refresh, 
               file = glue("{out_dir}/{out_base_name}_{y_var}_ordinal_model_summary_{Sys.Date()}"), 
               save_model = glue("{out_dir}/{out_base_name}_{y_var}_ordinal_model_summary_{Sys.Date()}"))
    print(fit)
    
    sink(glue("{out_dir}/{out_base_name}_{y_var}_ordinal_model_summary_{Sys.Date()}.txt"))
    print(fit, digits = 5)
    sink()
    
    return(fit)
}


#' generates a grand-intercept to "ground" the time 1 model and speed up estimate times
estimate_intercept_in_logits <- function(y_var, y_hit, data, cond_var=NULL, cond_value=NULL){
    if(is.null(cond_var) | is.null(cond_value)) {
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


#' Lightweight wrapper around brms categorical model
brms_categorical_wrapper <- function(y_var, rhs, ref_grp, data, iter, warmup, chains, control_list, priors, out_dir, 
                                     out_base_name, refresh=refresh) {
    require(glue)
    require(brms)
    form <- glue("{y_var} ~ {rhs}")
    
    model <- bf(form)+categorical(refcat = ref_grp)
    
    fit <- brm(model, data=data, prior=priors, iter=iter, warmup=warmup, chains=chains, cores=chains, 
               control=control_list, refresh=refresh, backend="cmdstanr",
               file = glue("{out_dir}/{out_base_name}_{y_var}_categorical_model_summary_{Sys.Date()}"), 
               save_model = glue("{out_dir}/{out_base_name}_{y_var}_categorical_model_summary_{Sys.Date()}"))
    print(fit)
    
    sink(glue("{out_dir}/{out_base_name}_{y_var}_categorical_model_summary_{Sys.Date()}.txt"))
    print(fit, digits = 5)
    sink()
    
    return(fit)
} 

#' generate single prior based on posterior for a particular parameter. Currently handles generation of normal and 
#' lognormal prior distributions from a model's posterior. 
make_single_prior <- function(model_obj, prior_dist = "normal", class, coef = NULL, group = NULL, 
                              resp = NULL, dpar = NULL) {
    require(brms)
    require(glue)
    if(!(prior_dist %in% c("normal", "lognormal"))){
        warning(glue("Warning! normal and lognormal only supported distributions. Received: {prior_dist} instead."))
    }
    
    # Should work so long as if dpar is null it is coerced into an empty string... 
    if(is.na(dpar)){
        if(class == "sd" & !is.null(group)) { 
            par_name <- glue("sd_{group}__{coef}") %>% as.character()
        } else if(class == "sd" & is.null(group)) {
            warning("Missing group value for an sd prior")
        } else if(class == "b" | class == "Intercept") {
            par_name <- glue("b_{coef}")
        } else {
            warning("unable to match settings with parameter name definition")
        }
    } else {
        if(class == "sd" & !is.null(group)) { 
            par_name <- glue("sd_{group}__{dpar}_{coef}") %>% as.character()
        } else if(class == "sd" & is.null(group)) {
            warning("Missing group value for an sd prior")
        } else if(class == "b" | class == "Intercept") {
            par_name <- glue("b_{dpar}_{coef}")
        } else {
            warning("unable to match settings with parameter name definition")
        }
    }
    
    tmp_dist <- posterior_samples(model_obj, pars = par_name)
    
    if(prior_dist == "lognormal") {
        tmp_dist <- log(tmp_dist)
    }
    
    if(is.null(group)){
        group <- ""
    }
    
    if(is.na(dpar)){
        dpar <- ""
    }
    
    tmp_mean <- mean(tmp_dist[[par_name]])
    tmp_sd <- sd(tmp_dist[[par_name]])
    
    # Note that the model-based prior is widened by a factor of 3 to ensure sufficient room to move upon new evidence
    prior_dist <- glue("{prior_dist}({tmp_mean}, {tmp_sd*3})") %>% as.character()
    
    # cleaning up intercept - can probably make this more robust...
    coef <- ifelse(class == "Intercept", "", coef)
    prior_settings <- set_prior(prior_dist, class = class, coef = coef, group = group, resp = resp, dpar = dpar)
    return(prior_settings)
}

#' Creates a set of priors from the posterior distributions of a model
create_priors_from_posterior <- function(data_path, model_obj, prior_list, cor = TRUE) {
    require(magrittr)

    load(data_path)
    model_obj <- eval(parse(text=model_obj))
    
    priors <- c(set_prior("normal(0,5)", class = "b"))
    
    if(cor) {
        priors <- c(priors, set_prior("lkj(2)", class = "cor"))
    }
    
    for(p in 1:length(prior_list)){
        term <- prior_list[[p]]

        if(term["class"] == "sd") {
            # Assume a lognormal distribution for standard deviation terms...
            tmp_prior <- make_single_prior(model_obj = model_obj, prior_dist = term["dist"], class = term["class"], 
                                           coef = term["name"], group = term["group"], resp = term["resp"], 
                                           dpar = term["dpar"])
            priors <- c(priors, tmp_prior)
        } else if(term["class"] == "b") { 
            # Assume a normal distribution for the population effects
            tmp_prior <- make_single_prior(model_obj = model_obj, prior_dist = term["dist"], class = term["class"], 
                                           coef = term["name"], resp = term["resp"], dpar = term["dpar"])
            priors <- c(priors, tmp_prior)
        } else if(term["class"] == "Intercept") {
            # Assume a normal distribution for the population Intercepts - in the link scale
            tmp_prior <- make_single_prior(model_obj = model_obj, prior_dist = term["dist"], class = term["class"], 
                                           coef = term["name"], resp = term["resp"], dpar = term["dpar"])
            priors <- c(priors, tmp_prior)
        } else {
            warning("Unable to identify class!")
        }
    }
    
    # cleaning - removing underscore in resp
    # may want to move elsewhere in the future
    priors$resp <- gsub("\\_", "", priors$resp)
    return(priors)
}


#' Lighweight wrapper around brms gaussian model
brms_gaussian_wrapper <- function(y_var, rhs, data, iter, warmup, chains, control_list, priors, out_dir, out_base_name, 
                                 refresh) {
    require(glue)
    require(brms)
    form <- glue("{y_var} ~ {rhs}")
    
    model <- bf(form)+gaussian()
    
    fit <- brm(model, data=data, prior=priors, iter=iter, warmup=warmup, chains=chains, cores=chains, 
               control=control_list, refresh=refresh, backend="cmdstanr", 
               file = glue("{out_dir}/{out_base_name}_{y_var}_gaussian_model_summary_{Sys.Date()}"), 
               save_model = glue("{out_dir}/{out_base_name}_{y_var}_gaussian_model_summary_{Sys.Date()}"))
    print(fit)
    
    sink(glue("{out_dir}/{out_base_name}_{y_var}_gaussian_model_summary_{Sys.Date()}.txt"))
    print(fit, digits = 5)
    sink()
    
    return(fit)
}