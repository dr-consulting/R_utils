#' Function performs a lavaan-based moderated-mediation model when there is a
#' sole ordinal/dichotomous y variable for an outcome. Takes the data.frame, 
#' y, x, moderator, and mediator variable names, centers the rhs variables, 
#' creates in an interaction term, and forms the model. Output includes several
#' sets of direct and indirect effects - focusing mainly on the mod-med portion 
#' of the model. 

lavaan_moderated_mediation_single_categorical_y <- function(data, y_var, x_var, mod_var,
                                                            med_var){
  require(glue)
  require(lavaan)
  vars_to_center <- c(x_var, mod_var, med_var)

  tmp_df <- dat_med_mod[vars_to_center]
  tmp_df <- data.frame(lapply(tmp_df, FUN=function(x){scale(x, scale = FALSE)}))
  tmp_df[paste0(x_var, "_x_", mod_var)] <- tmp_df[x_var]*tmp_df[mod_var]

  tmp_df <- cbind(dat_med_mod[y_var], tmp_df)

  model <- "
  {y_var} ~ x_slp_y*{x_var} + m_slp_y*{mod_var} + xm_slp_y*{x_var}_x_{mod_var} + med_slp_y*{med_var}
  {med_var} ~ x_slp_med*{x_var} + m_slp_med*{mod_var} + xm_slp_med*{x_var}_x_{mod_var}

  # Getting mean and variance for low/high moderation analysis
  {mod_var} ~ {mod_var}_mean*1
  {mod_var} ~~ {mod_var}_var*{mod_var}

  # Indirect Effects
  indirect_SD_low := (x_slp_med + xm_slp_med*({mod_var}_mean - sqrt({mod_var}_var))) * med_slp_y
  indirect_SD_high := (x_slp_med + xm_slp_med*({mod_var}_mean + sqrt({mod_var}_var))) * med_slp_y

  # Direct Effects
  direct_SD_low := x_slp_y + xm_slp_y*({mod_var}_mean - sqrt({mod_var}_var))
  direct_SD_high := x_slp_y + xm_slp_y*({mod_var}_mean + sqrt({mod_var}_var))

  # Total Effects
  total_SD_low := direct_SD_low + indirect_SD_low
  total_SD_high := direct_SD_high + indirect_SD_high

  # Proportion of Indirect
  prop_med_SD_low := indirect_SD_low / total_SD_low
  prop_med_SD_high := indirect_SD_high / total_SD_high

  # Indirect effects - statistical significance
  {x_var}_indirect := x_slp_med * med_slp_y
  {mod_var}_indirect := m_slp_med * med_slp_y
  moderator_indirect := xm_slp_med * med_slp_y
  "

  model <- glue(model)

  fit <- sem(model=model, data=tmp_df, ordered = y_var, conditional.x=FALSE, 
	     fixed.x=FALSE)

  summary(fit, rsquare=TRUE)
  return(fit)
}

