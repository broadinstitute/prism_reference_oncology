library(taigr)
library(tidyverse)
library(magrittr)
library(useful)
library(scam)
library(mixtools)
library(scales)
library(ggthemes)
library(dr4pl)
library(drc)

#----
# LOAD THE RAW DATA
#----

analyte_meta = data.table::fread("data/PRISM Oncology Reference 23Q2 - Analyte Meta.csv")
inst_meta <- data.table::fread("data/PRISM Oncology Reference 23Q2 - Inst Meta.csv")
LMFI.long <- data.table::fread("data/PRISM Oncology Reference 23Q2 - LMFI.csv") 


# ----
# COMPUTE THE REFERENCE VALUES FOR CONTROL BARCODES
# ----


REF <- analyte_meta %>% 
  dplyr::filter(pool_id == "CTLBC") %>% 
  dplyr::inner_join(LMFI.long) %>% 
  dplyr::inner_join(inst_meta %>% 
                     dplyr::filter(pert_type == "ctl_vehicle")) %>% 
  dplyr::group_by(prism_replicate, profile_id) %>% 
  dplyr::mutate(m = median(LMFI, na.rm = T)) %>% 
  dplyr::group_by(prism_replicate) %>%  
  dplyr::mutate(m = m - median(m, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(LMFI = LMFI - m) %>% 
  dplyr::group_by(prism_replicate, analyte_id) %>%
  dplyr::summarise(LMFI.ref = median(LMFI, na.rm = T)) %>% 
  dplyr::ungroup()

# ----
# NORMALIZE USING SPLINE FITS
# ----


LMFI.long <- inst_meta %>%
  dplyr::distinct(profile_id, prism_replicate) %>%
  dplyr::right_join(LMFI.long) %>%
  dplyr::left_join(REF) 

profiles = LMFI.long$profile_id %>% unique

LMFI.normalized.long = list()
ix = 1
for (profile in profiles) {
  print(ix)
  temp = LMFI.long %>%
    dplyr::filter(profile_id == profile) 
  
  
  g <- tryCatch(
    scam(
      LMFI.ref ~ s(LMFI, bs = "micx", k = 4),
      data = temp %>% dplyr::filter(is.finite(LMFI.ref)),
      optimizer = "nlm.fd"
    ),
    error = function(e)
      NULL
  )
  
  if (!is.null(g)) {
    p = predict(g, temp, se.fit = T)
    temp$LMFI.normalized = p$fit
    temp$LMFI.normalized.se = p$se.fit
  }
  
  LMFI.normalized.long[[ix]] = temp
  ix = ix + 1
}

LMFI.long <- LMFI.normalized.long %>%  dplyr::bind_rows()
rm(LMFI.normalized.long)


# ----
# QC TABLE
# ----




QC <- inst_meta %>%
  dplyr::filter(pert_type %in% c("ctl_vehicle", "trt_poscon")) %>%
  dplyr::inner_join(LMFI.long) %>%
  dplyr::inner_join(analyte_meta) %>% 
  dplyr::group_by(compound_plate, prism_replicate, culture, analyte_id, ccle_name, screen) %>%
  dplyr::filter(is.finite(LMFI.normalized)) %>% 
  dplyr::summarise(
    error_rate =  min(PRROC::roc.curve(scores.class0 = LMFI.normalized,
                                       weights.class0 = pert_type == "ctl_vehicle",
                                       curve = TRUE)$curve[,1] + 1 -
                        PRROC::roc.curve(scores.class0 = LMFI.normalized,
                                         weights.class0 = pert_type == "ctl_vehicle",
                                         curve = TRUE )$curve[,2])/2,
    NC.median = median(LMFI.normalized[pert_type == "ctl_vehicle"], na.rm = T),
    NC.mad = mad(LMFI.normalized[pert_type == "ctl_vehicle"], na.rm = T),
    PC.median = median(LMFI.normalized[pert_type == "trt_poscon"], na.rm = T),
    PC.mad = mad( LMFI.normalized[pert_type == "trt_poscon"], na.rm = T)) %>%
  dplyr::mutate(DR = NC.median - PC.median,
                SSMD = DR / sqrt(NC.mad ^ 2 + PC.mad ^ 2)) %>%
  dplyr::mutate(PASS = (error_rate <= 0.05) & (DR > 2)) %>%
  dplyr::group_by(compound_plate, analyte_id) %>%
  dplyr::mutate(n.PASS = sum(PASS, na.rm = T)) %>%
  dplyr::ungroup() %>% dplyr::distinct() 



QC %>% 
  write_csv("data/PRISM Oncology Reference 23Q2 - QC Table.csv")




# ----
# COMPUTE LOG-FOLD-CHANGES
# ----

LFC<- LMFI.long %>%  
  dplyr::left_join(QC %>% dplyr::select(analyte_id, prism_replicate, NC.median, PASS, n.PASS)) %>%
  dplyr::mutate(LFC = LMFI.normalized - NC.median) %>%
  dplyr::select(profile_id, prism_replicate, analyte_id, LFC, PASS, n.PASS, screen, compound_plate,culture) 

# -----
# APPLY COMBAT FOR THE POOL EFFECTS 
# -----

apply_combat <- function(Y) {
  # create "condition" column to be used as "batches"
  df <- Y %>%
    dplyr::distinct(analyte_id, profile_id, LFC, culture, pool_id, PASS) %>%
    tidyr::unite(cond, culture, pool_id, profile_id, sep = "::") %>%
    dplyr::filter(is.finite(LFC), PASS)
  
  # calculate means and sd's of each condition
  batch <- df$cond
  m <- rbind(df$LFC,
             rnorm(length(df$LFC),
                   mean =  mean(df$LFC, na.rm = TRUE),
                   sd = sd(df$LFC, na.rm = TRUE)))
  
  # use ComBat to align means and sd's of conditions
  combat <- sva::ComBat(dat = m, batch = batch) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(analyte_id = df$analyte_id, cond = df$cond) %>%
    dplyr::rename(LFC_cb = V1) %>%
    dplyr::mutate(culture = stringr::word(cond, 1, sep = stringr::fixed("::")),
                  pool_id = stringr::word(cond, 2, sep = stringr::fixed("::")),
                  profile_id = stringr::word(cond, 3, sep = stringr::fixed("::"))) %>%
    dplyr::select(-cond, -V2)
  
  combat_corrected <- Y %>%
    dplyr::left_join(combat, by = c("profile_id", "analyte_id", "pool_id", "culture")) %>%
    .$LFC_cb
  
  return(combat_corrected)
}

LFC %<>% 
  dplyr::left_join(analyte_meta %>% 
                     dplyr::distinct(pool_id, analyte_id, culture, compound_plate, screen)) %>% 
  dplyr::select(analyte_id, pool_id, culture, screen, profile_id, LFC, PASS, n.PASS) %>%
  dplyr::left_join(inst_meta %>% 
                     dplyr::distinct(profile_id, compound_plate, screen, pert_type, broad_id, pert_dose)) %>%
  dplyr::filter(is.finite(LFC), pool_id != "CTLBC", !is.na(pool_id), pert_type == "trt_cp") %>% 
  tidyr::unite(col = "condition",
               broad_id, pert_dose, compound_plate, sep = "::", remove = F) %>% 
  split(.$condition) %>% 
  purrr::map_dfr(~dplyr::mutate(.x, LFC_cb = apply_combat(.))) %>%
  dplyr::select(profile_id, analyte_id, LFC, LFC_cb, PASS, n.PASS, screen, compound_plate, culture) 


LFC %>% 
  write_csv("data/PRISM Oncology Reference - LFC.csv")

# -----
# COLLAPSE THE REPLICATES
# -----

LFC.Matrix <- LFC %>% 
  dplyr::left_join(analyte_meta) %>%
  dplyr::left_join(inst_meta) %>%
  dplyr::mutate(col_name = paste0(broad_id, "::", pert_dose, "::", compound_plate, "::", screen)) %>% 
  dplyr::filter(PASS, n.PASS > 1 ,!is.na(depmap_id), is.finite(LFC_cb)) %>%  
  reshape2::acast(depmap_id ~ col_name, value.var = "LFC_cb", fun.aggregate = median)
  
LFC.Matrix %>%
  write.csv("data/PRISM Oncology Reference 23Q2 - LFC Matrix.csv")


# ----
# FIT DOSE RESPONSE CURVES
# ----


compute_auc <- function(l, u, ec50, h, md, MD) {
  f1 = function(x) pmax(pmin((l + (u - l)/(1 + (2^x/ec50)^h)), 1, na.rm = T), 0, na.rm = T)
  return(tryCatch(integrate(f1, log2(md), log2(MD))$value/(log2(MD/md)),
                  error = function(e) {print(e); NA}))
}

compute_log_ic50 <- function(l, u, ec50, h, md, MD) {
  if((l >= 0.5) | (u <= 0.5)) {
    return(NA)
  } else {
    f1 = function(x) (l + (u - l)/(1 + (2^x/ec50)^h) - 0.5)
    return(tryCatch(uniroot(f1, c(log2(md), log2(MD)))$root,
                    error = function(x) NA))
  }
}

compute_MSE_MAD <- function(LFC,  UL, LL,  Slope, Inflection) {
  mse_compute <- LFC %>% 
    dplyr::filter(is.finite(FC)) %>% ## in case there are some nan values accidentally passed.
    dplyr::mutate(FC.pred = UL  + (LL -UL )/(1 + (pert_dose/Inflection)^Slope) ) %>% 
    dplyr::mutate(squared_deviation = (FC-FC.pred)^2, abs_deviation = abs(FC-FC.pred)) %>%
    dplyr::summarise(mse = mean(squared_deviation), mad= median(abs_deviation))
  return (mse_compute)
}

get_best_fit <- function(FC, dose, 
                         UL_low=0.8, UL_up=1.001, slope_decreasing=TRUE){
  ## get best fit among different attempts at fitting, and see if this fit works sufficiently well to be reported.
  
  #browser()
  LFC_filtered = tibble(pert_dose = dose, FC = FC) %>%
    tidyr::drop_na()
  
  var_data <- LFC_filtered$FC %>% var()
  
  all_fits.df <- fit_4param_drc(LFC_filtered, "pert_dose",  var_data, 
                                UL_low, UL_up, slope_decreasing)
  
  if (nrow(all_fits.df)>0){
    res.df <- all_fits.df %>%
      dplyr::top_n(1, frac_var_explained) %>%
      dplyr::mutate(successful_fit = frac_var_explained >= 0) 
    ## fit has to do somewhat better than predicting just the mean of the data to be called successful
  }else{
    res.df  <- data.frame(successful_fit=FALSE)
  }
  
  return (res.df)
}

fit_4param_drc <- function(LFC_filtered, dose_var = "pert_dose",  var_data, 
                           UL_low=0.8, UL_up=1.001, slope_decreasing=TRUE) {
  #fits a number of alternate models  to the DRC and passes the results to the calling function (which chooses the best fit.)
  
  results.df <- list(); ix = 1
  
  slope_bound <- ifelse(slope_decreasing, 0, Inf)  # bound the slopes by default unless passed another option
  
  # warning is in with Calling Handlers and error is in tryCatch
  drc_model <-  tryCatch(withCallingHandlers(
    drc::drm(as.formula(paste("FC ~ ", dose_var)), data=LFC_filtered,
             fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
             lowerl = c(-slope_bound,0.0, UL_low, -Inf),upperl = c(Inf,1.0,UL_up, Inf)),
    warning = function (w){
      invokeRestart("muffleWarning")
    }
  ),
  error = function(e)
  {
    return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))}
  )
  # "slope" in drc package is -ve of slope in dr4pl package
  
  if (drc_model$fit$convergence){
    
    mse_df <- compute_MSE_MAD(LFC_filtered, as.numeric(drc_model$coefficients [[3]]), as.numeric(drc_model$coefficients [[2]]),
                              -as.numeric(drc_model$coefficients [[1]]), as.numeric(drc_model$coefficients [[4]]))
    # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.
    
    results.df[[ix]] <- tibble( fit_name = "drc_drm", 
                                Lower_Limit = as.numeric(drc_model$coefficients [[2]]),
                                Upper_Limit = as.numeric(drc_model$coefficients [[3]]),
                                Slope = -as.numeric(drc_model$coefficients [[1]]),
                                Inflection = as.numeric(drc_model$coefficients [[4]]),
                                MSE = mse_df$mse, MAD = mse_df$mad, frac_var_explained = 1-mse_df$mse/var_data,
                                Input_Parameters = "constrained-drc")
    ix = ix + 1 
  }
  
  
  dr4pl_initL <- tryCatch(dr4pl(as.formula(paste("FC ~ ", dose_var)), data = LFC_filtered,
                                init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.005),
                                method.init = "logistic",
                                lowerl = c(UL_low, -Inf, -Inf, 0),
                                upperl = c(UL_up, Inf, slope_bound, 1.0)),
                          error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )
  
  if (dr4pl_initL$convergence){
    mse_df <- compute_MSE_MAD(LFC_filtered, dr4pl_initL$parameters[[1]], dr4pl_initL$parameters[[4]],
                              dr4pl_initL$parameters[[3]], dr4pl_initL$parameters [[2]])
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_initL", 
                                Lower_Limit = as.numeric(dr4pl_initL$parameters [[4]]),
                                Upper_Limit = as.numeric(dr4pl_initL$parameters [[1]]),
                                Slope = as.numeric(dr4pl_initL$parameters [[3]]),
                                Inflection = as.numeric(dr4pl_initL$parameters [[2]]),
                                MSE = mse_df$mse, MAD = mse_df$mad, frac_var_explained = 1-mse_df$mse/var_data,
                                Input_Parameters = "constrained|init_logistic")
    ix = ix + 1 
    
  }
  
  
  dr4pl_initM_optNM <- tryCatch(dr4pl(as.formula(paste("FC ~ ", dose_var)), data = LFC_filtered,
                                      init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.01),
                                      method.init = "Mead",
                                      lowerl = c(UL_low, -Inf, -Inf, 0),
                                      upperl = c(UL_up, Inf, slope_bound, 1.0),
                                      method.optim="Nelder-Mead"),
                                error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )
  if (dr4pl_initM_optNM$convergence){
    mse_df <- compute_MSE_MAD(LFC_filtered, dr4pl_initM_optNM$parameters[[1]], dr4pl_initM_optNM$parameters[[4]],
                              dr4pl_initM_optNM$parameters[[3]], dr4pl_initM_optNM$parameters [[2]])
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_initM_optNM", 
                                Lower_Limit = as.numeric(dr4pl_initM_optNM$parameters [[4]]),
                                Upper_Limit = as.numeric(dr4pl_initM_optNM$parameters [[1]]),
                                Slope = as.numeric(dr4pl_initM_optNM$parameters [[3]]),
                                Inflection = as.numeric(dr4pl_initM_optNM$parameters [[2]]),
                                MSE = mse_df$mse, MAD = mse_df$mad, frac_var_explained = 1-mse_df$mse/var_data,
                                Input_Parameters = "constrained|init_Mead|optim_Nelder-Mead")
    ix = ix + 1 
  }
  
  
  dr4pl_initL_optB <- tryCatch(dr4pl(as.formula(paste("FC ~ ", dose_var)), data = LFC_filtered,
                                     init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.01),
                                     method.init = "logistic",
                                     lowerl = c(UL_low, -Inf, -Inf, 0),
                                     upperl = c(UL_up, Inf, slope_bound, 1.0),
                                     method.optim="BFGS"),
                               error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )
  if (dr4pl_initL_optB$convergence){
    mse_df <- compute_MSE_MAD(LFC_filtered, dr4pl_initL_optB$parameters[[1]], dr4pl_initL_optB$parameters[[4]],
                              dr4pl_initL_optB$parameters[[3]], dr4pl_initL_optB$parameters [[2]])
    
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_optB", 
                                Lower_Limit = as.numeric(dr4pl_initL_optB$parameters [[4]]),
                                Upper_Limit = as.numeric(dr4pl_initL_optB$parameters [[1]]),
                                Slope = as.numeric(dr4pl_initL_optB$parameters [[3]]),
                                Inflection = as.numeric(dr4pl_initL_optB$parameters [[2]]),
                                MSE = mse_df$mse, MAD = mse_df$mad, frac_var_explained = 1-mse_df$mse/var_data,
                                Input_Parameters = "constrained|init_logistic|optim_BFGS")
    ix = ix + 1 
    
  }
  
  dr4pl_lossHuber<- tryCatch(dr4pl(as.formula(paste("FC ~ ", dose_var)), data = LFC_filtered,
                                   init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.01),
                                   method.robust="Huber",
                                   lowerl = c(UL_low, -Inf, -Inf, 0),
                                   upperl = c(UL_up, Inf, slope_bound, 1.0)),
                             error= function(e){return(list(convergence=FALSE, error=TRUE))}
  )
  if (dr4pl_lossHuber$convergence){
    mse_df <- compute_MSE_MAD(LFC_filtered, dr4pl_lossHuber$parameters[[1]], dr4pl_lossHuber$parameters[[4]],
                              dr4pl_lossHuber$parameters[[3]], dr4pl_lossHuber$parameters [[2]])
    
    
    results.df[[ix]] <- tibble( fit_name = "dr4pl_lossHuber", 
                                Lower_Limit = as.numeric(dr4pl_lossHuber$parameters [[4]]),
                                Upper_Limit = as.numeric(dr4pl_lossHuber$parameters [[1]]),
                                Slope = as.numeric(dr4pl_lossHuber$parameters [[3]]),
                                Inflection = as.numeric(dr4pl_lossHuber$parameters [[2]]),
                                MSE = mse_df$mse, MAD = mse_df$mad, frac_var_explained = 1-mse_df$mse/var_data,
                                Input_Parameters = "constrained|loss_Huber")
    ix = ix + 1 
    
  }
  
  return (dplyr::bind_rows(results.df))
}


DRC <-  LFC %>% 
  dplyr::left_join(inst_meta) %>%
  dplyr::filter(PASS, n.PASS > 1, is.finite(LFC_cb)) %>% 
  dplyr::group_by(broad_id, compound_plate, analyte_id, culture, screen) %>% 
  dplyr::summarise(get_best_fit(2^LFC_cb, pert_dose)) %>%
  dplyr::ungroup() 

DRC <- inst_meta %>% 
  dplyr::group_by(broad_id, compound_plate, screen) %>% 
  dplyr::summarise(md = min(pert_dose, na.rm = T),
                   MD = max(pert_dose, na.rm = T)) %>%
  dplyr::ungroup() %>% 
  dplyr::left_join(DRC) %>% 
  dplyr::filter(successful_fit) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(AUC = compute_auc(Lower_Limit, Upper_Limit, Inflection, -Slope, md, MD),
                log2.IC50 = compute_log_ic50(Lower_Limit, Upper_Limit, Inflection, -Slope, md, MD)) %>%
  dplyr::select(-successful_fit)


DRC %>% 
  write_csv("data/PRISM Oncology Reference 23Q2 - Dose Response Parameters.csv")





# ----
# PORTAL FILES 
# ----
# These files are produced to be consumed by the portal


M <- LFC %>% 
  dplyr::left_join(analyte_meta) %>%
  dplyr::left_join(inst_meta) %>% 
  dplyr::filter(PASS, n.PASS > 1, pert_type == "trt_cp", is.finite(LFC_cb)) %>%
  reshape2::acast(profile_id ~ depmap_id, value.var = "LFC_cb",
                  fun.aggregate = median)


Model.annotations <- tibble(cell_line_name = colnames(M),
                            index = 0:(dim(M)[2]-1)) 


Condition.annotations <- tibble(profile_id = rownames(M),
                                index = 0:(dim(M)[1]-1)) %>%
  dplyr::left_join(inst_meta) %>% 
  dplyr::distinct(index, broad_id, pert_dose, prism_replicate, compound_plate) %>%
  #dplyr::left_join(compound_list, by = c("pert_mfc_id" = "IDs")) %>%
  dplyr::rename(compound_name = broad_id,
                dose = pert_dose) %>%
  dplyr::mutate(masked = "F") %>%
  dplyr::group_by(compound_plate, compound_name, dose) %>% 
  dplyr::arrange(prism_replicate) %>% dplyr::mutate(replicate = 1:n()) %>% 
  dplyr::ungroup() %>% dplyr::select(-compound_plate, -prism_replicate) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(compound_name = paste0("BRD:", compound_name))

colnames(M) <- Model.annotations$index
rownames(M) <- Condition.annotations$index


Response.curves <- DRC %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::rename(cell_line_name = depmap_id,
                compound_name = broad_id,
                ec50 = Inflection,
                lower_asymptote = Lower_Limit,
                upper_asymptote = Upper_Limit,
                slope = Slope) %>% 
  dplyr::group_by(cell_line_name, compound_name) %>%
  dplyr::top_n(1, frac_var_explained) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(cell_line_name, compound_name, ec50, lower_asymptote, upper_asymptote, slope) %>%
  dplyr::mutate(compound_name = paste0("BRD:", compound_name))


AUC = DRC %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::mutate(compound_name = paste0("BRD:", broad_id)) %>% 
  dplyr::filter(is.finite(AUC)) %>% 
  reshape2::acast(compound_name ~ depmap_id, value.var = "AUC", fun.aggregate = median)

IC50 = DRC %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::mutate(compound_name = paste0("BRD:", broad_id)) %>% 
  dplyr::filter(is.finite(log2.IC50)) %>% 
  dplyr::mutate(IC50 = 2^log2.IC50) %>% 
  reshape2::acast(compound_name ~ depmap_id, value.var = "IC50", fun.aggregate = median)


Compound_List %>%
  tidyr::separate_rows(broad_id, sep = ";") %>% 
  dplyr::mutate(broad_id = paste0("BRD:", broad_id)) %>% 
  dplyr::group_by(Drug.Name, MOA, Synonyms, Target, screen) %>% 
  dplyr::summarise(IDs = paste0(broad_id, collapse = ";")) %>% 
  dplyr::rename(repurposing_target = Target) %>%
  write_csv("portal_files/compound_list.csv")

(2^M) %>% 
  write.csv("portal_files/Viability_matrix.csv")

Condition.annotations %>% 
  write_csv("portal_files/Condition_annotations.csv")

Model.annotations %>%  
  write_csv("portal_files/Model_annotations.csv")

Response.curves %>% 
  dplyr::distinct() %>% 
  write_csv("portal_files/Response_curves.csv")

AUC %>% 
  write.csv("portal_files/AUC_matrix.csv")

IC50 %>% 
  write.csv("portal_files/IC50_matrix.csv")






