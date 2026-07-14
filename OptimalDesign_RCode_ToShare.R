# Clear existing data and graphics
rm(list=ls())
#graphics.off()


# Load libraries ----
pacman::p_load(
  tidyverse, 
  rvinecopulib,  # for r-vine copula implementation
  e1071,  # for r-vine copula implementation
  caret,  # for r-vine copula implementation and to train ML models using CV
  EnvStats,  # for r-vine copula implementation - NOTE: predict OVERWRITES predict from base R stats package!!!
  truncnorm,  # for r-vine copula implementation
  gridExtra,
  mice,  # to perform multiple imputation by chained equations
  purrr,
  foreach,
  doParallel,
  doRNG,
  furrr,
  stats  # for t-test
)

# Functions to load ----

## Expit function ----
expit <- function(x) {exp(x)/(1+exp(x))}

## R-vine copula functions ----

# Function to transform uniform distribution to original scale according to empirical distribution
PseudoObsInverse <- function(DataBaseline, UniformData) {
  # DataBaseline is the matrix of covariates (original data)
  # UniformData is the matrix of uniform data to be transformed to original scale
  PsInverse <- list()
  for (j in 1:ncol(DataBaseline)){ 
    ecdfj <- ecdf(DataBaseline[, j])  # empirical cdf
    ECDFvar <- get("x", environment(ecdfj))
    ECDFjump <- get("y", environment(ecdfj))
    PsInverse[[j]] <- stepfun(ECDFjump[-length(ECDFjump)], ECDFvar)  # define step function
  } 
  ScaledData <- matrix(0, nrow(UniformData), ncol(UniformData))
  for (j in 1:ncol(UniformData)){ ScaledData[, j] <- PsInverse[[j]](UniformData[, j]) }
  ScaledData <- as.data.frame(ScaledData)
  # output
  return(ScaledData)
}

# Function to estimate the 'R-vine copula' model of the baseline data (covariates)
# Note: treating HADS score as continuous
Estimation_Copula <- function(DataBaseline)   {
  # DataBaseline is the matrix of covariates (original data)
  
  # Data  preparation
  # Transformation of continuous variables (in original scale) into uniform distribution variables  
  # Pseudo-observations compute using rvinecopulib package   
  U_cont <- pseudo_obs(DataBaseline[, 6:7])  # columns 6-7 are the continuous variables (age, hads_base)
  
  # Distribution of the discrete variables
  disc_1 <- as.integer(DataBaseline[, 1])  # binary variable should have levels 0, 1
  disc_2 <- as.integer(DataBaseline[, 2])  # categorical variable should have levels 0, 1, 2, etc.
  disc_3 <- as.integer(DataBaseline[, 3])
  disc_4 <- as.integer(DataBaseline[, 4])
  disc_5 <- as.integer(DataBaseline[, 5])
  # disc_7 <- as.integer(DataBaseline[, 7])  # HADS score
  freq_disc1 <- prop.table(table(DataBaseline[, 1]))
  freq_disc2 <- prop.table(table(DataBaseline[, 2]))
  freq_disc3 <- prop.table(table(DataBaseline[, 3]))
  freq_disc4 <- prop.table(table(DataBaseline[, 4]))
  freq_disc5 <- prop.table(table(DataBaseline[, 5]))
  # freq_disc7 <- prop.table(table(DataBaseline[, 7]))
  
  # Preparation of the discrete variables needed to use 'vinecop' function for mixed data (package rvinecopulib)
  Freq_disc_t1 <- cbind(pdiscrete(disc_1 + 1, freq_disc1), pdiscrete(disc_2 + 1, freq_disc2),
                        pdiscrete(disc_3 + 1, freq_disc3), pdiscrete(disc_4 + 1, freq_disc4),
                        pdiscrete(disc_5 + 1, freq_disc5))
  Freq_disc_t0 <- cbind(pdiscrete(disc_1, freq_disc1), pdiscrete(disc_2, freq_disc2),
                        pdiscrete(disc_3, freq_disc3), pdiscrete(disc_4, freq_disc4),
                        pdiscrete(disc_5, freq_disc5))
  U_mixte <- cbind(Freq_disc_t1, U_cont, Freq_disc_t0) # need Freq_disc_t0 to handle discrete obs (check details of rdocumentation)
  #density: ddiscrete(x+1, freq)
  #distribution function: pdiscrete(x+1, freq)
  #quantile function: qdiscrete(u[, 1], freq) - 1
  
  # Estimation of the R-vine model for mixed data using rvinecopulib package
  fit_DataDriven <- vinecop(U_mixte, var_types = c(rep("d", 5), rep("c", 2)))  # 5 discrete vars, 2 continuous vars
  #summary(fit_DataDriven)
  #plot(fit_DataDriven)
  #contour(fit_DataDriven)
  
  # Definition of the R-vine distribution 
  Fit_dist <- vinecop_dist(fit_DataDriven$pair_copulas, fit_DataDriven$structure, fit_DataDriven$var_types)
  
  ## Output
  return(Fit_dist)
} 

# Function to generate a sample according to the estimated R-vine model and baseline data (covariates)
# Note: treating HADS score as continuous
Simulation_Copula <- function(N, Fit_dist, DataBaseline)   {
  # N is number of observations to be generated (sample size)
  # Fit_dist is the R-vine model estimated on original data  
  # DataBaseline is the matrix of covariates (original data)
  
  # Generation of a uniform sample using the estimated R-vine copula distribution
  U_Simu <- rvinecop(N, Fit_dist)
  # Transform uniform distribution to original scale according to empirical distribution  
  # (reverse function for 'pseudo_obs' one)
  # This function is defined above
  VGenCop <- PseudoObsInverse(DataBaseline, U_Simu)
  # Data preparation
  for (i in 1:5){ VGenCop[,i] = as.factor(VGenCop[, i]) }  # discrete vars
  for (i in 6:7){ VGenCop[,i] = as.numeric(as.character(VGenCop[, i])) }  # continuous vars
  colnames(VGenCop) <- colnames(DataBaseline)
  levels(VGenCop[, 1]) <- c('0', '1')  # levels in original data 
  levels(VGenCop[, 2]) <- c('0', '1') 
  levels(VGenCop[, 3]) <- c('0', '1') 
  levels(VGenCop[, 4]) <- c('0', '1') 
  levels(VGenCop[, 5]) <- c('0', '1', '2') 
  
  # Output
  return(VGenCop)
}  

## Execution models ----
# Function to generate DT at baseline
# We will use all other baseline variables
# Note: we perform a data post-processing step of rounding the generated value to the nearest integer value
#       b/c this variable can only take on integer values between 0-10 (inclusive)
Simulation_DT_Base <- function(Model, Covariates) {
  # Model is the prediction model used to predict DT at baseline
  # Covariates are the synthetic covariates used to generate synthetic dt_base
  
  # Prediction
  dt_base_predict <- stats::predict(Model, newdata = Covariates)
  
  # Residuals
  dt_base_resid <- residuals(Model)
  
  # Initialize vector to store synthetic DT baseline values
  dt_base_synthetic = rep(NA, nrow(Covariates))
  
  for (i in 1:length(dt_base_predict)) {
    # Get all possible values of pred + resid for the i'th observation
    pred_resid_sums <- dt_base_predict[i] + dt_base_resid
    
    # All 0 <= pred + resid <= 10
    pred_resid_sums_pos <- pred_resid_sums[pred_resid_sums >= 0 & pred_resid_sums <= 10]
    
    # Randomly sample from sum values in range [0, 10]
    sample_val = sample(pred_resid_sums_pos, size = 1)
    
    # Save sample
    dt_base_synthetic[i] <- sample_val %>% round()
  }
  return(dt_base_synthetic)
}

# Function to generate random treatment assignment at Stage 1 - sequential stratified block randomization
# Note: this function currently returns the entire data set, not just tx
#       (in prior workflow, function to generate synthetic treatment returned only tx, not the entire data set)
Simulation_Treatment_1_HADS3levels <- function(data, seed) {
  
  set.seed(seed)
  
  ## Helper: generate one randomized block
  new_block <- function() {
    blk <- sample(c(2L, 4L), size = 1)
    
    tibble(
      treatment = sample(
        c(
          rep("CT+LG", blk / 2),   # At Stage 1, Intervention = Coping-Together + Lay Guidance
          rep("CT", blk / 2)  # At Stage 1, Control = Self-directed Coping-Together
        )
      ),
      blk_size = blk
    )
  }
  
  nested <- data %>%
    arrange(ID) %>%   # enrollment order
    mutate(
      anxiety_stratum = case_when(
        hads_base < 8 ~ "HADS < 8",
        hads_base >= 8 & hads_base <= 10 ~ "HADS 8-10",
        hads_base >= 11 ~ "HADS 11-21"
      )
    ) %>%
    group_by(anxiety_stratum) %>%
    nest() %>%
    mutate(
      randomization = map(data, function(df) {
        
        current_block <- new_block()
        block_number  <- 1L
        block_pos     <- 1L
        
        map_dfr(seq_len(nrow(df)), function(i) {
          
          ## Start a new block if needed
          if (block_pos > nrow(current_block)) {
            current_block <<- new_block()
            block_number  <<- block_number + 1L
            block_pos     <<- 1L
          }
          
          out <- tibble(
            treatment = current_block$treatment[block_pos],
            blk_size = current_block$blk_size[block_pos],
            block_number = block_number,
            position_in_block = block_pos
          )
          
          block_pos <<- block_pos + 1L
          out
        })
      })
    ) 
  
  final <- map2_dfr(
    nested$data,
    nested$randomization,
    ~ bind_cols(.x, .y)
  ) %>%
    mutate(anxiety_stratum = c(rep(nested$anxiety_stratum[1], nrow(nested[[2]][[1]])), 
                               rep(nested$anxiety_stratum[2], nrow(nested[[2]][[2]])),
                               rep(nested$anxiety_stratum[3], nrow(nested[[2]][[3]])))) %>%
    arrange(ID) %>%
    select(ID, treatment) %>%
    rename(tx_1 = treatment)
  
  return(final)
}

# Function to generate random treatment assignment at Stage 2
# If responded to treatment at Stage 1, then remain on same treatment
# If not, then re-randomize to stepped-up care or same treatment as Stage 1
Simulation_Treatment_2 <- function(tx_1_data_synth, responder_1_data_synth, seed) {
  
  set.seed(seed)
  
  data <- cbind(tx_1 = tx_1_data_synth, responder_1 = responder_1_data_synth)
  TxSimu_2 <- data %>%
    mutate(tx_2 = case_when(responder_1 == "1" & tx_1 == "CT+LG" ~ "CT+LG",
                            responder_1 == "1" & tx_1 == "CT" ~ "CT",
                            responder_1 == "0" & tx_1 == "CT+LG" ~ sample(c("CT+LG", "CT+MI"), n(), replace = TRUE, prob = c(0.5, 0.5)),
                            responder_1 == "0" & tx_1 == "CT" ~ sample(c("CT", "CT+LG"), n(), replace = TRUE, prob = c(0.5, 0.5)))) %>%
    select(tx_2)
  
  
  return(TxSimu_2)
}

## Function to calculate deltas for fixed betas ----
deltadist <- function(realdata_augmentedbaseline, realdata_SMART, random_seed, n_obs,
                      m, delta, S_hat, tx_1_haseffect) {
  # Input: realdata_SMART is the real data dataframe including all variables from the internal SMART
  #        random_seed is a number for the random seed,
  #        n_obs is number of observations to generate,
  #        m is the number of missing data imputations to perform
  #        delta = standardized effect size
  #        S_hat = estimated variance component, sqrt((Var(Y|A_1=1,A_2=strat1)+Var(Y_A_1=0,A_2=strat2))/2)
  #        tx_1_haseffect = indicator of whether we impose no treatment effect at stage 1 (1 = has effect)
  # Output: vector of delta values
  
  Cov_Discrete <- realdata_augmentedbaseline %>% 
    dplyr::select(c(sex, marital, language, birthcountry, educ)) %>%
    lapply(., as.factor) %>%
    as.data.frame()
  
  # Definition of the matrix of continuous covariates (at baseline)
  Cov_Cont <- realdata_augmentedbaseline %>%
    dplyr::select(c(age, hads_base))
  
  # Definition of the matrix of covariates (at baseline)
  Cov <- c(Cov_Discrete, Cov_Cont) %>% as.data.frame()
  
  # Marginal distribution of the covariates are estimated using empirical estimator 
  set.seed(random_seed)
  
  # Estimation of the R-vine model based on original data
  Rvine_dist <- Estimation_Copula(Cov)
  
  # Simulation of virtual patients based on the R-vine model and empirical distribution (of the original data)
  DataSimu <- Simulation_Copula(n_obs, Rvine_dist, Cov)  # n_obs = number of rows in original data
  
  # Next, generate synthetic DT at baseline (baseline var, but not shared by internal and external data sources)
  dt_base_mod <- lm(dt_base ~ age + as.factor(sex) + as.factor(marital) + as.factor(educ) +
                      as.factor(language) + as.factor(birthcountry) + as.numeric(hads_base),
                    data = realdata_SMART)
  
  DTSimu_base <- Simulation_DT_Base(Model = dt_base_mod, Covariates = DataSimu)
  
  TxSimu_1_withID <- Simulation_Treatment_1_HADS3levels(
    data = DataSimu %>% mutate(ID = 1:nrow(DataSimu), hads_base = as.numeric(hads_base)),
    seed = random_seed
  )
  
  TxSimu_1 <- TxSimu_1_withID %>% select(-c(ID))
  
  # Post-randomization variables have missingness, perform MI to get m complete data sets
  Imp_data <- realdata_SMART %>% 
    select(c(age, sex, marital, language, birthcountry, educ, dt_base, hads_base, 
             tx_1, dt_1, hads_1, responder_1, tx_2, hads_2, ltfu)) %>%
    mutate(age = as.numeric(age), sex = as.factor(sex), marital = as.factor(marital),
           language = as.factor(language), birthcountry = as.factor(birthcountry),
           educ = as.factor(educ), dt_base = as.numeric(dt_base), hads_base = as.numeric(hads_base),
           tx_1 = as.factor(tx_1), dt_1 = as.numeric(dt_1), hads_1 = as.numeric(hads_1),
           responder_1 = as.factor(responder_1), tx_2 = as.factor(tx_2), hads_2 = as.numeric(hads_2),
           ltfu = as.factor(ltfu))
  
  # Need to specify deterministic relationship between:
  # responder_1 and dt_base, dt_1
  # tx_2 and responder_1, tx_1
  
  ini <- mice(Imp_data, max = 0, print = FALSE)
  
  # Define imputation methods
  meth <- ini$meth
  meth["responder_1"] <- "~ ifelse(dt_base > 4, ifelse(dt_1 - dt_base < 0, '1', '0'), ifelse(dt_1 - dt_base < 2, '1', '0'))"
  meth["tx_2"] <- "~ ifelse(responder_1 == '1' & tx_1 == 'CT+LG', 'CT+LG',
       ifelse(responder_1 == '1' & tx_1 == 'CT', 'CT',
              ifelse(responder_1 == 0 & tx_1 == 'CT+LG', sample(c('CT+LG', 'CT+MI'), 1),
                     ifelse(responder_1 == 0 & tx_1 == 'CT', sample(c('CT', 'CT+LG'), 1), NA))))"
  
  # Define predictor matrix
  pred <- ini$pred
  # do NOT want to impute dt_base, dt_1 using responder_1 
  pred[c("dt_base", "dt_1"), "responder_1"] <- 0
  # do NOT want to impute tx_1, responder_1 using tx_2
  pred[c("tx_1", "responder_1"), "tx_2"] <- 0
  # ltfu is collinear with post-baseline missingness --> do NOT impute dt_1, hads_1, responder_1, tx_2, hads_2 using ltfu
  pred[c("dt_1", "hads_1", "responder_1", "tx_2", "hads_2"), "ltfu"] <- 0
  
  imp <- mice(data = Imp_data, 
              meth = meth,
              pred = pred,
              m = m, 
              maxit = 10, 
              seed = random_seed, 
              print = FALSE)
  
  # Simulation of post-randomization variables at T1 (dt_1, hads_1, responder_1)
  
  # dt_1:
  # First, fit execution model to each of the m imputed data sets
  dt_1_model_imp <- with(imp,
                         lm(dt_1 ~ as.factor(sex) + as.factor(marital) + as.factor(language) +
                              as.factor(birthcountry) + as.factor(educ) + age + as.numeric(hads_base) +
                              as.numeric(dt_base) + as.factor(tx_1))
  )
  
  # Predict dt_1 using each of the m model fits
  # Object to store predictions 
  # (rows represent synthetic patients/obs, columns represent predictions for the jth imputed data set)
  dt_1_preds <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m fitted models
    
    # Model parameters
    param <- dt_1_model_imp$analyses[[j]]$coefficients %>% as.matrix()
    
    # synthetic X,A with column of 1's for intercept
    data <- model.matrix(~ ., data = cbind(DataSimu, DTSimu_base, TxSimu_1))
    # data <- model.matrix(~ ., data = cbind(DTSimu_base, TxSimu_1))
    
    # Predicted synthetic Z_11, dt_1
    pred <- t(param) %*% t(data) %>% t()
    
    dt_1_preds[, j] <- pred
    
  }
  
  # Generate synthetic dt_1 for each of the m model fits
  # Object to store m sets of generated dt_1
  dt_1_mgens <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m model fits
    
    # Predicted values using j'th model
    dt_1_predict <- dt_1_preds[, j]
    
    # Residuals from j'th model
    dt_1_resid <- dt_1_model_imp$analyses[[j]]$residuals
    
    for (i in 1:nrow(dt_1_mgens)) {
      # Get all possible values of pred + resid for the i'th observation
      pred_resid_sums <- dt_1_predict[i] + dt_1_resid
      
      # All 0 <= pred + resid <= 10
      pred_resid_sums_pos <- pred_resid_sums[pred_resid_sums >= 0 & pred_resid_sums <= 10]
      
      # Randomly sample from sum values in range [0, 10]
      # It is possible that there are no sums in [0, 10] --> if this is the case, set sample_val = 0
      sample_val <- ifelse(length(pred_resid_sums_pos) != 0, sample(pred_resid_sums_pos, size = 1), 0) 
      # sample_val = sample(pred_resid_sums_pos, size = 1)
      
      # Save sample
      dt_1_mgens[i, j] <- sample_val
    }
    
  }
  
  # Generate final synthetic dt_1, Z^s_11 where _11 represents time T1 and variable 1 (dt_1) (all observed)
  DTSimu_1 <- apply(dt_1_mgens, 1, sample, size = 1) %>% round()
  
  
  # hads_1:
  # First, fit execution model to each of the m imputed data sets
  hads_1_model_imp <- with(imp,
                           lm(as.numeric(hads_1) ~ as.factor(sex) + as.factor(marital) + as.factor(language) +
                                as.factor(birthcountry) + as.factor(educ) + age + as.numeric(hads_base) + 
                                as.numeric(dt_base) + as.factor(tx_1))
  )
  
  # Predict hads_1 using each of the m model fits
  # Object to store predictions 
  # (rows represent synthetic patients/obs, columns represent predictions for the jth imputed data set)
  hads_1_preds <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m fitted models
    
    # Model parameters
    param <- hads_1_model_imp$analyses[[j]]$coefficients %>% as.matrix()
    
    # synthetic X,A with column of 1's for intercept
    data <- model.matrix(~ ., data = cbind(DataSimu, DTSimu_base, TxSimu_1))
    
    # Predicted synthetic Z_12, hads_1
    pred <- t(param) %*% t(data) %>% t()
    
    hads_1_preds[, j] <- pred
    
  }
  
  # Generate synthetic hads_1 for each of the m model fits
  # Object to store m sets of generated hads_1
  hads_1_mgens <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m model fits
    
    # Predicted values using j'th model
    hads_1_predict <- hads_1_preds[, j]
    
    # Residuals from j'th model
    hads_1_resid <- hads_1_model_imp$analyses[[j]]$residuals
    
    for (i in 1:nrow(hads_1_mgens)) {
      # Get all possible values of pred + resid for the i'th observation
      pred_resid_sums <- hads_1_predict[i] + hads_1_resid
      
      # All 0 <= pred + resid <= 21
      pred_resid_sums_pos <- pred_resid_sums[pred_resid_sums >= 0 & pred_resid_sums <= 21]
      
      # Randomly sample from sum values in range [0, 21]
      # It is possible that there are no sums in [0, 21] --> if this is the case, set sample_val = 0
      sample_val <- ifelse(length(pred_resid_sums_pos) != 0, sample(pred_resid_sums_pos, size = 1), 0) 
      # sample_val = sample(pred_resid_sums_pos, size = 1)
      
      # Save sample
      hads_1_mgens[i, j] <- sample_val
    }
    
  }
  
  # Generate final synthetic hads_1, Z^s_12 where _12 represents time T1 and variable 2 (hads_1) (all observed)
  HADSSimu_1 <- apply(hads_1_mgens, 1, sample, size = 1)
  
  
  # responder_1 (Z_13)
  # Note: we don't need to fit a model for this, as responder status is determined based on
  #       tailoring variables
  
  # DESIGN A ORIGINAL:
  # If dt_base >= 5, then responder = Yes if dt_1 - dt_base >= 1
  # If dt_base < 5, then responder = Yes if dt_1 < dt_base
  ResponderSimu_1_A <- cbind(DTSimu_base, DTSimu_1) %>%
    as.data.frame() %>%
    mutate(responder_1 = case_when(DTSimu_base >= 5 & DTSimu_1 - DTSimu_base < 0 ~ "1",
                                   DTSimu_base < 5 & DTSimu_1 - DTSimu_base < 2 ~ "1",
                                   TRUE ~ "0")) %>%
    select(responder_1)
  
  ResponderSimu_1_A_withtx <- cbind(TxSimu_1, ResponderSimu_1_A)
  p1_obs_designA <- (ResponderSimu_1_A_withtx %>% 
                       filter(tx_1 == "CT+LG") %>% 
                       select(responder_1) %>% 
                       unlist() %>% as.numeric() %>% 
                       sum)/nrow(ResponderSimu_1_A_withtx %>% filter(tx_1 == "CT+LG"))
  p0_obs_designA <- (ResponderSimu_1_A_withtx %>%
                       filter(tx_1 == "CT") %>% 
                       select(responder_1) %>% 
                       unlist() %>% as.numeric() %>% 
                       sum)/nrow(ResponderSimu_1_A_withtx %>% filter(tx_1 == "CT"))
  
  # DESIGN B:
  # If dt_1 <= median(dt_1) for tx_1 subgroup, then responder = Yes
  dt_1_med_1 <- cbind(DataSimu, DTSimu_base, TxSimu_1, DTSimu_1) %>% 
    as.data.frame() %>%
    filter(tx_1 == "CT+LG") %>%
    dplyr::select(DTSimu_1) %>% 
    unlist() %>%
    median()
  
  dt_1_med_0 <- cbind(DataSimu, DTSimu_base, TxSimu_1, DTSimu_1) %>% 
    as.data.frame() %>%
    filter(tx_1 == "CT") %>% 
    dplyr::select(DTSimu_1) %>% 
    unlist() %>%
    median()
  
  ResponderSimu_1_B <- cbind(TxSimu_1, DTSimu_1) %>%
    as.data.frame() %>%
    mutate(responder_1 = case_when(TxSimu_1 == "CT+LG" & DTSimu_1 <= dt_1_med_1 ~ "1",
                                   TxSimu_1 == "CT" & DTSimu_1 <= dt_1_med_0 ~ "1",
                                   TRUE ~ "0")) %>%
    select(responder_1)
  
  ResponderSimu_1_B_withtx <- cbind(TxSimu_1, ResponderSimu_1_B)
  (p1_obs_designB <- (ResponderSimu_1_B_withtx %>% 
                        filter(tx_1 == "CT+LG") %>% 
                        select(responder_1) %>% 
                        unlist() %>% as.numeric() %>% 
                        sum)/nrow(ResponderSimu_1_B_withtx %>% 
                                    filter(tx_1 == "CT+LG")))
  (p0_obs_designB <- (ResponderSimu_1_B_withtx %>% 
                        filter(tx_1 == "CT") %>% 
                        select(responder_1) %>% 
                        unlist() %>% as.numeric() %>% 
                        sum)/nrow(ResponderSimu_1_B_withtx %>% 
                                    filter(tx_1 == "CT")))
  
  # DESIGN C:
  # If dt_1 < dt_base - 1 AND hads_1 < hads_base - 1, then responder = Yes
  # i.e., both DT and HADS need to drop by 2 or more in order to be considered a responder
  ResponderSimu_1_C <- cbind(DTSimu_base, DTSimu_1, DataSimu %>% dplyr::select(hads_base), HADSSimu_1) %>%
    as.data.frame() %>%
    mutate(responder_1 = case_when(DTSimu_1 < DTSimu_base - 1 & HADSSimu_1 < hads_base - 1 ~ "1",
                                   TRUE ~ "0")) %>%
    select(responder_1)
  
  ResponderSimu_1_C_withtx <- cbind(TxSimu_1, ResponderSimu_1_C)
  (p1_obs_designC <- (ResponderSimu_1_C_withtx %>% 
                        filter(tx_1 == "CT+LG") %>% 
                        select(responder_1) %>% 
                        unlist() %>% as.numeric() %>% 
                        sum)/nrow(ResponderSimu_1_C_withtx %>% 
                                    filter(tx_1 == "CT+LG")))
  (p0_obs_designC <- (ResponderSimu_1_C_withtx %>% 
                        filter(tx_1 == "CT") %>% 
                        select(responder_1) %>% 
                        unlist() %>% as.numeric() %>% 
                        sum)/nrow(ResponderSimu_1_C_withtx %>% 
                                    filter(tx_1 == "CT")))
  
  
  # Simulation of synthetic treatment allocation - Stage 2
  
  # DESIGN A:
  TxSimu_2_A <- Simulation_Treatment_2(tx_1_data_synth = TxSimu_1,
                                       responder_1_data_synth = ResponderSimu_1_A,
                                       seed = random_seed)
  
  # Simulation of post-randomization variable at T2 (hads_2)
  
  # First, fit execution model to each of the m imputed data sets
  hads_2_model_imp <- with(imp,
                           lm(as.numeric(hads_2) ~ as.factor(sex) + as.factor(marital) + as.factor(language) +
                                as.factor(birthcountry) + as.factor(educ) + age + as.numeric(hads_base) + 
                                as.numeric(dt_base) + as.factor(tx_1) + as.numeric(dt_1) + as.numeric(hads_1) +
                                as.factor(tx_2))
  )
  
  # Calculate fixed beta4 values using each of the m model fits
  # Object to store beta values 
  # (rows represent values for given imputed data set, columns represent beta2, beta3, beta4)
  betas_fixed <- matrix(data = NA, nrow = m, ncol = 3) %>% as.data.frame()
  colnames(betas_fixed) <- c("beta2", "beta3", "beta4")
  
  for(j in 1:m) {  # iterate through m fitted models
    
    # Model parameters
    param <- hads_2_model_imp$analyses[[j]]$coefficients %>% as.matrix()
    
    # Fix parameter for tx_2_CT+MI = 1/(1 - p1)*(delta*S_hat - beta2 - beta3(p1-p0) - beta4(p1+p0-1)) 
    # tx_2_CT+MI is stepped-up care not shared by both strategies
    # delta = effect size, S_hat = estimated variance component
    # p1 = prop. responders at stage 1 to tx=CT+LG, p0 = prop. responders at stage 1 to tx=CT
    param_simdelta <- param 
    
    # If we want to impose that tx_1 has no effect, then set beta2 = 0
    if(tx_1_haseffect == 1){
      param_simdelta[rownames(param_simdelta) %in% c("as.factor(tx_1)CT+LG")] <- 0
    }
    
    beta2_hat <- param_simdelta[rownames(param_simdelta) %in% c("as.factor(tx_1)CT+LG")]
    beta3_hat <- param_simdelta[rownames(param_simdelta) %in% c("as.factor(tx_2)CT+LG")]
    
    # Fix beta4 based on estimated coef of beta2, beta3, delta, and S, then store beta4
    betas_fixed[j, 3] <- (1/(1 - p1_obs_designA))*((delta*S_hat) - beta2_hat - (beta3_hat*(p1_obs_designA + p0_obs_designA - 1)))
    
    # Store beta2, beta3
    betas_fixed[j, 1] <- beta2_hat
    betas_fixed[j, 2] <- beta3_hat
  }
  
  # Now calculate delta using p1_hat_alt, p0_hat_alt
  deltas_B <- (betas_fixed[, 1] + betas_fixed[, 2]*(p1_obs_designB + p0_obs_designB - 1) + betas_fixed[, 3]*(1 - p1_obs_designB))/S_hat
  deltas_C <- (betas_fixed[, 1] + betas_fixed[, 2]*(p1_obs_designC + p0_obs_designC - 1) + betas_fixed[, 3]*(1 - p1_obs_designC))/S_hat
  
  return(list(deltas_B, deltas_C))
  
}


## Function to generate 1 data set under fixed mechanism ----
generate1dataset_alldesigns_fixedmech <- function(realdata_augmentedbaseline, realdata_SMART,
                                                  random_seed, n_obs, delta, S_hat, tx_1_haseffect) {
  # Input: realdata_augmentedbaseline is the real data dataframe that includes only baseline data (augmented),
  #        realdata_SMART is the real data dataframe including all variables from the internal SMART
  #        random_seed is a number for the random seed,
  #        n_obs is the number of observations (i.e., rows) to generate
  #        delta = standardized effect size
  #        S_hat = estimated variance component, sqrt((Var(Y|A_1=1,A_2=strat1)+Var(Y_A_1=0,A_2=strat2))/2)
  #        tx_1_haseffect = indicator of whether we impose no treatment effect at stage 1 (1 = imposing tx_1 has no effect)
  # Output: dataframe of synthetic data under Design A, B, and C
  
  # Definition of the matrix of discrete covariates (at baseline) - SHARED by internal and external data
  # i.e., don't include dt_base (we will use execution model to generate this variable, not R-vine copula)
  Cov_Discrete <- realdata_augmentedbaseline %>% 
    dplyr::select(c(sex, marital, language, birthcountry, educ)) %>%
    lapply(., as.factor) %>%
    as.data.frame()
  
  # Definition of the matrix of continuous covariates (at baseline)
  Cov_Cont <- realdata_augmentedbaseline %>%
    dplyr::select(c(age, hads_base))
  
  # Definition of the matrix of covariates (at baseline)
  Cov <- c(Cov_Discrete, Cov_Cont) %>% as.data.frame()
  
  # Marginal distribution of the covariates are estimated using empirical estimator 
  set.seed(random_seed)
  
  # Estimation of the R-vine model based on original data
  Rvine_dist <- Estimation_Copula(Cov)
  
  # Simulation of virtual patients based on the R-vine model and empirical distribution (of the original data)
  DataSimu <- Simulation_Copula(n_obs, Rvine_dist, Cov)  # n_obs = number of rows in original data
  
  # Next, generate synthetic DT at baseline (baseline var, but not shared by internal and external data sources)
  dt_base_mod <- lm(dt_base ~ age + as.factor(sex) + as.factor(marital) + as.factor(educ) + 
                      as.factor(language) + as.factor(birthcountry) + as.numeric(hads_base),
                    data = realdata_SMART)
  
  DTSimu_base <- Simulation_DT_Base(Model = dt_base_mod, Covariates = DataSimu)
  
  TxSimu_1_withID <- Simulation_Treatment_1_HADS3levels(
    data = DataSimu %>% mutate(ID = 1:nrow(DataSimu), hads_base = as.numeric(hads_base)),
    seed = random_seed
  )
  
  TxSimu_1 <- TxSimu_1_withID %>% select(-c(ID))
  
  # Post-randomization variables have missingness, perform MI to get m complete data sets
  Imp_data <- realdata_SMART %>% 
    select(c(age, sex, marital, language, birthcountry, educ, dt_base, hads_base, 
             tx_1, dt_1, hads_1, responder_1, tx_2, hads_2, ltfu)) %>%
    mutate(age = as.numeric(age), sex = as.factor(sex), marital = as.factor(marital),
           language = as.factor(language), birthcountry = as.factor(birthcountry),
           educ = as.factor(educ), dt_base = as.numeric(dt_base), hads_base = as.numeric(hads_base),
           tx_1 = as.factor(tx_1), dt_1 = as.numeric(dt_1), hads_1 = as.numeric(hads_1),
           responder_1 = as.factor(responder_1), tx_2 = as.factor(tx_2), hads_2 = as.numeric(hads_2),
           ltfu = as.factor(ltfu))
  
  # Calculate m, where m ~ max(prop. of missingness)*100
  m <- round(colMeans(is.na(Imp_data)) %>% max()*100)  # m ~ prop. of missingness*100
  
  # Need to specify deterministic relationship between:
  # responder_1 and dt_base, dt_1
  # tx_2 and responder_1, tx_1
  
  ini <- mice(Imp_data, max = 0, print = FALSE)
  
  # Define imputation methods
  meth <- ini$meth
  meth["responder_1"] <- "~ ifelse(dt_base > 4, ifelse(dt_1 - dt_base < 0, '1', '0'), ifelse(dt_1 - dt_base < 2, '1', '0'))"
  meth["tx_2"] <- "~ ifelse(responder_1 == '1' & tx_1 == 'CT+LG', 'CT+LG',
       ifelse(responder_1 == '1' & tx_1 == 'CT', 'CT',
              ifelse(responder_1 == 0 & tx_1 == 'CT+LG', sample(c('CT+LG', 'CT+MI'), 1),
                     ifelse(responder_1 == 0 & tx_1 == 'CT', sample(c('CT', 'CT+LG'), 1), NA))))"
  
  # Define predictor matrix
  pred <- ini$pred
  # do NOT want to impute dt_base, dt_1 using responder_1 
  pred[c("dt_base", "dt_1"), "responder_1"] <- 0
  # do NOT want to impute tx_1, responder_1 using tx_2
  pred[c("tx_1", "responder_1"), "tx_2"] <- 0
  # ltfu is collinear with post-baseline missingness --> do NOT impute dt_1, hads_1, responder_1, tx_2, hads_2 using ltfu
  pred[c("dt_1", "hads_1", "responder_1", "tx_2", "hads_2"), "ltfu"] <- 0
  
  imp <- mice(data = Imp_data, 
              meth = meth,
              pred = pred,
              m = m, 
              maxit = 10, 
              seed = random_seed, 
              print = FALSE)
  
  # Simulation of post-randomization variables at T1 (dt_1, hads_1, responder_1)
  
  # dt_1:
  # First, fit execution model to each of the m imputed data sets
  dt_1_model_imp <- with(imp,
                         lm(dt_1 ~ as.factor(sex) + as.factor(marital) + as.factor(language) +
                              as.factor(birthcountry) + as.factor(educ) + age + as.numeric(hads_base) + 
                              as.numeric(dt_base) + as.factor(tx_1))
  )
  
  # Predict dt_1 using each of the m model fits
  # Object to store predictions 
  # (rows represent synthetic patients/obs, columns represent predictions for the jth imputed data set)
  dt_1_preds <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m fitted models
    
    # Model parameters
    param <- dt_1_model_imp$analyses[[j]]$coefficients %>% as.matrix()
    
    # synthetic X,A with column of 1's for intercept
    data <- model.matrix(~ ., data = cbind(DataSimu, DTSimu_base, TxSimu_1))
    
    # Predicted synthetic Z_11, dt_1
    pred <- t(param) %*% t(data) %>% t()
    
    dt_1_preds[, j] <- pred
    
  }
  
  # Generate synthetic dt_1 for each of the m model fits
  # Object to store m sets of generated dt_1
  dt_1_mgens <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m model fits
    
    # Predicted values using j'th model
    dt_1_predict <- dt_1_preds[, j]
    
    # Residuals from j'th model
    dt_1_resid <- dt_1_model_imp$analyses[[j]]$residuals
    
    for (i in 1:nrow(dt_1_mgens)) {
      # Get all possible values of pred + resid for the i'th observation
      pred_resid_sums <- dt_1_predict[i] + dt_1_resid
      
      # All 0 <= pred + resid <= 10
      pred_resid_sums_pos <- pred_resid_sums[pred_resid_sums >= 0 & pred_resid_sums <= 10]
      
      # Randomly sample from sum values in range [0, 10]
      # It is possible that there are no sums in [0, 10] --> if this is the case, set sample_val = 0
      sample_val <- ifelse(length(pred_resid_sums_pos) != 0, sample(pred_resid_sums_pos, size = 1), 0) 
      # sample_val = sample(pred_resid_sums_pos, size = 1)
      
      # Save sample
      dt_1_mgens[i, j] <- sample_val
    }
    
  }
  
  # Generate final synthetic dt_1, Z^s_11 where _11 represents time T1 and variable 1 (dt_1) (all observed)
  DTSimu_1 <- apply(dt_1_mgens, 1, sample, size = 1) %>% round()
  
  
  # hads_1:
  # First, fit execution model to each of the m imputed data sets
  hads_1_model_imp <- with(imp,
                           lm(as.numeric(hads_1) ~ as.factor(sex) + as.factor(marital) + as.factor(language) +
                                as.factor(birthcountry) + as.factor(educ) + age + as.numeric(hads_base) + 
                                as.numeric(dt_base) + as.factor(tx_1))
  )
  
  # Predict hads_1 using each of the m model fits
  # Object to store predictions 
  # (rows represent synthetic patients/obs, columns represent predictions for the jth imputed data set)
  hads_1_preds <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m fitted models
    
    # Model parameters
    param <- hads_1_model_imp$analyses[[j]]$coefficients %>% as.matrix()
    
    # synthetic X,A with column of 1's for intercept
    data <- model.matrix(~ ., data = cbind(DataSimu, DTSimu_base, TxSimu_1))
    
    # Predicted synthetic Z_12, hads_1
    pred <- t(param) %*% t(data) %>% t()
    
    hads_1_preds[, j] <- pred
    
  }
  
  # Generate synthetic hads_1 for each of the m model fits
  # Object to store m sets of generated hads_1
  hads_1_mgens <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m model fits
    
    # Predicted values using j'th model
    hads_1_predict <- hads_1_preds[, j]
    
    # Residuals from j'th model
    hads_1_resid <- hads_1_model_imp$analyses[[j]]$residuals
    
    for (i in 1:nrow(hads_1_mgens)) {
      # Get all possible values of pred + resid for the i'th observation
      pred_resid_sums <- hads_1_predict[i] + hads_1_resid
      
      # All 0 <= pred + resid <= 21
      pred_resid_sums_pos <- pred_resid_sums[pred_resid_sums >= 0 & pred_resid_sums <= 21]
      
      # Randomly sample from sum values in range [0, 21]
      # It is possible that there are no sums in [0, 21] --> if this is the case, set sample_val = 0
      sample_val <- ifelse(length(pred_resid_sums_pos) != 0, sample(pred_resid_sums_pos, size = 1), 0) 
      # sample_val = sample(pred_resid_sums_pos, size = 1)
      
      # Save sample
      hads_1_mgens[i, j] <- sample_val
    }
    
  }
  
  # Generate final synthetic hads_1, Z^s_12 where _12 represents time T1 and variable 2 (hads_1) (all observed)
  HADSSimu_1 <- apply(hads_1_mgens, 1, sample, size = 1)
  
  
  # responder_1 (Z_13)
  # Note: we don't need to fit a model for this, as responder status is determined based on
  #       the difference in dt_1 and dt_base
  
  # DESIGN A:
  # If dt_base >= 5, then responder = Yes if dt_1 - dt_base >= 1
  # If dt_base < 5, then responder = Yes if dt_1 < dt_base
  ResponderSimu_1_A <- cbind(DTSimu_base, DTSimu_1) %>%
    as.data.frame() %>%
    mutate(responder_1 = case_when(DTSimu_base >= 5 & DTSimu_1 - DTSimu_base < 0 ~ "1",
                                   DTSimu_base < 5 & DTSimu_1 - DTSimu_base < 2 ~ "1",
                                   TRUE ~ "0")) %>%
    select(responder_1)
  
  ResponderSimu_1_A_withtx <- cbind(TxSimu_1, ResponderSimu_1_A)
  p1_obs_designA <- (ResponderSimu_1_A_withtx %>% 
                       filter(tx_1 == "CT+LG") %>% 
                       select(responder_1) %>% 
                       unlist() %>% as.numeric() %>% 
                       sum)/nrow(ResponderSimu_1_A_withtx %>% filter(tx_1 == "CT+LG"))
  p0_obs_designA <- (ResponderSimu_1_A_withtx %>% 
                       filter(tx_1 == "CT") %>% 
                       select(responder_1) %>% 
                       unlist() %>% as.numeric() %>% 
                       sum)/nrow(ResponderSimu_1_A_withtx %>% 
                                   filter(tx_1 == "CT"))
  
  # DESIGN B:
  # If dt_1 <= median(dt_1) for tx_1 subgroup, then responder = Yes
  dt_1_med_1 <- cbind(DataSimu, DTSimu_base, TxSimu_1, DTSimu_1) %>% 
    as.data.frame() %>%
    filter(tx_1 == "CT+LG") %>%
    dplyr::select(DTSimu_1) %>% 
    unlist() %>%
    median()
  
  dt_1_med_0 <- cbind(DataSimu, DTSimu_base, TxSimu_1, DTSimu_1) %>% 
    as.data.frame() %>%
    filter(tx_1 == "CT") %>% 
    dplyr::select(DTSimu_1) %>% 
    unlist() %>%
    median()
  
  ResponderSimu_1_B <- cbind(TxSimu_1, DTSimu_1) %>%
    as.data.frame() %>%
    mutate(responder_1 = case_when(TxSimu_1 == "CT+LG" & DTSimu_1 <= dt_1_med_1 ~ "1",
                                   TxSimu_1 == "CT" & DTSimu_1 <= dt_1_med_0 ~ "1",
                                   TRUE ~ "0")) %>%
    select(responder_1)
  
  # DESIGN C:
  # If dt_1 < dt_base - 1 AND hads_1 < hads_base - 1, then responder = Yes
  # i.e., both DT and HADS need to drop by 2 or more in order to be considered a responder
  ResponderSimu_1_C <- cbind(DTSimu_base, DTSimu_1, DataSimu %>% dplyr::select(hads_base), HADSSimu_1) %>%
    as.data.frame() %>%
    mutate(responder_1 = case_when(DTSimu_1 < DTSimu_base - 1 & HADSSimu_1 < hads_base - 1 ~ "1",
                                   TRUE ~ "0")) %>%
    select(responder_1)
  
  
  # Simulation of synthetic treatment allocation - Stage 2
  
  # DESIGN A:
  TxSimu_2_A <- Simulation_Treatment_2(tx_1_data_synth = TxSimu_1,
                                       responder_1_data_synth = ResponderSimu_1_A,
                                       seed = random_seed)
  
  # DESIGN B:
  TxSimu_2_B <- Simulation_Treatment_2(tx_1_data_synth = TxSimu_1,
                                       responder_1_data_synth = ResponderSimu_1_B,
                                       seed = random_seed)
  
  # DESIGN C:
  TxSimu_2_C <- Simulation_Treatment_2(tx_1_data_synth = TxSimu_1,
                                       responder_1_data_synth = ResponderSimu_1_C,
                                       seed = random_seed)
  
  # Simulation of post-randomization variable at T2 (hads_2)
  
  # hads_2:
  # First, fit execution model to each of the m imputed data sets
  hads_2_model_imp <- with(imp,
                           lm(as.numeric(hads_2) ~ as.factor(sex) + as.factor(marital) + as.factor(language) +
                                as.factor(birthcountry) + as.factor(educ) + age + as.numeric(hads_base) + 
                                as.numeric(dt_base) + as.factor(tx_1) + as.numeric(dt_1) + as.numeric(hads_1) +
                                as.factor(tx_2))
  )
  
  # Predict hads_2 using each of the m model fits
  # Object to store predictions 
  # (rows represent synthetic patients/obs, columns represent predictions for the jth imputed data set)
  # Design A
  hads_2_preds_A <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  # Design B
  hads_2_preds_B <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  # Design C
  hads_2_preds_C <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m fitted models
    
    # Model parameters
    param <- hads_2_model_imp$analyses[[j]]$coefficients %>% as.matrix()
    
    # Fix parameter for tx_2_CT+MI = 1/(1 - p1)*(delta*S_hat - beta2 - beta3(p1+p0-1)) 
    # tx_2_CT+MI is stepped-up care not shared by both strategies
    # delta = effect size, S_hat = estimated variance component
    # p1 = prop. responders at stage 1 to tx=CT+LG, p0 = prop. responders at stage 1 to tx=CT
    param_simdelta <- param 
    
    # If we want to impose that tx_1 has no effect, then set beta2 = 0
    if(tx_1_haseffect == 1){
      param_simdelta[rownames(param_simdelta) %in% c("as.factor(tx_1)CT+LG")] <- 0
    }
    
    beta2_hat <- param_simdelta[rownames(param_simdelta) %in% c("as.factor(tx_1)CT+LG")]
    beta3_hat <- param_simdelta[rownames(param_simdelta) %in% c("as.factor(tx_2)CT+LG")]
    
    # Fix beta4 based on estimated coef of beta2, beta3
    param_simdelta[rownames(param_simdelta) %in% c("as.factor(tx_2)CT+MI")] <- (1/(1 - p1_obs_designA))*((delta*S_hat) - beta2_hat - (beta3_hat*(p1_obs_designA + p0_obs_designA - 1)))
    
    # synthetic X,A with column of 1's for intercept
    
    # Design A
    data_A <- model.matrix(~ ., data = cbind(DataSimu, DTSimu_base, TxSimu_1, DTSimu_1, 
                                             HADSSimu_1, TxSimu_2_A))
    
    # Design B
    data_B <- model.matrix(~ ., data = cbind(DataSimu, DTSimu_base, TxSimu_1, DTSimu_1, 
                                             HADSSimu_1, TxSimu_2_B))
    
    # Design C
    data_C <- model.matrix(~ ., data = cbind(DataSimu, DTSimu_base, TxSimu_1, DTSimu_1, 
                                             HADSSimu_1, TxSimu_2_C))
    
    # Predicted synthetic dt_2
    
    # Design A
    pred_A <- t(param_simdelta) %*% t(data_A) %>% t()
    
    # Design B
    pred_B <- t(param_simdelta) %*% t(data_B) %>% t()
    
    # Design C
    pred_C <- t(param_simdelta) %*% t(data_C) %>% t()
    
    # Store prediction for each design
    hads_2_preds_A[, j] <- pred_A
    hads_2_preds_B[, j] <- pred_B
    hads_2_preds_C[, j] <- pred_C
    
  }
  
  # Generate synthetic hads_2 for each of the m model fits
  # Object to store m sets of generated dt_1
  # Design A
  hads_2_mgens_A <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m model fits
    
    # Predicted values using j'th model
    hads_2_predict_A <- hads_2_preds_A[, j]
    
    # Residuals from j'th model
    hads_2_resid <- hads_2_model_imp$analyses[[j]]$residuals
    
    for (i in 1:nrow(hads_2_mgens_A)) {
      # Get all possible values of pred + resid for the i'th observation
      pred_resid_sums_A <- hads_2_predict_A[i] + hads_2_resid
      
      # All 0 <= pred + resid <= 21
      pred_resid_sums_pos_A <- pred_resid_sums_A[pred_resid_sums_A >= 0 & pred_resid_sums_A <= 21]
      
      # Randomly sample from sum values in range [0, 21]
      # It is possible that there are no sums in [0, 21] --> if this is the case, set sample_val = 0
      sample_val_A <- ifelse(length(pred_resid_sums_pos_A) != 0, sample(pred_resid_sums_pos_A, size = 1), 0) 
      
      # Save sample
      hads_2_mgens_A[i, j] <- sample_val_A
    }
    
  }
  
  # Generate final synthetic hads_2, Z^s_21 where _21 represents time T2 and variable 1 (hads_2) (all observed)
  HADSSimu_2_A <- apply(hads_2_mgens_A, 1, sample, size = 1)
  
  # Design B
  hads_2_mgens_B <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m model fits
    
    # Predicted values using j'th model
    hads_2_predict_B <- hads_2_preds_B[, j]
    
    # Residuals from j'th model
    hads_2_resid <- hads_2_model_imp$analyses[[j]]$residuals
    
    for (i in 1:nrow(hads_2_mgens_B)) {
      # Get all possible values of pred + resid for the i'th observation
      pred_resid_sums_B <- hads_2_predict_B[i] + hads_2_resid
      
      # All 0 <= pred + resid <= 21
      pred_resid_sums_pos_B <- pred_resid_sums_B[pred_resid_sums_B >= 0 & pred_resid_sums_B <= 21]
      
      # Randomly sample from sum values in range [0, 21]
      # It is possible that there are no sums in [0, 21] --> if this is the case, set sample_val = 0
      sample_val_B <- ifelse(length(pred_resid_sums_pos_B) != 0, sample(pred_resid_sums_pos_B, size = 1), 0) 
      
      # Save sample
      hads_2_mgens_B[i, j] <- sample_val_B
    }
    
  }
  
  # Generate final synthetic hads_2, Z^s_21 where _21 represents time T2 and variable 1 (hads_2) (all observed)
  HADSSimu_2_B <- apply(hads_2_mgens_B, 1, sample, size = 1)
  
  # Design C
  hads_2_mgens_C <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m model fits
    
    # Predicted values using j'th model
    hads_2_predict_C <- hads_2_preds_C[, j]
    
    # Residuals from j'th model
    hads_2_resid <- hads_2_model_imp$analyses[[j]]$residuals
    
    for (i in 1:nrow(hads_2_mgens_C)) {
      # Get all possible values of pred + resid for the i'th observation
      pred_resid_sums_C <- hads_2_predict_C[i] + hads_2_resid
      
      # All 0 <= pred + resid <= 21
      pred_resid_sums_pos_C <- pred_resid_sums_C[pred_resid_sums_C >= 0 & pred_resid_sums_C <= 21]
      
      # Randomly sample from sum values in range [0, 21]
      # It is possible that there are no sums in [0, 21] --> if this is the case, set sample_val = 0
      sample_val_C <- ifelse(length(pred_resid_sums_pos_C) != 0, sample(pred_resid_sums_pos_C, size = 1), 0) 
      
      # Save sample
      hads_2_mgens_C[i, j] <- sample_val_C
    }
    
  }
  
  # Generate final synthetic hads_2, Z^s_21 where _21 represents time T2 and variable 1 (hads_2) (all observed)
  HADSSimu_2_C <- apply(hads_2_mgens_C, 1, sample, size = 1)
  
  
  # Finally, generate synthetic missingness (ltfu)
  
  # Too few observations to fit execution model
  # Instead, randomly select participants to drop out with observed proportion of drop-out at each stage 
  
  prop_dropout_1 <- sum(Imp_data$ltfu == "T1")/nrow(Imp_data)
  prop_dropout_2 <- sum(Imp_data$ltfu == "T2")/nrow(Imp_data[Imp_data$ltfu != "T1", ])
  
  # Design A
  data_synth_A <- cbind(DataSimu, DTSimu_base, TxSimu_1_withID, DTSimu_1, HADSSimu_1, 
                        ResponderSimu_1_A, TxSimu_2_A, HADSSimu_2 = HADSSimu_2_A) %>%
    select(c(ID, colnames(DataSimu), DTSimu_base, tx_1, DTSimu_1, HADSSimu_1, responder_1,
             tx_2, HADSSimu_2))
  
  # Randomly select participants to drop-out (T1)
  ltfu_1_A <- rbinom(n = nrow(data_synth_A), size = 1, prob = prop_dropout_1)
  
  # Add ltfu at T1 to synthetic data
  data_synth_withltfu1_A <- cbind(data_synth_A, ltfu_1_A)
  
  # Participants who did not drop out at T1
  data_synth_notltfu1_A <- data_synth_withltfu1_A %>% filter(ltfu_1_A == 0)
  
  # Randomly select participants to drop-out (T2)
  ltfu_2_A <- rbinom(n = nrow(data_synth_notltfu1_A), size = 1, prob = prop_dropout_2)
  
  # Add ltfu at T2 to synthetic data
  data_synth_notltfu1_withltfu2_A <- cbind(data_synth_notltfu1_A, ltfu_2_A)
  
  # ltfu_1 and ltfu_2 in the same df
  data_synth_withltfu_A <- data_synth_withltfu1_A %>%
    left_join(data_synth_notltfu1_withltfu2_A %>% select(c(ID, ltfu_2_A)), by = "ID") %>%
    mutate(ltfu_2_A = ifelse(is.na(ltfu_2_A), 0, ltfu_2_A)) %>%
    mutate(ltfu_A = ifelse(ltfu_1_A == 1, "T1", ifelse(ltfu_2_A == 1, "T2", "N")))
  
  data_synth_A <- data_synth_withltfu_A %>%
    rename(dt_base = DTSimu_base,
           dt_1_complete = DTSimu_1,
           hads_1_complete = HADSSimu_1,
           responder_1_complete = responder_1,
           tx_2_complete = tx_2,
           hads_2_complete = HADSSimu_2) %>%
    mutate(dt_1 = ifelse(ltfu_A == "T1", NA, dt_1_complete),
           hads_1 = ifelse(ltfu_A == "T1", NA, hads_1_complete),
           responder_1 = ifelse(ltfu_A == "T1", NA, responder_1_complete),
           tx_2 = ifelse(ltfu_A == "T1", NA, tx_2_complete),
           hads_2 = ifelse(ltfu_A == "T2" | ltfu_A == "T1", NA, hads_2_complete)) %>%
    rename(ltfu = ltfu_A) %>%
    select(c(ID, sex, marital, language, birthcountry, educ, age, hads_base, dt_base, 
             tx_1, dt_1_complete, dt_1, hads_1_complete, hads_1, responder_1_complete,
             responder_1, tx_2_complete, tx_2, hads_2_complete, hads_2, ltfu))
  
  # Design B
  data_synth_B <- cbind(DataSimu, DTSimu_base, TxSimu_1_withID, DTSimu_1, HADSSimu_1, 
                        ResponderSimu_1_B, TxSimu_2_B, HADSSimu_2 = HADSSimu_2_B) %>%
    select(c(ID, colnames(DataSimu), DTSimu_base, tx_1, DTSimu_1, HADSSimu_1, responder_1,
             tx_2, HADSSimu_2))
  
  # Randomly select participants to drop-out (T1)
  ltfu_1_B <- rbinom(n = nrow(data_synth_B), size = 1, prob = prop_dropout_1)
  
  # Add ltfu at T1 to synthetic data
  data_synth_withltfu1_B <- cbind(data_synth_B, ltfu_1_B)
  
  # Participants who did not drop out at T1
  data_synth_notltfu1_B <- data_synth_withltfu1_B %>% filter(ltfu_1_B == 0)
  
  # Randomly select participants to drop-out (T2)
  ltfu_2_B <- rbinom(n = nrow(data_synth_notltfu1_B), size = 1, prob = prop_dropout_2)
  
  # Add ltfu at T2 to synthetic data
  data_synth_notltfu1_withltfu2_B <- cbind(data_synth_notltfu1_B, ltfu_2_B)
  
  # ltfu_1 and ltfu_2 in the same df
  data_synth_withltfu_B <- data_synth_withltfu1_B %>%
    left_join(data_synth_notltfu1_withltfu2_B %>% select(c(ID, ltfu_2_B)), by = "ID") %>%
    mutate(ltfu_2_B = ifelse(is.na(ltfu_2_B), 0, ltfu_2_B)) %>%
    mutate(ltfu_B = ifelse(ltfu_1_B == 1, "T1", ifelse(ltfu_2_B == 1, "T2", "N")))
  
  data_synth_B <- data_synth_withltfu_B %>%
    rename(dt_base = DTSimu_base,
           dt_1_complete = DTSimu_1,
           hads_1_complete = HADSSimu_1,
           responder_1_complete = responder_1,
           tx_2_complete = tx_2,
           hads_2_complete = HADSSimu_2) %>%
    mutate(dt_1 = ifelse(ltfu_B == "T1", NA, dt_1_complete),
           hads_1 = ifelse(ltfu_B == "T1", NA, hads_1_complete),
           responder_1 = ifelse(ltfu_B == "T1", NA, responder_1_complete),
           tx_2 = ifelse(ltfu_B == "T1", NA, tx_2_complete),
           hads_2 = ifelse(ltfu_B == "T2" | ltfu_B == "T1", NA, hads_2_complete)) %>%
    rename(ltfu = ltfu_B) %>%
    select(c(ID, sex, marital, language, birthcountry, educ, age, hads_base, dt_base, 
             tx_1, dt_1_complete, dt_1, hads_1_complete, hads_1, responder_1_complete,
             responder_1, tx_2_complete, tx_2, hads_2_complete, hads_2, ltfu))
  
  # Design C
  data_synth_C <- cbind(DataSimu, DTSimu_base, TxSimu_1_withID, DTSimu_1, HADSSimu_1, 
                        ResponderSimu_1_C, TxSimu_2_C, HADSSimu_2 = HADSSimu_2_C) %>%
    select(c(ID, colnames(DataSimu), DTSimu_base, tx_1, DTSimu_1, HADSSimu_1, responder_1,
             tx_2, HADSSimu_2))
  
  # Randomly select participants to drop-out (T1)
  ltfu_1_C <- rbinom(n = nrow(data_synth_C), size = 1, prob = prop_dropout_1)
  
  # Add ltfu at T1 to synthetic data
  data_synth_withltfu1_C <- cbind(data_synth_C, ltfu_1_C)
  
  # Participants who did not drop out at T1
  data_synth_notltfu1_C <- data_synth_withltfu1_C %>% filter(ltfu_1_C == 0)
  
  # Randomly select participants to drop-out (T2)
  ltfu_2_C <- rbinom(n = nrow(data_synth_notltfu1_C), size = 1, prob = prop_dropout_2)
  
  # Add ltfu at T2 to synthetic data
  data_synth_notltfu1_withltfu2_C <- cbind(data_synth_notltfu1_C, ltfu_2_C)
  
  # ltfu_1 and ltfu_2 in the same df
  data_synth_withltfu_C <- data_synth_withltfu1_C %>%
    left_join(data_synth_notltfu1_withltfu2_C %>% select(c(ID, ltfu_2_C)), by = "ID") %>%
    mutate(ltfu_2_C = ifelse(is.na(ltfu_2_C), 0, ltfu_2_C)) %>%
    mutate(ltfu_C = ifelse(ltfu_1_C == 1, "T1", ifelse(ltfu_2_C == 1, "T2", "N")))
  
  data_synth_C <- data_synth_withltfu_C %>%
    rename(dt_base = DTSimu_base,
           dt_1_complete = DTSimu_1,
           hads_1_complete = HADSSimu_1,
           responder_1_complete = responder_1,
           tx_2_complete = tx_2,
           hads_2_complete = HADSSimu_2) %>%
    mutate(dt_1 = ifelse(ltfu_C == "T1", NA, dt_1_complete),
           hads_1 = ifelse(ltfu_C == "T1", NA, hads_1_complete),
           responder_1 = ifelse(ltfu_C == "T1", NA, responder_1_complete),
           tx_2 = ifelse(ltfu_C == "T1", NA, tx_2_complete),
           hads_2 = ifelse(ltfu_C == "T2" | ltfu_C == "T1", NA, hads_2_complete)) %>%
    rename(ltfu = ltfu_C) %>%
    select(c(ID, sex, marital, language, birthcountry, educ, age, hads_base, dt_base, 
             tx_1, dt_1_complete, dt_1, hads_1_complete, hads_1, responder_1_complete,
             responder_1, tx_2_complete, tx_2, hads_2_complete, hads_2, ltfu))
  
  return(list(data_synth_A, data_synth_B, data_synth_C))
  
}

## Function to generate 1 data set under a given design (varying data generating mechanism) ----
# Function to generate 1 synthetic data set - ALL THREE DESIGNS (setting delta, re-defining beta4 for each design)
generate1dataset_designsABC <- function(realdata_augmentedbaseline, realdata_SMART,
                                        random_seed, n_obs, delta, S_hat, tx_1_haseffect) {
  # Input: realdata_augmentedbaseline is the real data dataframe that includes only baseline data (augmented),
  #        realdata_SMART is the real data dataframe including all variables from the internal SMART
  #        random_seed is a number for the random seed,
  #        n_obs is the number of observations (i.e., rows) to generate
  #        delta = standardized effect size
  #        S_hat = estimated variance component, sqrt((Var(Y|A_1=1,A_2=strat1)+Var(Y_A_1=0,A_2=strat2))/2)
  #        tx_1_haseffect = indicator of whether we impose no treatment effect at stage 1 (1 = has effect)
  # Output: dataframe of synthetic data under Design A, Design B, and Design C
  
  Cov_Discrete <- realdata_augmentedbaseline %>% 
    dplyr::select(c(sex, marital, language, birthcountry, educ)) %>%
    lapply(., as.factor) %>%
    as.data.frame()
  
  # Definition of the matrix of continuous covariates (at baseline)
  Cov_Cont <- realdata_augmentedbaseline %>%
    dplyr::select(c(age, hads_base))
  
  # Definition of the matrix of covariates (at baseline)
  Cov <- c(Cov_Discrete, Cov_Cont) %>% as.data.frame()
  
  # Marginal distribution of the covariates are estimated using empirical estimator 
  set.seed(random_seed)
  
  # Estimation of the R-vine model based on original data
  Rvine_dist <- Estimation_Copula(Cov)
  
  # Simulation of virtual patients based on the R-vine model and empirical distribution (of the original data)
  DataSimu <- Simulation_Copula(n_obs, Rvine_dist, Cov)  # n_obs = number of rows in original data
  
  # Next, generate synthetic DT at baseline (baseline var, but not shared by internal and external data sources)
  dt_base_mod <- lm(dt_base ~ age + as.factor(sex) + as.factor(marital) + as.factor(educ) +
                      as.factor(language) + as.factor(birthcountry) + as.numeric(hads_base),
                    data = realdata_SMART)
  
  DTSimu_base <- Simulation_DT_Base(Model = dt_base_mod, Covariates = DataSimu)
  
  TxSimu_1_withID <- Simulation_Treatment_1_HADS3levels(
    data = DataSimu %>% mutate(ID = 1:nrow(DataSimu), hads_base = as.numeric(hads_base)),
    seed = random_seed
  )
  
  TxSimu_1 <- TxSimu_1_withID %>% select(-c(ID))
  
  # Post-randomization variables have missingness, perform MI to get m complete data sets
  Imp_data <- realdata_SMART %>% 
    select(c(age, sex, marital, language, birthcountry, educ, dt_base, hads_base, 
             tx_1, dt_1, hads_1, responder_1, tx_2, hads_2, ltfu)) %>%
    mutate(age = as.numeric(age), sex = as.factor(sex), marital = as.factor(marital),
           language = as.factor(language), birthcountry = as.factor(birthcountry),
           educ = as.factor(educ), dt_base = as.numeric(dt_base), hads_base = as.numeric(hads_base),
           tx_1 = as.factor(tx_1), dt_1 = as.numeric(dt_1), hads_1 = as.numeric(hads_1),
           responder_1 = as.factor(responder_1), tx_2 = as.factor(tx_2), hads_2 = as.numeric(hads_2),
           ltfu = as.factor(ltfu))
  
  # Calculate m, where m ~ max(prop. of missingness)*100
  m <- round(colMeans(is.na(Imp_data)) %>% max()*100)  # m ~ prop. of missingness*100
  
  # Need to specify deterministic relationship between:
  # responder_1 and dt_base, dt_1
  # tx_2 and responder_1, tx_1
  
  ini <- mice(Imp_data, max = 0, print = FALSE)
  
  # Define imputation methods
  meth <- ini$meth
  meth["responder_1"] <- "~ ifelse(dt_base > 4, ifelse(dt_1 - dt_base < 0, '1', '0'), ifelse(dt_1 - dt_base < 2, '1', '0'))"
  meth["tx_2"] <- "~ ifelse(responder_1 == '1' & tx_1 == 'CT+LG', 'CT+LG',
       ifelse(responder_1 == '1' & tx_1 == 'CT', 'CT',
              ifelse(responder_1 == 0 & tx_1 == 'CT+LG', sample(c('CT+LG', 'CT+MI'), 1),
                     ifelse(responder_1 == 0 & tx_1 == 'CT', sample(c('CT', 'CT+LG'), 1), NA))))"
  
  # Define predictor matrix
  pred <- ini$pred
  # do NOT want to impute dt_base, dt_1 using responder_1 
  pred[c("dt_base", "dt_1"), "responder_1"] <- 0
  # do NOT want to impute tx_1, responder_1 using tx_2
  pred[c("tx_1", "responder_1"), "tx_2"] <- 0
  # ltfu is collinear with post-baseline missingness --> do NOT impute dt_1, hads_1, responder_1, tx_2, hads_2 using ltfu
  pred[c("dt_1", "hads_1", "responder_1", "tx_2", "hads_2"), "ltfu"] <- 0
  
  imp <- mice(data = Imp_data, 
              meth = meth,
              pred = pred,
              m = m, 
              maxit = 10, 
              seed = random_seed, 
              print = FALSE)
  
  # Simulation of post-randomization variables at T1 (dt_1, hads_1, responder_1)
  
  # dt_1:
  # First, fit execution model to each of the m imputed data sets
  dt_1_model_imp <- with(imp,
                         lm(dt_1 ~ as.factor(sex) + as.factor(marital) + as.factor(language) +
                              as.factor(birthcountry) + as.factor(educ) + age + as.numeric(hads_base) +
                              as.numeric(dt_base) + as.factor(tx_1))
  )
  
  # Predict dt_1 using each of the m model fits
  # Object to store predictions 
  # (rows represent synthetic patients/obs, columns represent predictions for the jth imputed data set)
  dt_1_preds <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m fitted models
    
    # Model parameters
    param <- dt_1_model_imp$analyses[[j]]$coefficients %>% as.matrix()
    
    # synthetic X,A with column of 1's for intercept
    data <- model.matrix(~ ., data = cbind(DataSimu, DTSimu_base, TxSimu_1))
    # data <- model.matrix(~ ., data = cbind(DTSimu_base, TxSimu_1))
    
    # Predicted synthetic Z_11, dt_1
    pred <- t(param) %*% t(data) %>% t()
    
    dt_1_preds[, j] <- pred
    
  }
  
  # Generate synthetic dt_1 for each of the m model fits
  # Object to store m sets of generated dt_1
  dt_1_mgens <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m model fits
    
    # Predicted values using j'th model
    dt_1_predict <- dt_1_preds[, j]
    
    # Residuals from j'th model
    dt_1_resid <- dt_1_model_imp$analyses[[j]]$residuals
    
    for (i in 1:nrow(dt_1_mgens)) {
      # Get all possible values of pred + resid for the i'th observation
      pred_resid_sums <- dt_1_predict[i] + dt_1_resid
      
      # All 0 <= pred + resid <= 10
      pred_resid_sums_pos <- pred_resid_sums[pred_resid_sums >= 0 & pred_resid_sums <= 10]
      
      # Randomly sample from sum values in range [0, 10]
      # It is possible that there are no sums in [0, 10] --> if this is the case, set sample_val = 0
      sample_val <- ifelse(length(pred_resid_sums_pos) != 0, sample(pred_resid_sums_pos, size = 1), 0) 
      # sample_val = sample(pred_resid_sums_pos, size = 1)
      
      # Save sample
      dt_1_mgens[i, j] <- sample_val
    }
    
  }
  
  # Generate final synthetic dt_1, Z^s_11 where _11 represents time T1 and variable 1 (dt_1) (all observed)
  DTSimu_1 <- apply(dt_1_mgens, 1, sample, size = 1) %>% round()
  
  
  # hads_1:
  # First, fit execution model to each of the m imputed data sets
  hads_1_model_imp <- with(imp,
                           lm(as.numeric(hads_1) ~ as.factor(sex) + as.factor(marital) + as.factor(language) +
                                as.factor(birthcountry) + as.factor(educ) + age + as.numeric(hads_base) + 
                                as.numeric(dt_base) + as.factor(tx_1))
  )
  
  # Predict hads_1 using each of the m model fits
  # Object to store predictions 
  # (rows represent synthetic patients/obs, columns represent predictions for the jth imputed data set)
  hads_1_preds <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m fitted models
    
    # Model parameters
    param <- hads_1_model_imp$analyses[[j]]$coefficients %>% as.matrix()
    
    # synthetic X,A with column of 1's for intercept
    data <- model.matrix(~ ., data = cbind(DataSimu, DTSimu_base, TxSimu_1))
    
    # Predicted synthetic Z_12, hads_1
    pred <- t(param) %*% t(data) %>% t()
    
    hads_1_preds[, j] <- pred
    
  }
  
  # Generate synthetic hads_1 for each of the m model fits
  # Object to store m sets of generated hads_1
  hads_1_mgens <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m model fits
    
    # Predicted values using j'th model
    hads_1_predict <- hads_1_preds[, j]
    
    # Residuals from j'th model
    hads_1_resid <- hads_1_model_imp$analyses[[j]]$residuals
    
    for (i in 1:nrow(hads_1_mgens)) {
      # Get all possible values of pred + resid for the i'th observation
      pred_resid_sums <- hads_1_predict[i] + hads_1_resid
      
      # All 0 <= pred + resid <= 21
      pred_resid_sums_pos <- pred_resid_sums[pred_resid_sums >= 0 & pred_resid_sums <= 21]
      
      # Randomly sample from sum values in range [0, 21]
      # It is possible that there are no sums in [0, 21] --> if this is the case, set sample_val = 0
      sample_val <- ifelse(length(pred_resid_sums_pos) != 0, sample(pred_resid_sums_pos, size = 1), 0) 
      # sample_val = sample(pred_resid_sums_pos, size = 1)
      
      # Save sample
      hads_1_mgens[i, j] <- sample_val
    }
    
  }
  
  # Generate final synthetic hads_1, Z^s_12 where _12 represents time T1 and variable 2 (hads_1) (all observed)
  HADSSimu_1 <- apply(hads_1_mgens, 1, sample, size = 1)
  
  
  # responder_1 (Z_13)
  # Note: we don't need to fit a model for this, as responder status is determined based on
  #       tailoring variables
  
  # DESIGN A ORIGINAL:
  # If dt_base >= 5, then responder = Yes if dt_1 - dt_base >= 1
  # If dt_base < 5, then responder = Yes if dt_1 < dt_base
  ResponderSimu_1_A <- cbind(DTSimu_base, DTSimu_1) %>%
    as.data.frame() %>%
    mutate(responder_1 = case_when(DTSimu_base >= 5 & DTSimu_1 - DTSimu_base < 0 ~ "1",
                                   DTSimu_base < 5 & DTSimu_1 - DTSimu_base < 2 ~ "1",
                                   TRUE ~ "0")) %>%
    select(responder_1)
  
  ResponderSimu_1_A_withtx <- cbind(TxSimu_1, ResponderSimu_1_A)
  p1_obs_designA <- (ResponderSimu_1_A_withtx %>% 
                       filter(tx_1 == "CT+LG") %>% 
                       select(responder_1) %>% 
                       unlist() %>% as.numeric() %>% 
                       sum)/nrow(ResponderSimu_1_A_withtx %>% filter(tx_1 == "CT+LG"))
  p0_obs_designA <- (ResponderSimu_1_A_withtx %>%
                       filter(tx_1 == "CT") %>% 
                       select(responder_1) %>% 
                       unlist() %>% as.numeric() %>% 
                       sum)/nrow(ResponderSimu_1_A_withtx %>% filter(tx_1 == "CT"))
  
  # DESIGN B:
  # If dt_1 <= median(dt_1) for tx_1 subgroup, then responder = Yes
  dt_1_med_1 <- cbind(DataSimu, DTSimu_base, TxSimu_1, DTSimu_1) %>% 
    as.data.frame() %>%
    filter(tx_1 == "CT+LG") %>%
    dplyr::select(DTSimu_1) %>% 
    unlist() %>%
    median()
  
  dt_1_med_0 <- cbind(DataSimu, DTSimu_base, TxSimu_1, DTSimu_1) %>% 
    as.data.frame() %>%
    filter(tx_1 == "CT") %>% 
    dplyr::select(DTSimu_1) %>% 
    unlist() %>%
    median()
  
  ResponderSimu_1_B <- cbind(TxSimu_1, DTSimu_1) %>%
    as.data.frame() %>%
    mutate(responder_1 = case_when(TxSimu_1 == "CT+LG" & DTSimu_1 <= dt_1_med_1 ~ "1",
                                   TxSimu_1 == "CT" & DTSimu_1 <= dt_1_med_0 ~ "1",
                                   TRUE ~ "0")) %>%
    select(responder_1)
  
  ResponderSimu_1_B_withtx <- cbind(TxSimu_1, ResponderSimu_1_B)
  (p1_obs_designB <- (ResponderSimu_1_B_withtx %>% 
                        filter(tx_1 == "CT+LG") %>% 
                        select(responder_1) %>% 
                        unlist() %>% as.numeric() %>% 
                        sum)/nrow(ResponderSimu_1_B_withtx %>% 
                                    filter(tx_1 == "CT+LG")))
  (p0_obs_designB <- (ResponderSimu_1_B_withtx %>% 
                        filter(tx_1 == "CT") %>% 
                        select(responder_1) %>% 
                        unlist() %>% as.numeric() %>% 
                        sum)/nrow(ResponderSimu_1_B_withtx %>% 
                                    filter(tx_1 == "CT")))
  
  # DESIGN C:
  # If dt_1 < dt_base - 1 AND hads_1 < hads_base - 1, then responder = Yes
  # i.e., both DT and HADS need to drop by 2 or more in order to be considered a responder
  ResponderSimu_1_C <- cbind(DTSimu_base, DTSimu_1, DataSimu %>% dplyr::select(hads_base), HADSSimu_1) %>%
    as.data.frame() %>%
    mutate(responder_1 = case_when(DTSimu_1 < DTSimu_base - 1 & HADSSimu_1 < hads_base - 1 ~ "1",
                                   TRUE ~ "0")) %>%
    select(responder_1)
  
  ResponderSimu_1_C_withtx <- cbind(TxSimu_1, ResponderSimu_1_C)
  (p1_obs_designC <- (ResponderSimu_1_C_withtx %>% 
                        filter(tx_1 == "CT+LG") %>% 
                        select(responder_1) %>% 
                        unlist() %>% as.numeric() %>% 
                        sum)/nrow(ResponderSimu_1_C_withtx %>% 
                                    filter(tx_1 == "CT+LG")))
  (p0_obs_designC <- (ResponderSimu_1_C_withtx %>% 
                        filter(tx_1 == "CT") %>% 
                        select(responder_1) %>% 
                        unlist() %>% as.numeric() %>% 
                        sum)/nrow(ResponderSimu_1_C_withtx %>% 
                                    filter(tx_1 == "CT")))
  
  
  # Simulation of synthetic treatment allocation - Stage 2
  
  # DESIGN A:
  TxSimu_2_A <- Simulation_Treatment_2(tx_1_data_synth = TxSimu_1,
                                       responder_1_data_synth = ResponderSimu_1_A,
                                       seed = random_seed)
  
  # DESIGN B:
  TxSimu_2_B <- Simulation_Treatment_2(tx_1_data_synth = TxSimu_1,
                                       responder_1_data_synth = ResponderSimu_1_B,
                                       seed = random_seed)
  
  # DESIGN C:
  TxSimu_2_C <- Simulation_Treatment_2(tx_1_data_synth = TxSimu_1,
                                       responder_1_data_synth = ResponderSimu_1_C,
                                       seed = random_seed)
  
  # Simulation of post-randomization variable at T2 (hads_2)
  
  # First, fit execution model to each of the m imputed data sets
  hads_2_model_imp <- with(imp,
                           lm(as.numeric(hads_2) ~ as.factor(sex) + as.factor(marital) + as.factor(language) +
                                as.factor(birthcountry) + as.factor(educ) + age + as.numeric(hads_base) + 
                                as.numeric(dt_base) + as.factor(tx_1) + as.numeric(dt_1) + as.numeric(hads_1) +
                                as.factor(tx_2))
  )
  
  # Predict hads_2 using each of the m model fits
  # Object to store predictions 
  # (rows represent synthetic patients/obs, columns represent predictions for the jth imputed data set)
  hads_2_preds_A <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  hads_2_preds_B <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  hads_2_preds_C <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m fitted models
    
    # Model parameters
    param <- hads_2_model_imp$analyses[[j]]$coefficients %>% as.matrix()
    
    # Fix parameter for tx_2_CT+MI = 1/(1 - p1)*(delta*S_hat - beta2 - beta3(p1+p0-1)) 
    # tx_2_CT+MI is stepped-up care not shared by both strategies
    # delta = effect size, S_hat = estimated variance component
    # p1 = prop. responders at stage 1 to tx=CT+LG, p0 = prop. responders at stage 1 to tx=CT
    param_simdelta_A <- param 
    param_simdelta_B <- param 
    param_simdelta_C <- param 
    
    beta2_hat <- param[rownames(param) %in% c("as.factor(tx_1)CT+LG")]
    beta3_hat <- param[rownames(param) %in% c("as.factor(tx_2)CT+LG")]
    
    # Fix beta4 based on estimated coef of beta2, beta3
    # Use observed p1, p0 from synthetic data
    param_simdelta_A[rownames(param_simdelta_A) %in% c("as.factor(tx_2)CT+MI")] <- (1/(1 - p1_obs_designA))*((delta*S_hat) - beta2_hat - (beta3_hat*(p1_obs_designA + p0_obs_designA - 1)))
    param_simdelta_B[rownames(param_simdelta_B) %in% c("as.factor(tx_2)CT+MI")] <- (1/(1 - p1_obs_designB))*((delta*S_hat) - beta2_hat - (beta3_hat*(p1_obs_designB + p0_obs_designB - 1)))
    param_simdelta_C[rownames(param_simdelta_C) %in% c("as.factor(tx_2)CT+MI")] <- (1/(1 - p1_obs_designC))*((delta*S_hat) - beta2_hat - (beta3_hat*(p1_obs_designC + p0_obs_designC - 1)))
    
    # synthetic X,A with column of 1's for intercept
    data_A <- model.matrix(~ ., data = cbind(DataSimu, DTSimu_base, TxSimu_1, DTSimu_1,
                                             HADSSimu_1, TxSimu_2_A))
    data_B <- model.matrix(~ ., data = cbind(DataSimu, DTSimu_base, TxSimu_1, DTSimu_1, 
                                             HADSSimu_1, TxSimu_2_B))
    data_C <- model.matrix(~ ., data = cbind(DataSimu, DTSimu_base, TxSimu_1, DTSimu_1, 
                                             HADSSimu_1, TxSimu_2_C))
    
    # Predicted synthetic dt_2
    pred_A <- t(param_simdelta_A) %*% t(data_A) %>% t()
    pred_B <- t(param_simdelta_B) %*% t(data_B) %>% t()
    pred_C <- t(param_simdelta_C) %*% t(data_C) %>% t()
    
    hads_2_preds_A[, j] <- pred_A
    hads_2_preds_B[, j] <- pred_B
    hads_2_preds_C[, j] <- pred_C
    
  }
  
  # Generate synthetic hads_2 for each of the m model fits
  # Object to store m sets of generated dt_1
  hads_2_mgens_A <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  hads_2_mgens_B <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  hads_2_mgens_C <- matrix(data = NA, nrow = n_obs, ncol = m) %>% as.data.frame()
  
  for(j in 1:m) {  # iterate through all columns i.e., all m model fits
    
    # Predicted values using j'th model
    hads_2_predict_A <- hads_2_preds_A[, j]
    hads_2_predict_B <- hads_2_preds_B[, j]
    hads_2_predict_C <- hads_2_preds_C[, j]
    
    # Residuals from j'th model
    hads_2_resid <- hads_2_model_imp$analyses[[j]]$residuals
    
    for (i in 1:n_obs) {
      # Get all possible values of pred + resid for the i'th observation
      pred_resid_sums_A <- hads_2_predict_A[i] + hads_2_resid
      pred_resid_sums_B <- hads_2_predict_B[i] + hads_2_resid
      pred_resid_sums_C <- hads_2_predict_C[i] + hads_2_resid
      
      # All 0 <= pred + resid <= 21
      pred_resid_sums_pos_A <- pred_resid_sums_A[pred_resid_sums_A >= 0 & pred_resid_sums_A <= 21]
      pred_resid_sums_pos_B <- pred_resid_sums_B[pred_resid_sums_B >= 0 & pred_resid_sums_B <= 21]
      pred_resid_sums_pos_C <- pred_resid_sums_C[pred_resid_sums_C >= 0 & pred_resid_sums_C <= 21]
      
      # Randomly sample from sum values in range [0, 21]
      # It is possible that there are no sums in [0, 21] --> if this is the case, set sample_val = 0
      sample_val_A <- ifelse(length(pred_resid_sums_pos_A) != 0, sample(pred_resid_sums_pos_A, size = 1), 0) 
      sample_val_B <- ifelse(length(pred_resid_sums_pos_B) != 0, sample(pred_resid_sums_pos_B, size = 1), 0) 
      sample_val_C <- ifelse(length(pred_resid_sums_pos_C) != 0, sample(pred_resid_sums_pos_C, size = 1), 0) 
      
      # Save sample
      hads_2_mgens_A[i, j] <- sample_val_A
      hads_2_mgens_B[i, j] <- sample_val_B
      hads_2_mgens_C[i, j] <- sample_val_C
    }
    
  }
  
  # Generate final synthetic hads_2, Z^s_21 where _21 represents time T2 and variable 1 (hads_2) (all observed)
  HADSSimu_2_A <- apply(hads_2_mgens_A, 1, sample, size = 1)
  HADSSimu_2_B <- apply(hads_2_mgens_B, 1, sample, size = 1)
  HADSSimu_2_C <- apply(hads_2_mgens_C, 1, sample, size = 1)
  
  # Finally, generate synthetic missingness (ltfu)
  
  # Too few observations to fit execution model
  # Instead, randomly select participants to drop out with observed proportion of drop-out at each stage 
  
  prop_dropout_1 <- sum(Imp_data$ltfu == "T1")/nrow(Imp_data)
  prop_dropout_2 <- sum(Imp_data$ltfu == "T2")/nrow(Imp_data[Imp_data$ltfu != "T1", ])
  
  # Design A
  data_synth_A <- cbind(DataSimu, DTSimu_base, TxSimu_1_withID, DTSimu_1, HADSSimu_1, 
                        ResponderSimu_1_A, TxSimu_2_A, HADSSimu_2 = HADSSimu_2_A) %>%
    select(c(ID, colnames(DataSimu), DTSimu_base, tx_1, DTSimu_1, HADSSimu_1, responder_1,
             tx_2, HADSSimu_2))
  
  # Randomly select participants to drop-out (T1)
  ltfu_1_A <- rbinom(n = nrow(data_synth_A), size = 1, prob = prop_dropout_1)
  
  # Add ltfu at T1 to synthetic data
  data_synth_withltfu1_A <- cbind(data_synth_A, ltfu_1_A)
  
  # Participants who did not drop out at T1
  data_synth_notltfu1_A <- data_synth_withltfu1_A %>% filter(ltfu_1_A == 0)
  
  # Randomly select participants to drop-out (T2)
  ltfu_2_A <- rbinom(n = nrow(data_synth_notltfu1_A), size = 1, prob = prop_dropout_2)
  
  # Add ltfu at T2 to synthetic data
  data_synth_notltfu1_withltfu2_A <- cbind(data_synth_notltfu1_A, ltfu_2_A)
  
  # ltfu_1 and ltfu_2 in the same df
  data_synth_withltfu_A <- data_synth_withltfu1_A %>%
    left_join(data_synth_notltfu1_withltfu2_A %>% select(c(ID, ltfu_2_A)), by = "ID") %>%
    mutate(ltfu_2_A = ifelse(is.na(ltfu_2_A), 0, ltfu_2_A)) %>%
    mutate(ltfu_A = ifelse(ltfu_1_A == 1, "T1", ifelse(ltfu_2_A == 1, "T2", "N")))
  
  data_synth_A <- data_synth_withltfu_A %>%
    rename(dt_base = DTSimu_base,
           dt_1_complete = DTSimu_1,
           hads_1_complete = HADSSimu_1,
           responder_1_complete = responder_1,
           tx_2_complete = tx_2,
           hads_2_complete = HADSSimu_2) %>%
    mutate(dt_1 = ifelse(ltfu_A == "T1", NA, dt_1_complete),
           hads_1 = ifelse(ltfu_A == "T1", NA, hads_1_complete),
           responder_1 = ifelse(ltfu_A == "T1", NA, responder_1_complete),
           tx_2 = ifelse(ltfu_A == "T1", NA, tx_2_complete),
           hads_2 = ifelse(ltfu_A == "T2" | ltfu_A == "T1", NA, hads_2_complete)) %>%
    rename(ltfu = ltfu_A) %>%
    select(c(ID, sex, marital, language, birthcountry, educ, age, hads_base, dt_base, 
             tx_1, dt_1_complete, dt_1, hads_1_complete, hads_1, responder_1_complete,
             responder_1, tx_2_complete, tx_2, hads_2_complete, hads_2, ltfu))
  
  # Design B
  data_synth_B <- cbind(DataSimu, DTSimu_base, TxSimu_1_withID, DTSimu_1, HADSSimu_1, 
                        ResponderSimu_1_B, TxSimu_2_B, HADSSimu_2 = HADSSimu_2_B) %>%
    select(c(ID, colnames(DataSimu), DTSimu_base, tx_1, DTSimu_1, HADSSimu_1, responder_1,
             tx_2, HADSSimu_2))
  
  # Randomly select participants to drop-out (T1)
  ltfu_1_B <- rbinom(n = nrow(data_synth_B), size = 1, prob = prop_dropout_1)
  
  # Add ltfu at T1 to synthetic data
  data_synth_withltfu1_B <- cbind(data_synth_B, ltfu_1_B)
  
  # Participants who did not drop out at T1
  data_synth_notltfu1_B <- data_synth_withltfu1_B %>% filter(ltfu_1_B == 0)
  
  # Randomly select participants to drop-out (T2)
  ltfu_2_B <- rbinom(n = nrow(data_synth_notltfu1_B), size = 1, prob = prop_dropout_2)
  
  # Add ltfu at T2 to synthetic data
  data_synth_notltfu1_withltfu2_B <- cbind(data_synth_notltfu1_B, ltfu_2_B)
  
  # ltfu_1 and ltfu_2 in the same df
  data_synth_withltfu_B <- data_synth_withltfu1_B %>%
    left_join(data_synth_notltfu1_withltfu2_B %>% select(c(ID, ltfu_2_B)), by = "ID") %>%
    mutate(ltfu_2_B = ifelse(is.na(ltfu_2_B), 0, ltfu_2_B)) %>%
    mutate(ltfu_B = ifelse(ltfu_1_B == 1, "T1", ifelse(ltfu_2_B == 1, "T2", "N")))
  
  data_synth_B <- data_synth_withltfu_B %>%
    rename(dt_base = DTSimu_base,
           dt_1_complete = DTSimu_1,
           hads_1_complete = HADSSimu_1,
           responder_1_complete = responder_1,
           tx_2_complete = tx_2,
           hads_2_complete = HADSSimu_2) %>%
    mutate(dt_1 = ifelse(ltfu_B == "T1", NA, dt_1_complete),
           hads_1 = ifelse(ltfu_B == "T1", NA, hads_1_complete),
           responder_1 = ifelse(ltfu_B == "T1", NA, responder_1_complete),
           tx_2 = ifelse(ltfu_B == "T1", NA, tx_2_complete),
           hads_2 = ifelse(ltfu_B == "T2" | ltfu_B == "T1", NA, hads_2_complete)) %>%
    rename(ltfu = ltfu_B) %>%
    select(c(ID, sex, marital, language, birthcountry, educ, age, hads_base, dt_base, 
             tx_1, dt_1_complete, dt_1, hads_1_complete, hads_1, responder_1_complete,
             responder_1, tx_2_complete, tx_2, hads_2_complete, hads_2, ltfu))
  
  # Design C
  data_synth_C <- cbind(DataSimu, DTSimu_base, TxSimu_1_withID, DTSimu_1, HADSSimu_1, 
                        ResponderSimu_1_C, TxSimu_2_C, HADSSimu_2 = HADSSimu_2_C) %>%
    select(c(ID, colnames(DataSimu), DTSimu_base, tx_1, DTSimu_1, HADSSimu_1, responder_1,
             tx_2, HADSSimu_2))
  
  # Randomly select participants to drop-out (T1)
  ltfu_1_C <- rbinom(n = nrow(data_synth_C), size = 1, prob = prop_dropout_1)
  
  # Add ltfu at T1 to synthetic data
  data_synth_withltfu1_C <- cbind(data_synth_C, ltfu_1_C)
  
  # Participants who did not drop out at T1
  data_synth_notltfu1_C <- data_synth_withltfu1_C %>% filter(ltfu_1_C == 0)
  
  # Randomly select participants to drop-out (T2)
  ltfu_2_C <- rbinom(n = nrow(data_synth_notltfu1_C), size = 1, prob = prop_dropout_2)
  
  # Add ltfu at T2 to synthetic data
  data_synth_notltfu1_withltfu2_C <- cbind(data_synth_notltfu1_C, ltfu_2_C)
  
  # ltfu_1 and ltfu_2 in the same df
  data_synth_withltfu_C <- data_synth_withltfu1_C %>%
    left_join(data_synth_notltfu1_withltfu2_C %>% select(c(ID, ltfu_2_C)), by = "ID") %>%
    mutate(ltfu_2_C = ifelse(is.na(ltfu_2_C), 0, ltfu_2_C)) %>%
    mutate(ltfu_C = ifelse(ltfu_1_C == 1, "T1", ifelse(ltfu_2_C == 1, "T2", "N")))
  
  data_synth_C <- data_synth_withltfu_C %>%
    rename(dt_base = DTSimu_base,
           dt_1_complete = DTSimu_1,
           hads_1_complete = HADSSimu_1,
           responder_1_complete = responder_1,
           tx_2_complete = tx_2,
           hads_2_complete = HADSSimu_2) %>%
    mutate(dt_1 = ifelse(ltfu_C == "T1", NA, dt_1_complete),
           hads_1 = ifelse(ltfu_C == "T1", NA, hads_1_complete),
           responder_1 = ifelse(ltfu_C == "T1", NA, responder_1_complete),
           tx_2 = ifelse(ltfu_C == "T1", NA, tx_2_complete),
           hads_2 = ifelse(ltfu_C == "T2" | ltfu_C == "T1", NA, hads_2_complete)) %>%
    rename(ltfu = ltfu_C) %>%
    select(c(ID, sex, marital, language, birthcountry, educ, age, hads_base, dt_base, 
             tx_1, dt_1_complete, dt_1, hads_1_complete, hads_1, responder_1_complete,
             responder_1, tx_2_complete, tx_2, hads_2_complete, hads_2, ltfu))
  
  return(list(data_synth_A, data_synth_B, data_synth_C))
  
}


# Delta distribution under competing designs ----

# Implementation:
# Under Design A, fix delta (e.g., 0.2) --> for each imputed data set, this gives a value of beta4
# --> under design A, there is a distribution of beta4 values ("nature")
# Fixing beta4 as values found in previous step, apply new definition of responder status
# --> apply Design B, Design C
# Determine value of delta under competing designs with fixed beta4
# --> recover distribution of delta over imputed data sets
#     i.e., each imputed data set has 1 beta4 value, which corresponds to 1 delta value
#           then repeat for all m imputed data sets to get distribution of delta values

# Calculate deltas for each competing design

# Designs B and C
deltas_designsBC <- deltadist(
  # NOTE: data not publicly-available, but this data set is the augmented baseline data (internal pilot + external)
  realdata_augmentedbaseline = data_cancer_augmentbase_nodyad,
  # NOTE: data not publicly-available, but this data set is the pilot SMART (all variables), 
  #       restricted to rows with no missingness at baseline
  realdata_SMART = data_cancer_smart_ccbase,
  random_seed = 20260323,
  n_obs = 10000,  # a sufficiently-large generated data set
  m = 1000,
  delta = 0.2,
  S_hat = S_hat_pilot,
  tx_1_haseffect = 0
)

# Power analysis for competing designs under fixed data generating mechanism ----

realdata_augmentedbaseline_fix <- data_cancer_augmentbase_nodyad
realdata_SMART_fix <- data_cancer_smart_ccbase
random_seed <- 20260302
delta_fix <- 0.2 # standardized effect size
tx_1_haseffect_fix <- 0  # do NOT impose that tx_1 has no effect

# Set up cluster
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
registerDoRNG(20260302)  # any fixed seed

start_time <- Sys.time()
# Parallel loop over sample sizes
sim.power.allN.fixedmechanism.alldesigns.StratAvStratC.parallel <- foreach(
  n_obs_current = seq(100, 3000, by = 100),
  .combine = rbind, 
  .packages = c("tidyverse", "rvinecopulib", "e1071", "caret", "EnvStats", "truncnorm", 
                "mice", "stats")
) %dopar% {
  # Define number of simulations per sample size
  n_sim <- 1000
  
  # Initialize object to store power estimates per simulation run
  # Design A
  pvals_currentn_A <- numeric(n_sim)
  # Design B
  pvals_currentn_B <- numeric(n_sim) 
  # Design C
  pvals_currentn_C <- numeric(n_sim) 
  
  # Calculate power across n_sim simulation runs
  for(j in 1:n_sim) {
    
    # Generate data set for each design
    data_synth_alldesigns <- generate1dataset_alldesigns_fixedmech(
      realdata_augmentedbaseline = realdata_augmentedbaseline_fix,
      realdata_SMART = realdata_SMART_fix,
      random_seed = sample.int(.Machine$integer.max, 1),
      n_obs = n_obs_current,  
      delta = delta_fix,  # standardized effect size
      S_hat = S_hat_pilot,  # empirical variance from pilot
      p1_hat_A = p1_hat_pilot_designA,  # empirical proportion of stage 1 responders, tx = 1 (CT+LG) (Design A)
      p0_hat_A = p0_hat_pilot_designA,  # empirical proportion of stage 1 responders, tx = 0 (CT) (Design A)
      p1_hat_B = p1_hat_pilot_designB,  # empirical proportion of stage 1 responders, tx = 1 (CT+LG) (Design B)
      p0_hat_B = p0_hat_pilot_designB,  # empirical proportion of stage 1 responders, tx = 0 (CT) (Design B)
      p1_hat_C = p1_hat_pilot_designC,  # empirical proportion of stage 1 responders, tx = 1 (CT+LG) (Design C)
      p0_hat_C = p0_hat_pilot_designC,  # empirical proportion of stage 1 responders, tx = 0 (CT) (Design C)
      tx_1_haseffect = tx_1_haseffect_fix  # whether tx_1 has no effect
    )
    
    # Select participants who received Strategy A or Strategy C
    
    # Design A
    data_synth_StratAvStratC_A <- data_synth_alldesigns[[1]] %>%
      mutate(strategyA = case_when(tx_1 == "CT+LG" & 
                                     (responder_1_complete == 1 | 
                                        responder_1_complete == 0 & tx_2_complete == "CT+MI") ~ 1,
                                   TRUE ~ 0),
             strategyB = case_when(tx_1 == "CT+LG" & 
                                     (responder_1_complete == 1 | 
                                        responder_1_complete == 0 & tx_2_complete == "CT+LG") ~ 1,
                                   TRUE ~ 0),
             strategyC = case_when(tx_1 == "CT" &
                                     (responder_1_complete == 1 |
                                        responder_1_complete == 0 & tx_2_complete == "CT+LG") ~ 1,
                                   TRUE ~ 0),
             strategyD = case_when(tx_1 == "CT" &
                                     (responder_1_complete == 1 |
                                        responder_1_complete == 0 & tx_2_complete == "CT") ~ 1,
                                   TRUE ~ 0)) %>%
      filter(strategyA == 1 | strategyC == 1) %>%
      mutate(strategy = case_when(strategyA == 1 ~ "A",
                                  strategyC == 1 ~ "C"))
    
    # Design B
    data_synth_StratAvStratC_B <- data_synth_alldesigns[[2]] %>%
      mutate(strategyA = case_when(tx_1 == "CT+LG" & 
                                     (responder_1_complete == 1 | 
                                        responder_1_complete == 0 & tx_2_complete == "CT+MI") ~ 1,
                                   TRUE ~ 0),
             strategyB = case_when(tx_1 == "CT+LG" & 
                                     (responder_1_complete == 1 | 
                                        responder_1_complete == 0 & tx_2_complete == "CT+LG") ~ 1,
                                   TRUE ~ 0),
             strategyC = case_when(tx_1 == "CT" &
                                     (responder_1_complete == 1 |
                                        responder_1_complete == 0 & tx_2_complete == "CT+LG") ~ 1,
                                   TRUE ~ 0),
             strategyD = case_when(tx_1 == "CT" &
                                     (responder_1_complete == 1 |
                                        responder_1_complete == 0 & tx_2_complete == "CT") ~ 1,
                                   TRUE ~ 0)) %>%
      filter(strategyA == 1 | strategyC == 1) %>%
      mutate(strategy = case_when(strategyA == 1 ~ "A",
                                  strategyC == 1 ~ "C"))
    
    # Design C
    data_synth_StratAvStratC_C <- data_synth_alldesigns[[3]] %>%
      mutate(strategyA = case_when(tx_1 == "CT+LG" & 
                                     (responder_1_complete == 1 | 
                                        responder_1_complete == 0 & tx_2_complete == "CT+MI") ~ 1,
                                   TRUE ~ 0),
             strategyB = case_when(tx_1 == "CT+LG" & 
                                     (responder_1_complete == 1 | 
                                        responder_1_complete == 0 & tx_2_complete == "CT+LG") ~ 1,
                                   TRUE ~ 0),
             strategyC = case_when(tx_1 == "CT" &
                                     (responder_1_complete == 1 |
                                        responder_1_complete == 0 & tx_2_complete == "CT+LG") ~ 1,
                                   TRUE ~ 0),
             strategyD = case_when(tx_1 == "CT" &
                                     (responder_1_complete == 1 |
                                        responder_1_complete == 0 & tx_2_complete == "CT") ~ 1,
                                   TRUE ~ 0)) %>%
      filter(strategyA == 1 | strategyC == 1) %>%
      mutate(strategy = case_when(strategyA == 1 ~ "A",
                                  strategyC == 1 ~ "C"))
    
    # Perform t-test and record p-value (using OBSERVED data - missing values are removed)
    
    # It is technically possible that there is no one who received 
    # Strategy A or C (has yet to happen in our sims)
    # To prevent this from crashing the code, include a fail-safe catch
    # for this scenario
    
    # Design A
    if (length(unique(data_synth_StratAvStratC_A$strategy)) == 2) {
      
      pvals_currentn_A[j] <- t.test(
        hads_2 ~ strategy,
        data = data_synth_StratAvStratC_A
      )$p.value
      
    } else {
      
      pvals_currentn_A[j] <- NA
      
    }
    
    # Design B
    if (length(unique(data_synth_StratAvStratC_B$strategy)) == 2) {
      
      pvals_currentn_B[j] <- t.test(
        hads_2 ~ strategy,
        data = data_synth_StratAvStratC_B
      )$p.value
      
    } else {
      
      pvals_currentn_B[j] <- NA
      
    }
    
    # Design C
    if (length(unique(data_synth_StratAvStratC_C$strategy)) == 2) {
      
      pvals_currentn_C[j] <- t.test(
        hads_2 ~ strategy,
        data = data_synth_StratAvStratC_C
      )$p.value
      
    } else {
      
      pvals_currentn_C[j] <- NA
      
    }
  }
  
  # Estimate power as proportion of null hypothesis rejections out of number of simulation runs
  # i.e., count(p-vals < 0.05)/n_sim
  # Design A
  power_est_A <- mean(pvals_currentn_A < 0.05, na.rm = TRUE)
  # Design B
  power_est_B <- mean(pvals_currentn_B < 0.05, na.rm = TRUE)
  # Design C
  power_est_C <- mean(pvals_currentn_C < 0.05, na.rm = TRUE)
  
  # Store the number of "valid" simulations (i.e., t-test was successfully performed)
  # Design A
  n_valid_A <- sum(!is.na(pvals_currentn_A))
  # Design B
  n_valid_B <- sum(!is.na(pvals_currentn_B))
  # Design C
  n_valid_C <- sum(!is.na(pvals_currentn_C))
  
  data.frame(
    n = rep(n_obs_current, 3),
    design = c("A", "B", "C"),
    power = c(power_est_A, power_est_B, power_est_C),
    n_sim = c(n_valid_A, n_valid_B, n_valid_C)
  )
  
}

stopCluster(cl)

end_time <- Sys.time()
(time_taken <- end_time - start_time)

# Estimate power over N simulations, per sample size (fix delta and vary data generating mechanism) ----


# All Designs, VARYING DELTA

realdata_augmentedbaseline_fix <- data_cancer_augmentbase_nodyad
realdata_SMART_fix <- data_cancer_smart_ccbase
random_seed <- 20260302
delta_fix <- 0.2 # standardized effect size (0.2, 0.5)
tx_1_haseffect_fix <- 0  # do NOT impose that tx_1 has no effect
realdata_augmentedbaseline = realdata_augmentedbaseline_fix
realdata_SMART = realdata_SMART_fix
random_seed = random_seed
delta = delta_fix  # standardized effect size
S_hat = S_hat_pilot  # empirical variance from pilot
tx_1_haseffect = tx_1_haseffect_fix  # whether tx_1 has no effect

# Set up cluster
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
registerDoRNG(20260302)  # any fixed seed

start_time <- Sys.time()
# Parallel loop over sample sizes
sim.power.allN.alldesigns.StratAvStratC.parallel <- foreach(
  n_obs_current = seq(100, 3000, by = 100),
  .combine = rbind, 
  .packages = c("tidyverse", "rvinecopulib", "e1071", "caret", "EnvStats", "truncnorm", 
                "mice", "stats")
) %dopar% {
  # Define number of simulations per sample size
  n_sim <- 1000
  
  # Initialize object to store power estimates per simulation run
  pvals_currentn_A <- numeric(n_sim)
  pvals_currentn_B <- numeric(n_sim)
  pvals_currentn_C <- numeric(n_sim)
  
  # Calculate power across n_sim simulation runs
  for(j in 1:n_sim) {
    
    # Generate data set 
    data_synth_ABC <- generate1dataset_designsABC(
      realdata_augmentedbaseline = realdata_augmentedbaseline_fix,
      realdata_SMART = realdata_SMART_fix,
      random_seed = sample.int(.Machine$integer.max, 1),
      n_obs = n_obs_current,  
      delta = delta_fix,  # standardized effect size
      S_hat = S_hat_pilot,  # empirical variance from pilot
      tx_1_haseffect = tx_1_haseffect_fix  # whether tx_1 has no effect
    )
    
    # Select participants who received Strategy A or Strategy C
    
    # Design A
    data_synth_A_StratAvStratC <- data_synth_ABC[[1]] %>%
      mutate(strategyA = case_when(tx_1 == "CT+LG" & 
                                     (responder_1_complete == 1 | 
                                        responder_1_complete == 0 & tx_2_complete == "CT+MI") ~ 1,
                                   TRUE ~ 0),
             strategyB = case_when(tx_1 == "CT+LG" & 
                                     (responder_1_complete == 1 | 
                                        responder_1_complete == 0 & tx_2_complete == "CT+LG") ~ 1,
                                   TRUE ~ 0),
             strategyC = case_when(tx_1 == "CT" &
                                     (responder_1_complete == 1 |
                                        responder_1_complete == 0 & tx_2_complete == "CT+LG") ~ 1,
                                   TRUE ~ 0),
             strategyD = case_when(tx_1 == "CT" &
                                     (responder_1_complete == 1 |
                                        responder_1_complete == 0 & tx_2_complete == "CT") ~ 1,
                                   TRUE ~ 0)) %>%
      filter(strategyA == 1 | strategyC == 1) %>%
      mutate(strategy = case_when(strategyA == 1 ~ "A",
                                  strategyC == 1 ~ "C"))
    
    # Design B
    data_synth_B_StratAvStratC <- data_synth_ABC[[2]] %>%
      mutate(strategyA = case_when(tx_1 == "CT+LG" & 
                                     (responder_1_complete == 1 | 
                                        responder_1_complete == 0 & tx_2_complete == "CT+MI") ~ 1,
                                   TRUE ~ 0),
             strategyB = case_when(tx_1 == "CT+LG" & 
                                     (responder_1_complete == 1 | 
                                        responder_1_complete == 0 & tx_2_complete == "CT+LG") ~ 1,
                                   TRUE ~ 0),
             strategyC = case_when(tx_1 == "CT" &
                                     (responder_1_complete == 1 |
                                        responder_1_complete == 0 & tx_2_complete == "CT+LG") ~ 1,
                                   TRUE ~ 0),
             strategyD = case_when(tx_1 == "CT" &
                                     (responder_1_complete == 1 |
                                        responder_1_complete == 0 & tx_2_complete == "CT") ~ 1,
                                   TRUE ~ 0)) %>%
      filter(strategyA == 1 | strategyC == 1) %>%
      mutate(strategy = case_when(strategyA == 1 ~ "A",
                                  strategyC == 1 ~ "C"))
    
    # Design C
    data_synth_C_StratAvStratC <- data_synth_ABC[[3]] %>%
      mutate(strategyA = case_when(tx_1 == "CT+LG" & 
                                     (responder_1_complete == 1 | 
                                        responder_1_complete == 0 & tx_2_complete == "CT+MI") ~ 1,
                                   TRUE ~ 0),
             strategyB = case_when(tx_1 == "CT+LG" & 
                                     (responder_1_complete == 1 | 
                                        responder_1_complete == 0 & tx_2_complete == "CT+LG") ~ 1,
                                   TRUE ~ 0),
             strategyC = case_when(tx_1 == "CT" &
                                     (responder_1_complete == 1 |
                                        responder_1_complete == 0 & tx_2_complete == "CT+LG") ~ 1,
                                   TRUE ~ 0),
             strategyD = case_when(tx_1 == "CT" &
                                     (responder_1_complete == 1 |
                                        responder_1_complete == 0 & tx_2_complete == "CT") ~ 1,
                                   TRUE ~ 0)) %>%
      filter(strategyA == 1 | strategyC == 1) %>%
      mutate(strategy = case_when(strategyA == 1 ~ "A",
                                  strategyC == 1 ~ "C"))
    
    # Perform t-test and record p-value (using OBSERVED data - missing values are removed)
    
    # It is technically possible that there is no one who received 
    # Strategy A or C (has yet to happen in our sims)
    # To prevent this from crashing the code, include a fail-safe catch
    # for this scenario
    
    # Design A
    if (length(unique(data_synth_A_StratAvStratC$strategy)) == 2) {
      
      pvals_currentn_A[j] <- t.test(
        hads_2 ~ strategy,
        data = data_synth_A_StratAvStratC
      )$p.value
      
    } else {
      
      pvals_currentn_A[j] <- NA
      
    }
    
    # Design B
    if (length(unique(data_synth_B_StratAvStratC$strategy)) == 2) {
      
      pvals_currentn_B[j] <- t.test(
        hads_2 ~ strategy,
        data = data_synth_B_StratAvStratC
      )$p.value
      
    } else {
      
      pvals_currentn_B[j] <- NA
      
    }
    
    # Design C
    if (length(unique(data_synth_C_StratAvStratC$strategy)) == 2) {
      
      pvals_currentn_C[j] <- t.test(
        hads_2 ~ strategy,
        data = data_synth_C_StratAvStratC
      )$p.value
      
    } else {
      
      pvals_currentn_C[j] <- NA
      
    }
    
  }
  
  # Estimate power as proportion of null hypothesis rejections out of number of simulation runs
  # i.e., count(p-vals < 0.05)/n_sim
  power_est_A <- mean(pvals_currentn_A < 0.05, na.rm = TRUE)
  power_est_B <- mean(pvals_currentn_B < 0.05, na.rm = TRUE)
  power_est_C <- mean(pvals_currentn_C < 0.05, na.rm = TRUE)
  
  # Store the number of "valid" simulations (i.e., t-test was successfully performed)
  n_valid_A <- sum(!is.na(pvals_currentn_A))
  n_valid_B <- sum(!is.na(pvals_currentn_B))
  n_valid_C <- sum(!is.na(pvals_currentn_C))
  
  # Record n, power, number of simulations that did not fail
  c(n = n_obs_current, power_A = power_est_A, n_sim_A = n_valid_A, 
    power_B = power_est_B, n_sim_B = n_valid_B, power_C = power_est_C, n_sim_C = n_valid_C)
}

stopCluster(cl)

end_time <- Sys.time()
(time_taken <- end_time - start_time)

power_persamplesize_ABC <- as.data.frame(sim.power.allN.alldesigns.StratAvStratC.parallel)
rownames(power_persamplesize_ABC) <- NULL

