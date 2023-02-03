#this script is an attempt at replicating Williamson et al 2019 FLGI-known simulations
#refer to "1. Normal ER.R", the present script takes the same structure only with the allocation rule changed

##################################################
#Simulate 2-armed trial in Williamson et al., 2019
##################################################

#Karrison et al. 2017 Phase II cancer trial
#two-armed, T = N = 72, mu_0 = 0.155, mu_1 = 0.529, true_sd = 0.64, trial replication R = 50,000
#limit to block size = 1, discount factor = 0.995

#the gittins index (known var) values
GI_table <- readxl::read_xlsx('GI Normal (Known Var).xlsx')

#outcome function, output a random reward following normal with given parameters
simulate_normal_outcome <- function(mean, sd){
  outcome <- rnorm(1, mean = mean, sd = sd)
  outcome
}

#trials function, output a trial outcome dataframe containing statistics of each simulated trials
simulate_2_armed_trials <- function(N, mu_0, mu_1, true_sd, R, implicit_prior_n){
  #create a dataframe to store each trial simulation's statistics
  trial.outcomes <- data.frame(Prop.Superior = rep(NA, R), Test.Stat = rep(NA, R), Bias = rep(NA, R))
  #set gittin index parameters
  discount_factor <- '0.995' #pre-set discount factor, can change later
  #in a for loop, simulate each trials and populate the frame row by row
  for(i in 1:R){
    #simulate one trial that sequentially allocate patients to the two arms at random with fixed equal probability
    patient.outcome <- data.frame(Arm = rep(NA, N), Outcome = rep(NA, N))
    #set initial index values of each arm; this is updated during the trial and decides on the allocation of the next patient
    GIs <- c(mu_0 + true_sd * GI_table[[discount_factor]][GI_table$n == signif(implicit_prior_n, 1)] / (implicit_prior_n * ((1 - as.double(discount_factor)) ** 0.5)), 
             mu_0 + true_sd * GI_table[[discount_factor]][GI_table$n == signif(implicit_prior_n, 1)] / (implicit_prior_n * ((1 - as.double(discount_factor)) ** 0.5)))
    for(t in 1:N){ #take in patients sequentially
      #patient allocation following GI-known rules; tie-break at random
      patient.outcome$Arm[t] <- ifelse(GIs[1] != GIs[2], which.max(GIs) - 1, ifelse(runif(1) >= 0.5, 1, 0))
      #simulate outcomes based on arms
      patient.outcome$Outcome[t] <- ifelse(patient.outcome$Arm[t] == 0, 
                                           simulate_normal_outcome(mean = mu_0, sd = true_sd), 
                                           simulate_normal_outcome(mean = mu_1, sd = true_sd))
      #update the GIs based on outcomes
      posterior_obs <- patient.outcome$Outcome[patient.outcome$Arm == patient.outcome$Arm[t] & !is.na(patient.outcome$Outcome)] #list of all obs in this arm so far
      GIs[patient.outcome$Arm[t] + 1] <- mean(posterior_obs) + 
        true_sd * GI_table[[discount_factor]][GI_table$n == signif(length(posterior_obs) + implicit_prior_n, 1)] / 
        ((length(posterior_obs) + implicit_prior_n) * ((1 - as.double(discount_factor)) ** 0.5))
    }
    #extract trial statistics from this outcome dataframe
    trial.outcomes$Prop.Superior[i] <- nrow(patient.outcome[patient.outcome$Arm == 1,]) / N
    trial.outcomes$Test.Stat[i] <- 
      (mean(patient.outcome$Outcome[patient.outcome$Arm == 1]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])) /
      sqrt(true_sd^2/nrow(patient.outcome[patient.outcome$Arm == 1,]) + true_sd^2/nrow(patient.outcome[patient.outcome$Arm == 0,]))
    trial.outcomes$Bias[i] <- (mean(patient.outcome$Outcome[patient.outcome$Arm == 1]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])) - (mu_1 - mu_0)
    #print for progress
    print(paste(i, 'out of', R, 'trials simulated'))
  }
  trial.outcomes #output
}

#run simulation for when null is true and when alt is true
null.true <- simulate_2_armed_trials(N = 72, mu_0 = 0.155, mu_1 = 0.155, true_sd = 0.64, R = 10000, implicit_prior_n = 1)
alte.true <- simulate_2_armed_trials(N = 72, mu_0 = 0.155, mu_1 = 0.529, true_sd = 0.64, R = 10000, implicit_prior_n = 1)

#compute trial properties
cutoff <- 1.991
alpha <- nrow(null.true[null.true$Test.Stat > cutoff,]) / nrow(null.true)
power <- 1 - nrow(alte.true[alte.true$Test.Stat < cutoff,]) / nrow(alte.true)
prop_superior_in_null <- mean(null.true$Prop.Superior)
prop_sup_in_null_sd <- sd(null.true$Prop.Superior)
prop_superior_in_alte <- mean(alte.true$Prop.Superior)
prop_sup_in_alte_sd <- sd(alte.true$Prop.Superior)
bias_in_null <- mean(null.true$Bias)
bias_in_null_sd <- sd(null.true$Bias)
bias_in_alte <- mean(alte.true$Bias)
bias_in_alte_sd <- sd(alte.true$Bias)
test_stat_mean_in_null <- mean(null.true$Test.Stat)
test_stat_mean_in_alte <- mean(alte.true$Test.Stat)
test_stat_sd_in_null <- sd(null.true$Test.Stat)
test_stat_sd_in_alte <- sd(alte.true$Test.Stat)

#report results
cat(paste(nrow(null.true), ' trials were simulated.\nAt a test statistics cut-off of ', cutoff, 
          ', Type I error = ', alpha, '; Power = ', power, 
          '\nH0 == T: P(superior) = ', round(prop_superior_in_null, 4), ' (SD = ', round(prop_sup_in_null_sd, 2), 
          '), Bias = ', round(bias_in_null, 5), ' (SD = ', round(bias_in_null_sd, 2),
          '), Mean Test Statistics = ', round(test_stat_mean_in_null, 3), ' (SD = ', round(test_stat_sd_in_null, 2), ');',
          '\nH1 == T: P(superior) = ', round(prop_superior_in_alte, 4), ' (SD = ', round(prop_sup_in_alte_sd, 2), 
          '), Bias = ', round(bias_in_alte, 5), ' (SD = ', round(bias_in_alte_sd, 2), 
          '), Mean Test Statistics = ', round(test_stat_mean_in_alte, 3), ' (SD = ', round(test_stat_sd_in_alte, 2), ').', sep = ''))


