#initial attempt at simulating simple trials with exponential endpoints 
#and evaluate the performance of equal randomization designs (F statistics)

#try to extract measures of trial design performance: 
#1) proportion of patients in superior arm "prop_superior" 
#2) power: 1 - beta or 1 - P(false negative)
#3) expected total outcome
#4) s.d. of estimate as an accuracy measure
#5) Bias

##################################################
#Simulate 2-armed trial with Exp endpoints
##################################################

#Arbitrarily deciding:
#T = N = 100, mu_0 = 0.5, mu_1 = 0.7

#outcome function, output a random reward following exponential with given mean
simulate_exponential_outcome <- function(mean){
  outcome <- rexp(1, rate = 1/mean)
  outcome
}

#trials function, output a trial outcome dataframe containing statistics of each simulated trials
simulate_trials_ER <- function(N, mu_0, mu_1, R){
  #create a dataframe to store each trial simulation's statistics
  trial.outcomes <- data.frame(Prop.Superior = rep(NA, R), Test.Stat = rep(NA, R), Bias = rep(NA, R), Estimate = rep(NA, R))
  #in a for loop, simulate each trials and populate the frame row by row
  for(i in 1:R){
    #simulate one trial that sequentially allocate patients to the two arms at random with fixed equal probability
    patient.outcome <- data.frame(Arm = rep(NA, N), Outcome = rep(NA, N))
    for(t in 1:N){ #take in patients sequentially
      patient.outcome$Arm[t] <- ifelse(runif(1) >= 0.5, 1, 0) #allocation to arms at random - this should be the thing that changes when we run other trial designs
      patient.outcome$Outcome[t] <- ifelse(patient.outcome$Arm[t] == 0, 
                                           simulate_exponential_outcome(mu_0), 
                                           simulate_exponential_outcome(mu_1)) #simulate outcomes based on arms
    }
    #extract trial statistics from this outcome dataframe
    trial.outcomes$Prop.Superior[i] <- nrow(patient.outcome[patient.outcome$Arm == 1,]) / N
    #test statistics here is the p-value generated by a F distribution F(2n_1, 2n_2)
    trial.outcomes$Test.Stat[i] <- pf(mean(patient.outcome$Outcome[patient.outcome$Arm == 0])/mean(patient.outcome$Outcome[patient.outcome$Arm == 1]),
                                      df1 = 2*nrow(patient.outcome[patient.outcome$Arm == 0,]), df2 = 2*nrow(patient.outcome[patient.outcome$Arm == 1,]))
    trial.outcomes$Bias[i] <- (mean(patient.outcome$Outcome[patient.outcome$Arm == 1]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])) - (mu_1 - mu_0)
    #estimate here is the unbiased ratio of treatment mean / control mean
    trial.outcomes$Estimate[i] <- mean(patient.outcome$Outcome[patient.outcome$Arm == 1]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])
    #print for progress
    print(paste(i, 'out of', R, 'trials simulated'))
  }
  trial.outcomes #output
}

#run simulation for when null is true and when alt is true
mu_0 <- 0.5
mu_1 <- 0.85
null.true <- simulate_trials_ER(N = 100, mu_0, mu_0, R = 1000)
alte.true <- simulate_trials_ER(N = 100, mu_0, mu_1, R = 1000)

#compute trial properties
#since the stats here is p-value, we set cutoff to 0.05, and reverse the inequalities 
#again, I am not sure if this is the right test for exponential endpoints
cutoff <- 0.05
alpha <- nrow(null.true[null.true$Test.Stat < cutoff,]) / nrow(null.true)
power <- 1 - nrow(alte.true[alte.true$Test.Stat > cutoff,]) / nrow(alte.true)
prop_superior_in_null <- mean(null.true$Prop.Superior)
prop_sup_in_null_sd <- sd(null.true$Prop.Superior)
prop_superior_in_alte <- mean(alte.true$Prop.Superior)
prop_sup_in_alte_sd <- sd(alte.true$Prop.Superior)
bias_in_null <- mean(null.true$Bias)
bias_in_null_sd <- sd(null.true$Bias)
bias_in_alte <- mean(alte.true$Bias)
bias_in_alte_sd <- sd(alte.true$Bias)
est_mean_in_null <- mean(null.true$Estimate)
est_mean_in_alte <- mean(alte.true$Estimate)
est_sd_in_null <- sd(null.true$Estimate)
est_sd_in_alte <- sd(alte.true$Estimate)

#true effect size
true_effect <- mu_1 - mu_0

#report results
cat(paste(nrow(null.true), ' trials were simulated; true effect size = ', round(true_effect, 3),
          '\nAt a test statistics cut-off of ', cutoff, ', Type I error = ', alpha, '; Power = ', power, 
          '\nH0 == T: P(superior) = ', round(prop_superior_in_null, 4), ' (SD = ', round(prop_sup_in_null_sd, 2), 
          '), Bias = ', round(bias_in_null, 5), ' (SD = ', round(bias_in_null_sd, 2),
          '), Mean Estimate = ', round(est_mean_in_null, 3), ' (SD = ', round(est_sd_in_null, 2), ');',
          '\nH1 == T: P(superior) = ', round(prop_superior_in_alte, 4), ' (SD = ', round(prop_sup_in_alte_sd, 2), 
          '), Bias = ', round(bias_in_alte, 5), ' (SD = ', round(bias_in_alte_sd, 2), 
          '), Mean Estimate = ', round(est_mean_in_alte, 3), ' (SD = ', round(est_sd_in_alte, 2), ').', sep = ''))

