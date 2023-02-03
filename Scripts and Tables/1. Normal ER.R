#initial attempt at simulating simple trials with normal endpoints 
#and evaluate the performance of equal randomization designs

#plan: have two-armed first, then extend to three-arms

#try to extract measures of trial design performance: 
#1) proportion of patients in superior arm "prop_superior" 
#2) power: 1 - beta or 1 - P(false negative)
#3) alpha or type I error rate (false positives)
#4) bias
#5) expected total outcome change???

##################################################
#Simulate 2-armed trial in Williamson et al., 2019
##################################################

#Karrison et al. 2017 Phase II cancer trial
#two-armed, T = N = 72, mu_0 = 0.155, mu_1 = 0.529, true_sd = 0.64, trial replication R = 50,000

#outcome function, output a random reward following normal with given parameters
simulate_normal_outcome <- function(mean, sd){
  outcome <- rnorm(1, mean = mean, sd = sd)
  outcome
}

#trials function, output a trial outcome dataframe containing statistics of each simulated trials
simulate_2_armed_trials <- function(N, mu_0, mu_1, true_sd, R){
  #create a dataframe to store each trial simulation's statistics
  trial.outcomes <- data.frame(Prop.Superior = rep(NA, R), Test.Stat = rep(NA, R), Bias = rep(NA, R))
  #in a for loop, simulate each trials and populate the frame row by row
  for(i in 1:R){
    #simulate one trial that sequentially allocate patients to the two arms at random with fixed equal probability
    patient.outcome <- data.frame(Arm = rep(NA, N), Outcome = rep(NA, N))
    for(t in 1:N){ #take in patients sequentially
      patient.outcome$Arm[t] <- ifelse(runif(1) >= 0.5, 1, 0) #allocation to arms at random - this should be the thing that changes when we run other trial designs
      patient.outcome$Outcome[t] <- ifelse(patient.outcome$Arm[t] == 0, 
                                           simulate_normal_outcome(mean = mu_0, sd = true_sd), 
                                           simulate_normal_outcome(mean = mu_1, sd = true_sd)) #simulate outcomes based on arms
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
null.true <- simulate_2_armed_trials(N = 72, mu_0 = 0.155, mu_1 = 0.155, true_sd = 0.64, R = 1000)
alte.true <- simulate_2_armed_trials(N = 72, mu_0 = 0.155, mu_1 = 0.529, true_sd = 0.64, R = 1000)

#compute trial properties
cutoff <- 1.654
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

#report results
cat(paste(nrow(null.true), ' trials were simulated.\nAt a test statistics cut-off of ', cutoff, 
          ', this design has a Type I error rate of ', alpha, ' and a power of ', power, 
          '\nWhen H0 is true, the design assigned ', round(prop_superior_in_null, 4), ' (SD = ', round(prop_sup_in_null_sd, 2), 
          ') of patients to the superior arm, Bias = ', round(bias_in_null, 5), ' (SD = ', round(bias_in_null_sd, 2), ');',
          '\nWhen H1 is true, the design assigned ', round(prop_superior_in_alte, 4), ' (SD = ', round(prop_sup_in_alte_sd, 2), 
          ') of patients to the superior arm, Bias = ', round(bias_in_alte, 5), ' (SD = ', round(bias_in_alte_sd, 2), ').', sep = ''))


##################################################
#Simulate 3-armed trial in Williamson et al., 2019
##################################################

#3-armed, T = N = 120, mu_0 = -0.05, mu_1 = 0.07, mu_2 = 0.13, true_sd = 0.346, R = 50,000

#the "simulate_normal_outcome" function from the previous section is used here as well.
simulate_normal_outcome <- function(mean, sd){
  outcome <- rnorm(1, mean = mean, sd = sd)
  outcome
}
#trials function, output a trial outcome dataframe containing statistics of each simulated trials
simulate_3_armed_trials <- function(N, mu_0, mu_1, mu_2, true_sd, R){
  #create a dataframe to store each trial simulation's statistics
  trial.outcomes <- data.frame(Prop.Superior = rep(NA, R), Test.Stat.1 = rep(NA, R), Test.Stat.2 = rep(NA, R), Bias.1 = rep(NA, R), Bias.2 = rep(NA, R))
  #in a for loop, simulate each trials and populate the frame row by row
  for(i in 1:R){
    #simulate one trial that sequentially allocate patients to the three arms at random with fixed equal probability
    patient.outcome <- data.frame(Arm = rep(NA, N), Outcome = rep(NA, N))
    for(t in 1:N){ #take in patients sequentially
      patient.outcome$Arm[t] <- ifelse(runif(1) <= 1/3, 2, ifelse(runif(1) <= 0.5, 1, 0)) #allocation to arms at random
      #simulate outcomes based on arms
      if(patient.outcome$Arm[t] == 0){
        patient.outcome$Outcome[t] <- simulate_normal_outcome(mean = mu_0, sd = true_sd)
      } else if(patient.outcome$Arm[t] == 1){
        patient.outcome$Outcome[t] <- simulate_normal_outcome(mean = mu_1, sd = true_sd)
      } else{
        patient.outcome$Outcome[t] <- simulate_normal_outcome(mean = mu_2, sd = true_sd)
      }
    }
    #extract trial statistics from this outcome dataframe
    trial.outcomes$Prop.Superior[i] <- nrow(patient.outcome[patient.outcome$Arm == 2,]) / N
    trial.outcomes$Test.Stat.1[i] <- 
      (mean(patient.outcome$Outcome[patient.outcome$Arm == 1]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])) /
      sqrt(true_sd^2/nrow(patient.outcome[patient.outcome$Arm == 1,]) + true_sd^2/nrow(patient.outcome[patient.outcome$Arm == 0,]))
    trial.outcomes$Test.Stat.2[i] <- 
      (mean(patient.outcome$Outcome[patient.outcome$Arm == 2]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])) /
      sqrt(true_sd^2/nrow(patient.outcome[patient.outcome$Arm == 2,]) + true_sd^2/nrow(patient.outcome[patient.outcome$Arm == 0,]))
    trial.outcomes$Bias.1[i] <- (mean(patient.outcome$Outcome[patient.outcome$Arm == 1]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])) - (mu_1 - mu_0)
    trial.outcomes$Bias.2[i] <- (mean(patient.outcome$Outcome[patient.outcome$Arm == 2]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])) - (mu_2 - mu_0)
    #print for progress
    print(paste(i, 'out of', R, 'trials simulated'))
  }
  trial.outcomes #output
}

#run simulation for when null is true and when alt is true
null.true <- simulate_3_armed_trials(N = 120, mu_0 = -0.05, mu_1 = -0.05, mu_2 = -0.05, true_sd = 0.346, R = 50000)
alte.true <- simulate_3_armed_trials(N = 120, mu_0 = -0.05, mu_1 = 0.07, mu_2 = 0.13, true_sd = 0.346, R = 50000)

#compute trial properties
cutoff <- 1.595
alpha <- nrow(null.true[null.true$Test.Stat.1 > cutoff | null.true$Test.Stat.2 > cutoff,]) / nrow(null.true)
power <- 1 - nrow(alte.true[alte.true$Test.Stat.1 < cutoff & alte.true$Test.Stat.2 < cutoff,]) / nrow(alte.true)
prop_superior_in_null <- mean(null.true$Prop.Superior)
prop_sup_in_null_sd <- sd(null.true$Prop.Superior)
prop_superior_in_alte <- mean(alte.true$Prop.Superior)
prop_sup_in_alte_sd <- sd(alte.true$Prop.Superior)
#not sure if it is right to report pair-wise bias rather than family-wise bias
#bias_in_null <- mean(null.true$Bias)
#bias_in_null_sd <- sd(null.true$Bias)
#bias_in_alte <- mean(alte.true$Bias)
#bias_in_alte_sd <- sd(alte.true$Bias)

#report results
cat(paste(nrow(null.true), ' trials were simulated.\nAt a test statistics cut-off of ', cutoff, 
          ', this design has a Type I error rate of ', alpha, ' and a power of ', power, 
          '\nWhen H0 is true, the design assigned ', round(prop_superior_in_null, 4), ' (SD = ', round(prop_sup_in_null_sd, 2), 
          ') of patients to the superior arm', #'Bias = ', round(bias_in_null, 5), ' (SD = ', round(bias_in_null_sd, 2), ');',
          '\nWhen H1 is true, the design assigned ', round(prop_superior_in_alte, 4), ' (SD = ', round(prop_sup_in_alte_sd, 2), 
          ') of patients to the superior arm', #'Bias = ', round(bias_in_alte, 5), ' (SD = ', round(bias_in_alte_sd, 2), ').',
          sep = ''))




