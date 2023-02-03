#this script combines the previously written ER and GI simulations of 2-armed trials
#and simulate the two designs in different trial conditions to produce visualisations

##########################
#FILES, LIBS, and PREVIOUS FUNCTIONS
#from previous scripts attempts
##########################

#the gittins index (known var) values
GI_table <- readxl::read_xlsx('GI Normal (Known Var).xlsx')

#outcome function, output a random reward following normal with given parameters
simulate_normal_outcome <- function(mean, sd){
  outcome <- rnorm(1, mean = mean, sd = sd)
  outcome
}

#trials function, output a 2-armed trial outcome dataframe containing statistics of each simulated trials, under ER allocation rules
simulate_trials_ER <- function(N, mu_0, mu_1, true_sd, R){
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
    trial.outcomes$Prop.Superior[i] <- nrow(patient.outcome[patient.outcome$Arm == which.max(c(mu_0, mu_1)) - 1,]) / N
    trial.outcomes$Test.Stat[i] <- 
      (mean(patient.outcome$Outcome[patient.outcome$Arm == 1]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])) /
      sqrt(true_sd^2/nrow(patient.outcome[patient.outcome$Arm == 1,]) + true_sd^2/nrow(patient.outcome[patient.outcome$Arm == 0,]))
    trial.outcomes$Bias[i] <- (mean(patient.outcome$Outcome[patient.outcome$Arm == 1]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])) - (mu_1 - mu_0)
    #print for progress
    print(paste(i, 'out of', R, 'trials simulated'))
  }
  trial.outcomes #output
}

#trials function, output a 2-armed trial outcome dataframe containing statistics of each simulated trials, under GI allocation rules
simulate_trials_GI <- function(N, mu_0, mu_1, true_sd, R, implicit_prior_n){
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
      #patient allocation following GI-known rules; tie-break in favor of mu_0 (otherwise we have cases in large eff.size where all are assigned to the superior arm)
      patient.outcome$Arm[t] <- ifelse(GIs[1] != GIs[2], which.max(GIs) - 1, 0)
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
    trial.outcomes$Prop.Superior[i] <- nrow(patient.outcome[patient.outcome$Arm == which.max(c(mu_0, mu_1)) - 1,]) / N
    trial.outcomes$Test.Stat[i] <- 
      (mean(patient.outcome$Outcome[patient.outcome$Arm == 1]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])) /
      sqrt(true_sd^2/nrow(patient.outcome[patient.outcome$Arm == 1,]) + true_sd^2/nrow(patient.outcome[patient.outcome$Arm == 0,]))
    trial.outcomes$Bias[i] <- (mean(patient.outcome$Outcome[patient.outcome$Arm == 1]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])) - (mu_1 - mu_0)
    #print for progress
    print(paste(i, 'out of', R, 'trials simulated'))
  }
  trial.outcomes #output
}


##########################
#SIMULATE & COMPARE
#in different trial conditions
##########################

#mimicking figure 3 in Mavrogonatou et al, 2022, see power curve in different effect sizes
#set arbiturary mu_0 = 10, true_sd = 5.5, vary mu_1_lower = 0 and mu_1_upper = 20, mu_1_step = 1, implicit_prior_n = 1, N = 120, R = 100

#a function that outputs a data frame of x = mu_1 and y = power, accuracy, prop superior; for ER trials
simulate_mu_range_ER <- function(N, mu_0, mu_1_lower, mu_1_upper, mu_1_step, true_sd, R){
  #create a data frame to store outcome
  output <- data.frame(mu_1 = seq(mu_1_lower, mu_1_upper, by = mu_1_step), Alpha = NA, Power = NA, H1.Test.Stat.SD = NA, H1.Prop.Superior = NA)
  #go through each row to simulate
  for(row in 1:nrow(output)){
    #simulate by GI
    null.true <- simulate_trials_ER(N, mu_0, mu_0, true_sd, R)
    alte.true <- simulate_trials_ER(N, mu_0, mu_1 = output$mu_1[row], true_sd, R)
    #extract trial properties
    output$Alpha[row] <- nrow(null.true[abs(null.true$Test.Stat) > qnorm(0.95),]) / nrow(null.true)
    output$Power[row] <- 1 - nrow(alte.true[abs(alte.true$Test.Stat) < qnorm(0.95),]) / nrow(alte.true)
    output$H1.Test.Stat.SD[row] <- sd(alte.true$Test.Stat)
    output$H1.Prop.Superior[row] <- mean(alte.true$Prop.Superior)
  }
  output #output
}

#a function that outputs a data frame of x = mu_1 and y = power, accuracy, prop superior; for GI trials
simulate_mu_range_GI <- function(N, mu_0, mu_1_lower, mu_1_upper, mu_1_step, true_sd, R, implicit_prior_n){
  #create a data frame to store outcome
  output <- data.frame(mu_1 = seq(mu_1_lower, mu_1_upper, by = mu_1_step), Alpha = NA, Power = NA, H1.Test.Stat.SD = NA, H1.Prop.Superior = NA)
  #go through each row to simulate
  for(row in 1:nrow(output)){
    #simulate by GI
    null.true <- simulate_trials_GI(N, mu_0, mu_0, true_sd, R, implicit_prior_n)
    alte.true <- simulate_trials_GI(N, mu_0, mu_1 = output$mu_1[row], true_sd, R, implicit_prior_n)
    #extract trial properties
    output$Alpha[row] <- nrow(null.true[abs(null.true$Test.Stat) > qnorm(0.95),]) / nrow(null.true)
    output$Power[row] <- 1 - nrow(alte.true[abs(alte.true$Test.Stat) < qnorm(0.95),]) / nrow(alte.true)
    output$H1.Test.Stat.SD[row] <- sd(alte.true$Test.Stat)
    output$H1.Prop.Superior[row] <- mean(alte.true$Prop.Superior)
  }
  output #output
}

#get trial property outcomes
mu_range_outcome_ER <- simulate_mu_range_ER(N = 72, mu_0 = 1, mu_1_lower = 0.5, mu_1_upper = 1.5, mu_1_step = 0.01, true_sd = 0.6, R = 10000)
mu_range_outcome_GI <- simulate_mu_range_GI(N = 72, mu_0 = 1, mu_1_lower = 0.5, mu_1_upper = 1.5, mu_1_step = 0.01, true_sd = 0.6, R = 10000, implicit_prior_n = 1)

#save because it takes AGES to simulate
write.csv(mu_range_outcome_ER, 'mu_range_outcome_ER.csv')
write.csv(mu_range_outcome_GI, 'mu_range_outcome_GI.csv')


##########################
#VISUALISE COMPARISON
#and save to .pdf format
##########################

#get library and files
library(ggplot2)
library(gridExtra)
mu_range_outcome_ER <- read.csv('mu_range_outcome_ER.csv')[,-1]
mu_range_outcome_GI <- read.csv('mu_range_outcome_GI.csv')[,-1]

#first we combine the trial property outcome
mu_range_outcome_combined <- data.frame(rbind(mu_range_outcome_ER, mu_range_outcome_GI),
                                        Design = c(rep('ER', nrow(mu_range_outcome_ER)), rep('GI', nrow(mu_range_outcome_GI))))

#start plotting!
Power_Curve <- ggplot(data = mu_range_outcome_combined) +
  geom_point(aes(x = mu_1, y = Power, color = Design)) + 
  scale_color_manual(values = c('black', 'darkred')) +
  xlab(expression(paste('\U00B5'[1], ' (\U00B5'[0], ' = 1.00)'))) + ylab('') + 
  labs(title = '1 - \U03B2') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

Alpha_Curve <- ggplot(data = mu_range_outcome_combined) +
  geom_point(aes(x = mu_1, y = Alpha, color = Design)) + 
  scale_color_manual(values = c('black', 'darkred')) +
  xlab(expression(paste('\U00B5'[1], ' (\U00B5'[0], ' = 1.00)'))) + ylab('') + 
  labs(title = '\U03B1') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

Accuracy_Curve <- ggplot(data = mu_range_outcome_combined) +
  geom_point(aes(x = mu_1, y = H1.Test.Stat.SD, color = Design)) + 
  scale_color_manual(values = c('black', 'darkred')) +
  xlab(expression(paste('\U00B5'[1], ' (\U00B5'[0], ' = 1.00)'))) + ylab('') + 
  labs(title = expression(paste('\U03C3'['Test Statistic'], ' | H'[1]))) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

Prop_Sup_Curve <- ggplot(data = mu_range_outcome_combined) +
  geom_point(aes(x = mu_1, y = H1.Prop.Superior, color = Design)) + 
  scale_color_manual(values = c('black', 'darkred')) +
  xlab(expression(paste('\U00B5'[1], ' (\U00B5'[0], ' = 1.00)'))) + ylab('') + 
  labs(title = expression(paste('\U03C1'['superior'], ' | H'[1]))) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

#combine all 4 into one plot
quartz(type = 'pdf', file = 'Trial Properties Comparison (Normal GI v. ER).pdf', width = 9, height = 9)
grid.arrange(Alpha_Curve, Power_Curve, Accuracy_Curve, Prop_Sup_Curve, nrow = 2)
dev.off()

