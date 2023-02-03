#this script is to interpolate the exponential GI and try to use it in adaptive trials
#try to extract measures of trial design performance: 
#1) proportion of patients in superior arm "prop_superior" 
#2) power: 1 - beta or 1 - P(false negative)
#3) expected total outcome
#4) s.d. of estimate as an accuracy measure
#5) Bias

###########################
#FILE & Interpolation
###########################

#import the table copied from Gittins et al 2011 (Table 8.13 trimmed)
GI_table <- readxl::read_xlsx('GI Exponential.xlsx')
#make a new empty data frame to fill data for each n
GI_interpolated <- data.frame(2:900, NA, NA, NA, NA)
colnames(GI_interpolated) <- c('n', '0.8', '0.9', '0.95', '0.99')
#if the n is recorded in the original table, copy it directly
GI_interpolated[GI_interpolated$n %in% GI_table$n, -1] <- GI_table[,-1]
#interpolate:
for(row in 1:nrow(GI_interpolated)){
  #for readability we have a short-hand n
  n <- GI_interpolated$n[row]
  if(!n %in% GI_table$n){ # if this n needs interpolating
    #we retrieve the available n before it and after it 
    n_priv <- max(GI_table$n[GI_table$n < n])
    n_next <- min(GI_table$n[GI_table$n > n])
    #get the step assuming even intervals
    step <- (GI_table[GI_table$n == n_next,-1] - GI_table[GI_table$n == n_priv,-1]) / (n_next - n_priv)
    #interpolate at the position of n
    GI_interpolated[GI_interpolated$n == n,-1] <- GI_table[GI_table$n == n_priv,-1] + step * (n - n_priv)
  }
  #print for progress
  print(paste(row, 'out of', nrow(GI_interpolated), 'interpolated'))
}


###########################
#SIMULATION
###########################

#outcome function, output a random reward following exponential with given mean
simulate_exponential_outcome <- function(mean){
  outcome <- rexp(1, rate = 1/mean)
  outcome
}

#a function for the GI allocation rule, 
#with constraint that ensures at least 1/k patient get allocated to each arm
#k should be < N to make sure that at least some patients are in the control, otherwise index gives treatment to everyone
simulate_trials_GI_constraint <- function(N, mu_0, mu_1, R, implicit_prior_n = 2, k = 5){
  #create a dataframe to store each trial simulation's statistics
  trial.outcomes <- data.frame(Prop.Superior = rep(NA, R), Test.Stat = rep(NA, R), Estimate = rep(NA, R), Total.Outcome = rep(NA, R))
  #set gittin index parameters
  discount_factor <- '0.99' #pre-set discount factor, can change later
  #in a for loop, simulate each trials and populate the frame row by row
  for(i in 1:R){
    
    #simulate one trial that sequentially allocate patients to the two arms at random with fixed equal probability
    patient.outcome <- data.frame(Arm = rep(NA, N), Outcome = rep(NA, N))
    #set initial index values of each arm; this is updated during the trial and decides on the allocation of the next patient
    GIs <- c(mu_0 * GI_interpolated[[discount_factor]][GI_interpolated$n == implicit_prior_n],
             mu_1 * GI_interpolated[[discount_factor]][GI_interpolated$n == implicit_prior_n])
    
    for(t in 1:N){ #take in patients sequentially
      
      #first, calculate the ratio of existing obs in each arms
      arm_0_ratio <- nrow(patient.outcome[!is.na(patient.outcome$Arm) & patient.outcome$Arm == 0,]) / 
        nrow(patient.outcome[!is.na(patient.outcome$Outcome),])
      arm_1_ratio <- 1 - arm_0_ratio
      
      #patient allocation following GI rules; ensuring at least 1/k of patients are in each arm; tie-break at random
      patient.outcome$Arm[t] <- 
        ifelse(nrow(patient.outcome[!is.na(patient.outcome$Outcome),]) == 0, #if there are no obs yet
               ifelse(GIs[1] != GIs[2], which.max(GIs) - 1, ifelse(runif(1) >= 0.5, 1, 0)), #assign by prior GI, tiebreak at random
               #otherwise, if there are observations,
               ifelse(arm_0_ratio < 1/k | arm_1_ratio < 1/k, #and if less than 1/k of existing obs are in each arm
                      which.min(c(arm_0_ratio,arm_1_ratio))-1, #assign to the lower one to get it up to 1/k
                      #otherwise, assign by prior GI, tiebreak at random
                      ifelse(GIs[1] != GIs[2], which.max(GIs) - 1, ifelse(runif(1) >= 0.5, 1, 0))))
      
      #simulate outcomes based on arms
      patient.outcome$Outcome[t] <- ifelse(patient.outcome$Arm[t] == 0, 
                                           simulate_exponential_outcome(mu_0), 
                                           simulate_exponential_outcome(mu_1)) 
      #update the GIs based on outcomes
      posterior_obs <- patient.outcome$Outcome[patient.outcome$Arm == patient.outcome$Arm[t] & !is.na(patient.outcome$Outcome)] #list of all obs in this arm so far
      GIs[patient.outcome$Arm[t] + 1] <- 
        (mu_0 * implicit_prior_n + sum(posterior_obs)) / (implicit_prior_n + length(posterior_obs)) * #posterior mean
        GI_interpolated[[discount_factor]][GI_interpolated$n == implicit_prior_n + length(posterior_obs)] #table value of v(∑,n) where lambda = 1
    }
    
    #extract trial statistics from this outcome dataframe
    trial.outcomes$Prop.Superior[i] <- nrow(patient.outcome[patient.outcome$Arm == which.max(c(mu_0,mu_1))-1,]) / N
    #test statistics here is the p-value generated by a F distribution F(2n_1, 2n_2)
    trial.outcomes$Test.Stat[i] <- pf(mean(patient.outcome$Outcome[patient.outcome$Arm == 0])/mean(patient.outcome$Outcome[patient.outcome$Arm == 1]),
                                      df1 = 2*nrow(patient.outcome[patient.outcome$Arm == 0,]), df2 = 2*nrow(patient.outcome[patient.outcome$Arm == 1,]))
    trial.outcomes$Estimate[i] <- mean(patient.outcome$Outcome[patient.outcome$Arm == 1])-mean(patient.outcome$Outcome[patient.outcome$Arm == 0])
    trial.outcomes$Total.Outcome[i] <- sum(patient.outcome$Outcome)
    #print for progress
    print(paste(i, 'out of', R, 'trials simulated'))
  }
  trial.outcomes #output
}


#run simulation for when null is true and when alt is true
mu_0 <- 0.5
mu_1 <- 0.2
N <- 100
null.true <- simulate_trials_GI_constraint(N, mu_0, mu_0, R = 100, k = 5)
alte.true <- simulate_trials_GI_constraint(N, mu_0, mu_1, R = 100, k = 5)

#compute trial properties
#since the stats here is p-value, we set cutoff to 0.05, and reverse the inequalities 
#again, I am not sure if this is the right test for exponential endpoints
cutoff <- 0.05
alpha <- ifelse(mu_1>mu_0, #split in cases for one-tailed tests
                nrow(null.true[null.true$Test.Stat < cutoff,]) / nrow(null.true),
                nrow(null.true[null.true$Test.Stat > 1-cutoff,]) / nrow(null.true))
power <- ifelse(mu_1>mu_0,
                1 - nrow(alte.true[alte.true$Test.Stat > cutoff,]) / nrow(alte.true),
                1 - nrow(alte.true[alte.true$Test.Stat < 1-cutoff,]) / nrow(alte.true))
prop_superior_in_null <- mean(null.true$Prop.Superior)
prop_sup_in_null_sd <- sd(null.true$Prop.Superior)
prop_superior_in_alte <- mean(alte.true$Prop.Superior)
prop_sup_in_alte_sd <- sd(alte.true$Prop.Superior)
est_mean_in_null <- mean(null.true$Estimate)
est_mean_in_alte <- mean(alte.true$Estimate)
est_sd_in_null <- sd(null.true$Estimate)
est_sd_in_alte <- sd(alte.true$Estimate)
#added later (August 22nd, not reported below)
ETO_change_percentage_null <- mean((null.true$Total.Outcome - mu_0*N)/(mu_0*N)*100)
ETO_change_percentage_alte <- mean((alte.true$Total.Outcome - (mu_0*N/2+mu_1*N/2))/(mu_0*N/2+mu_1*N/2)*100)

#true effect size
true_effect <- mu_1 - mu_0

#report results
cat(paste(nrow(null.true), ' trials were simulated; true effect size = ', round(true_effect, 3),
          '\nAt a test statistics cut-off of ', cutoff, ', Type I error = ', alpha, '; Power = ', round(power,4), 
          '\nH0 == T: P(superior) = ', round(prop_superior_in_null, 4), ' (SD = ', round(prop_sup_in_null_sd, 2), 
          '), Bias = ', round(bias_in_null, 5), ' (SD = ', round(bias_in_null_sd, 2),
          '), Mean Estimate = ', round(est_mean_in_null, 3), ' (SD = ', round(est_sd_in_null, 2), ');',
          '\nH1 == T: P(superior) = ', round(prop_superior_in_alte, 4), ' (SD = ', round(prop_sup_in_alte_sd, 2), 
          '), Bias = ', round(bias_in_alte, 5), ' (SD = ', round(bias_in_alte_sd, 2), 
          '), Mean Estimate = ', round(est_mean_in_alte, 3), ' (SD = ', round(est_sd_in_alte, 2), ').', sep = ''))

