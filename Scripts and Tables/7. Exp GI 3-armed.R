#this script is to extend script 6 to 3-armed;

##########################
#FILES, LIBS, and PREVIOUS FUNCTIONS
#from previous scripts attempts
##########################

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


#outcome function, output a random reward following exp. with given parameters
simulate_exponential_outcome <- function(mean){
  outcome <- rexp(1, rate = 1/mean)
  outcome
}

#trials function, output a trial outcome dataframe containing statistics of each simulated trials
#also include a previous message input for progress printing
simulate_trials_ER <- function(N, mu_0, mu_1, mu_2, R, previous_message = ''){
  #create a dataframe to store each trial simulation's statistics
  trial.outcomes <- data.frame(Prop.Superior = rep(NA, R), 
                               Test.Stat.1 = rep(NA, R), Test.Stat.2 = rep(NA, R), 
                               Estimate.1 = rep(NA, R), Estimate.2 = rep(NA, R), 
                               Total.Outcome = rep(NA, R))
  
  #in a for loop, simulate each trials and populate the frame row by row
  for(i in 1:R){
    
    #simulate one trial that sequentially allocate patients to the two arms at random with fixed equal probability
    patient.outcome <- data.frame(Arm = rep(NA, N), Outcome = rep(NA, N))
    for(t in 1:N){ 
      #take in patients sequentially to allocate
      patient.outcome$Arm[t] <- ifelse(runif(1) >= 2/3, #at 1/3 of the time
                                       0, #allocate to arm 0
                                       ifelse(runif(1) >= 1/2, #then at half of the 2/3 (1/3)
                                              1, 2)) #allocate to arm 1 and 2
      #simulate outcomes based on arms
      patient.outcome$Outcome[t] <- ifelse(patient.outcome$Arm[t] == 0, 
                                           simulate_exponential_outcome(mu_0), 
                                           ifelse(patient.outcome$Arm[t] == 1,
                                                  simulate_exponential_outcome(mu_1),
                                                  simulate_exponential_outcome(mu_2)))
    }
    
    #extract trial statistics from this outcome dataframe
    trial.outcomes$Prop.Superior[i] <- 
      nrow(patient.outcome[patient.outcome$Arm == which.max(c(mu_0, mu_1, mu_2)) - 1,]) / N
    #test statistics here is the p-value generated by a f test - test mu_1-0 and mu_2-0
    trial.outcomes$Test.Stat.1[i] <- pf(mean(patient.outcome$Outcome[patient.outcome$Arm == 0])/mean(patient.outcome$Outcome[patient.outcome$Arm == 1]),
                                      df1 = 2*nrow(patient.outcome[patient.outcome$Arm == 0,]), df2 = 2*nrow(patient.outcome[patient.outcome$Arm == 1,]))
    trial.outcomes$Test.Stat.2[i] <- pf(mean(patient.outcome$Outcome[patient.outcome$Arm == 0])/mean(patient.outcome$Outcome[patient.outcome$Arm == 2]),
                                        df1 = 2*nrow(patient.outcome[patient.outcome$Arm == 0,]), df2 = 2*nrow(patient.outcome[patient.outcome$Arm == 2,]))
    #record the estimates
    trial.outcomes$Estimate.1[i] <- mean(patient.outcome$Outcome[patient.outcome$Arm == 1]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])
    trial.outcomes$Estimate.2[i] <- mean(patient.outcome$Outcome[patient.outcome$Arm == 2]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])
    #and outcomes
    trial.outcomes$Total.Outcome[i] <- sum(patient.outcome$Outcome)
    #print for progress
    print(paste(previous_message, i, 'out of', R, 'trials simulated'))
  }
  
  trial.outcomes #output
}

#a function for the GI allocation rule, 
#with constraint that ensures at least 1/k patient get allocated to each arm
#k should be < N to make sure that at least some patients are in the control, otherwise index gives treatment to everyone
simulate_trials_GI_constraint <- function(N, mu_0, mu_1, mu_2, R, implicit_prior_n = 2, k = 5, previous_message = ''){
  
  #create a dataframe to store each trial simulation's statistics
  trial.outcomes <- data.frame(Prop.Superior = rep(NA, R), Test.Stat = rep(NA, R), Estimate = rep(NA, R), Total.Outcome = rep(NA, R))
  #set gittin index parameters
  discount_factor <- '0.99' #pre-set discount factor, can change later
  
  #in a for loop, simulate each trials and populate the frame row by row
  for(i in 1:R){
    
    #simulate one trial that sequentially allocate patients to the two arms at random with fixed equal probability
    patient.outcome <- data.frame(Arm = rep(NA, N), Outcome = rep(NA, N))
    #set initial index values of each arm; this is updated during the trial and decides on the allocation of the next patient
    #here we test the performance of the normal GI, approximating mu_0 as the true sd
    GIs <- c(mu_0 * GI_interpolated[[discount_factor]][GI_interpolated$n == implicit_prior_n],
             mu_0 * GI_interpolated[[discount_factor]][GI_interpolated$n == implicit_prior_n],
             mu_0 * GI_interpolated[[discount_factor]][GI_interpolated$n == implicit_prior_n])    
    
    for(t in 1:N){ #take in patients sequentially
      
      #first, calculate the ratio of existing obs in each arms
      arm_0_ratio <- nrow(patient.outcome[!is.na(patient.outcome$Arm) & patient.outcome$Arm == 0,]) / 
        nrow(patient.outcome[!is.na(patient.outcome$Outcome),])
      arm_1_ratio <- nrow(patient.outcome[!is.na(patient.outcome$Arm) & patient.outcome$Arm == 1,]) / 
        nrow(patient.outcome[!is.na(patient.outcome$Outcome),])
      arm_2_ratio <- 1 - arm_0_ratio - arm_1_ratio
      
      #patient allocation following GI rules; ensuring at least 1/k of patients are in each arm; tie-break at random
      patient.outcome$Arm[t] <- 
        ifelse(nrow(patient.outcome[!is.na(patient.outcome$Outcome),]) == 0, #if there are no obs yet
               ifelse(GIs[1] != GIs[2] & GIs[1] != GIs[3], which.max(GIs) - 1, #assign by prior GI
                      ifelse(runif(1) >= 2/3, 0, ifelse(runif(1) >= 1/2, 1, 2))), #tiebreak at random
               #otherwise, if there are observations,
               ifelse(arm_0_ratio < 1/k | arm_1_ratio < 1/k | arm_2_ratio < 1/k, #and if less than 1/k of existing obs are in each arm
                      which.min(c(arm_0_ratio,arm_1_ratio,arm_2_ratio))-1, #assign to the lower one to get it up to 1/k
                      #otherwise
                      ifelse(GIs[1] != GIs[2] & GIs[1] != GIs[3], which.max(GIs) - 1, #assign by prior GI
                             ifelse(runif(1) >= 2/3, 0, ifelse(runif(1) >= 1/2, 1, 2))))) #tie-break at random
      
      #simulate outcomes based on arms
      patient.outcome$Outcome[t] <- ifelse(patient.outcome$Arm[t] == 0, 
                                           simulate_exponential_outcome(mu_0), 
                                           ifelse(patient.outcome$Arm[t] == 1,
                                                  simulate_exponential_outcome(mu_1),
                                                  simulate_exponential_outcome(mu_2)))
      #update the GIs based on outcomes
      posterior_obs <- patient.outcome$Outcome[patient.outcome$Arm == patient.outcome$Arm[t] & !is.na(patient.outcome$Outcome)] #list of all obs in this arm so far
      GIs[patient.outcome$Arm[t] + 1] <- 
        (mu_0 * implicit_prior_n + sum(posterior_obs)) / (implicit_prior_n + length(posterior_obs)) * #posterior mean
        GI_interpolated[[discount_factor]][GI_interpolated$n == implicit_prior_n + length(posterior_obs)] #table value of v(∑,n) where lambda = 1
    }
    
    #extract trial statistics from this outcome dataframe
    trial.outcomes$Prop.Superior[i] <- 
      nrow(patient.outcome[patient.outcome$Arm == which.max(c(mu_0, mu_1, mu_2)) - 1,]) / N
    #test statistics here is the p-value generated by a f test - test mu_1-0 and mu_2-0
    trial.outcomes$Test.Stat.1[i] <- pf(mean(patient.outcome$Outcome[patient.outcome$Arm == 0])/mean(patient.outcome$Outcome[patient.outcome$Arm == 1]),
                                        df1 = 2*nrow(patient.outcome[patient.outcome$Arm == 0,]), df2 = 2*nrow(patient.outcome[patient.outcome$Arm == 1,]))
    trial.outcomes$Test.Stat.2[i] <- pf(mean(patient.outcome$Outcome[patient.outcome$Arm == 0])/mean(patient.outcome$Outcome[patient.outcome$Arm == 2]),
                                        df1 = 2*nrow(patient.outcome[patient.outcome$Arm == 0,]), df2 = 2*nrow(patient.outcome[patient.outcome$Arm == 2,]))
    #record the estimates
    trial.outcomes$Estimate.1[i] <- mean(patient.outcome$Outcome[patient.outcome$Arm == 1]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])
    trial.outcomes$Estimate.2[i] <- mean(patient.outcome$Outcome[patient.outcome$Arm == 2]) - mean(patient.outcome$Outcome[patient.outcome$Arm == 0])
    #and outcomes
    trial.outcomes$Total.Outcome[i] <- sum(patient.outcome$Outcome)
    #print for progress
    print(paste(previous_message, i, 'out of', R, 'trials simulated'))
  }
  
  trial.outcomes #output
}


#run simulation for when null is true and when alt is true
mu_0 <- 0.5
mu_1 <- 0.7
mu_2 <- 0.8
N <- 100
null.true <- simulate_trials_GI_constraint(N, mu_0, mu_0, mu_0, R = 1000, k = 5)
alte.true <- simulate_trials_GI_constraint(N, mu_0, mu_1, mu_2, R = 1000, k = 5)

#compute trial properties
#since the stats here is p-value, we set cutoff to 0.05, and reverse the inequalities 
#again, I am not sure if this is the right test for exponential endpoints
cutoff <- 0.05 / 2 #Bonferroni correction 
#get family-wise alpha and power in different cases
alpha <- 
  if(mu_1>mu_0 & mu_2>mu_0){
    nrow(null.true[null.true$Test.Stat.1 < cutoff | null.true$Test.Stat.2 < cutoff,]) / nrow(null.true)
  }else if(mu_1>mu_0 & mu_2<mu_0){
    nrow(null.true[null.true$Test.Stat.1 < cutoff | null.true$Test.Stat.2 > 1-cutoff,]) / nrow(null.true)
  }else if(mu_1<mu_0 & mu_2>mu_0){
    nrow(null.true[null.true$Test.Stat.1 > 1-cutoff | null.true$Test.Stat.2 < cutoff,]) / nrow(null.true)
  }else{
    nrow(null.true[null.true$Test.Stat.1 > 1-cutoff | null.true$Test.Stat.2 > 1-cutoff,]) / nrow(null.true)
  }

power <- 
  if(mu_1>mu_0 & mu_2>mu_0){
    1 - nrow(alte.true[alte.true$Test.Stat.1 > cutoff | alte.true$Test.Stat.2 > cutoff,]) / nrow(alte.true)
  }else if(mu_1>mu_0 & mu_2<mu_0){
    1 - nrow(alte.true[alte.true$Test.Stat.1 > cutoff | alte.true$Test.Stat.2 < 1-cutoff,]) / nrow(alte.true)
  }else if(mu_1<mu_0 & mu_2>mu_0){
    1 - nrow(alte.true[alte.true$Test.Stat.1 < 1-cutoff | alte.true$Test.Stat.2 > cutoff,]) / nrow(alte.true)
  }else{
    1 - nrow(alte.true[alte.true$Test.Stat.1 < 1-cutoff | alte.true$Test.Stat.2 < 1-cutoff,]) / nrow(alte.true)
  }
  
prop_superior_in_null <- mean(null.true$Prop.Superior)
prop_sup_in_null_sd <- sd(null.true$Prop.Superior)
prop_superior_in_alte <- mean(alte.true$Prop.Superior)
prop_sup_in_alte_sd <- sd(alte.true$Prop.Superior)
#added later (August 22nd, not reported below)
ETO_change_percentage_null <- mean((null.true$Total.Outcome - mu_0*N)/(mu_0*N)*100)
ETO_change_percentage_alte <- mean((alte.true$Total.Outcome - (mu_0*N/3+mu_1*N/3+mu_2*N/3))/(mu_0*N/3+mu_1*N/3+mu_2*N/3)*100)

#true effect size
true_effect_1 <- mu_1 - mu_0
true_effect_2 <- mu_2 - mu_0

#report results
cat(paste(nrow(null.true), ' trials were simulated; true effect size = ', 
          round(true_effect_1, 3), ', ', round(true_effect_2, 3),
          '\nAt a test statistics cut-off of ', cutoff, 
          ', Type I error = ', alpha, '; Power = ', round(power,4), 
          '\nH0 == T: P(superior) = ', round(prop_superior_in_null, 4), 
          ' (SD = ', round(prop_sup_in_null_sd, 2), '), ETO change = ', 
          round(ETO_change_percentage_null,2), '%;',
          '\nH1 == T: P(superior) = ', round(prop_superior_in_alte, 4), 
          ' (SD = ', round(prop_sup_in_alte_sd, 2), '), ETO change = ', 
          round(ETO_change_percentage_alte,2), '%.', sep = ''))

