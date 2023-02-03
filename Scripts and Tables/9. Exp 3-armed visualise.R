#this script visualise the simulation results from script 8

##########################
#FILES and LIBS
##########################

library(plotly) #for 3d plots
library(ggplot2) #for 2d plots

mu_range_outcome_ER <- read.csv('3-armed EXP/EXP_mu_range_outcome_ER.csv')[,-1]
mu_range_outcome_GI_const_5 <- read.csv('3-armed EXP/EXP_mu_range_outcome_GI_const_5.csv')[,-1]
mu_range_outcome_GI_const_9 <- read.csv('3-armed EXP/EXP_mu_range_outcome_GI_const_9.csv')[,-1]
mu_range_outcome_GI_near_op <- read.csv('3-armed EXP/EXP_mu_range_outcome_GI_near_op.csv')[,-1]

#combine the trial property outcome
mu_range_outcome_combined <- data.frame(rbind(mu_range_outcome_ER, mu_range_outcome_GI_const_5, 
                                              mu_range_outcome_GI_const_9, mu_range_outcome_GI_near_op),
                                        Design = c(rep('ER', nrow(mu_range_outcome_ER)), 
                                                   rep('GI-const-5', nrow(mu_range_outcome_GI_const_5)), 
                                                   rep('GI-const-9', nrow(mu_range_outcome_GI_const_9)),
                                                   rep('GI-near_op', nrow(mu_range_outcome_GI_near_op))))

##########################
#Visualise
##########################

#3d plots - save as html using built-in
plot_ly(data = mu_range_outcome_combined,
        x = ~mu_1, y = ~mu_2, z = ~Alpha, color = ~Design,
        type = "scatter3d", mode = "markers")

plot_ly(data = mu_range_outcome_combined,
        x = ~mu_1, y = ~mu_2, z = ~Power, color = ~Design,
        type = "scatter3d", mode = "markers")

plot_ly(data = mu_range_outcome_combined,
        x = ~mu_1, y = ~mu_2, z = ~H1.Estimate.1.SD, color = ~Design,
        type = "scatter3d", mode = "markers")

plot_ly(data = mu_range_outcome_combined,
        x = ~mu_1, y = ~mu_2, z = ~H1.Estimate.2.SD, color = ~Design,
        type = "scatter3d", mode = "markers")

plot_ly(data = mu_range_outcome_combined,
        x = ~mu_1, y = ~mu_2, z = ~H1.Prop.Superior, color = ~Design,
        type = "scatter3d", mode = "markers")

plot_ly(data = mu_range_outcome_combined,
        x = ~mu_1, y = ~mu_2, z = ~H1.ETO.Change.Percent, color = ~Design,
        type = "scatter3d", mode = "markers")


