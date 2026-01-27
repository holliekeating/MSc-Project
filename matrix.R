library(data.table)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(plotly)
library(fda)
library(gridExtra)

setwd("C:/Users/holli/OneDrive - National University of Ireland, Galway/Modules/Research Project")

# Loop which reads in separate excel file for each patient and compiles all data into a dataframe
# Merges time data and demogrpahic data which contains whether patients are controls or cases


times = fread("archive/time.csv")
demographic <- fread("archive/demographic.csv")
  
filenames = list.files("archive/patients/")

p_list <- list()

for(i in 1:length(filenames)){
  
  pathtomyfile = paste0("archive/patients/", filenames[i])
  temp <- fread(pathtomyfile, sep = ",", header = "auto")
  temp <- temp[ ,1:12]
  colnames(temp)<- c("subject", "trial", "condition", "sample", "fp1", "af7", "af3", "f1", "f3", "f5", "f7", "ft7")
  temp <- merge(temp, times, by ="sample")
  temp <- merge(temp, demographic, by = "subject")
  temp<- temp%>%
    filter (time_ms > -100)%>%   # time from -100ms to remove unnecessary data while including some pre beep curve info (100ms is where the beep occurs)
    filter(trial==1) # Only one trial out of 100 choosen. Trial 1 so patient hasn't learnt the test/timing/ keep up concealment
  temp<- rename(temp, schizophrenia = group)
  
  
  
  p_list[[i]] <- temp
  
}

p_df <- do.call(rbind, p_list)




p_df<- p_df%>%
  mutate(
    subject = as.factor(subject),
    schizophrenia = as.factor(schizophrenia)
  )


# Plot showing fp1 over time for conditions 1,2 and 3 colored by cases/controls

str(p_df$subject)
plot1 <- ggplot(p_df, aes(x = time_ms, y = fp1, 
                          group = subject,
                          colour = as.factor(condition)))+geom_line(size = 0.1)+
  facet_wrap(~schizophrenia, ncol = 3)

ggplotly(plot1)

plot2 <- ggplot(p_df, aes(x = time_ms, y = fp1, 
                          group = subject,
                          colour = as.factor(schizophrenia)))+geom_line(size = 0.1, se=F)+
                          facet_wrap(~ condition, ncol = 3)

ggplotly(plot2)

p_df <- p_df %>%
  arrange(subject, condition, time_ms)



# =============================================================================
## Smooth the data - functional principle component analysis


## Create a B-spline basis (mybasis) - building blocks to represent the curves
K = 100 # basis functions. More = more flexible
#?create.bspline.basis

d = 3 # basis degree e.g. 3 = cubic
mybasis = create.bspline.basis(rangeval = c(0,1638), # x-axis range for curves 0-1638 timepoints
                               nbasis = K, 
                               norder = d+1) # cubic is an order 4 polynomial 



## Smoothing parameter object

## create an fd object - functional data object

penalty_order = 2 # 2nd derivative ie penalizes wigglier
penalty = 0.05 #amount of smoothing to be applied. Larger penalty = smoother.
fdParObj <- fdPar(fd(matrix(0, K, 1), mybasis), # fd object with k=100 coeff, 1 curve
                  Lfdobj = penalty_order, 
                  lambda = penalty)
#?fdPar





# ============================================================================
# Comparison of cases vs controls for each condition (1,2,3)
#=============================================================================

# Condition 1 - cases vs controls

# Control
cond1_s0 <- p_df %>% 
  filter(condition == 1) %>% 
  filter(schizophrenia == 0)

# Schizophrenia
cond1_s1 <- p_df %>% 
  filter(condition == 1) %>%
  filter(schizophrenia == 1)


# Matrix to compare condition 1 - cases vs controls

  # 1638 - no of rows for patient 1 (timepoints)
  # 25 controls (ncols)
  # 11 patients with schizophrenia (ncols)

matrix_cond1_s0 <- matrix(cond1_s0$fp1, nrow = 1638, ncol = 25, byrow= FALSE)
matrix_cond1_s0


matrix_cond1_s1 <- matrix(cond1_s1$fp1, nrow = 1638, ncol = 11, byrow= FALSE)
matrix_cond1_s1


## smooth the data - cond 1 S 0 (control)
sm_cond1_s0 <- smooth.basis(seq(1, 1638, length.out = 1638),
                                  matrix_cond1_s0, 
                                  fdParObj)$fd
plot(sm_cond1_s0)

## smooth the data - cond 1 s 1
# smooth.basis(x, y, fdParObj)
sm_cond1_s1 <- smooth.basis(seq(1, 1638, length.out = 1638), # x
                            matrix_cond1_s1, # y
                            fdParObj)$fd     # extracts the smoothed functional data (fd) object ie the curves
plot(sm_cond1_s1)


#=============================================================================
# Condition 1 - Plot one mean line for each group cases and controls for comparison
# Expecting cond 1 to be most statistically different 
# Plot cases mean function and then add the controls mean function

mean_cond1_s0 <- mean.fd(sm_cond1_s0)
mean_cond1_s1 <- mean.fd(sm_cond1_s1)

# Plot both mean curves together
mean.fd <- plot(mean_cond1_s1, 
     col = "blue", lwd = 2,
     xlab = "Time", ylab = "Value",
     main = "Mean Smoothed Curves: Cases vs Controls",
     ylim = c(-20, 40))

lines(mean_cond1_s0, col = "green", lwd = 2)
abline(v = 100, lty = 2, col = "red", lwd = 2) # threshold line at 100ms when stimulus occurs

legend("topright",
       legend = c("Controls", "Cases"),
       col = c("blue", "green"),
       lwd = 2)

mean.fd



# ============================================================================
  
# Condition 2 - cases vs controls
  
# Control
  cond2_s0 <- p_df %>% 
  filter(condition == 2) %>% 
  filter(schizophrenia == 0)

# Schizophrenia
cond2_s1 <- p_df %>% 
  filter(condition == 2) %>%
  filter(schizophrenia == 1)



# Matrix to compare condition 2 - cases vs controls

matrix_cond2_s0 <- matrix(cond2_s0$fp1, nrow = 1638, ncol = 25, byrow= FALSE)
matrix_cond2_s0


matrix_cond2_s1 <- matrix(cond2_s1$fp1, nrow = 1638, ncol = 11, byrow= FALSE)
matrix_cond2_s1



## smooth the data - cond 2 s 0
sm_cond2_s0 <- smooth.basis(seq(1, 1638, length.out = 1638),
                            matrix_cond2_s0, 
                            fdParObj)$fd
plot(sm_cond2_s0)

## smooth the data - cond 2 s 1
sm_cond2_s1 <- smooth.basis(seq(1, 1638, length.out = 1638),
                            matrix_cond2_s1, 
                            fdParObj)$fd
plot(sm_cond2_s1)










#===========================================================================
# Condition 3 - cases vs controls

# Control
cond3_s0 <- p_df %>% 
  filter(condition == 3) %>% 
  filter(schizophrenia == 0)

# Schizophrenia
cond3_s1 <- p_df %>% 
  filter(condition == 3) %>%
  filter(schizophrenia == 1)


# Matrix to compare condition 3 - cases vs controls

matrix_cond3_s0 <- matrix(cond3_s0$fp1, nrow = 1638, ncol = 25, byrow= FALSE)
matrix_cond3_s0


matrix_cond3_s1 <- matrix(cond3_s1$fp1, nrow = 1638, ncol = 11, byrow= FALSE)
matrix_cond3_s1



## smooth the data - cond 3 s 0
sm_cond3_s0 <- smooth.basis(seq(1, 1638, length.out = 1638),
                            matrix_cond3_s0, 
                            fdParObj)$fd
plot(sm_cond3_s0)

## smooth the data - cond 3 s 1
sm_cond3_s1 <- smooth.basis(seq(1, 1638, length.out = 1638),
                            matrix_cond3_s1, 
                            fdParObj)$fd
plot(sm_cond3_s1)


# ===========================================================================
# Permutation T tests

# ===========================================================================
# Permutation t tests compares entire curves between two groups as oppose to single data points
# n=200 permutations/resamplings
# 0.05 significance level


# Permutation t test for comparing condition 1 control vs case 

tperm_cond1 <- tperm.fd(sm_cond1_s0, sm_cond1_s1, nperm=200, q=0.05, argvals=NULL, plotres= TRUE)
tperm_cond1

# Permutation t test for comparing condition 2 control vs case 

tperm_cond2 <- tperm.fd(sm_cond2_s0, sm_cond2_s1, nperm=200, q=0.05, argvals=NULL, plotres= TRUE)
tperm_cond2

# Permutation t test for comparing condition 3 control vs case 

tperm_cond3 <- tperm.fd(sm_cond3_s0, sm_cond3_s1, nperm=200, q=0.05, argvals=NULL, plotres= TRUE)
tperm_cond3




