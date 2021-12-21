# Simulations to compare bootstrapped 95%CI metrics for highly skewed populations (count data) and small sample size (N = 8)
# Metrics compared: standard-t, bootstrap-t interval, bootstrap percentile method using quantile type 6 vs. boot.ci default settings, boot.ci bca
# If a simulated sample had no variability (e.g., all 0's or all 1's), the sample was not used

# UPSHOT IS THAT BOOT-T AND STANDARD-T HAD BEST COVERAGE, THOUGH BOTH COULD SOMETIMES LEAD TO NEGATIVE VALUES FOR LOWER 95%CI. COVERAGE WAS TOO LOW FOR SMALL (N=8) VERY SKEWED POPS NO MATTER WHICH OF THESE METRICS WAS USED. UNDER THESE CONDITIONS, BOOT-T SOMETIMES HAD CRAZY LARGE UPPER 95% CI VALUES AND STANDARD-T OFTEN HAD MORE NEGATIVE LOWER 95% CI VALUES COMPARED TO BOOT-T. FOR VERY SKEWED, BOOT-T COVERAGE WAS REASONABLY GOOD ONCE SAMPLE SIZE REACHED N = 16.
# FOR BOTH STANDARD-T AND BOOT-T, I'M USING STANDARD ERROR CALCULATED FROM THE SAMPLE

simfun <- function(pop, samp_n) {

  repeat{ # pick a sample from the pop
    samp <- sample(pop, samp_n)
    if (length(unique(samp))>1) break # do not calculate for samples that are all the same number
  }
  
  # for calculating normal-approximation 95%CI with t-interval
  t_se = sd(samp)/sqrt(samp_n)
  t_margin = qt(0.975,df=samp_n-1)*t_se
  
  # bootstrapped 95%CI's
  bs <- boot::boot(samp,statistic=function(xx,index) {
    bs_mean <- mean(xx[index])
    bs_sd <- sd(xx[index])
    bs_t <- (bs_mean - mean(samp))/(bs_sd/sqrt(samp_n)) # for calculating bootstrap-t interval
    c(bs_mean, bs_t) # these are the bootstrap metrics output
  }, R=1000) # bootstrap from the sample to get a bootstrap distribution
  
  boots_ci = boot::boot.ci(boot.out = bs,
                           type = c("perc", "bca")) # using boot.ci for comparison. Calculate percentile and bca intervals
  
  out <- rbind(
    c("standard_t", mean(samp) - t_margin, mean(samp) + t_margin), # standard t
    c("boot_t", mean(samp)-unname(quantile(bs$t[,2][!is.infinite(bs$t[,2])], probs=c(0.975, 0.025), type = 6, na.rm = TRUE))*t_se), # for bootstrap-t interval need to swap the probs
    c("boot_p", unname(quantile(bs$t[,1], probs=c(0.025, 0.975), type = 6, na.rm = TRUE))), # boostrap percentile interval with type 6 quantile calc
    c("bootci_p", boots_ci$percent[4], boots_ci$percent[5]), # bootstrap percentile interval, using boot.ci
    c("bootci_bca", boots_ci$bca[4], boots_ci$bca[5]) # bca interval
  )
}

# Run the simulations
pop <- rnbinom(n=10000, mu=1.2, size=1.2) # <<<<<< create a negative binomial distribution population
# pop <- rnorm(n=10000, mean=1.2, sd=0.4)  # create a normal distribution population
CI_list <- list() # create list to hold CI estimates for the different methods

for(i in 1:10000) {
  CI_list[[i]] <- simfun(pop, samp_n = 8) # <<<<<<<<<<<<<<<< change the sample size here
}

CI_df <- do.call(rbind.data.frame, CI_list)
names(CI_df) <- c("CI_metric", "low95", "high95")
CI_df$low95 <- round(as.numeric(CI_df$low95), 3)
CI_df$high95 <- round(as.numeric(CI_df$high95), 3)
CI_df$includes_truth <- CI_df$low95 <= mean(pop) & CI_df$high95 >= mean(pop) # create column indicating if CI included the true pop mean

# Calculate summaries of % coverage AND mean and median 95%CI size by metric
CI_df %>% 
  group_by(CI_metric) %>% 
  summarize(CI_coverage = mean(includes_truth, na.rm = T),
            mean_CI = mean(high95 - low95),
            median_CI = median(high95 - low95),
            perc_truth_below_low = mean(mean(pop) < low95, na.rm = T),
            perc_truth_above_high = mean(mean(pop) > high95, na.rm = T))%>% 
  knitr::kable(digits = 3)
  
# Save data frame
write_csv(CI_df, "CIsims_nbinom_mu1pt2_size1pt2_n8.csv") # <<<<<<

### 95% CI COVERAGE WITH DIFFERENT SIMULATION PARAMETERS ----
### REMEMBER THAT WHEN A BOOTSTRAP SAMPLE WAS ALL THE SAME NUMBER (SO VARIABILITY COULD NOT BE CALCULATED), IT WAS DISCARDED. WITH THE SKEWED POPS, THESE CASES WERE ALMOST CERTAINLY WHEN THE ENTIRE SAMPLE WAS ZERO'S. WE COULD HAVE ALTERNATIVELY CONSIDERED THAT AS 'DID NOT INCLUDE TRUTH' BECAUSE THE ESTIMATED MEAN AND CI WOULD HAVE BEEN ALL ZERO'S
# RESULTS FOR HIGHLY SKEWED POPS -------------------------------

### HIGHLY SKEWED, LARGE SAMPLE SIZE (NEGBINOM, MU=1.2, SIZE=0.4, N = 80) >> BOOT-T IS BEST
### 95%CI data frame saved as "CIsims_nbinom_mu1pt2_size0pt4_n80.csv"
# |CI_metric  | CI_coverage| mean_CI| median_CI| perc_truth_below_low| perc_truth_above_high|
# |:----------|-----------:|-------:|---------:|--------------------:|---------------------:|
# |boot_p     |       0.930|   0.904|     0.875|                0.009|                 0.070|
# |boot_t     |       0.950|   1.027|     0.970|                0.016|                 0.035|
# |bootci_bca |       0.933|   0.945|     0.902|                0.018|                 0.052|
# |bootci_p   |       0.930|   0.904|     0.875|                0.009|                 0.070|
# |standard_t |       0.929|   0.924|     0.895|                0.005|                 0.075|

### HIGHLY SKEWED, LARGE SAMPLE SIZE (NEGBINOM, MU=1.2, SIZE=0.4, N = 16) >> BOOT-T IS BEST
### 95%CI data frame saved as "CIsims_nbinom_mu1pt2_size0pt4_n16.csv"
# |CI_metric  | CI_coverage| mean_CI| median_CI| perc_truth_below_low| perc_truth_above_high|
# |:----------|-----------:|-------:|---------:|--------------------:|---------------------:|
# |boot_p     |       0.847|   1.889|     1.688|                0.008|                 0.145|
# |boot_t     |       0.927|   3.597|     2.696|                0.009|                 0.065|
# |bootci_bca |       0.858|   2.124|     1.858|                0.018|                 0.124|
# |bootci_p   |       0.848|   1.881|     1.688|                0.008|                 0.145|
# |standard_t |       0.858|   2.162|     1.946|                0.001|                 0.141|


### HIGHLY SKEWED, PRETTY SMALL SAMPLE SIZE (NEGBINOM, MU=1.2, SIZE=0.4, N = 12) >> BOOT-T IS BEST BUT ALL ARE BAD
### 95%CI data frame saved as "CIsims_nbinom_mu1pt2_size0pt4_n12.csv"
# |CI_metric  | CI_coverage| mean_CI| median_CI| perc_truth_below_low| perc_truth_above_high|
# |:----------|-----------:|-------:|---------:|--------------------:|---------------------:|
# |boot_p     |       0.832|   2.001|     1.833|                0.011|                 0.156|
# |boot_t     |       0.899|   4.021|     2.941|                0.009|                 0.092|
# |bootci_bca |       0.835|   2.215|     1.917|                0.019|                 0.146|
# |bootci_p   |       0.833|   1.977|     1.751|                0.011|                 0.155|
# |standard_t |       0.863|   2.401|     2.176|                0.002|                 0.135|

### HIGHLY SKEWED, SMALL SAMPLE SIZE (NEGBINOM, MU=1.2, SIZE=0.4, N = 8) >> STANDARD T IS BEST BUT ALL ARE BAD
### 95%CI data frame saved as "CIsims_nbinom_mu1pt2_size0pt4_n8.csv". I DID THIS ONE W/10K BOOTSTRAP REPS ON EACH OF 10K BOOSTRAP SIMS, JUST TO CONFIRM. RESULTS WERE VERY SIMILAR TO JUST USING 1K BOOTSTRAP REPS ON EACH OF 10K BOOTSTRAP SIMS
# > boot_t has larger CI's but lower coverage--these conditions can lead to crazy high upper 95% CI for bootstrap-t. On the flip side, standard-t leads to more negative lower 95% CI than does bootstrap-t
# |CI_metric  | CI_coverage| mean_CI| median_CI| perc_truth_below_low| perc_truth_above_high|
# |:----------|-----------:|-------:|---------:|--------------------:|---------------------:|
# |boot_p     |       0.809|   2.279|     2.000|                0.012|                 0.179|
# |boot_t     |       0.835|   4.936|     3.158|                0.003|                 0.162|
# |bootci_bca |       0.791|   2.503|     2.000|                0.015|                 0.194|
# |bootci_p   |       0.808|   2.208|     1.875|                0.014|                 0.179|
# |standard_t |       0.845|   3.102|     2.672|                0.001|                 0.154|
  
# TRY A COUPLE OTHER LEVELS OF SKEW, WITH N = 8 -------------------------------

### HIGHLY HIGHLY SKEWED, SMALL SAMPLE SIZE (NEGBINOM, MU=1.2, SIZE=0.1, N = 8) >> ALL TERRIBLE! BUT STANDARD-T IS BEST

# |CI_metric  | CI_coverage| mean_CI| median_CI| perc_truth_below_low| perc_truth_above_high|
# |:----------|-----------:|-------:|---------:|--------------------:|---------------------:|
# |boot_p     |       0.714|   3.318|     2.250|                0.004|                 0.282|
# |boot_t     |       0.626|  10.598|     2.562|                0.001|                 0.372|
# |bootci_bca |       0.694|   3.722|     2.375|                0.010|                 0.296|
# |bootci_p   |       0.654|   3.004|     2.000|                0.069|                 0.277|
# |standard_t |       0.753|   4.777|     3.096|                0.000|                 0.247|

### LESS SKEWED, SMALL SAMPLE SIZE (NEGBINOM, MU=1.2, SIZE = 0.8, N = 8) >> ALL TERRIBLE! BUT BOOT-T IS BEST BY A HAIR

# |CI_metric  | CI_coverage| mean_CI| median_CI| perc_truth_below_low| perc_truth_above_high|
# |:----------|-----------:|-------:|---------:|--------------------:|---------------------:|
# |boot_p     |       0.841|   1.949|     1.750|                0.017|                 0.142|
# |boot_t     |       0.894|   3.566|     2.790|                0.005|                 0.101|
# |bootci_bca |       0.816|   2.055|     1.875|                0.019|                 0.165|
# |bootci_p   |       0.841|   1.908|     1.750|                0.017|                 0.141|
# |standard_t |       0.884|   2.586|     2.364|                0.002|                 0.114|

### EVEN LESS SKEWED, SMALL SAMPLE SIZE (NEGBINOM, MU=1.2, SIZE = 1.2, N = 8) >> ALL TERRIBLE! BUT BOOT-T IS BEST BY A HAIR
# |CI_metric  | CI_coverage| mean_CI| median_CI| perc_truth_below_low| perc_truth_above_high|
# |:----------|-----------:|-------:|---------:|--------------------:|---------------------:|
# |boot_p     |       0.847|   1.774|     1.625|                0.020|                 0.133|
# |boot_t     |       0.908|   3.025|     2.497|                0.007|                 0.085|
# |bootci_bca |       0.820|   1.848|     1.625|                0.020|                 0.159|
# |bootci_p   |       0.847|   1.746|     1.625|                0.021|                 0.132|
# |standard_t |       0.893|   2.342|     2.178|                0.003|                 0.104|
  
### SKEWED AND WITH LARGER MEAN, SMALL SAMPLE SIZE (NEGBINOM, MU=1.8, SIZE = 0.4, N = 8) >> ALL TERRIBLE! BUT BOOT-T IS BEST
  
# |CI_metric  | CI_coverage| mean_CI| median_CI| perc_truth_below_low| perc_truth_above_high|
# |:----------|-----------:|-------:|---------:|--------------------:|---------------------:|
# |boot_p     |       0.797|   3.453|     3.000|                0.013|                 0.190|
# |boot_t     |       0.859|   9.176|     5.435|                0.003|                 0.138|
# |bootci_bca |       0.796|   3.853|     3.203|                0.020|                 0.184|
# |bootci_p   |       0.798|   3.390|     2.875|                0.014|                 0.188|
# |standard_t |       0.830|   4.658|     3.990|                0.000|                 0.169|

# NOW FOR NORMALLY DISTRIBUTED POPS -------------------------------

### NORMAL, LARGE SAMPLE SIZE (NORMAL, MEAN=1.2, SD=0.4, N = 80) >> BOOT-T AND STANDARD-T ESSENTIALLY SPOT ON, ALL ARE GOOD ON COVERAGE BUT CI IS OFFSET
### 95%CI data frame saved as "CIsims_norm_mean1pt2_sd0pt4_n80.csv"
# |CI_metric  | CI_coverage| mean_CI| median_CI| perc_truth_below_low| perc_truth_above_high|
# |:----------|-----------:|-------:|---------:|--------------------:|---------------------:|
# |boot_p     |       0.946|   0.174|     0.173|                0.021|                 0.034|
# |boot_t     |       0.951|   0.178|     0.177|                0.020|                 0.031|
# |bootci_bca |       0.946|   0.174|     0.173|                0.021|                 0.033|
# |bootci_p   |       0.946|   0.174|     0.173|                0.021|                 0.034|
# |standard_t |       0.952|   0.177|     0.177|                0.020|                 0.030|

### NORMAL, PRETTY SMALL SAMPLE SIZE (NORMAL, MEAN=1.2, SD=0.4, N = 12) >> BOOT-T AND STANDARD-T ESSENTIALLY SPOT ON WITH COVERAGE, STANDARD-T HAS SMALLER INTERVALS
### 95%CI data frame saved as "CIsims_norm_mean1pt2_sd0pt4_n12.csv"
# |CI_metric  | CI_coverage| mean_CI| median_CI| perc_truth_below_low| perc_truth_above_high|
# |:----------|-----------:|-------:|---------:|--------------------:|---------------------:|
# |boot_p     |       0.911|   0.425|     0.422|                0.045|                 0.046|
# |boot_t     |       0.950|   0.525|     0.519|                0.025|                 0.025|
# |bootci_bca |       0.911|   0.430|     0.427|                0.044|                 0.046|
# |bootci_p   |       0.911|   0.425|     0.422|                0.045|                 0.046|
# |standard_t |       0.948|   0.498|     0.494|                0.026|                 0.025|

### NORMAL, SMALL SAMPLE SIZE (NORMAL, MEAN=1.2, SD=0.4, N = 8) >> BOOT-T AND STANDARD-T ESSENTIALLY SPOT ON WITH COVERAGE, STANDARD-T HAS SMALLER INTERVALS
### 95%CI data frame saved as "CIsims_norm_mean1pt2_sd0pt4_n8.csv"
# |CI_metric  | CI_coverage| mean_CI| median_CI| perc_truth_below_low| perc_truth_above_high|
# |:----------|-----------:|-------:|---------:|--------------------:|---------------------:|
# |boot_p     |       0.888|   0.505|     0.500|                0.056|                 0.055|
# |boot_t     |       0.951|   0.742|     0.715|                0.024|                 0.026|
# |bootci_bca |       0.886|   0.517|     0.512|                0.057|                 0.057|
# |bootci_p   |       0.888|   0.505|     0.500|                0.056|                 0.055|
# |standard_t |       0.952|   0.654|     0.648|                0.024|                 0.024|