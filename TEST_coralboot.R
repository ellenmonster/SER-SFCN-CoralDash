### CORAL BOOTSTRAP CODE FOR TOM TO LAUGH AT. CLEARLY NOT THE 8-LINE CODE TOM HAD IN MIND. 

# USER-SPECIFIED FOR TESTING CODE ----
boot_master_dat <- readRDS("boot_master_dat.RDS") # import the test data
TEST_n_bootsamples <- 10000 # number of parametric bootstrap samples per transect-survey (default value in function is 1M)
TEST_n_bootreps <- 1000 # number of bootstrap samples per site-survey, for calculating bootstrapped confidence intervals (default value in function is 10K)

### Load packages ----
pkgList_pre <- c("magrittr", 
                 "plyr",
                 "lubridate",
                 "purrr", # to map functions to elements
                 "tidyverse") # broom, units
inst_pre <- pkgList_pre %in% installed.packages()
if (length(pkgList_pre[!inst_pre]) > 0) install.packages(pkgList_pre[!inst_pre],dep=TRUE)
lapply(pkgList_pre, library, character.only = TRUE)

### Bootstrap functions ----

# Function to calculate bootstrap quantiles from long-format data and format output with one col per quantile
FuncPullQuant <- purrr::map(
  c(0.25, 0.75, 0.125, 0.875, 0.05, 0.95), 
  ~purrr::partial(quantile, probs = .x, na.rm = TRUE)) %>% # apply quantile function to each specified probability, resulting format has a separate col for each quantile (instead of a separate row for each quantile)
  set_names(nm = c("50%CI_low", "50%CI_high", "75%CI_low", "75%CI_high", "90%CI_low", "90%CI_high"))

FuncBootSamples <- function(dat, n_bootsamples = 1000000) { 
  # Function to generate parametric bootstrap samples for each transect of a site-survey
  #
  # Args:
  #   dat = data frame of the original data to bootstrap, with the taxa as cols. Includes all transects for a single site-survey
  #
  # Returns:
  #   Data frame of samples from parametric bootstrapping
  
  # Convert to % cover by transect-survey
  dat[, 5:ncol(dat)] <- dat[, 5:ncol(dat)]/dat$AdjTot 
  
  # Create list of parametric bootstrap samples, each list element is a data frame of samples for one transect-survey
  boot_samples_list <- apply(dat, 1, FUN = function(z) { # for each transect-survey
    prob_vec <- z[5:length(z)][!is.na(z[5:length(z)])] # probability vector w/o NA's
    boot_out <- rmultinom(n = n_bootsamples, size = as.numeric(z["AdjTot"]), prob = prob_vec) # generate parametric bootstrap samples of "hits" for all the taxa
    boot_df <- as.data.frame(cbind("TransectSurveyID" = z["TransectSurveyID"], t(boot_out))) # each list element is a data frame with TransectSurveyID and the taxa as cols (values are # of hits in transect-survey). Don't calculate % cover yet because for relative cover by category will need the actual number of hits
  })
  
  boot_samples_df <- plyr::rbind.fill(boot_samples_list) # combine list elements in a single data frame
  boot_samples_df[2:ncol(boot_samples_df)] <- sapply(boot_samples_df[2:ncol(boot_samples_df)], as.numeric) # convert cols to numeric
  return(boot_samples_df)
}

FuncBootDraws <- function(boot_dat, actual_dat, groups_df, n_bootreps = 10000) {
  # Function to calculate bootstrap CI's for % cover and relative % cover, by site-survey (for intensive sites)
  #
  # Args:
  #   boot_dat = data frame of the parametric bootstrap samples for all transects of a single site-survey. Requires cols for "TransectSurveyID" and each taxon. 
  #   actual_dat = data frame of the original count data used to generate parametric bootstrap samples. Same cols as boot_dat, but includes an AdjTot col.
  #   groups_df = data frame assigning each taxon to higher level groups
  #   n_bootreps = number of bootstrap reps to calculate CI's over
  #
  # Returns:
  #   A data frame with bootstrapped CI's for % cover and relative % cover. Has these cols: DenomGroup (e.g., AdjTot, ALGAE, CORAL, GORGO, SPONGE), NumerGroup (e.g., Category, FunctionalGroup, Taxon), NumerLevel(e.g., ALGAE). Final columns have the quantiles for bootstrapped 50%, 75%, 95%CI.
  #
  # NOTE: Not estimating CI's for trend -- just individually by site-survey.
  
  n_transects <- nrow(actual_dat) # number of transects surveyed in that site-survey
  
  samps <- boot_dat %>%
    dplyr::slice_sample(n = n_draws * n_transects, replace = TRUE) %>% # randomly sample from the parametric bootstrap samples
    as.data.frame() # can't add repeating vector column to tibble
  
  samps[, 2:ncol(samps)][is.na(samps[, 2:ncol(samps)])] <- 0 # replace NA's with zero
  samps$BootRep = rep(1:n_draws, each = n_transects) # add col of bootstrap replicate number; each replicate has same number of transects as original data for the site-survey
  
  # Add the actual data so actual % cover and relative % cover can be simultaneously calculated for each combination of numerator and denominator groups/levels. BootRep 0 is the actual data and will be excluded for CI estimates
  samps <- plyr::rbind.fill(samps, actual_dat %>% dplyr::select(-AdjTot) %>% dplyr::mutate(BootRep = 0)) %>%
    tibble::rowid_to_column(var = "RowID") # after converting to long format, still need to be able to group data by row at times
  
  # This is the master data frame from which to calculate % cover at different grouping levels.
  samps_master <- samps %>%
    tidyr::gather(key = "Taxon", value = "Count", -BootRep, -RowID, -TransectSurveyID) %>% # convert to long format
    dplyr::left_join(groups_df, by = "Taxon")
  
  # Calculate % cover and RELATIVE % cover for certain combinations of numerator and denominator groups/levels:
  # > % cover by Category. Denominator is adjusted total count for the sampled transect. Transects have equal weight in a BootRep.
  # > Drill down in certain categories--Calculate % cover and RELATIVE % cover by Functional Group for ALGAE and CORAL (NumerGroup = FunctionalGroup, DenomGroup = ALGAE or CORAL); and by Taxon for ALGAE, CORAL, GORGO, SPONGE (e.g., NumerGroup = Taxon, DenomGroup = GORGO). Weighted by the # of hits of the category in each transect.
  
  # Calculate bootstrapped CI's for these combinations (can add to this as requested by SFCN)
  combos_df <- data.frame(rbind(
    c("Category", "AdjTot"),
    c("FunctionalGroup", "ALGAE"),
    c("FunctionalGroup", "CORAL"),
    c("Taxon", "ALGAE"),
    c("Taxon", "CORAL"),
    c("Taxon", "GORGO"),
    c("Taxon", "SPONGE")
  ))
  names(combos_df) <- c("NumerGroup", "DenomGroup")
  
  boot_CIs_list <- apply(combos_df, 1, FUN = function(x) {
    
    tmp <- samps_master %>%
      dplyr::rename(NumerLevel = x[["NumerGroup"]]) %>%
      dplyr::filter(if(x[["DenomGroup"]] != "AdjTot") Category == x[["DenomGroup"]] else TRUE) %>% # for drill down (i.e., DenomGroup is a Category instead of AdjTot), restrict to one category; otherwise keep all categories
      dplyr::group_by(BootRep, TransectSurveyID, RowID, NumerLevel) %>%
      dplyr::summarize(GroupCount = sum(Count), .groups = "drop")
    
    # For % cover, calculate per transect (RowID) and then take mean across all transects in the BootRep. The adjusted total count (denominator value) for any bootstrapped transect is the same as for the original data used to generate parametric bootstrap samples
    out_df <- tmp %>%
      dplyr::left_join(actual_dat[c("TransectSurveyID", "AdjTot")], by = "TransectSurveyID") %>%
      dplyr::mutate(PercCov = GroupCount / AdjTot) %>%
      dplyr::group_by(BootRep, NumerLevel) %>%
      dplyr::summarize(EstimCov = mean(PercCov), .groups = "drop")
    
    # Extract the means calculated for the original data (BootRep = 0) so can be added to the final data frame
    actual_df <- subset(out_df, BootRep == 0) %>% 
      dplyr::select(-BootRep) 
    
    # Calculate % cover CI's
    out_df %<>% 
      dplyr::filter(BootRep != 0) %>% # remove the estimates from the actual data
      dplyr::group_by(NumerLevel) %>%
      dplyr::summarize_at(vars(EstimCov), FuncPullQuant, .groups = "drop") %>%
      dplyr::mutate(NumerGroup = x[["NumerGroup"]], DenomGroup = "AdjTot") %>%
      dplyr::left_join(actual_df, by = "NumerLevel") # add in the mean % cov from actual data
    
    # For relative % cover, sum(numerator across all transects)/sum(denom across all transects)
    if(x[["DenomGroup"]] != "AdjTot") {
      categtot_df <- tmp %>%
        dplyr::group_by(BootRep) %>%
        dplyr::summarize(CategTot = sum(GroupCount), .groups = "drop")
      
      relcov_df <- tmp %>%
        dplyr::group_by(BootRep, NumerLevel) %>%
        dplyr::summarize(SampGroupCount = sum(GroupCount), .groups = "drop") %>%
        dplyr::left_join(categtot_df, by = "BootRep") %>%
        dplyr::mutate(EstimCov = SampGroupCount/CategTot)
      
      # Extract the means calculated for the original data (BootRep = 0) so can be added to the final data frame
      actual_relcov_df <- subset(relcov_df, BootRep == 0) %>% 
        dplyr::select(NumerLevel, EstimCov) 
      
      # Calculate % cover CI's
      relcov_df %<>% 
        dplyr::filter(BootRep != 0) %>% # remove the estimates from the actual data
        dplyr::group_by(NumerLevel) %>%
        dplyr::summarize_at(vars(EstimCov), FuncPullQuant, .groups = "drop") %>%
        dplyr::mutate(NumerGroup = x[["NumerGroup"]], DenomGroup = x[["DenomGroup"]]) %>%
        dplyr::left_join(actual_relcov_df, by = "NumerLevel") # add in the mean % cov from actual data
      
      out_df <- rbind(out_df, relcov_df)
    }
    
    out_df %<>%
      dplyr::mutate_if(is.numeric, round, 3) %>%  # calculate quantiles and output in nice format
      dplyr::mutate(Category = ifelse(x[["DenomGroup"]] == "AdjTot", "All", x[["DenomGroup"]])) %>%
      dplyr::select(Category, NumerGroup, NumerLevel, DenomGroup, EstimCov, everything())
  })
  
  boot_CIs_df <- as.data.frame(do.call("rbind", boot_CIs_list)) # combine list elements in a single data frame
  return(boot_CIs_df)
}

### MAIN CODE ----

# Make sure each taxon belongs to only one Functional Group and only one Category
grps_df <- unique(boot_master_dat[c("Category", "FunctionalGroup", "Taxon")])
if(length(unique(grps_df$Taxon)) != nrow(grps_df)) stop("EACH TAXON CAN ONLY BE ASSIGNED TO ONE FUNCTIONAL GROUP AND TO ONE CATEGORY")

# For intensive sites (typically 20 transects), create data frame with one row per TransectSurveyID and with Taxa as columns. Values are the % cover for each taxon for that TransectSurveyID 
boot_site_dat <- boot_master_dat %>%
  dplyr::filter(ReportingSite == Site) %>% # only bootstrapping intensive sites. For the extensive sites, will permute and do hierarchical sampling.
  dplyr::group_by(Site, SurvDate, TransectSurveyID, Taxon) %>%
  dplyr::summarize(Count = sum(CountOfTaxon), .groups = "drop") %>% # for each transect-survey, count hits per taxon
  tidyr::spread(key = Taxon, value = Count) %>% # each taxon becomes a column
  dplyr::ungroup() %>%
  dplyr::mutate(AdjTot = rowSums(across(4:ncol(.)), na.rm = T)) %>% # calculate total hits per transect-survey (already excludes shadow and equipment points)
  dplyr::select(Site, SurvDate, TransectSurveyID, AdjTot, everything()) # rearrange cols

# For each transect-survey, generate parametric bootstrap samples (default is N = 1M) using the probabilities from the collected data. Run one site-survey at a time to avoid computer memory problems ---- 
site_CIs_list <- apply(unique(boot_site_dat[c("Site", "SurvDate")]), 1, FUN = function(x) { 
  subdat <- boot_site_dat %>%
    dplyr::filter(Site == x[["Site"]] & SurvDate == x[["SurvDate"]]) %>% # pull data for the specific site and survey date
    dplyr::select(where(~!all(is.na(.x)))) # get rid of cols with no data
  subdat[, 5:ncol(subdat)][is.na(subdat[, 5:ncol(subdat)])] <- 0 # and then replace the remaining NA's with 0
  
  # Generate parametric bootstrap samples
  boot_samples <- FuncBootSamples(dat = subdat, n_bootsamples = TEST_n_bootsamples)
  
  # For various combinations of numerator groups and denominator groups, calculate quantiles from bootstrap distribution
  boot_site_estimates <- FuncBootDraws(boot_dat = boot_samples, actual_dat = subdat %>% dplyr::select(-Site, -SurvDate), groups_df = grps_df, n_bootreps = TEST_n_bootreps)
  
  cbind("Site" = x[["Site"]], "SurvDate" = x[["SurvDate"]], boot_site_estimates)
})

site_CIs_df <- as.data.frame(do.call("rbind", site_CIs_list))