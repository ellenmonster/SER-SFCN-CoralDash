############################################################
### CORAL DASHBOARD DATA IMPORT, FORMAT AND CALCULATIONS ###
############################################################

### LOAD PACKAGES ----
# Look for packages on local machine, install if necessary, then load all
pkgList_pre <- c("magrittr", 
                 "vroom",
                 "plyr",
                 "lubridate",
                 "ggpubr", # to get legends
                 "tidyselect", # for where() function
                 "purrr", # to map functions to elements
                 "tidyverse") # broom, units
inst_pre <- pkgList_pre %in% installed.packages()
if (length(pkgList_pre[!inst_pre]) > 0) install.packages(pkgList_pre[!inst_pre],dep=TRUE)
lapply(pkgList_pre, library, character.only = TRUE)

### FUNCTIONS ----
# Function to calculate bootstrap quantiles from long-format data and format output with one col per quantile
FuncPullQuant <- function() {
  purrr::map(
  c(0.25, 0.75, 0.16, 0.84, 0.05, 0.95), 
  ~purrr::partial(quantile, probs = .x, na.rm = TRUE)) %>% # apply quantile function to each specified probability, resulting format has a separate col for each quantile (instead of a separate row for each quantile)
  set_names(nm = c("50%CI_low", "50%CI_high", "68%CI_low(1SE)", "68%CI_high(1SE)", "90%CI_low", "90%CI_high"))
}

FuncBootSamples <- function(dat, start_col, n_bootsamples = 100000) { 
  # Function to generate parametric bootstrap samples for each row of dat
  #
  # Args:
  #   dat = data frame of the original data to bootstrap, with the counts spread across cols
  #   start_col = starting column number for the count data
  #   n_bootsamples = number of parametric bootstrap samples to generate per row of dat
  #
  # Returns:
  #   Data frame of samples from parametric bootstrapping
  
  # Convert to % cover by row
  dat[, start_col:ncol(dat)] <- dat[, start_col:ncol(dat)]/dat$AdjTot 
  
  # Create list of parametric bootstrap samples, each list element is a data frame of samples for one row of dat
  boot_samples_list <- apply(data.frame(dat), 1, FUN = function(z) { # for each row
    prob_vec <- z[start_col:length(z)][!is.na(z[start_col:length(z)])] # probability vector w/o NA's
    boot_out <- rmultinom(n = n_bootsamples, size = as.numeric(z["AdjTot"]), prob = prob_vec) # generate parametric bootstrap samples of "hits" for all the taxa
    
    boot_df <- as.data.frame(cbind(data.frame(t(z[1:(start_col-1)]), row.names = NULL), t(boot_out)))
    
    boot_df
  })
  
  boot_samples_df <- plyr::rbind.fill(boot_samples_list) # combine list elements in a single data frame
  
  boot_samples_df[start_col:ncol(boot_samples_df)] <- sapply(boot_samples_df[start_col:ncol(boot_samples_df)], as.numeric) # convert cols to numeric
  return(boot_samples_df)
}

FuncBootDraws <- function(boot_dat, actual_dat, start_col, groups_df, combos_df, key_col, n_bootreps = 5000) {
  # Function to calculate bootstrap CI's for % cover and relative % cover, by site-survey (for intensive sites)
  #
  # Args:
  #   boot_dat = data frame of the parametric bootstrap samples for all transects of a single site-survey. 
  #   actual_dat = data frame of the original count data used to generate parametric bootstrap samples.
  #   start_col = starting column number for the count data
  #   groups_df = data frame assigning each items to higher level groups
  #   combos_df = data frame of the numerator and denominator group combinations to calculate CI's for
  #   key_col = name of column of smallest grouping level
  #   n_bootreps = number of bootstrap reps to calculate CI's over
  #
  # Returns:
  #   A data frame with bootstrapped CI's for % cover and relative % cover. 
  #
  # NOTE: Not estimating CI's for trend -- just individually by site-survey.
  
  n_transects <- nrow(actual_dat) # number of transects surveyed in that site-survey
  
  samps <- boot_dat %>%
    dplyr::slice_sample(n = n_bootreps * n_transects, replace = TRUE) %>% # randomly sample from the parametric bootstrap samples
    as.data.frame() # can't add repeating vector column to tibble
  
  samps[, start_col:ncol(samps)][is.na(samps[, start_col:ncol(samps)])] <- 0 # replace NA's with zero
  samps$BootRep = rep(1:n_bootreps, each = n_transects) # add col of bootstrap replicate number; each replicate has same number of transects as original data for the site-survey
  
  # Add the actual data so actual % cover and relative % cover can be simultaneously calculated for each combination of numerator and denominator groups/levels. BootRep 0 is the actual data and will be excluded for CI estimates
  samps_df <- plyr::rbind.fill(samps, actual_dat %>% dplyr::mutate(BootRep = 0)) %>%
    tibble::rowid_to_column(var = "RowID") # after converting to long format, still need to be able to group data by row at times
  
  # This is the master data frame from which to calculate % cover at different grouping levels.
  samps_master <- samps_df %>%
    tidyr::gather(key = !!key_col, value = "Count", -Site, -SurvDate, -AdjTot, -BootRep, -RowID, -TransectSurveyID) %>% # convert to long format
    dplyr::left_join(groups_df, by = key_col)
  
  # Calculate % cover and RELATIVE % cover for certain combinations of numerator and denominator groups/levels
  
  boot_CIs_list <- apply(combos_df, 1, FUN = function(x) {
    
    tmp <- samps_master %>%
      dplyr::mutate(NumerLevel = get(x[["NumerGroup"]])) %>%
      dplyr::filter(if(any(unlist(x[["DenomGroup"]]) %in% unique(samps_master$Category))) Category %in% unlist(x[["DenomGroup"]]) else TRUE) %>% # filter for specified Categories unless DenomGroup == "AdjTot"
      dplyr::group_by(BootRep, TransectSurveyID, RowID, Category, NumerLevel) %>%
      dplyr::summarize(GroupCount = sum(Count), .groups = "drop")
    
    # For % cover, calculate per transect (RowID) and then take mean across all transects in the BootRep. The adjusted total count (denominator value) for any bootstrapped transect is the same as for the original data used to generate parametric bootstrap samples
    tmp2 <- tmp %>%
      dplyr::left_join(actual_dat[c("TransectSurveyID", "AdjTot")], by = "TransectSurveyID") %>%
      dplyr::mutate(PercCov = 100*(GroupCount / AdjTot)) %>%
      dplyr::group_by(BootRep, Category, NumerLevel) %>% # grouping by Category just to keep the column in the output
      dplyr::summarize(EstimCov = mean(PercCov), .groups = "drop")
    
    # Extract the means calculated for the original data (BootRep = 0) so can be added to the final data frame
    actual_df <- subset(tmp2, BootRep == 0) %>% 
      dplyr::select(-BootRep) 
    
    # Calculate % cover CI's
    out_df <- tmp2 %>% 
      dplyr::filter(BootRep != 0) %>% # remove the estimates from the actual data
      dplyr::group_by(Category, NumerLevel) %>%
      dplyr::summarize_at(vars(EstimCov), FuncPullQuant(), .groups = "drop") %>%
      dplyr::mutate(NumerGroup = x[["NumerGroup"]], DenomGroup = "TransectCount") %>%
      dplyr::left_join(actual_df, by = c("Category","NumerLevel")) # add in the mean % cov from actual data
    
    # # If denominator includes categories, add any missing zero-items
    # if(any(unlist(x[["DenomGroup"]]) %in% unique(samps_master$Category))) {
    #   perccov_template <- groups_df %>%
    #     dplyr::select(NumerLevel = x[["NumerGroup"]], Category) %>%
    #     dplyr::filter(Category %in% x[["DenomGroup"]]) %>%
    #     dplyr::mutate(NumerGroup = x[["NumerGroup"]], DenomGroup = "TransectCount") %>%
    #     dplyr::distinct()
    # 
    #   out_df %<>%
    #     right_join(perccov_template, by = c("Category", "NumerLevel", "NumerGroup", "DenomGroup"))
    #   
    # }
    
    # For relative % cover, sum(numerator across all transects)/sum(denom across all transects)
    if(any(unlist(x[["DenomGroup"]]) %in% unique(samps_master$Category))) {
      categtot_df <- tmp %>%
        dplyr::group_by(BootRep, Category) %>%
        dplyr::summarize(CategTot = sum(GroupCount), .groups = "drop")
      
      relcov_df <- tmp %>%
        dplyr::group_by(BootRep, Category, NumerLevel) %>%
        dplyr::summarize(SampGroupCount = sum(GroupCount), .groups = "drop") %>%
        dplyr::left_join(categtot_df, by = c("BootRep", "Category")) %>%
        dplyr::mutate(EstimCov = 100*(SampGroupCount/CategTot))
      
      # Extract the means calculated for the original data (BootRep = 0) so can be added to the final data frame
      actual_relcov_df <- subset(relcov_df, BootRep == 0) %>% 
        dplyr::select(Category, NumerLevel, EstimCov) 
      
      # Calculate % cover CI's
      relcov_df %<>% 
        dplyr::filter(BootRep != 0) %>% # remove the estimates from the actual data
        dplyr::group_by(Category, NumerLevel) %>%
        dplyr::summarize_at(vars(EstimCov), FuncPullQuant(), .groups = "drop") %>%
        dplyr::mutate(NumerGroup = x[["NumerGroup"]], DenomGroup = "CategoryCount") %>%
        dplyr::left_join(actual_relcov_df, by = c("Category", "NumerLevel")) # add in the mean % cov from actual data
      
      # # Add any missing zero-items
      # relcov_template <- groups_df %>% 
      #   dplyr::select(NumerLevel = x[["NumerGroup"]], Category) %>% 
      #   dplyr::filter(Category %in% x[["DenomGroup"]]) %>%
      #   dplyr::mutate(NumerGroup = x[["NumerGroup"]], DenomGroup = "CategoryCount") %>% 
      #   dplyr::distinct()
      # 
      # relcov_df %<>% 
      #   right_join(relcov_template, by = c("Category", "NumerLevel", "NumerGroup", "DenomGroup"))
      
      out_df <- rbind(out_df, relcov_df)
    }
    
    out_df %<>%
      dplyr::mutate_if(is.numeric, round, 2) %>%  # calculate quantiles and output in nice format
      dplyr::mutate(N = n_transects) %>%
      dplyr::select(Category, NumerGroup, NumerLevel, DenomGroup, N, EstimCov, everything())
    
    out_df$EstimCov[is.na(out_df$EstimCov)] <- 0
    
    # out_df[, 5:ncol(out_df)][is.na(out_df[, 5:ncol(out_df)])] <- 0
    
    out_df
    })
  
  boot_CIs_df <- as.data.frame(do.call("rbind", boot_CIs_list)) # combine list elements in a single data frame
  return(boot_CIs_df)
  }

FuncBootDrawsRS <- function(boot_dat, actual_dat, start_col, groups_df, combos_df, key_col, n_bootreps = 5000) {
  # Function to calculate bootstrap CI's for % cover and relative % cover, by site-survey (for RS sites). Some repetitiveness with FuncBootDraws
  #
  # Args:
  #   boot_list = data frame of the parametric bootstrap samples, each list element is a subsite with all its transects
  #   actual_dat = data frame of the original count data used to generate parametric bootstrap samples.
  #   start_col = starting column number for the count data
  #   groups_df = data frame assigning each items to higher level groups
  #   combos_df = data frame of the numerator and denominator group combinations to calculate CI's for
  #   key_col = name of column of smallest grouping level
  #   n_bootreps = number of bootstrap reps to calculate CI's over
  #
  # Returns:
  #   A data frame with bootstrapped CI's for RS % cover and relative % cover. 
  #
  # NOTE: Not estimating CI's for trend -- just individually by site-survey.
  
  n_transects <- 4 # number of transects surveyed in each subsite-survey
  subsites <- names(boot_dat)
  
  shuffle_samples <- lapply(boot_dat, FUN = function(x) {
    shuffle_dat <- x[sample(nrow(x), replace = F),] # shuffle the bootstrap samples
    shuffle_dat$BootSetID = rep(1:(nrow(shuffle_dat)/n_transects), each = n_transects)
    grouped_shuffle <- shuffle_dat %>% tidyr::nest(BootTransects = c(-Site, -SurvDate, -BootSetID)) %>%
      dplyr::select(-BootSetID) # the list-col now has 4 random parametric bootstraps of transects for that site-survey
  }) %>% 
    plyr::rbind.fill()
  
  samps <- shuffle_samples %>%
    dplyr::slice_sample(n = n_bootreps * length(subsites), replace = TRUE) %>% # randomly sample from the parametric bootstrap samples
    as.data.frame() # can't add repeating vector column to tibble
  
  samps$BootRep = rep(1:n_bootreps, each = length(subsites)) # add col of bootstrap replicate number; each replicate has same number of subsites as original data for the RSS-survey
  samps %<>% 
    tidyr::unnest(cols = BootTransects) %>%
    dplyr::select(Site, SurvDate, BootRep, TransectSurveyID, AdjTot, everything())
  samps[, (start_col+1):ncol(samps)][is.na(samps[, (start_col+1):ncol(samps)])] <- 0 # replace NA's with zero
  
  # Add the actual data so actual % cover and relative % cover can be simultaneously calculated for each combination of numerator and denominator groups/levels. BootRep 0 is the actual data and will be excluded for CI estimates
  samps_df <- plyr::rbind.fill(samps, actual_dat  %>% dplyr::mutate(BootRep = 0)) %>%
    tibble::rowid_to_column(var = "RowID") # after converting to long format, still need to be able to group data by row at times
  
  # This is the master data frame from which to calculate % cover at different grouping levels.
  samps_master <- samps_df %>%
    tidyr::gather(key = !!key_col, value = "Count", -Site, -SurvDate, -AdjTot, -BootRep, -RowID, -TransectSurveyID) %>% # convert to long format
    dplyr::left_join(groups_df, by = key_col)
  
  # Calculate % cover and RELATIVE % cover for certain combinations of numerator and denominator groups/levels
  
  boot_CIs_list <- apply(combos_df, 1, FUN = function(x) {
    
    tmp <- samps_master %>%
      dplyr::mutate(NumerLevel = get(x[["NumerGroup"]])) %>%
      dplyr::filter(if(any(unlist(x[["DenomGroup"]]) %in% unique(samps_master$Category))) Category %in% unlist(x[["DenomGroup"]]) else TRUE) %>% # filter for specified Categories unless DenomGroup == "AdjTot"
      dplyr::group_by(BootRep, Site, SurvDate, TransectSurveyID, RowID, Category, NumerLevel) %>%
      dplyr::summarize(GroupCount = sum(Count), .groups = "drop")
    
    # For % cover, calculate per RowID and then take mean across all RowID in the BootRep. The adjusted total count (denominator value) for any bootstrapped transect is the same as for the original data used to generate parametric bootstrap samples
    tmp2 <- tmp %>%
      dplyr::left_join(actual_dat[c("TransectSurveyID", "AdjTot")], by = "TransectSurveyID") %>%
      dplyr::mutate(PercCov = 100*(GroupCount / AdjTot)) %>%
      dplyr::group_by(BootRep, Category, NumerLevel) %>% # grouping by Category just to keep the column in the output
      dplyr::summarize(EstimCov = mean(PercCov), .groups = "drop")
    
    # Extract the means calculated for the original data (BootRep = 0) so can be added to the final data frame
    actual_df <- subset(tmp2, BootRep == 0) %>% 
      dplyr::select(-BootRep) 
    
    # Calculate % cover CI's
    out_df <- tmp2 %>% 
      dplyr::filter(BootRep != 0) %>% # remove the estimates from the actual data
      dplyr::group_by(Category, NumerLevel) %>%
      dplyr::summarize_at(vars(EstimCov), FuncPullQuant(), .groups = "drop") %>%
      dplyr::mutate(NumerGroup = x[["NumerGroup"]], DenomGroup = "TransectCount") %>%
      dplyr::left_join(actual_df, by = c("Category","NumerLevel")) # add in the mean % cov from actual data
    
    # # If denominator includes categories, add any missing zero-items
    # if(any(unlist(x[["DenomGroup"]]) %in% unique(samps_master$Category))) {
    #   perccov_template <- groups_df %>% 
    #     dplyr::select(NumerLevel = x[["NumerGroup"]], Category) %>% 
    #     dplyr::filter(Category %in% x[["DenomGroup"]]) %>%
    #     dplyr::mutate(NumerGroup = x[["NumerGroup"]], DenomGroup = "TransectCount") %>%
    #     dplyr::distinct()
    #   
    #   out_df %<>% 
    #     right_join(perccov_template, by = c("Category", "NumerLevel", "NumerGroup", "DenomGroup"))
    # }
    
    # For relative % cover, sum(numerator across all transects)/sum(denom across all transects)
    if(any(unlist(x[["DenomGroup"]]) %in% unique(samps_master$Category))) {
      categtot_df <- tmp %>%
        dplyr::group_by(BootRep, Category) %>%
        dplyr::summarize(CategTot = sum(GroupCount), .groups = "drop")
      
      relcov_df <- tmp %>%
        dplyr::group_by(BootRep, Category, NumerLevel) %>%
        dplyr::summarize(SampGroupCount = sum(GroupCount), .groups = "drop") %>%
        dplyr::left_join(categtot_df, by = c("BootRep", "Category")) %>%
        dplyr::mutate(EstimCov = SampGroupCount/CategTot)
      
      # Extract the means calculated for the original data (BootRep = 0) so can be added to the final data frame
      actual_relcov_df <- subset(relcov_df, BootRep == 0) %>% 
        dplyr::select(Category, NumerLevel, EstimCov) 
      
      # Calculate % cover CI's
      relcov_df %<>% 
        dplyr::filter(BootRep != 0) %>% # remove the estimates from the actual data
        dplyr::group_by(Category, NumerLevel) %>%
        dplyr::summarize_at(vars(EstimCov), FuncPullQuant(), .groups = "drop") %>%
        dplyr::mutate(NumerGroup = x[["NumerGroup"]], DenomGroup = "CategoryCount") %>%
        dplyr::left_join(actual_relcov_df, by = c("Category", "NumerLevel")) # add in the mean % cov from actual data
      
      # # Add any missing zero-items
      # relcov_template <- groups_df %>% 
      #   dplyr::select(NumerLevel = x[["NumerGroup"]], Category) %>% 
      #   dplyr::filter(Category %in% x[["DenomGroup"]]) %>%
      #   dplyr::mutate(NumerGroup = x[["NumerGroup"]], DenomGroup = "CategoryCount") %>% 
      #   dplyr::distinct()
      # 
      # relcov_df %<>% 
      #   right_join(relcov_template, by = c("Category", "NumerLevel", "NumerGroup", "DenomGroup"))
      
      out_df <- rbind(out_df, relcov_df)
    }
    
    out_df %<>%
      dplyr::mutate_if(is.numeric, round, 2) %>%  # calculate quantiles and output in nice format
      dplyr::mutate(N = length(subsites)) %>%
      dplyr::select(Category, NumerGroup, NumerLevel, DenomGroup, N, EstimCov, everything())
    out_df$EstimCov[is.na(out_df$EstimCov)] <- 0
    
    # out_df[, 5:ncol(out_df)][is.na(out_df[, 5:ncol(out_df)])] <- 0
    
    out_df
  })
  
  boot_CIs_df <- as.data.frame(do.call("rbind", boot_CIs_list)) # combine list elements in a single data frame
  return(boot_CIs_df)
}

FuncCorals <- function(filenam, sitesfilenam = NULL, out_prefix) {

# ### IMPORT AND FORMAT DATA ----
coral <- vroom::vroom(filenam)
if(!is.null(sitesfilenam)) tbl_link <- vroom::vroom(sitesfilenam)

# #  <<<<<<<< TESTING >>>>>>>>>>>>>> ----
# coral <- vroom::vroom("Data_LOCAL_ONLY/BUIS_CoralVideo Summary by Transect.csv") # <- vroom::vroom("Data_LOCAL_ONLY/demo_CoralDat.csv")
# out_prefix="test"
# tbl_link <- vroom::vroom("Data_LOCAL_ONLY/SFCN_CoralSites.csv")

# Format columns ----
coral$Date <- lubridate::mdy(coral$Date)
coral %<>%
  dplyr::rename(SurvDate = Date, Subcategory = SubCategory, Taxon = TaxonCode, CountOfTaxon = CountOfTaxonCode) %>%
  dplyr::mutate(TransectSurveyID = paste0(Site, "_", Transect, "_", SurvDate))  # create unique transect-specific EventID

# If the file doesn't already have a Reporting Site column, get the information from the link table
if(!"ReportingSite" %in% colnames(coral) & exists("tbl_link")) {
  coral %<>% dplyr::left_join(subset(tbl_link, select=-c(IsActive)), by = c("ParkCode", "Site"))
}

# If the file doesn't already have an Activity Status column, get the information from the link table
if(!"IsActive" %in% colnames(coral) & exists("tbl_link")) {
  coral %<>% dplyr::left_join(subset(tbl_link, select=c(ParkCode, Site, IsActive)), by = c("ParkCode", "Site"))
}

# Rename unknown group classifications to "UNK" or (for Functional Group) "UNGROUPED"--add to list as necessary
coral$Category[coral$Category %in% c("UNKNOWN") | is.na(coral$Category)] <- "UNK"
coral$Taxon[coral$Taxon %in% c("No Taxon") | is.na(coral$Taxon)] <- "UNK"
coral$FunctionalGroup[is.na(coral$FunctionalGroup)] <- paste(coral$Category[is.na(coral$FunctionalGroup)], "UNGROUPED", sep = "_")

### POPULATE LIST OF WARNINGS ----
warn_list <- sapply(c("AltPurp", "UNKTaxon", "BleachCode", "TransCount", "TaxonCountDiv", "TaxonCount", "EquipShadowCount", "UNKCount"), function(x) NULL)

warn_list$AltPurp <- coral %>%
  dplyr::arrange(SurvDate) %>%
  dplyr::select(TripName, Purpose) %>%
  dplyr::filter(!Purpose %in% c("Annual", "Episodic")) %>%
  distinct() %>%
  dplyr::rename("Trip" = TripName)

warn_list$UNKTaxon <- coral %>%
  dplyr::select(TransectSurveyID, Taxon) %>%
  dplyr::filter(Taxon == "UNK") %>%
  dplyr::rename("Site_Transect_SurveyDate" = "TransectSurveyID")

warn_list$BleachCode <- coral %>%
  dplyr::select(TransectSurveyID, Category, BleachingCode) %>%
  dplyr::filter(Category != "CORAL" & !is.na(BleachingCode)) %>%
  dplyr::rename("Site_Transect_SurveyDate" = "TransectSurveyID")

warn_list$TransCount <- coral %>%
  dplyr::select(Site, Transect, SurvDate, Purpose) %>%
  distinct() %>%
  dplyr::group_by(Site, SurvDate, Purpose) %>%
  dplyr::arrange(Site, SurvDate, Purpose) %>%
  dplyr::summarize(NumTransects = n(), .groups = "drop") %>%
  dplyr::filter(Purpose %in% c("Annual", "Episodic") & !NumTransects %in% c(4, 20)) %>% # 'off' transect counts only matters for annual or episodic surveys
  dplyr::mutate(SurvDate = lubridate::ymd(SurvDate)) %>% # renderTable does not play nice with dates--need to format as character
  dplyr::rename("Survey Date" = "SurvDate")

Err_counts <- coral %>%
  dplyr::group_by(TransectSurveyID) %>%
  dplyr::summarize(TotPoints = sum(CountOfTaxon, na.rm=TRUE),
            EquipShadowPoints = sum(CountOfTaxon[Category %in% c("EQUIP", "SHADOW")]),
            UNKPoints = sum(CountOfTaxon[Category %in% c("UNK", "UNKNOWN") | is.na(Category)]), .groups = "drop")

warn_list$TaxonCountDiv <- Err_counts %>%
  dplyr::mutate(ModOut = TotPoints %% 10) %>%
  dplyr::filter(ModOut != 0) %>%
  dplyr::select(TransectSurveyID, TotPoints) %>%
  dplyr::rename("Site_Transect_SurveyDate" = "TransectSurveyID", "Total Points" = TotPoints)

warn_list$TaxonCount <- Err_counts %>%
  dplyr::filter(TotPoints < 200 | TotPoints > 480) %>%
  dplyr::select(TransectSurveyID, TotPoints) %>%
  dplyr::rename("Site_Transect_SurveyDate" = "TransectSurveyID", "Total Points" = TotPoints)

warn_list$EquipShadowCount <- Err_counts %>%
  dplyr::mutate(PercEquipShadow = (EquipShadowPoints / TotPoints)*100) %>%
  dplyr::filter(PercEquipShadow > 5.0) %>%
  dplyr::select(TransectSurveyID, PercEquipShadow) %>%
  dplyr::arrange(desc(PercEquipShadow)) %>%
  dplyr::rename("Site_Transect_SurveyDate" = "TransectSurveyID", "% Equip or Shadow" = PercEquipShadow)

warn_list$UNKCount <- Err_counts %>%
  dplyr::mutate(PercUNK = (UNKPoints / TotPoints)*100) %>%
  dplyr::filter(PercUNK > 5.0) %>%
  dplyr::select(TransectSurveyID, PercUNK) %>%
  dplyr::arrange(desc(PercUNK)) %>%
  dplyr::rename("Site_Transect_SurveyDate" = "TransectSurveyID", "% Unknown" = PercUNK)


# Final formatting steps ----

# Account for different classification of bleaching prior to 10/1/05
coral$BleachingCode[is.na(coral$BleachingCode) & coral$Category =="CORAL" & coral$SurvDate < "2005-10-01"] <- "NoData"
coral$BleachingCode[is.na(coral$BleachingCode) & coral$Category =="CORAL" & coral$SurvDate >= "2005-10-01"] <- "UNBL"
coral$BleachingCode[coral$BleachingCode == "BL" & coral$Category =="CORAL"] <- "BL1"
coral$BleachingCode[coral$Category != "CORAL"] <- NA

coral%<>%
  dplyr::filter(Purpose %in% c("Annual", "Episodic")) %>% # only include annual and episodic data
  dplyr::select(ParkCode, ReportingSite, ReportingSiteName, Site, Transect, Latitude, Longitude, Year, SurvDate, TransectSurveyID, TripName, Purpose, IsActive, Category, Subcategory, FunctionalGroup, Taxon, BleachingCode, CountOfTaxon) %>%
  dplyr::arrange(Site, SurvDate, Transect)

# STOP-Error check--make sure all records have associated Reporting Site and Activity Status 
missing_RS <- coral[is.na(coral$ReportingSite), ] %>%
  distinct()
if(nrow(missing_RS) > 0) stop("These data records do not have reporting site information:\n", paste(capture.output(print(missing_RS)), collapse = "\n")) else cat("OK >>> All records have reporting site info")
missing_IsActive <- coral[is.na(coral$IsActive), ] %>%
  distinct()
if(nrow(missing_IsActive) > 0) stop("These data records do not have activity status information:\n", paste(capture.output(print(missing_IsActive)), collapse = "\n")) else cat("OK >>> All records have activity status info")

# Uncomment (i.e., remove the '#') the line below to output a CSV of the cleaned data
# write.csv(coral, paste0(out_prefix, "_cleandat.csv"), row.names = FALSE) 

### DATA FOR MAP ----
mapdat <- coral %>%
  dplyr::select(ReportingSiteName, ReportingSite, Site, Transect, Latitude, Longitude, Year, IsActive) %>% # pull only the data relevant for mapping by Site
  dplyr::group_by(Site) %>%
  dplyr::summarize( # site-level summary stats for mapping
    ReportingSiteName = unique(ReportingSiteName),
    ReportingSite = unique(ReportingSite),
    IsActive = unique(IsActive),
    NumTransects = length(unique(Transect)),
    MinYr = min(Year, na.rm = TRUE),
    MaxYr = max(Year, na.rm = TRUE),
    MedLat = median(Latitude, na.rm = TRUE),
    MedLong = median(Longitude, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    PopText = paste0(Site, " (", MinYr, " - ", MaxYr, ")"))

### BOOTSTRAPPING ESTIMATES ----
# Format data for bootstrapping
boot_master_dat <- dplyr::filter(coral, !Category %in% c("EQUIP", "SHADOW")) %>%  # remove EQUIP & SHADOW data, for all analyses
  dplyr::select(ReportingSite, Site, Transect, SurvDate, TransectSurveyID, Category, FunctionalGroup, Taxon, CountOfTaxon, BleachingCode)

# Make sure each taxon belongs to only one Functional Group and only one Category
grps_df <- unique(boot_master_dat[c("Category", "FunctionalGroup", "Taxon")])
if(length(unique(grps_df$Taxon)) != nrow(grps_df)) stop("EACH TAXON CAN ONLY BE ASSIGNED TO ONE FUNCTIONAL GROUP AND TO ONE CATEGORY")

# Create data frame with one row per TransectSurveyID and with Taxa as columns. Values are the % cover for each taxon for that TransectSurveyID 
boot_site_dat <- boot_master_dat %>%
  dplyr::group_by(Site, SurvDate, TransectSurveyID, Taxon) %>%
  dplyr::summarize(Count = sum(CountOfTaxon), .groups = "drop") %>% # for each transect-survey, count hits per taxon
  tidyr::spread(key = Taxon, value = Count) %>% # each taxon becomes a column
  dplyr::mutate(AdjTot = rowSums(across(4:ncol(.)), na.rm = T)) %>% # calculate total hits per transect-survey (already excludes shadow and equipment points)
  dplyr::select(Site, SurvDate, TransectSurveyID, AdjTot, everything()) # rearrange cols

# Calculate bootstrapped cover CI's for these combinations. 
combos_df <- tibble(
  NumerGroup = c("Category", "FunctionalGroup", "Taxon"), 
  DenomGroup = list(c("AdjTot"), c("ALGAE", "CORAL"), c("ALGAE", "CORAL", "GORGO", "SPONGE")))

# Calculate bootstrapped bleaching CI's for these combinations.

bleach_grps_df <- data.frame (Category = c("CORAL", "CORAL", "CORAL", "CORAL", "CORAL", "CORAL", "NONCORAL"), BleachingCode = c("UNBL", "BL1", "BL2", "BL3", "BL4", "NoData", "NONCORAL"))

bleach_combos_df <- tibble(
  NumerGroup = c("BleachingCode"), 
  DenomGroup = c("CORAL"))

## SITE-LEVEL ESTIMATES ----

# Site-level cover estimates ----
incProgress(1/5, detail = paste0("...site-level % cover estimates"))
# First, work with taxon count data
# For each transect-survey, generate parametric bootstrap samples (default is N = 1M) using the probabilities from the collected data. Run one site-survey at a time to avoid computer memory problems
site_CIs_list <- apply(unique(boot_site_dat[c( "Site", "SurvDate")]), 1, FUN = function(x) {
  
incProgress(1/5, detail = paste0("...site-level % cover estimates for ", x[["Site"]], x[["SurvDate"]]))
  cat("site-level % cover estimates")
  cat(x[["Site"]], x[["SurvDate"]])
  subdat <- boot_site_dat %>%
    dplyr::filter(Site == x[["Site"]] & SurvDate == x[["SurvDate"]]) %>% # pull data for the specific site and survey date
    dplyr::select(where(~!all(is.na(.x)))) # get rid of cols with no data
  start_col = 5
  
  # Generate parametric bootstrap samples
  boot_samples <- FuncBootSamples(dat = subdat, start_col = start_col)
  boot_samples %<>% dplyr::mutate(across(where(is.numeric), as.numeric)) # this addresses an unusual problem of columns being nested matrices
  
  # Grab the actual data to send to boot draws function
  actual_dat <- subdat
  
  actual_dat[, start_col:ncol(actual_dat)][is.na(actual_dat[,  start_col:ncol(actual_dat)])] <- 0 # and then replace the remaining NA's with 0

  # For various combinations of numerator groups and denominator groups, calculate quantiles from bootstrap distribution
  boot_site_estimates <- FuncBootDraws(boot_dat = boot_samples, start_col = start_col, actual_dat = actual_dat, groups_df = grps_df, combos_df = combos_df, key_col = "Taxon")
  
  cbind("SiteScale" = "Site", "ReportingSite" = mapdat$ReportingSite[mapdat$Site == x[["Site"]]], "Site" = x[["Site"]], "RSS" = x[["Site"]], "SurvDate" = x[["SurvDate"]], boot_site_estimates)
  })
site_CIs_df <- as.data.frame(do.call("rbind", site_CIs_list))
site_CIs_df$SurvDate <- lubridate::ymd(site_CIs_df$SurvDate)

# Site-level coral bleaching estimates ----
incProgress(2/5, detail = paste0("...site-level coral bleaching estimates"))
# Seems like most of the bleaching is in the Orbicella functional group. Within that, OFAV is the one least affected, there are 3 taxa highly affected. Raw data plots will show this.

# Now work with coral bleach data
boot_bleach_dat <- boot_master_dat %>%
  dplyr::group_by(Site, SurvDate, TransectSurveyID, BleachingCode) %>%
  dplyr::summarize(Count = sum(CountOfTaxon), .groups = "drop") %>% 
  tidyr::replace_na(list(BleachingCode = "NONCORAL")) %>%
  tidyr::spread(key = BleachingCode, value = Count) %>% 
  dplyr::mutate(AdjTot = rowSums(across(4:ncol(.)), na.rm = T)) %>%
  dplyr::select(Site, SurvDate, TransectSurveyID, AdjTot, everything())

# For each transect-survey, generate parametric bootstrap samples 
site_bleach_CIs_list <- apply(unique(boot_bleach_dat[c("Site", "SurvDate")]), 1, FUN = function(x) {
  
incProgress(2/5, detail = paste0("...site-level coral bleaching estimates for ", x[["Site"]], x[["SurvDate"]]))
  cat("site-level coral bleaching estimates")
  cat(x[["Site"]], x[["SurvDate"]])
  bleach_subdat <- boot_bleach_dat %>%
    dplyr::filter(Site == x[["Site"]] & SurvDate == x[["SurvDate"]]) %>% # pull data for the specific site and survey date
    dplyr::select(where(~!all(is.na(.x)))) # get rid of cols with no data
  
  bleach_start_col = 5
  
  # Generate parametric bootstrap samples
  boot_bleach_samples <- FuncBootSamples(dat = bleach_subdat, start_col = bleach_start_col)
  boot_bleach_samples %<>% mutate(across(where(is.numeric), as.numeric)) # this addresses an unusual problem of columns being nested matrices
  
  # Grab the actual data to send to boot draws function
  actual_bleach_dat <- bleach_subdat
  actual_bleach_dat[, bleach_start_col:ncol(actual_bleach_dat)][is.na(actual_bleach_dat[,  bleach_start_col:ncol(actual_bleach_dat)])] <- 0 # and then replace the remaining NA's with 0
  
  # For various combinations of numerator groups and denominator groups, calculate quantiles from bootstrap distribution
  boot_bleach_estimates <- FuncBootDraws(boot_dat = boot_bleach_samples, start_col = bleach_start_col, actual_dat = actual_bleach_dat, groups_df = bleach_grps_df, combos_df = bleach_combos_df, key_col = "BleachingCode")
  cbind("SiteScale" = "Site", "ReportingSite" = mapdat$ReportingSite[mapdat$Site == x[["Site"]]], "Site" = x[["Site"]], "RSS" = x[["Site"]], "SurvDate" = x[["SurvDate"]], boot_bleach_estimates)
  })
site_bleach_CIs_df <- as.data.frame(do.call("rbind", site_bleach_CIs_list))
site_bleach_CIs_df$SurvDate <- lubridate::ymd(site_bleach_CIs_df$SurvDate)

## REPORTING SITE-LEVEL ESTIMATES ----
# For each reporting site, these are the only survey dates to use for reporting site estimates (because all sites were surveyed in these years)
RS_estim_dates <- boot_master_dat %>%
  dplyr::select(ReportingSite, Site, SurvDate) %>% 
  dplyr::filter(ReportingSite != Site) %>% 
  dplyr::distinct() %>%
  dplyr::group_by(ReportingSite) %>% 
  dplyr::mutate(NumSub = length(unique(Site))) %>%    
  dplyr::group_by(ReportingSite, SurvDate) %>% 
  dplyr::mutate(NumSubSurveyed = length(unique(Site))) %>%
  dplyr::filter(NumSub == NumSubSurveyed) %>%
  dplyr::select(-NumSub, -NumSubSurveyed)

# Reporting site-level cover estimates ----
incProgress(3/5, detail = paste0("...`reporting site`-level % cover estimates"))
if(nrow(RS_estim_dates) > 0) {
# For each combination of reporting site and survey date
RS_CIs_list <- apply(unique(RS_estim_dates[c( "ReportingSite", "SurvDate")]), 1, FUN = function(x) {
  
incProgress(3/5, detail = paste0("...`reporting site`-level % cover estimates for ", x[["ReportingSite"]], x[["SurvDate"]]))
  cat("reporting site-level % cover estimates")
  cat(x[["ReportingSite"]], x[["SurvDate"]])
  subsites <- RS_estim_dates %>% dplyr::filter(ReportingSite == x[["ReportingSite"]] & SurvDate == x[["SurvDate"]]) %>% dplyr::pull(Site) # these are the sites in that RS with that survey date
  RS_subdat <- boot_site_dat %>%
    dplyr::filter(Site %in% subsites & SurvDate == x[["SurvDate"]]) %>% # pull data for the specific site and survey date
    dplyr::select(where(~!all(is.na(.x)))) # get rid of cols with no data
  RS_start_col = 5
  
  # Generate parametric bootstrap samples
  boot_RS_samples <- FuncBootSamples(dat = RS_subdat, start_col = RS_start_col)
  boot_RS_samples %<>% mutate(across(where(is.numeric), as.numeric)) # this addresses an unusual problem of columns being nested matrices
  
  boot_RS_list <- split(boot_RS_samples, f = boot_RS_samples$Site) # split samples so each list element is a site
  
  # Grab the actual data to send to boot draws function
  actual_RS_dat <- RS_subdat
  
  actual_RS_dat[, RS_start_col:ncol(actual_RS_dat)][is.na(actual_RS_dat[,  RS_start_col:ncol(actual_RS_dat)])] <- 0 # and then replace the remaining NA's with 0
  
  # For various combinations of numerator groups and denominator groups, calculate RS quantiles from bootstrap distribution
  boot_RS_estimates <- FuncBootDrawsRS(boot_dat = boot_RS_list, start_col = RS_start_col, actual_dat = actual_RS_dat, groups_df = grps_df, combos_df = combos_df, key_col = "Taxon")
  
  cbind("SiteScale" = "ReportingSite", "ReportingSite" = x[["ReportingSite"]], "Site" = NA, "RSS" = x[["ReportingSite"]], "SurvDate" = x[["SurvDate"]], boot_RS_estimates)
})

RS_CIs_df <- as.data.frame(do.call("rbind", RS_CIs_list))
RS_CIs_df$SurvDate <- lubridate::ymd(RS_CIs_df$SurvDate)

# Reporting site-level coral bleaching estimates ----
incProgress(4/5, detail = paste0("...`reporting site`-level coral bleaching estimates"))
# For each reporting site-survey, generate parametric bootstrap samples 
RS_bleach_CIs_list <- apply(unique(RS_estim_dates[c( "ReportingSite", "SurvDate")]), 1, FUN = function(x) {
  
incProgress(4/5, detail = paste0("...`reporting site`-level coral bleaching estimates for ", x[["ReportingSite"]], x[["SurvDate"]]))
  cat("reporting site-level coral bleaching estimates")
  cat(x[["ReportingSite"]], x[["SurvDate"]])
  subsites <- RS_estim_dates %>% dplyr::filter(ReportingSite == x[["ReportingSite"]] & SurvDate == x[["SurvDate"]]) %>% dplyr::pull(Site) # these are the sites in that RS with that survey date
  
  RS_bleach_subdat <- boot_bleach_dat %>%
    dplyr::filter(Site %in% subsites & SurvDate == x[["SurvDate"]]) %>% # pull data for the specific site and survey date
    dplyr::select(where(~!all(is.na(.x)))) # get rid of cols with no data
  
  RS_bleach_start_col = 5
  
  # Generate parametric bootstrap samples
  boot_RS_bleach_samples <- FuncBootSamples(dat = RS_bleach_subdat, start_col = RS_bleach_start_col)
  boot_RS_bleach_samples %<>% mutate(across(where(is.numeric), as.numeric)) # this addresses an unusual problem of columns being nested matrices
  
  boot_RS_bleach_list <- split(boot_RS_bleach_samples, f = boot_RS_bleach_samples$Site) # split samples so each list element is a site
  
  # Grab the actual data to send to boot draws function
  actual_RS_bleach_dat <- RS_bleach_subdat
  actual_RS_bleach_dat[, RS_bleach_start_col:ncol(actual_RS_bleach_dat)][is.na(actual_RS_bleach_dat[,  RS_bleach_start_col:ncol(actual_RS_bleach_dat)])] <- 0 # and then replace the remaining NA's with 0
  
  # For various combinations of numerator groups and denominator groups, calculate quantiles from bootstrap distribution
  boot_RS_bleach_estimates <- FuncBootDrawsRS(boot_dat = boot_RS_bleach_list, start_col = RS_bleach_start_col, actual_dat = actual_RS_bleach_dat, groups_df = bleach_grps_df, combos_df = bleach_combos_df, key_col = "BleachingCode") 
  
  cbind("SiteScale" = "ReportingSite", "ReportingSite" = x[["ReportingSite"]], "Site" = NA, "RSS" = x[["ReportingSite"]], "SurvDate" = x[["SurvDate"]], boot_RS_bleach_estimates)
})

RS_bleach_CIs_df <- as.data.frame(do.call("rbind", RS_bleach_CIs_list))
RS_bleach_CIs_df$SurvDate <- lubridate::ymd(RS_bleach_CIs_df$SurvDate)
} else {
  RS_CIs_df <- RS_bleach_CIs_df <- NULL
}

### ADD CORAL NUMBERS TO MAPDAT ----

# For each site and survey, the mean % coral
CoralPercCov <- site_CIs_df %>% 
  dplyr::filter(NumerGroup == "Category" & NumerLevel == "CORAL" & DenomGroup == "TransectCount") %>% 
  dplyr::select(Site, SurvDate, EstimCov) %>%
  dplyr::mutate('%Coral' = EstimCov,
                '%NonCoral' = 100 - EstimCov) %>%
  dplyr::select(-EstimCov)

BleachRelCov <- site_bleach_CIs_df %>%
  dplyr::filter(NumerGroup == "BleachingCode" & DenomGroup == "CategoryCount") %>%
  dplyr::select(Site, SurvDate, NumerLevel, EstimCov) %>%
  # dplyr::mutate(EstimCov = round(EstimCov * 100)) %>%
  tidyr::spread(key = NumerLevel, value = EstimCov, fill = 0) %>%
  full_join(CoralPercCov, by = c("Site", "SurvDate"))

mapdat2 <- mapdat %>%
  full_join(BleachRelCov, by = "Site") %>%
  left_join(coral[c("Site", "SurvDate", "Purpose")] %>% distinct(), by = c("Site", "SurvDate")) %>% # add Purpose info
  mutate(PopText = paste0(PopText, ": ", SurvDate, "(", Purpose, ")")) %>%
  dplyr::select("ReportingSite", "ReportingSiteName", "Site", "IsActive", "NumTransects", "MinYr", "MaxYr", "MedLat", "MedLong", "SurvDate", "Purpose", "PopText", "%Coral", "%NonCoral", "UNBL", "BL1", "BL2", "BL3", "BL4", "NoData") # order columns

cover_CIs_df <- site_CIs_df %>% 
    left_join(coral[c("Site", "SurvDate", "IsActive", "Purpose")] %>% distinct(), by = c("Site", "SurvDate")) %>%
  dplyr::select(SiteScale, ReportingSite, Site, RSS, IsActive, SurvDate, Purpose, N, everything())

if (!is.null(RS_CIs_df)) {
cover_CIs_df <- rbind(
  cover_CIs_df,
  RS_CIs_df %>%
    left_join(coral[c("ReportingSite", "SurvDate", "IsActive", "Purpose")] %>% distinct(), by = c("ReportingSite", "SurvDate")) %>%
  dplyr::select(SiteScale, ReportingSite, Site, RSS, IsActive, SurvDate, Purpose, N, everything())
)
}

bleach_CIs_df <- site_bleach_CIs_df %>% 
    left_join(coral[c("Site", "SurvDate", "IsActive", "Purpose")] %>% distinct(), by = c("Site", "SurvDate")) %>%
  dplyr::select(SiteScale, ReportingSite, Site, RSS, IsActive, SurvDate, Purpose, N, everything())

if(!is.null(RS_bleach_CIs_df)) {
  bleach_CIs_df <- rbind(
    bleach_CIs_df,
    RS_bleach_CIs_df %>%
    left_join(coral[c("ReportingSite", "SurvDate", "IsActive", "Purpose")] %>% distinct(), by = c("ReportingSite", "SurvDate")) %>%
  dplyr::select(SiteScale, ReportingSite, Site, RSS, IsActive, SurvDate, Purpose, N, everything())
  )
}

# Calculate % cover by category for each transect
transect_counts_df <- boot_master_dat %>%
  tidyr::expand(TransectSurveyID, Category) %>%
  left_join(boot_master_dat[c("TransectSurveyID", "Category", "CountOfTaxon")], by = c("TransectSurveyID", "Category")) %>%
  dplyr::mutate(CountOfTaxon = replace_na(CountOfTaxon, 0)) %>%
  group_by(TransectSurveyID) %>%
  dplyr::mutate(AdjTot = sum(CountOfTaxon)) %>%
  group_by(TransectSurveyID, Category) %>%
  dplyr::mutate(TransectCount = sum(CountOfTaxon),
                PercCov = round(100*(TransectCount/AdjTot), 1)) %>%
  dplyr::select(-CountOfTaxon) %>%
  dplyr::distinct()

out_list <- list(raw_dat = coral, # the raw data include equipment and shadow counts; only annual and episodic surveys included
                 map_dat = mapdat2,
                 groups_df = grps_df,
                 warn_list = warn_list,
                 transect_counts_df = transect_counts_df,
                 cover_CIs_df = cover_CIs_df,
                 bleach_CIs_df = bleach_CIs_df)
saveRDS(out_list, paste0(out_prefix, "_coralsummary.RDS"))
}

