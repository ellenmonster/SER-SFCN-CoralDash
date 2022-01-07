############################################################
### CORAL DASHBOARD DATA IMPORT, FORMAT AND CALCULATIONS ###
############################################################

### LOAD PACKAGES ----
# Look for packages on local machine, install if necessary, then load all
pkgList_pre <- c("magrittr", 
                 "plyr",
                 "lubridate",
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

FuncBootDraws <- function(boot_dat, actual_dat, start_col, groups_df, combos_df, key_col, n_bootreps = 10000) {
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
  samps <- plyr::rbind.fill(samps, actual_dat %>% dplyr::mutate(BootRep = 0)) %>%
    tibble::rowid_to_column(var = "RowID") # after converting to long format, still need to be able to group data by row at times
  
  # This is the master data frame from which to calculate % cover at different grouping levels.
  samps_master <- samps %>%
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
      dplyr::mutate(PercCov = GroupCount / AdjTot) %>%
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
    
    # If denominator includes categories, add any missing zero-items
    if(any(unlist(x[["DenomGroup"]]) %in% unique(samps_master$Category))) {
      perccov_template <- groups_df %>% 
        dplyr::select(NumerLevel = x[["NumerGroup"]], Category) %>% 
        dplyr::filter(Category %in% x[["DenomGroup"]]) %>%
        dplyr::mutate(NumerGroup = x[["NumerGroup"]], DenomGroup = "TransectCount")
      
      out_df %<>% 
        right_join(perccov_template, by = c("Category", "NumerLevel", "NumerGroup", "DenomGroup"))
    }
    
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
      
      # Add any missing zero-items
      relcov_template <- groups_df %>% 
        dplyr::select(NumerLevel = x[["NumerGroup"]], Category) %>% 
        dplyr::filter(Category %in% x[["DenomGroup"]]) %>%
        dplyr::mutate(NumerGroup = x[["NumerGroup"]], DenomGroup = "CategoryCount")
      
      relcov_df %<>% 
        right_join(relcov_template, by = c("Category", "NumerLevel", "NumerGroup", "DenomGroup"))
      
      out_df <- rbind(out_df, relcov_df)
    }
    
    out_df %<>%
      dplyr::mutate_if(is.numeric, round, 3) %>%  # calculate quantiles and output in nice format
      dplyr::select(Category, NumerGroup, NumerLevel, DenomGroup, EstimCov, everything())
    out_df$EstimCov[is.na(out_df$EstimCov)] <- 0
    
    # out_df[, 5:ncol(out_df)][is.na(out_df[, 5:ncol(out_df)])] <- 0
    
    out_df
  })
  
  boot_CIs_df <- as.data.frame(do.call("rbind", boot_CIs_list)) # combine list elements in a single data frame
  return(boot_CIs_df)
}

# FuncCorals <- function(filenam, sitesfilenam = NULL, out_prefix, run_parallel, set_cores = 1) {

### IMPORT AND FORMAT DATA ----
# coral <- read_csv(filenam)
# if(!is.null(sitesfilenam)) tbl_link <- read_csv(sitesfilenam)

#  <<<<<<<< TESTING >>>>>>>>>>>>>> ----
coral <- read_csv("Data_LOCAL_ONLY/BUIS_CoralVideo Summary by Transect.csv") # <- read_csv("Data_LOCAL_ONLY/demo_CoralDat.csv")
out_prefix="test"
# run_parallel=TRUE
# set_cores=11
tbl_link <- read_csv("Data_LOCAL_ONLY/SFCN_CoralSites.csv")

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

# ### POPULATE LIST OF WARNINGS ----
# Warn_list <- sapply(c("AltPurp", "UNKTaxon", "BleachCode", "TransCount", "TaxonCountDiv", "TaxonCount", "EquipShadowCount", "UNKCount"), function(x) NULL)
# 
# Warn_list$AltPurp <- coral %>%
#   dplyr::arrange(SurvDate) %>%
#   dplyr::select(TripName, Purpose) %>%
#   dplyr::filter(!Purpose %in% c("Annual", "Episodic")) %>%
#   distinct() %>%
#   dplyr::rename("Trip" = TripName)
# 
# Warn_list$UNKTaxon <- coral %>%
#   dplyr::select(TransectSurveyID, Taxon) %>%
#   dplyr::filter(Taxon == "UNK") %>%
#   dplyr::rename("Site_Transect_SurveyDate" = "TransectSurveyID")
# 
# Warn_list$BleachCode <- coral %>%
#   dplyr::select(TransectSurveyID, Category, BleachingCode) %>%
#   dplyr::filter(Category != "CORAL" & !is.na(BleachingCode)) %>%
#   dplyr::rename("Site_Transect_SurveyDate" = "TransectSurveyID")
# 
# Warn_list$TransCount <- coral %>%
#   dplyr::select(Site, Transect, SurvDate, Purpose) %>%
#   distinct() %>%
#   dplyr::group_by(Site, SurvDate, Purpose) %>%
#   dplyr::arrange(Site, SurvDate, Purpose) %>%
#   dplyr::summarize(NumTransects = n(), .groups = "drop") %>%
#   dplyr::filter(Purpose %in% c("Annual", "Episodic") & !NumTransects %in% c(4, 20)) %>% # 'off' transect counts only matters for annual or episodic surveys
#   dplyr::mutate(SurvDate = lubridate::ymd(SurvDate)) %>% # renderTable does not play nice with dates--need to format as character
#   dplyr::rename("Survey Date" = "SurvDate")
# 
# Err_counts <- coral %>%
#   dplyr::group_by(TransectSurveyID) %>%
#   dplyr::summarize(TotPoints = sum(CountOfTaxon, na.rm=TRUE),
#             EquipShadowPoints = sum(CountOfTaxon[Category %in% c("EQUIP", "SHADOW")]),
#             UNKPoints = sum(CountOfTaxon[Category %in% c("UNK", "UNKNOWN") | is.na(Category)]), .groups = "drop")
# 
# Warn_list$TaxonCountDiv <- Err_counts %>%
#   dplyr::mutate(ModOut = TotPoints %% 10) %>%
#   dplyr::filter(ModOut != 0) %>%
#   dplyr::select(TransectSurveyID, TotPoints) %>%
#   dplyr::rename("Site_Transect_SurveyDate" = "TransectSurveyID", "Total Points" = TotPoints)
# 
# Warn_list$TaxonCount <- Err_counts %>%
#   dplyr::filter(TotPoints < 200 | TotPoints > 480) %>%
#   dplyr::select(TransectSurveyID, TotPoints) %>%
#   dplyr::rename("Site_Transect_SurveyDate" = "TransectSurveyID", "Total Points" = TotPoints)
# 
# Warn_list$EquipShadowCount <- Err_counts %>%
#   dplyr::mutate(PercEquipShadow = (EquipShadowPoints / TotPoints)*100) %>%
#   dplyr::filter(PercEquipShadow > 5.0) %>%
#   dplyr::select(TransectSurveyID, PercEquipShadow) %>%
#   dplyr::arrange(desc(PercEquipShadow)) %>%
#   dplyr::rename("Site_Transect_SurveyDate" = "TransectSurveyID", "% Equip or Shadow" = PercEquipShadow)
# 
# Warn_list$UNKCount <- Err_counts %>%
#   dplyr::mutate(PercUNK = (UNKPoints / TotPoints)*100) %>%
#   dplyr::filter(PercUNK > 5.0) %>%
#   dplyr::select(TransectSurveyID, PercUNK) %>%
#   dplyr::arrange(desc(PercUNK)) %>%
#   dplyr::rename("Site_Transect_SurveyDate" = "TransectSurveyID", "% Unknown" = PercUNK)
# 

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
# %>%  mutate_if(is.character, as.factor)

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
# For each active reporting site, these are the start and end years for reporting site plots (because all subsites were surveyed in these years). So estimates at RS level should be limited to these years.
RS_keepyrs <- mapdat %>%
  dplyr::select(ReportingSite, IsActive, MinYr, MaxYr) %>%
  filter(IsActive == TRUE) %>%
  dplyr::group_by(ReportingSite) %>%
  dplyr::summarize(StartYr = max(MinYr),
                   EndYr = min(MaxYr), .groups = "drop") %>%
  dplyr::left_join(mapdat[c("ReportingSite", "Site")], by = "ReportingSite") %>%
  dplyr::select(ReportingSite, Site, StartYr, EndYr)

# Format data for bootstrapping
boot_master_dat <- dplyr::filter(coral, !Category %in% c("EQUIP", "SHADOW")) %>%  # remove EQUIP & SHADOW data, for all analyses
  dplyr::select(Site, Transect, SurvDate, TransectSurveyID, Category, FunctionalGroup, Taxon, CountOfTaxon, BleachingCode)

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

# COVER ESTIMATES ----
# First, work with taxon count data
# For each transect-survey, generate parametric bootstrap samples (default is N = 1M) using the probabilities from the collected data. Run one site-survey at a time to avoid computer memory problems
site_CIs_list <- apply(unique(boot_site_dat[c( "Site", "SurvDate")]), 1, FUN = function(x) {
  cat(x)
  subdat <- boot_site_dat %>%
    dplyr::filter(Site == x[["Site"]] & SurvDate == x[["SurvDate"]]) %>% # pull data for the specific site and survey date
    dplyr::select(where(~!all(is.na(.x)))) # get rid of cols with no data
  start_col = 5
  
  # Generate parametric bootstrap samples
  boot_samples <- FuncBootSamples(dat = subdat, start_col = start_col)
  boot_samples %<>% mutate(across(where(is.numeric), as.numeric)) # this addresses an unusual problem of columns being nested matrices
  
  # Grab the actual data to send to boot draws function
  actual_dat <- subdat
  
  actual_dat[, start_col:ncol(actual_dat)][is.na(actual_dat[,  start_col:ncol(actual_dat)])] <- 0 # and then replace the remaining NA's with 0
  
  # Calculate bootstrapped CI's for these combinations. 
  combos_df <- tibble(
    NumerGroup = c("Category", "FunctionalGroup", "Taxon"), 
    DenomGroup = list(c("AdjTot"), c("ALGAE", "CORAL"), c("ALGAE", "CORAL", "GORGO", "SPONGE")))

  # For various combinations of numerator groups and denominator groups, calculate quantiles from bootstrap distribution
  boot_site_estimates <- FuncBootDraws(boot_dat = boot_samples, start_col = start_col, actual_dat = actual_dat, groups_df = grps_df, combos_df = combos_df, key_col = "Taxon", n_bootreps = 100) # <<<<<<<<<<<<<<< SMALL NUMBER FOR TESTING
  
  boot_dat = boot_samples; start_col = start_col; actual_dat = actual_dat; groups_df = grps_df; combos_df = combos_df; key_col = "Taxon"; n_bootreps = 100 # <<<<<<<<<<<<<<<<<<<<<<< TESTING FUNCBOOTDRAWS
  
  final_df <- cbind("Site" = x[["Site"]], "SurvDate" = x[["SurvDate"]], boot_site_estimates)
  })

site_CIs_df <- as.data.frame(do.call("rbind", site_CIs_list))

# CORAL BLEACHING ESTIMATES ---
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
  bleach_subdat <- boot_bleach_dat %>%
    dplyr::filter(Site == x[["Site"]] & SurvDate == x[["SurvDate"]]) %>% # pull data for the specific site and survey date
    dplyr::select(where(~!all(is.na(.x)))) # get rid of cols with no data
  
  bleach_start_col = 5
  
  # Generate parametric bootstrap samples
  boot_bleach_samples <- FuncBootSamples(dat = bleach_subdat, start_col = bleach_start_col, n_bootsamples = 100000)
  boot_bleach_samples %<>% mutate(across(where(is.numeric), as.numeric)) # this addresses an unusual problem of columns being nested matrices
  
  # Grab the actual data to send to boot draws function
  actual_bleach_dat <- bleach_subdat
  actual_bleach_dat[, bleach_start_col:ncol(actual_bleach_dat)][is.na(actual_bleach_dat[,  bleach_start_col:ncol(actual_bleach_dat)])] <- 0 # and then replace the remaining NA's with 0
  
  bleach_grps_df <- data.frame (Category = c("CORAL", "CORAL", "CORAL", "CORAL", "CORAL", "CORAL", "NONCORAL"), BleachingCode = c("UNBL", "BL1", "BL2", "BL3", "BL4", "NoData", "NONCORAL"))
  
  bleach_combos_df <- tibble(
    NumerGroup = c("BleachingCode"), 
    DenomGroup = c("CORAL"))
  
  # For various combinations of numerator groups and denominator groups, calculate quantiles from bootstrap distribution
  boot_bleach_estimates <- FuncBootDraws(boot_dat = boot_bleach_samples, start_col = bleach_start_col, actual_dat = actual_bleach_dat, groups_df = bleach_grps_df, combos_df = bleach_combos_df, key_col = "BleachingCode", n_bootreps = 1000) #<<<<<<<<<<< SMALL NUMBER JUST FOR TESTING
  cbind("Site" = x[["Site"]], "SurvDate" = x[["SurvDate"]], boot_bleach_estimates)
})

site_bleach_CIs_df <- as.data.frame(do.call("rbind", site_bleach_CIs_list))



  # Summarize by REPORTING SITE - INCLUDE ACTIVE SITES ONLY. START YEAR IS THE YEAR THAT FIRST ACCOMMODATES ALL ACTIVE SITES!!
  if(CalcRepSites) {
    sub.RSdat <- FuncRecalcRS(trans.dat, RS_keepyrs)

    RS_temp <- FuncBootGroup(GroupID = "RepSiteSurvID", GroupName = "ReportingSite", CalcRC, trans.dat = sub.RSdat[[2]], coral_sub, run_parallel, set_cores)
    
    # Merge the data from PC_Site[[x]]
    PC.add.dat <- PC_Site[[x]] %>%
      filter(RSS %in% sub.RSdat[[1]]$Site) %>%
      mutate(SiteLev = "RepSiteSurvID")
    PC_Site[[x]] <- rbind(PC_Site[[x]], RS_temp[["PC"]], PC.add.dat)
    
    # Merge the data from RC_Site[[x]]
    if(CalcRC) {
      RC.add.dat <- RC_Site[[x]] %>%
        filter(RSS %in% sub.RSdat[[1]]$Site) %>%
        mutate(SiteLev = "RepSiteSurvID")
      RC_Site[[x]] <- rbind(RC_Site[[x]], RS_temp[["RC"]], RC.add.dat)
    }
  } else {
    PC.add.dat <- PC_Site[[x]] %>%
      mutate(SiteLev = "RepSiteSurvID")
    PC_Site[[x]] <- rbind(PC_Site[[x]], PC.add.dat)
    if(CalcRC) {
      RC.add.dat <- RC_Site[[x]] %>%
        mutate(SiteLev = "RepSiteSurvID")
      RC_Site[[x]] <- rbind(RC_Site[[x]], RC.add.dat)
    }
  }
}

# SURVEY-SITE tables by Subcategory levels ----
SubCat <- Trans_Subcat %>%
  left_join(unique(coral_sub[c("EventID", "SiteSurvID", "Transect")]), by = "EventID")
SubCat$N[is.na(SubCat$N)] <- 0

SS_PointCount <- SS_PercCov <- vector("list", length = length(unique(SubCat$SiteSurvID))) 
names(SS_PointCount) <- names(SS_PercCov) <- unique(SubCat$SiteSurvID)

for(s in unique(SubCat$SiteSurvID)) {
  sub.dat <- as.data.frame(subset(SubCat, SiteSurvID == s))
  sub.dat$Transect <- as.numeric(sub.dat$Transect)
  point.tab <- sub.dat %>%
    dplyr::select(Sublevel, N, Transect) %>%
    arrange(Sublevel, Transect) %>%
    spread(key = Transect, value = N) %>%
    set_rownames(.$Sublevel) %>%  # this comes from packages 'magrittr', because column_to_rownames wasn't working
    dplyr::select(-Sublevel) %>%
    data.matrix(rownames.force = TRUE)
  point.tab <- addmargins(as.matrix(point.tab), FUN = list(Total = sum), quiet = TRUE)
  SS_PointCount[[s]] <- point.tab
  
  perccov.tab <- sub.dat %>%
    dplyr::select(Sublevel, PC, Transect) %>%
    spread(key = Transect, value = PC) %>%
    set_rownames(.$Sublevel) %>%  # this comes from packages 'magrittr', because column_to_rownames wasn't working
    dplyr::select(-Sublevel) %>%
    data.matrix(rownames.force = TRUE)
  
  perccov.summary <- subset(PC_Site$Subcategory, RSS_SurvID == s & SiteLev == "SiteSurvID", select = c("Sublevel", "mean", "mean_low95", "mean_high95"))
  
  perccov.tab <- cbind(perccov.tab, perccov.summary[match(rownames(perccov.tab), perccov.summary$Sublevel), c("mean", "mean_low95", "mean_high95")])
  
  SS_PercCov[[s]] <- perccov.tab
}

# Bleaching Trends ----
bl_template <- merge(unique(coral_sub[c("EventID", "ReportingSite", "Site", "Year", "SurvDate", "Purpose", "IsActive")]), unique(coral_sub["BleachingCode"])) 
bl_template <- bl_template[complete.cases(bl_template),]

bl_temp <- coral_sub[c(names(bl_template), "CountOfTaxon")] %>%
  right_join(bl_template, by = names(bl_template)) %>%  # add the adjusted points column
  dplyr::group_by(.dots = names(bl_template)) %>%
  dplyr::summarize(N = sum(CountOfTaxon, na.rm=TRUE), .groups = "drop") %>%
  ungroup() %>%
  left_join(N_event, by = "EventID") 

# By site
bl_site <- bl_temp %>%
  dplyr::select(-EventID, -ReportingSite) %>%
  dplyr::group_by(Site, Year, SurvDate, Purpose, IsActive, BleachingCode) %>%
  dplyr::summarize(N = sum(N),
            AdjPoints = sum(ADJ_TOT),
            CoralPoints = sum(CORAL), .groups = "drop") %>%
  mutate(PC = round((N/ADJ_TOT)*100, 2),
         RC = round((N/CORAL)*100, 2),
         SiteLev = "SiteSurvID") %>% 
  ungroup() %>%
  dplyr::rename(RSS = "Site") %>%
  dplyr::select(SiteLev, everything())

# By reporting site
bl_sub.RSdat <- FuncRecalcRS(bl_temp, RS_keepyrs)
bl_RSdat <- bl_sub.RSdat[[2]] %>%
  dplyr::select(-EventID, -Site) %>%
  dplyr::group_by(ReportingSite, Year, SurvDate, Purpose, IsActive, BleachingCode) %>%
  dplyr::summarize(N = sum(N),
            AdjPoints = sum(ADJ_TOT),
            CoralPoints = sum(CORAL), .groups = "drop") %>%
  mutate(PC = round((N/ADJ_TOT)*100, 2),
         RC = round((N/CORAL)*100, 2),
         SiteLev = "RepSiteSurvID") %>% 
  ungroup() %>%
  dplyr::rename(RSS = "ReportingSite") %>%
  dplyr::select(SiteLev, everything())

# Merge the data from bl_site
bl_add.dat <- bl_site %>%
  filter(RSS %in% bl_sub.RSdat[[1]]$Site) %>%
  mutate(SiteLev = "RepSiteSurvID")

Bleach <- rbind(bl_site, bl_RSdat, bl_add.dat) %>%
  arrange(SiteLev, RSS, SurvDate, BleachingCode)
Bleach$BleachingCode <- factor(Bleach$BleachingCode, levels = c("NoData", "UNBL", "BL1", "BL2", "BL3", "BL4"))

# Add rest of data to mapdat ----
PC_temp <- subset(PC_Site$Category, SiteLev == "SiteSurvID" & Sublevel == "CORAL", select = c("RSS", "SurvDate", "mean")) %>% # for each site and survey, the mean % coral
  dplyr::rename("%Coral" = "mean") %>%
  mutate(`%NonCoral` = 100 - `%Coral`,
         RSS= as.character(RSS))

BleachRC <- Bleach %>%
  filter(SiteLev == "SiteSurvID") %>%
  dplyr::select(RSS, SurvDate, BleachingCode, RC) %>%
  mutate(RSS = as.character(RSS)) %>%
  spread(key = BleachingCode, value = RC, drop = FALSE, fill = 0) %>%
  full_join(PC_temp, by = c("RSS", "SurvDate"))

mapdat2 <- mapdat %>%
  mutate(Site = as.character(Site)) %>%
  full_join(BleachRC, by = c("Site" = "RSS")) %>%
  mutate(PopText = paste0(PopText, ": ", SurvDate))

# Create options tables for flexdashboard user input ---
PC_options <- rbind(
  data.frame(Level = "Category", Sublevel = unique(coral_sub$Category), stringsAsFactors = FALSE),
  data.frame(Level = "Subcategory", Sublevel = unique(coral_sub$Subcategory), stringsAsFactors = FALSE),
  data.frame(Level = "FunctionalGroup", Sublevel = unique(coral_sub$FunctionalGroup), stringsAsFactors = FALSE),
  data.frame(Level = "Taxon", Sublevel = unique(coral_sub$Taxon), stringsAsFactors = FALSE))
PC_options <- PC_options[complete.cases(PC_options),]
PC_options <- arrange(PC_options, Level, Sublevel) 
# %>% mutate_if(is.character, as.factor)

out.list <- list(rawdat = coral, # the raw data include equipment and shadow counts; only annual and episodic surveys included
                 CalcRepSites = CalcRepSites,
                 Warn_list = Warn_list,
                 PC_options = PC_options,
                 PC_Site = PC_Site,
                 SS_PointCount = SS_PointCount,
                 SS_PercCov = SS_PercCov,
                 RC_Site = RC_Site,
                 Bleach = Bleach,
                 mapdat = mapdat2,
                 coraldat = coral_sub)
saveRDS(out.list, paste0(out_prefix, "_coralsummary.RDS"))
}

