# Load packages ----
# Look for packages on local machine, install if necessary, then load all
pkgList_pre <- c("magrittr", 
                 "plyr",
                 "lubridate",
                 "purrr", # to map functions to elements
                 "tidyverse") # broom, units
inst_pre <- pkgList_pre %in% installed.packages()
if (length(pkgList_pre[!inst_pre]) > 0) install.packages(pkgList_pre[!inst_pre],dep=TRUE)
lapply(pkgList_pre, library, character.only = TRUE)

# FuncCorals <- function(filenam, sitesfilenam = NULL, out_prefix, run_parallel, set_cores = 1) {

# Clean data ----
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

# # Populate list of warnings ----
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
# Format data for analyses ----

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

# STOP-Error check--make sure all records have associated Reporting Site and Activity Status ----
missing_RS <- coral[is.na(coral$ReportingSite), ] %>%
  distinct()
if(nrow(missing_RS) > 0) stop("These data records do not have reporting site information:\n", paste(capture.output(print(missing_RS)), collapse = "\n")) else cat("OK >>> All records have reporting site info")
missing_IsActive <- coral[is.na(coral$IsActive), ] %>%
  distinct()
if(nrow(missing_IsActive) > 0) stop("These data records do not have activity status information:\n", paste(capture.output(print(missing_IsActive)), collapse = "\n")) else cat("OK >>> All records have activity status info")

# Uncomment (i.e., remove the '#') the line below to output a CSV of the cleaned data
# write.csv(coral, paste0(out_prefix, "_cleandat.csv"), row.names = FALSE) 

# Map data ----
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

# Bootstrap functions ----

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

FuncBootDraws <- function(dat, tot_df, groups_df, n_bootreps = 10000) {
  # Function to calculate bootstrap CI's for % cover and relative % cover, by site-survey (for intensive sites)
  #
  # Args:
  #   dat = data frame of the parametric bootstrap samples for all transects of a single site-survey. Requires cols for "TransectSurveyID" and each taxon. 
  #   tot_df = data frame of the adjusted total (count) for each transect-survey for that site-survey. Use as denominator for % cover estimates
  #   groups_df = data frame assigning each taxon to higher level groups
  #   n_bootreps = number of bootstrap reps to calculate CI's over
  #
  # Returns:
  #   A data frame with bootstrapped CI's for % cover and relative % cover. Has these cols: DenomGroup (e.g., AdjTot, ALGAE, CORAL, GORGO, SPONGE), NumerGroup (e.g., Category, FunctionalGroup, Taxon), NumerLevel(e.g., ALGAE). Final columns have the quantiles for bootstrapped 50%, 75%, 95%CI.
  #
  # NOTE: Not estimating CI's for trend -- just individually by site-survey.
  
  n_transects <- length(unique(dat$TransectSurveyID)) # number of transects surveyed in that site-survey
  
  samps <- dat %>%
    dplyr::slice_sample(n = n_draws * n_transects, replace = TRUE) %>% # randomly sample from the parametric bootstrap samples
    tibble::rowid_to_column(var = "RowID") %>%
    as.data.frame() # can't add repeating vector column to tibble
  
  samps[, 3:ncol(samps)][is.na(samps[, 3:ncol(samps)])] <- 0 # replace NA's with zero
  samps$BootRep = rep(1:n_draws, each = n_transects) # add col of bootstrap replicate number; each replicate has same number of transects as original data for the site-survey
    
  # This is the master data frame from which to calculate % cover at different grouping levels
  samps_master <- samps %>%
    tidyr::gather(key = "Taxon", value = "Count", -BootRep, -RowID, -TransectSurveyID) %>% # convert to long format
    dplyr::left_join(groups_df, by = "Taxon")
  
  # Calculate % cover and RELATIVE % cover for certain groupings:
  # > % cover by Category. Denominator is adjusted total count for the sampled transect. Transects have equal weight in a BootRep.
  # > Drill down in certain categories--Calculate % cover and RELATIVE % cover by Functional Group for ALGAE and CORAL (NumerGroup = FunctionalGroup, DenomGroup = ALGAE or CORAL); and by Taxon for ALGAE, CORAL, GORGO, SPONGE (e.g., NumerGroup = Taxon, DenomGroup = GORGO). Weighted by the # of hits of the category in each transect.
  
  # Calculate bootstrapped CI's for these groupings
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
    dplyr::filter(if(x[["DenomGroup"]] != "AdjTot") Category == x[["DenomGroup"]] else TRUE) %>% # for drill down, restrict to one category
    dplyr::group_by(BootRep, TransectSurveyID, RowID, NumerLevel) %>%

    dplyr::summarize(GroupCount = sum(Count), .groups = "drop")

  # For % cover, calculate per transect and then take mean across all transects in the BootRep
  out_df <- tmp %>%
    dplyr::left_join(tot_df, by = "TransectSurveyID") %>%
    dplyr::mutate(PercCov = GroupCount / AdjTot) %>%
    dplyr::group_by(BootRep, NumerLevel) %>%
    dplyr::summarize(MeanPercCov = mean(PercCov), .groups = "drop") %>%
    dplyr::group_by(NumerLevel) %>%
    dplyr::summarize_at(vars(MeanPercCov), FuncPullQuant, .groups = "drop") %>%
    dplyr::mutate(NumerGroup = x[["NumerGroup"]], DenomGroup = "AdjTot")
  
  # For relative % cover, sum(numerator across all transects)/sum(denom across all transects)
  if(x[["DenomGroup"]] != "AdjTot") {
    categtot_df <- tmp %>%
      dplyr::group_by(BootRep) %>%
      dplyr::summarize(CategTot = sum(GroupCount), .groups = "drop")

    add_df <- tmp %>%
      dplyr::group_by(BootRep, NumerLevel) %>%
      dplyr::summarize(SampGroupCount = sum(GroupCount), .groups = "drop") %>%
      dplyr::left_join(categtot_df, by = "BootRep") %>%
      dplyr::mutate(SiteRelCov = SampGroupCount/CategTot) %>%
      dplyr::group_by(NumerLevel) %>%
      dplyr::summarize_at(vars(SiteRelCov), FuncPullQuant, .groups = "drop") %>%
      dplyr::mutate(NumerGroup = x[["NumerGroup"]], DenomGroup = x[["DenomGroup"]])

    out_df <- rbind(out_df, add_df)
  }
  
  out_df %<>%
    dplyr::mutate_if(is.numeric, round, 3) %>%  # calculate quantiles and output in nice format
    dplyr::mutate(Category = x[["DenomGroup"]]) %>%
    dplyr::select(Category, NumerGroup, NumerLevel, DenomGroup, everything())
  return(out_df)
  })
  
  boot_CIs_df <- as.data.frame(do.call("rbind", boot_CIs_list)) # combine list elements in a single data frame
  return(boot_CIs_df)
}
  
# Create master data frame for bootstrapping ----
boot_master_dat <- dplyr::filter(coral, !Category %in% c("EQUIP", "SHADOW")) %>%  # remove EQUIP & SHADOW data, for all analyses
  dplyr::select(ReportingSite, Site, Transect, SurvDate, TransectSurveyID, Category, FunctionalGroup, Taxon, CountOfTaxon)

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

# Convert to % cover by transect-survey
boot_site_dat[, 5:ncol(boot_site_dat)] <- boot_site_dat[, 5:ncol(boot_site_dat)]/boot_site_dat$AdjTot

# For each transect-survey, generate parametric bootstrap samples (default is N = 1M) using the probabilities from the collected data. Run one site-survey at a time to avoid computer memory problems ---- 
site_CIs_list <- apply(unique(boot_site_dat[c("Site", "SurvDate")]), 1, FUN = function(x) { 
  subdat <- subset(boot_site_dat, Site == x[["Site"]] & SurvDate == x[["SurvDate"]])
  boot_samples <- FuncBootSamples(dat = subdat, n_bootsamples = 100000 ) #<<<<<<<<<<< SMALL NUMBER JUST FOR TESTING

  # For each taxon, calculate quantiles from bootstrap distribution
  boot_site_estimates <- FuncBootDraws(dat = boot_samples, tot_df = subdat[, c("TransectSurveyID", "AdjTot")], groups_df = grps_df, n_bootreps = 1000) #<<<<<<<<<<< SMALL NUMBER JUST FOR TESTING
  cbind("Site" = x[["Site"]], "SurvDate" = x[["SurvDate"]], boot_site_estimates)
  })   
site_CIs_df <- as.data.frame(do.call("rbind", site_CIs_list))


# >>>>> PICK UP FROM HERE. ADD THE MEANS FROM ACTUAL DATA.









# # For each active reporting site, these are the start and end years for reporting site plots (because all subsites were surveyed in these years) ----
# RS_keepyrs <- mapdat %>%
#   dplyr::select(ReportingSite, IsActive, MinYr, MaxYr) %>%
#   filter(IsActive == TRUE) %>%
#   dplyr::group_by(ReportingSite) %>%
#   dplyr::summarize(startyr = max(MinYr),
#             endyr = min(MaxYr), .groups = "drop")

# (NOTE: BAR02_01_2018-05-02 has AdjTot = 10--this is a data error that is flagged in the data errors page of the dashboard and is supposed to be corrected by the user before proceeding)






  # if(CalcRC) {
  #   trans.dat$RC[is.na(trans.dat$RC)] <- 0 # if there are no hits, then RC is 0
  #   trans.dat$RC[trans.dat$CoralPoints == 0] <- NA # but if there are no coral hits in that event, then RC is NA
  #   trans.dat$RC[!trans.dat$Sublevel %in% coral_categ$Sublevel] <- NA # if it's not a coral category, then RC is NA
  #   } else { 
  #     trans.dat$RC <- NA # if it's not a grouping for which relative cover is calculated, then convert all values to NA
  #   }
  # trans.dat$PC[is.na(trans.dat$PC)] <- 0


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

