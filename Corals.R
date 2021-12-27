
# NOTE: Remember that for the following functions, need to add 'na.rm = TRUE': sum, mean

# Load packages ----
# Look for packages on local machine, install if necessary, then load all
pkgList_pre <- c("units", "tidyverse", "lubridate", "magrittr", "boot", "parallel", "multidplyr")
inst_pre <- pkgList_pre %in% installed.packages()
if (length(pkgList_pre[!inst_pre]) > 0) install.packages(pkgList_pre[!inst_pre],dep=TRUE)
lapply(pkgList_pre, library, character.only = TRUE)

FuncBootGroup <- function(GroupID, GroupName, CalcRC, trans.dat, coral.sub, run_parallel, set_cores) {
  
  FuncCalcVar <- function(boot_temp3, run_parallel, set_cores) {
    if(run_parallel) { # run in parallel. Boot.ci has a parallel processing function, but multidplyr does faster parallel processing here
      no_cores <- ifelse(set_cores > detectCores()-1, detectCores()-1, set_cores)
      
      clust <- new_cluster(no_cores)
      boot_temp4 <- boot_temp3 %>% partition(clust)
      cluster_library(cluster = clust, packages = c("dplyr", "purrr", "boot"))
      
      boot_temp5 <- boot_temp4 %>%
        mutate(
          simout_mean = purrr::map(.x = data, ~ boot(data = .x$cov, statistic = function(d,i) mean(d[i], na.rm = TRUE), R = 1000, stype = "i"))) %>%
        collect() %>% # Special collect() function to recombine partitions
        as_tibble() %>%
        ungroup()
    } else { # run serially
      boot_temp5 <- boot_temp3 %>%
        mutate(
          simout_mean = purrr::map(.x = data, ~ boot(data = .x$cov, statistic = function(d,i) mean(d[i], na.rm = TRUE), R = 1000, stype = "i"))) %>%
        ungroup()
    }
    
    for(i in 1:nrow(boot_temp5)) {
      boot_temp5[i, "mean_SE"] = sd(boot_temp5[i, "simout_mean"]$simout_mean[[1]]$t)
    }
    
    boot_temp5 %<>%
      mutate(
        boot_meanci = tryCatch(purrr::map(.x = simout_mean, ~ boot.ci(.x, conf = 0.95, type = "perc")), error = function(e) NA),
        mean_low95 = tryCatch(purrr::map(.x = boot_meanci, ~ .x$perc[[4]]), error = function(e) NA),
        mean_high95 = tryCatch(purrr::map(.x = boot_meanci, ~ .x$perc[[5]]), error = function(e) NA)) %>%
      dplyr::select(-data, -simout_mean, -boot_meanci)

    boot_temp5$mean_low95 <- sapply(boot_temp5$mean_low95, function(x) ifelse(is.null(x), NA, x))
    boot_temp5$mean_high95 <- sapply(boot_temp5$mean_high95, function(x) ifelse(is.null(x), NA, x))
    return(boot_temp5)
  }
  
  FuncCleanup <- function(boot_temp6, GroupID, GroupName, coral.sub) {
    boot_temp6 %<>%
      left_join(unique(coral.sub[c(GroupID, GroupName)]), by = c("RSS_SurvID" = GroupID)) %>%
      rename(RSS = GroupName) %>%
      mutate_at(c("mean", "mean_low95", "mean_high95", "mean_SE"), funs(round(., 2))) %>%
      dplyr::select("SiteLev", "RSS", "RSS_SurvID", everything()) %>%
      arrange(RSS_SurvID, Sublevel)
    return(boot_temp6)
  }
  
  # Summarize by grouping level
  boot_temp <- trans.dat %>%
    left_join(unique(coral.sub[, c("EventID", GroupID)]), by = "EventID") %>%
    rename("RSS_SurvID" = GroupID) %>%
    dplyr::select(RSS_SurvID, IsActive, Year, TripName, SurvDate, Purpose, Sublevel, PC, RC) %>%
    group_by(RSS_SurvID, IsActive, Year, TripName, SurvDate, Purpose, Sublevel)
  
  boot_temp2 <- boot_temp %>%
    add_tally() %>%  # number of transects per survey-site
    mutate(
      NumCoralTrans = sum(!is.na(RC)),
      PCtot = sum(PC, na.rm = TRUE),
      PCmean = PCtot/n, # mean % cover, averaged across transects
      RCtot = ifelse(is.na(sum(RC)), NA, sum(RC, na.rm = TRUE)), # if it's not a coral, then NA <<<<< MAKE SURE THIS IS THE WAY TO CALC
      RCmean = RCtot/NumCoralTrans) %>%
    dplyr::select(-PC, -RC) %>%
    distinct() %>%
    arrange(RSS_SurvID, Sublevel, Year) %>%
    ungroup() %>% # hold this summary information (doesn't require bootstrapping) for later merging...
    rename(NumTransect = n) %>%
    mutate(EstimPCVar = NumTransect >= 8 & PCtot > 0,
           EstimRCVar = NumTransect >= 8 & RCtot > 0)# if at least 8 transects and non-zero, OK to estimate variability
  boot_temp2$EstimRCVar[is.na(boot_temp2$EstimRCVar)] <- FALSE
  
  # Calculate variability for PC
  PCboot_temp3 <- boot_temp %>% # these are the ones to estimate PC variability for
    dplyr::select(-RC) %>%
    rename(cov = PC) %>%
    tidyr::nest() %>% # nests the PC data for each site
    left_join(unique(boot_temp2[c("RSS_SurvID", "Sublevel", "EstimPCVar")]), by = c("RSS_SurvID", "Sublevel")) %>%
    filter(EstimPCVar == TRUE) %>%
    dplyr::select(-EstimPCVar) 
  
  PCboot_temp5 <- FuncCalcVar(PCboot_temp3, run_parallel, set_cores)
         
  PCboot_temp6 <- boot_temp2 %>%
    dplyr::select(-RCtot, -RCmean, -EstimRCVar) %>%
    full_join(PCboot_temp5, by = c("RSS_SurvID", "IsActive", "Year", "TripName", "SurvDate", "Purpose", "Sublevel")) %>%
    rename(mean = PCmean)
  PCboot_temp6$SiteLev <- GroupID

  PCboot_temp7 <- FuncCleanup(PCboot_temp6, GroupID, GroupName, coral.sub)
  
  # Calculate variability for RC 
  if(CalcRC) {
    RCboot_temp3 <- boot_temp %>% # these are the ones to estimate RC variability for
      dplyr::select(-PC) %>%
      filter(!is.na(RC)) %>%
      rename(cov = RC) %>%
      tidyr::nest() %>% # nests the PC data for each site
      left_join(unique(boot_temp2[c("RSS_SurvID", "Sublevel", "EstimRCVar")]), by = c("RSS_SurvID", "Sublevel")) %>%
      filter(EstimRCVar == TRUE) %>%
      dplyr::select(-EstimRCVar) 
    
    RCboot_temp5 <- FuncCalcVar(RCboot_temp3, run_parallel, set_cores)
    RCboot_temp6 <- boot_temp2 %>%
      dplyr::select(-PCtot, -PCmean, -EstimPCVar) %>%
      filter(!is.na(RCtot)) %>%
      full_join(RCboot_temp5, by = c("RSS_SurvID", "IsActive", "Year", "TripName", "SurvDate", "Purpose", "Sublevel")) %>%
      rename(mean = RCmean)
    RCboot_temp6$SiteLev <- GroupID

    RCboot_temp7 <- FuncCleanup(RCboot_temp6, GroupID, GroupName, coral.sub)
    } else {
      RCboot_temp7 <- NA
      }
  
  return.list <- list(PC = PCboot_temp7, RC = RCboot_temp7)
  return(return.list)
}

FuncRecalcRS <- function(dat, RSyrs) {
  sub.dat <- dat
  sub.dat[c("ReportingSite", "Site")] <- sapply(sub.dat[c("ReportingSite", "Site")], as.character)
  add.sites <- unique(subset(sub.dat, ReportingSite == Site, select = Site)) # these are the sites to add back after calculations
  recalc.RSdat <- subset(sub.dat, ReportingSite != Site & IsActive==1) # for sites summarized at RS level, only include the Active sites
  if(nrow(recalc.RSdat) > 0) {
    recalc.RSdat %<>%
      rowwise() %>%
      mutate(TooEarly = Year < RSyrs$startyr[RSyrs$ReportingSite == ReportingSite]) %>%
      filter(TooEarly == FALSE) %>%
      dplyr::select(-TooEarly) %>%
      ungroup()
  }
  
  return.list <- list(add.sites, recalc.RSdat)
  return(return.list)
}

FuncCorals <- function(filenam, sitesfilenam = NULL, out_prefix, run_parallel, set_cores = 1) {

# Clean data ----
coral <- read_csv(filenam)
if(!is.null(sitesfilenam)) tbl_link <- read_csv(sitesfilenam)

#  <<<<<<<< TESTING >>>>>>>>>>>>>> ----
coral <- read_csv("Data_LOCAL_ONLY/BUIS_CoralVideo Summary by Transect.csv") # <- read_csv("Data_LOCAL_ONLY/demo_CoralDat.csv")
out_prefix="test"
run_parallel=TRUE
set_cores=11
tbl_link <- read_csv("Data_LOCAL_ONLY/SFCN_CoralSites.csv")

# Format columns ----
coral$Date <- mdy(coral$Date)
coral %<>%
  rename(SurvDate = Date, Subcategory = SubCategory, Taxon = TaxonCode, CountOfTaxon = CountOfTaxonCode) %>%
  mutate(EventID = paste0(ParkCode, "_", Site, "_", Transect, "_", TripName, "_", Purpose),  # create unique transect-specific EventID
         SiteSurvID = paste0(ParkCode, "_", Site, "_", TripName, "_", Purpose)) %>%
  arrange(Site, Transect, SurvDate, Purpose)

# Populate list of warnings ----
Warn_list <- sapply(c("AltPurp", "NoTaxon", "BleachCode", "TransCount", "TaxonCountDiv", "TaxonCount", "ShadowCount", "UNKCount"), function(x) NULL)

Warn_list$AltPurp <- coral %>%
  arrange(SurvDate) %>%
  dplyr::select(TripName, Purpose) %>%
  filter(!Purpose %in% c("Annual", "Episodic")) %>%
  distinct() %>%
  rename("Trip" = TripName)

Warn_list$NoTaxon <- coral %>%
  dplyr::select(EventID, Taxon) %>%
  filter(is.na(Taxon) | Taxon == "No Taxon") %>% # 
  rename("Park_Site_Transect_Trip_Purpose" = "EventID")

Warn_list$BleachCode <- coral %>%
  dplyr::select(EventID, Category, BleachingCode) %>%
  filter(Category != "CORAL" & !is.na(BleachingCode)) %>%
  rename("Park_Site_Transect_Trip_Purpose" = "EventID")
      
Warn_list$TransCount <- coral %>%
  dplyr::select(Site, Transect, SurvDate, Purpose) %>%
  distinct() %>%
  group_by(Site, SurvDate, Purpose) %>%
  arrange(Site, SurvDate, Purpose) %>%
  summarize(NumTransects = n()) %>%
  filter(Purpose %in% c("Annual", "Episodic") & !NumTransects %in% c(4, 20)) %>% # 'off' transect counts only matters for annual or episodic surveys
  dplyr::mutate(SurvDate = lubridate::ymd(SurvDate)) %>% # renderTable does not play nice with dates--need to format as character
  dplyr::rename("Survey Date" = "SurvDate")

Err_counts <- coral %>%
  group_by(EventID) %>%
  summarize(TotPoints = sum(CountOfTaxon, na.rm=TRUE),
            ShadowPoints = sum(CountOfTaxon[Category == "SHADOW"]),
            UNKPoints = sum(CountOfTaxon[Category %in% c("UNK", "UNKNOWN") | is.na(Category)]))

Warn_list$TaxonCountDiv <- Err_counts %>%
  mutate(ModOut = TotPoints %% 10) %>%
  filter(ModOut != 0) %>%
  dplyr::select(EventID, TotPoints) %>%
  rename("Park_Site_Transect_Trip_Purpose" = "EventID", "Total Points" = TotPoints)
  
Warn_list$TaxonCount <- Err_counts %>%
  filter(TotPoints < 200 | TotPoints > 480) %>%
  dplyr::select(EventID, TotPoints) %>%
  rename("Park_Site_Transect_Trip_Purpose" = "EventID", "Total Points" = TotPoints)

Warn_list$ShadowCount <- Err_counts %>%
  mutate(PercShadow = (ShadowPoints / TotPoints)*100) %>%
  filter(PercShadow > 5.0) %>%
  dplyr::select(EventID, PercShadow) %>%
  arrange(desc(PercShadow)) %>%
  rename("Park_Site_Transect_Trip_Purpose" = "EventID", "% Shadow" = PercShadow)
  
Warn_list$UNKCount <- Err_counts %>%
  mutate(PercUNK = (UNKPoints / TotPoints)*100) %>%
  filter(PercUNK > 5.0) %>%
  dplyr::select(EventID, PercUNK) %>%
  arrange(desc(PercUNK)) %>%
  rename("Park_Site_Transect_Trip_Purpose" = "EventID", "% Unknown" = PercUNK)

# Format data for analyses ----

# Account for different classification of bleaching prior to 10/1/05
coral$BleachingCode[is.na(coral$BleachingCode) & coral$Category =="CORAL" & coral$SurvDate < "2005-10-01"] <- "NoData"
coral$BleachingCode[is.na(coral$BleachingCode) & coral$Category =="CORAL" & coral$SurvDate >= "2005-10-01"] <- "UNBL"
coral$BleachingCode[coral$BleachingCode == "BL" & coral$Category =="CORAL"] <- "BL1"
coral$BleachingCode[coral$Category != "CORAL"] <- NA

coral.sub <- coral %>%
  filter(!Category %in% c("EQUIP", "SHADOW")) %>%  # remove EQUIP & SHADOW data, to calculate the adjusted points. 
  filter(Purpose %in% c("Annual", "Episodic")) %>% # only include annual and episodic data
  arrange(Site, SurvDate, Transect)

# If the file doesn't already have a Reporting Site column, get the information from the link table
if(!"ReportingSite" %in% colnames(coral.sub) & exists("tbl_link")) {
  coral.sub %<>% left_join(subset(tbl_link, select=-c(IsActive)), by = c("ParkCode", "Site"))
}

# If the file doesn't already have an Activity Status column, get the information from the link table
if(!"IsActive" %in% colnames(coral.sub) & exists("tbl_link")) {
  coral.sub %<>% left_join(subset(tbl_link, select=c(ParkCode, Site, IsActive)), by = c("ParkCode", "Site"))
}

coral.sub %<>%
  mutate(RepSiteSurvID = paste0(ParkCode, "_", ReportingSite, "_", TripName)) %>% # create a unique ID for each reporting site-survey event
  dplyr::select(Year, SurvDate, ParkCode, Site, Transect, Latitude, Longitude, Taxon, BleachingCode, Category, Subcategory, FunctionalGroup, CountOfTaxon, EventID, RepSiteSurvID, SiteSurvID, TripName, Purpose, ReportingSite, ReportingSiteName, IsActive) 
# %>%  mutate_if(is.character, as.factor)

# STOP-Error check--make sure all records have associated Reporting Site and Activity Status ----
missing_RS <- coral.sub[is.na(coral.sub$ReportingSite), ] %>%
  distinct()
if(nrow(missing_RS) > 0) stop("These data records do not have reporting site information:\n", paste(capture.output(print(missing_RS)), collapse = "\n")) else cat("OK >>> All records have reporting site info")
missing_IsActive <- coral.sub[is.na(coral.sub$IsActive), ] %>%
  distinct()
if(nrow(missing_IsActive) > 0) stop("These data records do not have activity status information:\n", paste(capture.output(print(missing_IsActive)), collapse = "\n")) else cat("OK >>> All records have activity status info")

CalcRepSites <- !identical(coral.sub$Site, coral.sub$ReportingSite) # if TRUE it means there are reporting sites that are different from sites, so need to separately summarize by ReportingSite (for those that are different)

# Uncomment (i.e., remove the '#') the line below to output a CSV of the cleaned data
# write.csv(coral.sub, paste0(out_prefix, "_cleandat.csv"), row.names = FALSE) 


# Map data ----
mapdat <- coral.sub %>%
  dplyr::select(ReportingSiteName, ReportingSite, Site, Latitude, Longitude, Year, IsActive) %>% # pull only the data relevant for mapping by Site
  group_by(Site) %>%
  summarize( # site-level summary stats for mapping
    ReportingSiteName = unique(ReportingSiteName),
    ReportingSite = unique(ReportingSite),
    IsActive = unique(IsActive),
    minyr = min(Year, na.rm = TRUE),
    maxyr = max(Year, na.rm = TRUE),
    medLat = median(Latitude, na.rm = TRUE),
    medLong = median(Longitude, na.rm = TRUE)) %>%
  mutate(
    poptext = paste0(Site, " (", minyr, " - ", maxyr, ")"))

# For each active reporting site, this is the year to start trend plots (because all active sites started by this year) ----
RS_startyr <- mapdat %>%
  dplyr::select(ReportingSite, IsActive, minyr) %>%
  filter(IsActive == TRUE) %>%
  group_by(ReportingSite) %>%
  summarize(startyr = max(minyr))

# Calculate adjusted total number of points and number of coral points per survey-transect ----
# >>>>>>>>>>>>>> PICK UP FROM HERE--ADD THE ALGAL FUNCTIONAL GROUPS FROM SOP 21
N.event <- coral.sub %>%
  group_by(EventID) %>%
  summarize(AdjPoints = sum(CountOfTaxon, na.rm=TRUE))  # use for denominator. This AdjPoints is the number of hits of "something" for that particular transect survey
N.event$AdjPoints[is.na(N.event$AdjPoints)] <- 0
N.categ <- coral.sub %>% # for categories sponge, gorgo, coral, calculate number of hits to use as denominator for relative cover
  filter(Category %in% c("CORAL", "GORGO", "SPONGE")) %>%
  group_by(Category, EventID) %>% 
  summarize(CategPoints = sum(CountOfTaxon, na.rm=TRUE))   # This is the number of hits of "something" in that category and particular transect survey
N.categ$CategPoints[is.na(N.categ$CategPoints)] <- 0

# Calculate percent cover and relative cover by different grouping variables ----
group.var <- c("Taxon", "Category", "Subcategory", "FunctionalGroup", "BleachingCode")
PC_Site <- vector("list", length = length(group.var)) # each list element is a data frame with percent cover information by that grouping variable
names(PC_Site) <- group.var 
RC_Site <- vector("list", length = 2)
names(RC_Site) <- c("Taxon", "FunctionalGroup")
  
for(x in group.var) {
  CalcRC <- x %in% names(RC_Site) # for these groupings, also calculate relative cover
  
  incProgress(1/(length(group.var) + 1987), detail = paste0(" Bootstrapping 95% CI for each level of ", x)) 
  
  cov_template <- merge(unique(coral.sub[c("EventID", "Site", "ReportingSite", "Year", "TripName", "SurvDate", "Purpose", "IsActive")]), unique(coral.sub[x]))  # create template to add zeros when species are not detected in a survey-transect
  cov_template <- cov_template[complete.cases(cov_template),]
  
  # Summarize by SURVEY-TRANSECT
  
  if(CalcRC) {
    # Calculate relative cover only for coral taxa and coral functional groups
    coral_categ <- coral.sub %>%
      dplyr::select(c(x, "Category")) %>%
      dplyr::rename("Sublevel" = x) %>%
      filter(Category == "CORAL") %>%
      dplyr::select(-Category) %>%
      filter(!is.na(Sublevel)) %>% # get rid of the NA functional groups
      distinct()
  }
  
  trans.dat <- coral.sub[c(names(cov_template), "CountOfTaxon")] %>%
    group_by(.dots = names(cov_template)) %>%
    summarize(N = sum(CountOfTaxon, na.rm=TRUE)) %>%
    right_join(cov_template, by = names(cov_template)) %>%  # add the adjusted points column
    left_join(N.event, by = "EventID") %>%
    mutate(PC = round((N/AdjPoints)*100, 2), # PC is the proportion of usable data points in each transect that is X
           RC = round((N/CoralPoints)*100, 2)) %>% # RC is relative cover (denominator is the number of CORAL hits), calculated only for coral
    arrange(Site, SurvDate) %>%
    ungroup() %>%
    dplyr::rename("Sublevel" = x)
  
  if(CalcRC) {
    trans.dat$RC[is.na(trans.dat$RC)] <- 0 # if there are no hits, then RC is 0
    trans.dat$RC[trans.dat$CoralPoints == 0] <- NA # but if there are no coral hits in that event, then RC is NA
    trans.dat$RC[!trans.dat$Sublevel %in% coral_categ$Sublevel] <- NA # if it's not a coral category, then RC is NA
    } else { 
      trans.dat$RC <- NA # if it's not a grouping for which relative cover is calculated, then convert all values to NA
    }
  trans.dat$PC[is.na(trans.dat$PC)] <- 0
  
  if(x == "Subcategory") Trans_Subcat <- trans.dat
  
  # Summarize by SITE
  SS_temp <- FuncBootGroup(GroupID = "SiteSurvID", GroupName = "Site", CalcRC, trans.dat, coral.sub, run_parallel, set_cores)
  
  PC_Site[[x]] <- SS_temp[["PC"]]
  
  if(CalcRC) {
    RC_Site[[x]] <- SS_temp[["RC"]]
  }
  
  # Summarize by REPORTING SITE - INCLUDE ACTIVE SITES ONLY. START YEAR IS THE YEAR THAT FIRST ACCOMMODATES ALL ACTIVE SITES!!
  if(CalcRepSites) {
    sub.RSdat <- FuncRecalcRS(trans.dat, RS_startyr)

    RS_temp <- FuncBootGroup(GroupID = "RepSiteSurvID", GroupName = "ReportingSite", CalcRC, trans.dat = sub.RSdat[[2]], coral.sub, run_parallel, set_cores)
    
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
  left_join(unique(coral.sub[c("EventID", "SiteSurvID", "Transect")]), by = "EventID")
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
bl_template <- merge(unique(coral.sub[c("EventID", "ReportingSite", "Site", "Year", "SurvDate", "Purpose", "IsActive")]), unique(coral.sub["BleachingCode"])) 
bl_template <- bl_template[complete.cases(bl_template),]

bl_temp <- coral.sub[c(names(bl_template), "CountOfTaxon")] %>%
  right_join(bl_template, by = names(bl_template)) %>%  # add the adjusted points column
  group_by(.dots = names(bl_template)) %>%
  summarize(N = sum(CountOfTaxon, na.rm=TRUE)) %>%
  ungroup() %>%
  left_join(N.event, by = "EventID") 

# By site
bl_site <- bl_temp %>%
  dplyr::select(-EventID, -ReportingSite) %>%
  group_by(Site, Year, SurvDate, Purpose, IsActive, BleachingCode) %>%
  summarize(N = sum(N),
            AdjPoints = sum(AdjPoints),
            CoralPoints = sum(CoralPoints)) %>%
  mutate(PC = round((N/AdjPoints)*100, 2),
         RC = round((N/CoralPoints)*100, 2),
         SiteLev = "SiteSurvID") %>% 
  ungroup() %>%
  dplyr::rename(RSS = "Site") %>%
  dplyr::select(SiteLev, everything())

# By reporting site
bl_sub.RSdat <- FuncRecalcRS(bl_temp, RS_startyr)
bl_RSdat <- bl_sub.RSdat[[2]] %>%
  dplyr::select(-EventID, -Site) %>%
  group_by(ReportingSite, Year, SurvDate, Purpose, IsActive, BleachingCode) %>%
  summarize(N = sum(N),
            AdjPoints = sum(AdjPoints),
            CoralPoints = sum(CoralPoints)) %>%
  mutate(PC = round((N/AdjPoints)*100, 2),
         RC = round((N/CoralPoints)*100, 2),
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
  mutate(poptext = paste0(poptext, ": ", SurvDate))

# Create options tables for flexdashboard user input ---
PC_options <- rbind(
  data.frame(Level = "Category", Sublevel = unique(coral.sub$Category), stringsAsFactors = FALSE),
  data.frame(Level = "Subcategory", Sublevel = unique(coral.sub$Subcategory), stringsAsFactors = FALSE),
  data.frame(Level = "FunctionalGroup", Sublevel = unique(coral.sub$FunctionalGroup), stringsAsFactors = FALSE),
  data.frame(Level = "Taxon", Sublevel = unique(coral.sub$Taxon), stringsAsFactors = FALSE))
PC_options <- PC_options[complete.cases(PC_options),]
PC_options <- arrange(PC_options, Level, Sublevel) 
# %>% mutate_if(is.character, as.factor)

out.list <- list(CalcRepSites = CalcRepSites,
                 Warn_list = Warn_list,
                 PC_options = PC_options,
                 PC_Site = PC_Site,
                 SS_PointCount = SS_PointCount,
                 SS_PercCov = SS_PercCov,
                 RC_Site = RC_Site,
                 Bleach = Bleach,
                 mapdat = mapdat2,
                 coraldat = coral.sub)
saveRDS(out.list, paste0(out_prefix, "_coralsummary.RDS"))
}

