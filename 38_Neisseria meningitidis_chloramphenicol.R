#===============================================================================
# Title		: Summarising NMFI-AMR systematic literature review database
# Data version	: Data received from Poojan Shrestha on 15-07-2019
# Manuscript	: NMFI-Antimicrobial resistance analysis
# URL			: https://www.iddo.org/non-malarial-febrile-illness-map
# Script date	: 14-04-2021
# Script		: Prabin Dahal; prabin.dahal@iddo.org
# R Version		: R version 4.0.4 (2021-02-15)
#===============================================================================
#rm(list=ls())
require(tableone)
require(dplyr)
require(plyr) # for revalue

# meta-analysis packages
require(meta)
require(metafor)
require(metasens)

#-----------------------
# Read analysis dataset 
#-----------------------
dat0<-read.csv("D:/_IDDO/_nmfi_amr/NMFI AMR/Data/nmfi_amr_analysis_set_final.csv")
dat<- dat0

#-----------------------------------
# Look at pathogen-drug combination
#-----------------------------------
drugbug <- dat[which(dat$pathogen_2=="Neisseria meningitidis" & dat$antibiotic %in% c("Chloramphenicol")),]
drugbug <- droplevels(drugbug)

# number of unique articles
nrow(unique<-drugbug[which(!duplicated(drugbug$article_id)),])

drugbug_africa<- drugbug[which(drugbug$region=="Africa"),]
drugbug_asia<- drugbug[which(drugbug$region=="Asia"),]

#===========================================================================================
# Part I: Summarise the number of events and total isolates tested by region and time-period
#===========================================================================================

#--------------------------------------------------------------
# Exclude studies with either missing numerator or denominator
#--------------------------------------------------------------
drugbug1 		<- drugbug[which(!is.na(drugbug$number_resistant) & !is.na(drugbug$total_invasive_isolates)),]
drugbug_africa1 	<- drugbug_africa[which(!is.na(drugbug_africa$number_resistant) & !is.na(drugbug_africa$total_invasive_isolates)),]
drugbug_asia1 	<- drugbug_asia[which(!is.na(drugbug_asia$number_resistant) & !is.na(drugbug_asia$total_invasive_isolates)),]

# number of unique articles
nrow(unique_1<-drugbug1[which(!duplicated(drugbug1$article_id)),])

#----------------------------------
# Total isolates and total events
#----------------------------------
drugbug1 %>% 
	dplyr::summarise(
		n_pub= length(unique(article_id)),
		n_site= length(unique(site_code)),
		n_resis = sum(number_resistant, na.rm=T),
		n_total= sum(total_invasive_isolates, na.rm=T),
		paste = paste(sum(number_resistant, na.rm=T), sum(total_invasive_isolates, na.rm=T), sep="/")
	)
#-------------------------------------------
# Total isolates and total events for Africa
#-------------------------------------------
drugbug_africa1 %>% 
	dplyr::summarise(
		n_site= length(unique(site_code)),
		n_resis = sum(number_resistant, na.rm=T),
		n_total= sum(total_invasive_isolates, na.rm=T),
		paste = paste(sum(number_resistant, na.rm=T), sum(total_invasive_isolates, na.rm=T), sep="/")
	)
drugbug_africa1 %>% 
	dplyr::group_by(year_cat) %>%
	dplyr::summarise(
		n_site= length(unique(site_code)),
		n_resis = sum(number_resistant, na.rm=T),
		n_total= sum(total_invasive_isolates, na.rm=T),
		paste = paste(sum(number_resistant, na.rm=T), sum(total_invasive_isolates, na.rm=T), sep="/")
	)
#-------------------------------------------
# Total isolates and total events for Asia
#-------------------------------------------
drugbug_asia1 %>% 
	dplyr::summarise(
		n_site= length(unique(site_code)),
		n_resis = sum(number_resistant, na.rm=T),
		n_total= sum(total_invasive_isolates, na.rm=T),
		paste = paste(sum(number_resistant, na.rm=T), sum(total_invasive_isolates, na.rm=T), sep="/")
	)
drugbug_asia1 %>% 
	dplyr::group_by(year_cat) %>%
	dplyr::summarise(
		n_site= length(unique(site_code)),
		n_resis = sum(number_resistant, na.rm=T),
		n_total= sum(total_invasive_isolates, na.rm=T),
		paste = paste(sum(number_resistant, na.rm=T), sum(total_invasive_isolates, na.rm=T), sep="/")
	)
#=======================================================
# Part II: Carry out meta-analysis of single proportions
#=======================================================
#--------------------------------------------------
# Some studies have split data into multiple rows
# Aggregate by study sites for meta-analysis
#--------------------------------------------------

drugbug1 <- drugbug1 %>% 
		dplyr::group_by(site_code,year_cat,region) %>%
		dplyr::summarise(
			number_resistant  = sum(number_resistant), 
			total_invasive_isolates =sum(total_invasive_isolates)
		)
drugbug_africa1 <- drugbug_africa1 %>% 
		dplyr::group_by(site_code,year_cat) %>%
		dplyr::summarise(
			number_resistant  = sum(number_resistant), 
			total_invasive_isolates =sum(total_invasive_isolates)
		)
drugbug_asia1 <- drugbug_asia1 %>% 
		dplyr::group_by(site_code,year_cat) %>%
		dplyr::summarise(
			number_resistant = sum(number_resistant), 
			total_invasive_isolates =sum(total_invasive_isolates)
		)

#----------------
# Meta-analysis
#----------------

(ma1 <- metaprop(
			data = drugbug1,
			number_resistant, 
			total_invasive_isolates,
			studlab = site_code
			)
		)

# Model failed to converge when applying the default logit transformation

# Use untransformed proportion for derivation
(ma1 <- metaprop(
			data = drugbug1,
			number_resistant, 
			total_invasive_isolates,
			studlab = site_code,
			sm="PRAW"
			)
		)

update.meta(ma1, byvar=region, comb.random = TRUE, comb.fixed = F)

# Examine the influence of any single study
metainf(ma1, pooled = "random")

# Repeat analysis by removing studies that only tested <10 isolates
update(ma1, subset = total_invasive_isolates > 1)
update(ma1, subset = total_invasive_isolates > 2)
update(ma1, subset = total_invasive_isolates > 5)
update(ma1, subset = total_invasive_isolates > 10)

#--------------------------------------------------------------------------------------
# Asses funnel plot asymmetry and derive adjusted estimated using trim-and-fill method
#--------------------------------------------------------------------------------------
metabias(ma1)
trimfill(ma1)

#----------------------
# Look at Africa only
#----------------------
(afr1	<- metaprop(
			data = drugbug_africa1,
			number_resistant, 
			total_invasive_isolates,
			studlab = site_code
			)
		)
update.meta(afr1, byvar=year_cat, comb.random = TRUE, comb.fixed = F)

require(binom)
binom.confint(c(0),c(20), method="wilson")


# Warning of potential unstable estimates were returned
# Time period sub-group (2010s) from 1 study site had no events

drugbug_africa2 <- drugbug_africa1[which(drugbug_africa1$year_cat!="2010-2016"),]

# Carry out overall meta-analysis of proportion
(afr2	<- metaprop(
			data = drugbug_africa2,
			number_resistant, 
			total_invasive_isolates,
			studlab = site_code
			)
		)
update.meta(afr2, byvar=year_cat, comb.random = TRUE, comb.fixed = F)

#----------------------
# Look at Asia only
#----------------------

# Model failed to converge when applying the default logit transformation
# Use untransformed proportion for derivation

(asia1	<- metaprop(
			data = drugbug_asia1 ,
			number_resistant, 
			total_invasive_isolates,
			studlab = site_code,
			sm="PRAW"
			)
		)
update.meta(asia1, byvar=year_cat, comb.random = TRUE, comb.fixed = F)

## End Code