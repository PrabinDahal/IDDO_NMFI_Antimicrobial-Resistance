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
require(maps)
require(rworldmap)
require(rworldxtra)
require(maptools)
require(classInt)
require(RColorBrewer)
require(doBy)
require(dplyr)
require(plyr) # for revalue
require(eply)
require(stringr)
require(reshape)
require(gtools) # for smartbind
require(ggpubr)
require(car)
require(patchwork)
require(utf8)
require(ggalt)  
require(viridis)
require(hrbrthemes)
require(tableone)
require(cowplot)

# meta-analysis packages
require(meta)
require(metafor)
require(metasens) 

#-----------------------
# Read analysis dataset 
#-----------------------
dat0<-read.csv("D:/_IDDO/_nmfi_amr/NMFI AMR/Data/nmfi_amr_analysis_set_final.csv")
dat<- dat0

#-------------------------------------------------
# Number of unique studies by study design
#-------------------------------------------------
unique<-dat[which(!duplicated(dat$article_id)),]
nrow(unique)

#----------------------------------------------------------
# Total number of unique articles by country
#----------------------------------------------------------
(country <- dat %>% 
	dplyr::group_by(country_name) %>%
	dplyr::summarise(number_of_articles = length(unique(article_id)))
)

CreateTableOne(
		vars 	= c("antibiotic_group"), 
		factorVars 	= c("sub_region"), 
		strata 	= "sub_region", 
		data 		= dat[!duplicated(dat$article_id),], 
		test		= FALSE
		)

CreateTableOne(
		vars 	= c("country_name"), 
		factorVars 	= c("sub_region"), 
		#strata 	= "sub_region", 
		data 		= dat[!duplicated(dat$article_id),], 
		test		= FALSE
		)

CreateTableOne(
		vars 	= c("sub_region"), 
		factorVars 	= c("sub_region"), 
		#strata 	= "sub_region", 
		data 		= dat[!duplicated(dat$article_id),], 
		test		= FALSE
		)

# Time-period of sample collection
dat %>% 
	dplyr::group_by(sub_region) %>%
	dplyr::summarise(
		min = min(study_start_year),
		max = max(study_end_year)
	)

#----------------------------------------------------------
# Number of articles from each country
#----------------------------------------------------------
ggplot(country, aes(x=reorder(country_name, number_of_articles), y=number_of_articles)) + 
	geom_point(size=3,stroke =2,col="#00AFBB")+
  	#geom_pointrange(aes(ymin=number_of_articles, ymax=number_of_articles),size=1,stroke =1,col="#00AFBB")+
  	#geom_pointrange(aes(ymin=re_l95, ymax=re_u95),size=1,stroke =1,col="#00AFBB")+
	ylab("Number of articles") +
	ylim(0,150)+
	xlab("")+
	ggtitle("Number of unique articles per country")+
	coord_flip()

#------------------------------
# Identify multi-centre studies
#------------------------------
n_sites <- dat %>% 
	dplyr::group_by(article_id) %>%
	dplyr::summarise(n_sites = length(unique(site_code))
	)

# Merge information on if the study was single or multi-centre
dat <- merge(dat, n_sites, by="article_id")

CreateTableOne(
		vars 	= c("n_sites"), 
		factorVars 	= c("n_sites"), 
		#strata 	= "sub_region", 
		data 		= n_sites,
		test		= FALSE
		)

CreateTableOne(
		vars 	= c("sub_region"), 
		factorVars 	= c("sub_region"), 
		#strata 	= "sub_region", 
		data 		= dat[!duplicated(dat$site_code),], 
		test		= FALSE
		)

(country <- dat %>% 
	dplyr::group_by(country_name) %>%
	dplyr::summarise(number_of_articles = length(unique(article_id)))
)

#----------------------------------------------------------
# Total number of unique articles by publication year
#----------------------------------------------------------
(year <- dat %>% 
	dplyr::group_by(year_published, sub_region ) %>%
	dplyr::summarise(number_of_articles = length(unique(article_id)))
)

(a <- ggbarplot(year, "year_published", "number_of_articles",
	fill = "sub_region", 
	#palette = "Paired",
	label = F, 
	lab.col = "white", lab.pos = "in")+
	ggtitle("A: Region")+
	ylab("Number of articles") +
	xlab("") +
	labs(fill = "")+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
	theme(legend.title = element_text(size = 8), 
               legend.text = element_text(size = 8))
)

#----------------------
# By AST Guidelines 
#----------------------
(year_guideline <- dat %>% 
	dplyr::group_by(year_published, amr_guideline) %>%
	dplyr::summarise(number_of_articles = length(unique(article_id)))
)

(b <- ggbarplot(year_guideline , "year_published", "number_of_articles",
	fill = "amr_guideline", 
	#palette = "Paired",
	label = F, 
	lab.col = "white", lab.pos = "in")+
	ggtitle("B: AST guidelines")+
	ylab("") +
	xlab("") +
	labs(fill= "")+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
	theme(legend.title = element_text(size = 6), 
               legend.text = element_text(size = 6))
)

#----------------------
# By AST testing 
#----------------------
(year_testing <- dat %>% 
	dplyr::group_by(year_published, amr_testing) %>%
	dplyr::summarise(number_of_articles = length(unique(article_id)))
)

(c <- ggbarplot(year_testing , "year_published", "number_of_articles",
	fill = "amr_testing", 
	#palette = "Paired",
	label = F, 
	lab.col = "white", lab.pos = "in")+
	ggtitle("C: AST method")+
	ylab("Number of articles") +
	xlab("Year published") +
	labs(fill= "")+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +	
	theme(legend.title = element_text(size = 8), 
      legend.text = element_text(size = 8))
)

#--------------------------------
# By study quality scores over time
#--------------------------------
year_micro <- dat[!duplicated(dat$article_id),]
summary(year_micro$Final_result)

(d <- ggplot(year_micro , aes(x=year_published, y=Final_result)) + 
	geom_jitter(col="salmon",shape=1, size=2.5)+
  	geom_smooth()+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
	#scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
	xlab("Year published") +
	ylab("Total scores") +
	ylim(0,15)+
	ggtitle("D: MICRO grading")+
	theme_minimal()+
	theme(axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=12, color="black"))
	+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black"))
	)

plot_grid(a, b,c,d)

#--------------------------------
# By Glass Priority Pathogens 
#--------------------------------
glass_pathogens <- dat[which(dat$glass=="glass"),]
table(glass_pathogens$pathogen_2)

(year_glass <- glass_pathogens %>% 
	dplyr::group_by(year_published, pathogen_2) %>%
	dplyr::summarise(number_of_articles = length(unique(article_id)))
)

(d <- ggbarplot(year_glass, "year_published", "number_of_articles",
	fill = "pathogen_2", 
	#palette = "Paired",
	label = F, 
	lab.col = "white", lab.pos = "in")+
	ggtitle("GLASS priority pathogens")+
	ylab("Number of articles") +
	xlab("Year published") +
	#labs(fill= "")+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +	
	theme(legend.title = element_text(size = 6), 
      legend.text = element_text(size = 6))
 )

#---------------------
# Age-group table
#---------------------
CreateTableOne(
		vars 	= c("age_label1"), 
		factorVars 	= c("age_label1"), 
		strata 	= "sub_region", 
		data 		= dat[!duplicated(dat$article_id),], 
		test		= FALSE
		)

#----------------------------------------------------------
# sample source used for confirmed presence of pathogens 
#----------------------------------------------------------
CreateTableOne(
		vars 	= c("sample_cat"), 
		factorVars 	= c("sample_cat"), 
		strata 	= "sub_region", 
		data 		= dat[!duplicated(dat$article_id),], 
		test		= FALSE
		)

#---------------------
# AMR Guidelines used
#---------------------
CreateTableOne(
		vars 	= c("amr_guideline"), 
		factorVars 	= c("amr_guideline"), 
		#strata 	= "sub_region", 
		data 		= dat[!duplicated(dat$article_id),], 
		test		= FALSE
		)

#---------------------
# AMR Testing 
#---------------------
CreateTableOne(
		vars 	= c("amr_testing"), 
		factorVars 	= c("amr_testing"), 
		strata 	= "sub_region", 
		data 		= dat[!duplicated(dat$article_id),], 
		test		= FALSE
		)

#----------------------------------------------------------
# Total number of unique articles by country and sub-region
#----------------------------------------------------------
(country_region <- dat %>% 
	dplyr::group_by(sub_region,country_name) %>%
	dplyr::summarise(number_of_articles = length(unique(article_id)))
)

ggdotchart(	dat		=	country_region, 
		x		=	"country_name", 
		y		=	"number_of_articles",
           	color 	= 	"sub_region",        
           	sorting 	= "descending",          
           	add 		= "segments",                
           	rotate 	= TRUE,                   
           	group 	= "sub_region", 
                      	dot.size 	= 5,                    
           	label 	= country_region$number_of_articles,                        
           	font.label 	= list(
				color = "white", size = 8, vjust = 0.5),               
           			ggtheme = theme_pubr())+
		ylab("Number of articles") +
		xlab("")+
		ggtitle("")+
		theme(axis.text.x = element_text(size=10, angle=0),
		          axis.text.y = element_text(size=10),
			 axis.title=element_text(size=14,face="bold"),
		plot.title = element_text(size = 20, face = "bold")
		)+
theme(legend.title = element_blank()) 

#--------------------------------------------
# Plot all the study sites in a single map 
#--------------------------------------------
unique_lat_lon<-dat[!duplicated(paste(dat$google_lon,dat$google_lat)),] 

asia_latlon <- unique_lat_lon[which(unique_lat_lon$sub_region %in% c("South-eastern Asia","Southern Asia")),]
africa_latlon <- unique_lat_lon[which(unique_lat_lon$sub_region %in% c("Northern Africa","Sub-Saharan Africa")),]

worldmap <- getMap(resolution = "high")
asia 	<- worldmap[which(worldmap$REGION=="Asia"),]      
africa <- worldmap[which(worldmap$REGION=="Africa"),]      

#tiff(file="Figure 1_map.tiff", 
#            width=17, height=16, units="cm", 
#           pointsize="7", compression = "lzw+p", 
#            bg="white", res=500, antialias = "none" )

par(mfrow=c(1,2))
plot(africa, lwd=2, main="", cex.main=2)
points(africa_latlon$google_lon,africa_latlon$google_lat,col=rgb(0.4,0.4,0.8,0.6), pch=19 , cex=1.3)

plot(asia, lwd=2, main="", cex.main=2)
points(asia_latlon$google_lon,asia_latlon$google_lat,col=rgb(0.4,0.4,0.8,0.6), pch=19 , cex=1.3)

#dev.off()

#=======================
# Drug-bug network graph
#=======================
bug<-as.data.frame(table(dat$pathogen))

drug_bug_model = list()
for (r in 1:nrow(bug)){
		k1<-dat[which(dat$pathogen==bug$Var1[[r]]),]
		k1<-droplevels(k1)
		k2<-as.data.frame(table(k1$antibiotic_group))
		OUT<-NULL;
		for (i in 1:nrow(k2)){
			k3<-k1[which(k1$antibiotic_group==k2$Var1[[i]]),]
			k3<-k3[which(!duplicated(k3$article_id)),]
			k3<-droplevels(k3)
			tmp<- merge(k2$Var1[[i]], nrow(k3))
			tmp$bug<-bug$Var1[[r]]
			OUT<-rbind(tmp,OUT)
		}
	drug_bug_model[[r]] <-OUT
	}
drug_bug_model_articles = do.call(rbind, drug_bug_model)	
colnames(drug_bug_model_articles)<-c("drug","total_articles","bug")

drug_bug_model_articles <- drug_bug_model_articles[order(drug_bug_model_articles$total_articles),]

#-----------------------------
# create bipartite network
#-----------------------------
library(reshape2)
library(multigraph)

test 	<- drug_bug_model_articles 
test  <- test [which(test$drug!="NULL"),] 
test  <- test [which(test$bug!="NULL"),] 
test 	<- droplevels(test )

# Keep combination with min of 10 articles for graphics
test <- test [order(test$total_articles),]
test1 <- test [which(test$total_articles >=10),] 
test1 <- droplevels(test1 )

test2 <- dcast(test1 , bug ~ drug, value.var="total_articles", fun.aggregate=sum)
row.names(test2) <- test2$bug
test2$bug <-NULL

test2<-as.matrix(test2)

setwd("D:/_IDDO/_nmfi_amr/NMFI AMR/Results")

#tiff(file="network_graph.tiff", 
#            width=18, 
#		height=12, 
#		units="cm", 
#           pointsize="8", 
#		compression = "lzw+p", 
#            bg="white", 
#		res=600, 
#		antialias = "none" 
#		)

bmgraph(test2, 
		layout ="bipc",
		asp = NA,cex=2.5,ecol="#00AFBB", vcol="#00AFBB",lwd=1.5,  pch=21,
		weighted =TRUE,
		fstyle=c("plain")
		)

#dev.off()

### End Script