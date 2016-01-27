#Written by Catrina Loucks, 2016-01-08

#download and load required libraries

if(!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}

if(!require("tidyr")) {
  install.packages("tidyr")
library(tidyr)
}

if(!require("plyr")) {
  install.packages("plyr")
  library(plyr)
}

if(!require("binom")) {
  install.packages("binom")
  library(binom)
}

args <- commandArgs(trailingOnly = TRUE)
buzz <- args[1]
press <- args[2]
basal <- args[3]

#load in data
buzz.data <- read.table(buzz, sep="\t", header = TRUE)
press.data <- read.table(press, sep="\t", header = TRUE)
basal.data <- read.table(basal, sep="\t", header = TRUE)

#tidy data (make a wide data format long)
buzz.tidy <- gather(buzz.data, strain, peak.ratio.change) #(data,key,value) 
press.tidy <- gather(press.data, strain, peak.ratio.change) #(data,key,value) 
basal.tidy <- gather(basal.data, strain, basal.fluorescence.ratio) #(data,key,value) 

#add a column with type of stimulus
buzz.tidy$stimulus <- "buzz"
press.tidy$stimulus <- "press"
  
#remove rows with NA values
buzz.tidy <- buzz.tidy[complete.cases(buzz.tidy),]
press.tidy <- press.tidy[complete.cases(press.tidy),]
basal.tidy <- basal.tidy[complete.cases(basal.tidy),]
calcium.tidy <- rbind(buzz.tidy,press.tidy)

#optional - remove values for a strain 
#calcium.tidy <- calcium.tidy[calcium.tidy$strain != "efhc.1",]


#Plot proportion of worms reponding to each stimulus
#calculate N and number of worms responding
calcium.summary <- ddply(calcium.tidy,
                         .(strain,stimulus),
                         summarise,
                         N=length(peak.ratio.change), #number of animals for each genotype
                         respond=sum(peak.ratio.change != 0.000000), #number of worms responding
                         proportion=respond/N, #number responding divided by number assessed
                         conf.int.lower = binom.confint(respond, N, methods = "exact")$lower,
                         conf.int.upper = binom.confint(respond, N, methods = "exact")$upper)

#plot data as proportion responding
barplot <- ggplot(calcium.summary, aes(x=stimulus, y=proportion, fill=strain)) +
  geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat="identity") +
  geom_errorbar(position=position_dodge(width=0.9), aes(ymin=conf.int.lower, ymax=conf.int.upper), width=0.4) +
  scale_y_continuous(limits=c(0.0, 1.0), breaks=c(0.0,0.2,0.4,0.6,0.8,1.0), expand =c(0,0)) + ##Set the y-axis limits to a range from 0 to 40, y-values every 10, and remove extra space above and below
  theme_classic() + ## remove background to make it white
  theme(legend.position="none") +
  labs(x=NULL, y="Proportion responding")

barplot

#stats
#BUZZ
## fit a logistic regression with proportion modeled against strain
my.glm.buzz <- glm(proportion ~ strain, weights = N, family = binomial(link = "logit"), 
                data = calcium.summary[calcium.summary$stimulus =="buzz",])
  
  ## Get the logistic regression summary statistics
  summary(my.glm.buzz)
  
  ## load the library to perform multiple comparisons
  library(multcomp)
  
  ## perform a Tukey's HSD multiple comparison posthoc test
  glht.my.glm.buzz <- glht(my.glm.buzz, mcp(strain="Tukey"))
  summary(glht.my.glm.buzz) 
  
#PRESS
  ## fit a logistic regression with proportion modeled against strain
  my.glm.press <- glm(proportion ~ strain, weights = N, family = binomial(link = "logit"), 
                     data = calcium.summary[calcium.summary$stimulus =="press",])
  
  ## Get the logistic regression summary statistics
  summary(my.glm.press)
  
  ## load the library to perform multiple comparisons
  library(multcomp)
  
  ## perform a Tukey's HSD multiple comparison posthoc test
  glht.my.glm.press <- glht(my.glm.press, mcp(strain="Tukey"))
  
  summary(glht.my.glm.press) 
  
  #save barplot
  ggsave(barplot, file="proportion_responding.pdf", useDingbats=FALSE, height=4, width=6, units="in", dpi=300)


  #make a data.frame including only values that are non-zero (responders)
  calcium.non.zero <- calcium.tidy[calcium.tidy$peak.ratio.change>0.000000,]
  calcium.non.zero <- calcium.non.zero[calcium.non.zero$re.order.x !="efhc.1-buzz",] #remove efhc-1 buzz because it only responded twice (don't include in boxplot)
  calcium.non.zero <- calcium.non.zero[calcium.non.zero$re.order.x !="trp.4-buzz",] #remove trp-4 because only responded once like above

  
#scatterplot  
scatterplot <- ggplot(calcium.tidy, aes(x=as.factor(stimulus), y=peak.ratio.change, fill=strain, colour=strain, shape=strain)) +
  geom_point(size=2,alpha=0.6, position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0.5, dodge.width = 0.9)) +
  geom_boxplot(alpha=0.2, outlier.size = 0, outlier.colour=NA, data=calcium.tidy, position=position_dodge(width=0.9))  +
  scale_y_continuous(limits=c(-1, 120), breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120)) + ##Set the y-axis limits to a range from 0 to 40, y-values every 10, and remove extra space above and below
  theme_classic() + ## remove background to make it white
  theme(legend.position="none") +
  labs(x=NULL, y="Peak ratio change (%)")

scatterplot

ggsave(scatterplot, file="scatterplot_peak_ratio_change.pdf", useDingbats=FALSE, height=4, width=6, units="in", dpi=300)

#scatterplot for basal fluorescence




