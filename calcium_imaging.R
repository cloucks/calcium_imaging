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

#load in data
buzz.data <- read.table(buzz, sep="\t", header = TRUE)
press.data <- read.table(press, sep="\t", header = TRUE)

#tidy data (make a wide data format long)
buzz.tidy <- buzz.data %>%
  gather(strain, peak.ratio.change) #(key,value) 

#add a column with type of stimulus
buzz.tidy$stimulus <- "buzz"
  
#remove rows with NA values
calcium.tidy <- calcium.tidy[complete.cases(calcium.tidy),]
calcium.non.zero <- calcium.tidy[calcium.tidy$peak.ratio.change>0.000000,]


#calculate N and number of worms responding
calcium.summary <- ddply(calcium.tidy,
                         .(strain),
                         summarise,
                         N=length(peak.ratio.change), #number of animals for each genotype
                         respond=sum(peak.ratio.change != 0.000000), #number of worms responding
                         proportion=respond/N, #number responding divided by number assessed
                         confintupper <- binom.confint(respond, N, methods = "exact")$upper,                                            confintlower <- binom.confint(respond, N, methods = "exact")$lower)

#plot data as a scatterplot
scatterplot <- ggplot(calcium.tidy, aes(x=as.factor(strain), y=peak.ratio.change, colour=strain)) +
  geom_point(alpha=0.6, position=position_jitter(width=0.3, height=0.5)) +
  geom_boxplot(alpha=0.2, outlier.size = 0, outlier.colour=NA, data=calcium.non.zero)  +
  scale_y_continuous(limits=c(-1, 40), breaks=c(0,10,20,30,40)) + ##Set the y-axis limits to a range from 0 to 40, y-values every 10, and remove extra space above and below
  theme_classic() + ## remove background to make it white
  theme(legend.position = "none") + ##remove legend
  labs(x="Strain", y="Peak ratio change (%)")

scatterplot

ggsave(scatterplot, file="scatterplot_peak_ratio_change.pdf", useDingbats=FALSE, height=4, width=6, units="in", dpi=300)

#plot data as proportion responding
barplot <- ggplot(calcium.summary, aes(x=strain, y=proportion)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=..5, ymax=..4), width=0.5) +
  scale_y_continuous(limits=c(0.0, 1.0), breaks=c(0.0,0.2,0.4,0.6,0.8,1.0)) + ##Set the y-axis limits to a range from 0 to 40, y-values every 10, and remove extra space above and below
  scale_color_discrete(name="") + ## remove legend title
  theme(legend.key = element_blank()) + ## remove boxes around legend values
  theme_classic() + ## remove background to make it white
  labs(x="Strain", y="Peak ratio change (%)")

barplot

ggsave(scatterplot, file="scatterplot_peak_ratio_change.pdf", useDingbats=FALSE, height=4, width=6, units="in", dpi=300)




