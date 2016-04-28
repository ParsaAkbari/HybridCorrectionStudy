setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis/')
fivemer_names <- read.table('6mer_names.csv', header=FALSE, sep=",")
fivemer_names <- as.vector(unlist(fivemer_names))
overall <- read.table('6mers.csv', header = T)
#overall <- overall[grep("channel_", overall$read_name),] # ebi output
#overall<- overall[grep("^[^channel_]", overall$read_name),] # tomas output

before <- numeric(0)
after <- numeric(0)
for (fivemer in fivemer_names) {
  before_point <- paste("before_corr_", fivemer, sep="")
  after_point <- paste("after_corr_", fivemer, sep="")
  before <- c(before, sum(overall[,before_point]))
  after <- c(after, sum(overall[,after_point]))
}
fivemer.df <- data.frame(df=fivemer_names, before=before, after=after)
rm(overall, after, before, before_point, after_point, fivemer, fivemer_names)
write.csv(fivemer.df, "sixmer.csv")
fivemer.df <- read.csv("sixmer.csv")
fivemer.df$change <- fivemer.df$after - fivemer.df$before
fivemer.df$percent_change <- (fivemer.df$change/fivemer.df$before)*100

# top four infinite values removed
percent.change.descending <- sort(fivemer.df$percent_change, decreasing=TRUE)[-1:-13]
removed_values <- sort(fivemer.df$percent_change, decreasing=TRUE)[1:13]
x <- 1:length(percent.change.descending)
data <- data.frame("x"=x, "percent.change.descending"=percent.change.descending)
p <- ggplot(data, aes(x, percent.change.descending))
p + geom_point() +
  geom_abline(intercept = 0, slope = 0, colour="gray") +
  labs(y="Change in percentage enrichment") +
  theme(axis.ticks.x=element_blank()) +
  theme_bw(base_size = 19, base_family = "Helvetica")

library('ggplot2')
library('Hmisc')
library(ggrepel)
fivemer.df$residuals <- abs(log(fivemer.df$before) - log(fivemer.df$after))/sqrt(1)
p <- ggplot(fivemer.df, aes(log(before), log(after)))
p + geom_point(alpha=0.2) +
  geom_abline(intercept = 0, slope = 1, colour="green") +
  geom_text(aes(label=ifelse(abs(residuals) > 1.5,as.character(df), '')),
            hjust=0, vjust=0,colour="black", fontface="bold") +
  labs(x="log(6-mer count before correction)", y="log(6-mer count after correction)") +
  xlim(3, 11) + ylim(3,11) +
  theme_bw(base_size = 19, base_family = "Helvetica")

top.residuals <- fivemer.df[order(fivemer.df$residuals,decreasing=TRUE)[1:5],] 
rm(fivemer.df)
write.csv(top.residuals, "top6mer.csv")
setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis/')
top.residuals <- read.csv("top6mer.csv")
library(ggplot2)
library(plyr)
library(reshape2)
top.residuals <- rename(top.residuals, c("before"="Before Correction", "after"="After Correction"))

top.residuals <- melt(top.residuals, id.vars=c("X", "df", "residuals"))
ggplot(top.residuals, aes(df, value)) +   
  geom_bar(aes(fill = variable), position = "dodge", stat="identity")+
  labs(x="k-mer", y="Count", fill="") +
  theme_bw(base_size = 19, base_family = "Helvetica")
