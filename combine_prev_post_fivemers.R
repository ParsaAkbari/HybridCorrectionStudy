setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis/prev_post_csv_output/')
# loop through all the files here adding them to one massive csv file which we then save
overall <- data.frame(fivemer_prev=character(), fivemer_post=character(), read_name=character())
files <- list.files(path="./", pattern="*.csv", full.names=T, recursive=FALSE)
for (file in files) {
  currfile <- read.table(file, header=T, sep=",")
  overall <- rbind(overall, currfile)
}
setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis')

library(plyr)
fivemer_prev_counts <- count(overall, "fivemer_prev")
fivemer_post_counts <- count(overall, "fivemer_post")
xaxis = numeric(0)
yaxis = numeric(0)
data_point_mer = character(0)
for (mer in fivemer_prev_counts$fivemer_prev) {
  xaxis <- c(xaxis, fivemer_post_counts[which(fivemer_post_counts$fivemer_post==mer),]$freq)
  yaxis <- c(yaxis, fivemer_prev_counts[which(fivemer_prev_counts$fivemer_prev==mer),]$freq)
  data_point_mer <- rbind(data_point_mer, mer)
}
library('ggplot2')
library('ggrepel')
beforeafter.df <- data.frame(x = xaxis, y=yaxis, currmer=as.vector(as.character(data_point_mer[,1])))
p <- ggplot(beforeafter.df, aes(x, y))
p + 
  geom_point(alpha=0.4, colour="black") +
  geom_text_repel(
    data=subset(beforeafter.df, x > 500 | y > 500),
    aes(label=currmer),
    fontface = 'bold',
    color="black",
    size=5
  ) +
  labs(x="Proceeding 5mer", y="Preceeding 5mer") +
  theme_bw(base_size = 19, base_family = "Helvetica")


## MEASURE ENRICHMENT OF 5MERS PROCEEDING/PRECEEDING AND NOT RAW VALUES
fivemer_numbers <- read.csv('before5mers.csv')
proceeding <- double(0)
preceeding <- double(0)
fivemer <- character(0)
for (kmer in fivemer_numbers$fivemer_names) {
  getrow <- beforeafter.df[beforeafter.df$currmer==kmer,]
  gettotal <- fivemer_numbers[fivemer_numbers$fivemer_names == kmer,]$before
  proceeding <- c(proceeding, (getrow$x*100/gettotal))
  preceeding <- c(preceeding, (getrow$y*100/gettotal))
}
total.df <- data.frame("proceeding"=proceeding, "preceeding"=preceeding, kmer=fivemer_numbers$fivemer_names)

# enriched proceeding
top.proceeding <- total.df[order(total.df$proceeding,decreasing=TRUE)[1:6],] 
# enriched preceeding
top.preceeding <- total.df[order(total.df$preceeding,decreasing=TRUE)[1:6],] 

p <- ggplot(total.df, aes(proceeding, preceeding))
p + 
  geom_point(alpha=0.4, colour="black") +
  geom_text_repel(
    data=subset(total.df, proceeding > 2.5 | preceeding > 2.5),
    aes(label=kmer),
    fontface = 'bold',
    color="red",
    size=5
  ) +
  labs(x="Percentage of occurances proceeding trim site", y="Percentage of occurances preceeding trim site") +
  theme_bw(base_size = 19, base_family = "Helvetica")

quantile(total.df$proceeding, c(0.99307136))

