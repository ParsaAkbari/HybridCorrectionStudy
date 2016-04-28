setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis/5mer_csv_output/')
# loop through all the files here adding them to one massive csv file which we then save
overall <- read.table('./5mer_0001.csv', header=T, sep=",")
files <- list.files(path="./", pattern="*.csv", full.names=T, recursive=FALSE)
for (file in files) {
  currfile <- read.table(file, header=T, sep=",")
  overall <- rbind(overall, currfile)
}
setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis')
write.table(overall, '5mers.csv')
# sum across all columns
overall_colsum <- overall[,grep("^[before|after]", names(overall), value = TRUE)]
overall_colsum$read_name <- NULL
overall_colsum <- colSums(overall_colsum)
# plot each before & after combination as a datapoint
# we have a list of 5mers saved in 5mer_names.csv
# we just loop through these and pull columns from overall_colsum
fivemer_names <- read.csv('5mer_names.csv', header=F)
data_points <- as.vector(unlist(fivemer_names[1,]))
before = numeric(0)
after = numeric(0)
for (point in data_points) {
  before_colname <- paste("before_corr_", point, sep="")
  after_colname <- paste("after_corr_", point, sep="")
  before <- c(before,overall_colsum[before_colname])
  after <- c(after,overall_colsum[after_colname])
}
fivemer.df <- data.frame(dp=data_points, before=before, after=after)
fivemer.df$offset <- abs(1-(fivemer.df$before / fivemer.df$after))
library(ggplot2)
p <- ggplot(fivemer.df, aes(before,after))
p + geom_point(alpha=0.2)+ geom_text(aes(label=ifelse(offset>0.35 & (before > 25000 | after > 25000),as.character(dp),'')),hjust=0.4,vjust=-1.2,colour="red")


p + geom_point(alpha=0.09) + geom_text(aes(label=ifelse(offset>0.35, "@", "")), hjust=-0.01, vjust=0, color="red")
