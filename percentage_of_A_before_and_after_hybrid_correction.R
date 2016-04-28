setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis/')
# load the csv of 2mers
overall <- read.table('1mer_compare.csv', header = T)
# change overall to just include ebi output or tomas output?
#overall <- overall[grep("channel_", overall$read_name),] # ebi output
overall<- overall[grep("^[^channel_]", overall$read_name),] # tomas o
# find totals of A T G C and N
overall$before_corr_total <- overall$before_corr_C + overall$before_corr_A +
  overall$before_corr_T + overall$before_corr_G + overall$before_corr_N
overall$after_corr_total <- overall$after_corr_C + overall$after_corr_A +
  overall$after_corr_T + overall$after_corr_G + overall$after_corr_N
# get percentages for each base pair before and after correction
percentages = data.frame(overall$read_name, overall$partition)
percentages$before_A <- overall$before_corr_A / overall$before_corr_total
percentages$after_A<- overall$after_corr_A / overall$after_corr_total
percentages$before_G<- overall$before_corr_G / overall$before_corr_total
percentages$after_G<- overall$after_corr_G / overall$after_corr_total
percentages$before_T <- overall$before_corr_T / overall$before_corr_total
percentages$after_T<- overall$after_corr_T / overall$after_corr_total
percentages$before_C <- overall$before_corr_C / overall$before_corr_total
percentages$after_C<- overall$after_corr_C / overall$after_corr_total
percentages$before_N <- overall$before_corr_N / overall$before_corr_total
percentages$after_N<- overall$after_corr_N / overall$after_corr_total

setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis/figures/')


# boxplot of A% before and after correction
after_col <- cbind(percentages$after_G*100, "After Correction")
before_col <- cbind(percentages$before_G*100, "Before Correction")
boxplot.df <- data.frame(rbind(after_col, before_col))
colnames(boxplot.df) <- c("percentage", "corrected")
boxplot.df$percentage <- as.numeric(as.character(boxplot.df$percentage))
library(ggplot2)
bp <- ggplot(data=boxplot.df, aes(x=corrected, y=percentage, fill=corrected)) +
  geom_boxplot() +
  scale_fill_discrete(name = "Nanocorr Correction") +
  labs(x="", y="Percentage of basepairs in read A %")
# remove the legend
bp + theme_bw(base_size = 19, base_family = "Helvetica")

###### frequency read length distribution ####
# y axis frequency
# x axis ONT read length
overall <- overall[overall$before_corr_total < 15000,]
before.col <- data.frame(total=overall$before_corr_total, corr="No")
after.col <- data.frame(total=overall$after_corr_total, corr="Yes")
comb <- rbind(before.col, after.col)
hist <- ggplot(comb, aes(total, fill=corr)) + 
  geom_histogram(alpha=0.2)
hist + 
  labs(x="Read Length", y="Frequency") +
  theme_bw(base_size = 19, base_family = "Helvetica")
library(devtools)
library(easyGgplot2)
ggplot2.histogram(data=comb, xName='total',
                  groupName='corr', legendPosition="top",
                  alpha=0.5, position="identity",
                  legendTitle="Corrected")+
  labs(x='Read Length', y='Frequency') +
  theme_bw(base_size = 19, base_family = "Helvetica") + theme(legend.position=c(0.9, 0.9))
