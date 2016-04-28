setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis/')
overall <- read.table('gcbias_illumina.csv', header=T)
# we have one NA value for some reason
overall <- overall[overall$read_name!='channel_59_read_9_twodirections',]
# boxplot of A% before and after correction
kept_gc <- cbind(overall$kept_gc, "Kept region")
chucked_gc <- cbind(overall$chucked_gc, "Trimmed region")
boxplot.df <- data.frame(rbind(kept_gc, chucked_gc))
colnames(boxplot.df) <- c("percentage", "kept")
boxplot.df$percentage <- as.numeric(as.character(boxplot.df$percentage))
library(ggplot2)
bp <- ggplot(data=boxplot.df, aes(x=kept, y=percentage, fill=kept)) +
  geom_boxplot() +
  scale_fill_discrete(name = "") +
  labs(x="", y="Percentage of GC %")
# remove the legend
bp + theme_bw(base_size = 19, base_family = "Helvetica")

t.test(overall$kept_gc, overall$chucked_gc)
mean(overall$chucked_gc)
mean(overall$kept_gc)
