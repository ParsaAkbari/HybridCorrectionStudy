setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis/gcbias/')
I <- read.table('I.txt')
II <- read.table('II.txt')
III <- read.table('III.txt')
IV <- read.table('IV.txt')
V <- read.table('V.txt')
X <- read.table('X.txt')
coverage <- rbind(I, II, III, IV, V, X)
names(coverage) <- c("cov")
rm(I, II, III, IV, V, X)
I <- read.csv('I_gc_percentage.csv', header = T)
II <- read.csv('II_gc_percentage.csv', header = T)
III <- read.csv('III_gc_percentage.csv', header = T)
IV <- read.csv('IV_gc_percentage.csv', header = T)
V <- read.csv('V_gc_percentage.csv', header = T)
X <- read.csv('X_gc_percentage.csv', header = T)
gcpercentage <- rbind(I, II, III, IV, V, X)
gcpercentage$gc_percentage = gcpercentage$gc_percentage
overall <- cbind(gcpercentage, coverage)
# 38652 obvs
rm(I, II, III, IV, V, X, gcpercentage, coverage)
library(ggplot2)

overall.lm <- lm(gc_percentage~cov, data=overall)
#20 = quantile(overall$cov, 0.99979948772)
overall <- overall[overall$cov < 11,]
# 38641 obvs (8 windows gone)
p <- ggplot(overall, aes(cov, gc_percentage))
p + geom_point(alpha=0.2) +
  labs(x="Coverage", y="GC Percentage") +
  geom_abline(intercept=35.3197, slope=0.067730, color="green") +
  theme_bw(base_size = 19, base_family = "Helvetica")
