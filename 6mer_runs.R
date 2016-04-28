setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis/')
library("ggplot2")
sixmer.df <- read.csv("sixmer.csv")
sixmer.df$change <- sixmer.df$after - sixmer.df$before
sixmer.df$percent_change <- (sixmer.df$change/sixmer.df$before)*100
sixmer.df <- sixmer.df[is.finite(sixmer.df$percent_change),]
three.t <- sixmer.df[grep("TTT", sixmer.df$df),]
three.t$bp <- "TTT"
three.a <- sixmer.df[grep("AAA", sixmer.df$df),]
three.a$bp <- "TTT"
three.g <- sixmer.df[grep("GGG", sixmer.df$df),]
three.g$bp <- "GGG"
three.c <- sixmer.df[grep("CCC", sixmer.df$df),]
three.c$bp <- "GGG"
overall <- rbind(three.t, three.a, three.g, three.c)
overall$bp = factor(overall$bp)
p <- ggplot(overall, aes(percent_change)) +
  geom_histogram(data=three.g, fill="green", color="green", alpha=0.7) +
  geom_histogram(data=three.t, fill="red", color="red", alpha=0.7)
p + 
  labs(x='6-mer', y='Percentage Enrichment') +
  theme_bw(base_size = 19, base_family = "Helvetica")