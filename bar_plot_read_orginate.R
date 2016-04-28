setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis/')
overall <- read.table('1mer_compare.csv', header = T)
tomas <- overall[grep("^[^channel_]", overall$read_name),]$read_name # tomas o
ebi <- overall[grep("channel_", overall$read_name),]$read_name
tomas <- data.frame(read_name=tomas, so ="Miska Lab")
ebi <- data.frame(read_name=ebi, so ="EBI")

total <- rbind(tomas,ebi)

c <- ggplot(total, aes(x=factor(so), fill=factor(so)))
c + geom_bar()+
  labs(x='Read Source', y='Number of Reads', fill="") +
  theme_bw(base_size = 19, base_family = "Helvetica")

