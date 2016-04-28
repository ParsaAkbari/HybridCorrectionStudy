setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis/')
overall <- read.table('accuracy_compare.csv', header = T)
overalln <- overall[overall$identity_after<61,]
rm(overall)
write.csv(overalln$read_name, "low_accuracy.csv")
