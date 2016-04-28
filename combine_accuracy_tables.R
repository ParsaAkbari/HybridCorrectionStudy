setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis/accuracy_compare/')
# loop through all the files here adding them to one massive csv file which we then save
#partition,read_name,identity_after,alignment_length_after,identity_before,alignment_length_before
overall <- data.frame(partition=numeric(0), read_name=character(0), identity_after=double(0),
                      alignment_length_after=numeric(0), read_length_after=numeric(0),
                      identity_before=double(0), alignment_length_before=numeric(0),
                      read_length_before=numeric(0))
files <- list.files(path="./", pattern="*.csv", full.names=T, recursive=FALSE)
for (file in files) {
  currfile <- read.csv(file, header=T)
  overall <- rbind(overall, currfile)
}
setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis')
write.table(overall, 'accuracy_compare.csv')

ab <- read.table('accuracy_compare.csv', header = T)
