setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis/')
overall <- read.csv('accuracy_channel_1.csv', header=T)
overall$isebi <- grepl("channel_", overall$read_name)
#overall <- overall[overall$isebi == TRUE,]
overall$read_identity_after <- (floor((overall$identity_after/100)*overall$alignment_length_after)/overall$read_length_after)*100
overall$read_identity_before <- (floor((overall$identity_before/100)*overall$alignment_length_before)/overall$read_length_before)*100

channel <- character(0)
mean_read_identity_before <- double(0)
isebi <- logical(0)
var_read_identity_before <- double(0)
# loop through all the unique values in the channel_number column
for (chan in unique(overall$channel_number)) {
  these.channels <- overall[overall$channel_number==chan,]
  channel <- c(channel, chan)
  mean_read_identity_before <- c(mean_read_identity_before,mean(these.channels$read_identity_before))
  var_read_identity_before <- c(var_read_identity_before, var(these.channels$read_identity_before))
  isebi <- c(isebi, grepl(".*_ebi", chan))
}
results <- data.frame(
  "channel"=channel,
  "mean_read_identity_before"=mean_read_identity_before,
  "var_read_identity_before"=var_read_identity_before,
  "isebi"=isebi)

library('ggplot2')
p <- ggplot(results, aes(mean_read_identity_before, 
                    fill=factor(isebi, labels=c("MISKA LAB", "EBI")), col=factor(isebi)))
p + geom_histogram(alpha=0.5) +
  labs(x='Channel Mean Read Identity (%)', y='Frequency', fill="Channel Origin:") +
  theme_bw(base_size = 19, base_family = "Helvetica") + theme(legend.position=c(0.854, 0.9))

p <- ggplot(results, aes(var_read_identity_before, 
                         fill=factor(isebi, labels=c("MISKA LAB", "EBI"))))
p + geom_histogram(alpha=0.5) +
  labs(x='Channel Variance Read Identity', y='Frequency', fill="Source of Data:") +
  theme_bw(base_size = 19, base_family = "Helvetica") + theme(legend.position=c(0.854, 0.9))


