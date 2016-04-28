setwd('/home/lm17/Documents/masters_project/c_elegans_data_analysis/')
overall <- read.table('accuracy_compare.csv', header = T)
overall$isebi <- grepl("channel_", overall$read_name)
library('ggplot2')
overall$read_identity_after <- (floor((overall$identity_after/100)*overall$alignment_length_after)/overall$read_length_after)*100
overall$read_identity_before <- (floor((overall$identity_before/100)*overall$alignment_length_before)/overall$read_length_before)*100
overall.ebi <- overall[overall$isebi == FALSE,]
overalln <- overall[overall$identity_after<61,]



bp <- ggplot(data=overall, aes(x=isebi, y=read_identity_before, fill=isebi)) +
  geom_boxplot() +
  scale_fill_discrete(name = "") +
  labs(x="", y="Read Identity before correction %")
# remove the legend
bp + theme_bw(base_size = 19, base_family = "Helvetica")

bp <- ggplot(data=overall, aes(x=isebi, y=read_identity_after, fill=isebi)) +
  geom_boxplot() +
  scale_fill_discrete(name = "") +
  labs(x="", y="Read Identity before correction %")
# remove the legend
bp + theme_bw(base_size = 19, base_family = "Helvetica")



p <- ggplot(overall, aes(y=read_identity_before, x=read_identity_after,
                         color=factor(isebi, labels=c("MISKA LAB", "EBI"))))+
   scale_color_manual(values=c("black", "red")) +
   guides(colour = guide_legend(override.aes = list(alpha = 1)))
p + geom_point(alpha=0.1) +
  labs(x="Read Identity after correction %", y="Read Identity before correction %",
       color="Source of Data") +
  geom_abline(slope = 1, intercept = 0, color="green")+
  xlim(0, 1) + ylim(0,1) +
  theme_bw(base_size = 19, base_family = "Helvetica") + theme(legend.position=c(0.09, 0.89))
  

# before_corr_read_length_accuracy
overall.lm <- lm(read_identity_before~read_length_before, overall.ebi)
p <- ggplot(overall, aes(y=read_identity_before, x=read_length_before,
                         color=factor(isebi, labels=c("MISKA LAB", "EBI"))))+
  scale_color_manual(values=c("black", "red")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
p + geom_point(alpha=0.1) +
  labs(x="Read length before correction (basepairs)", y="Read Identity before correction %",
       color="Source of Data") +
  #geom_abline(intercept=86.5365070, slope=-0.0001196, color="green") +
  theme_bw(base_size = 19, base_family = "Helvetica") + theme(legend.position=c(0.89, 0.9))


overall.lm <- lm(identity_after~read_length_after, overall)
p <- ggplot(overall, aes(y=read_identity_after, x=read_length_after,
                         color=factor(isebi, labels=c("MISKA LAB", "EBI"))))+
  scale_color_manual(values=c("black", "red")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
p + geom_point(alpha=0.1) +
  labs(x="Read length after correction (basepairs)", y="Read Identity after correction %",
       color="Source of Data") +
  stat_function(fun=function(x)100/x, geom="line", color="green")+
  #geom_abline(intercept = 60.8, slope=0) +
  #geom_abline(intercept=94.4312047, slope=-0.0008194, color="green") +
  theme_bw(base_size = 19, base_family = "Helvetica") + theme(legend.position=c(0.89, 0.9))


#Call:
#  lm(formula = identity_after ~ read_length_after, data = overall)
#
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-38.241  -2.466   2.852   5.567  19.569 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)        9.443e+01  6.515e-02 1449.44   <2e-16 ***
#  read_length_after -8.194e-04  1.676e-05  -48.88   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 8.383 on 41015 degrees of freedom
#Multiple R-squared:  0.05504,	Adjusted R-squared:  0.05502 
#F-statistic:  2389 on 1 and 41015 DF,  p-value: < 2.2e-16
