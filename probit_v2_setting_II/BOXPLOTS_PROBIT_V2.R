###########################
##### BOXPLOTS PROBIT #####
###########################
#####################
##### BOTH NNTs #####
#####################
library(reshape2)
library(ggplot2)

### READING FROM FILE ###
bbb_both2 <- read.csv("NNTm1000_PROBIT_5_3.csv", header = T)
bbb_both2 <- bbb_both2[bbb_both2$TYPE == "IV" & 
                         bbb_both2$NNT > 1  & 
                         bbb_both2$NNT < 10 &
                         !is.na(bbb_both2$n), ]
                 
ggplot(data = bbb_both2, 
       aes(x = factor(n), y = NNT, fill = factor(TYPE))) + 
  geom_boxplot() +
  geom_hline(yintercept = 5.3, linetype = "dashed",
             col        = "red", size = 1)     +
  xlab("Sample size")        + 
  ylab("NNT")           +
  ggtitle("IV-based Estimators of the NNT as a 
         Function of the Sample Size")  +
  scale_y_continuous(breaks = 1:10)          +
  theme_minimal() + 
  scale_fill_discrete(name = "NNT Type") + 
  theme(plot.title  = element_text(hjust = 0.5,     size = 25),
        axis.text.x = element_text(  size = 20),
        axis.text.y = element_text(  size = 20), 
        axis.title  = element_text(size  = 20), 
        legend.key.size = unit(1, 'cm'), 
        legend.text  = element_text(size = 20), 
        legend.title = element_text(size = 20))

ggsave("NNT_5_3_PROBITm1000R1.jpeg")
       
###############################
### EIN BOXPLOTS - ONLY ADJ ###
###############################
### READING FROM FILE ###
bbb2 <- read.csv("EINm1000_PROBIT_4_496.csv")

bbb2$EIN <- ifelse(bbb2$EIN > 10 | bbb2$EIN < 1, NA, bbb2$EIN)
bbb2 <- bbb2[bbb2$TYPE == "IV",]

ggplot(data = bbb2, 
       aes(x = factor(n), y = EIN, fill = factor(TYPE))) +
  geom_boxplot(fill = c("#F8766D")) +
  geom_hline(yintercept = 4.496, linetype = "dashed", col = "red", size = 1)     +
  xlab("Sample size")        + 
  ylab("EIN")           +
  ggtitle("IV-based Estimators of the EIN as a \n Function of the Sample Size")  +
  scale_y_continuous(breaks = 1:10)          +
  theme_minimal() + 
  theme(plot.title  = element_text(hjust = 0.5,     size = 25),
        axis.text.x = element_text(  size = 20),
        axis.text.y = element_text(  size = 20), 
        axis.title  = element_text(size  = 20), 
        legend.key.size = unit(1, 'cm'), 
        legend.text  = element_text(size = 20), 
        legend.title = element_text(size = 20))

ggsave("EIN_4_496_PROBITm1000R1.jpeg")

#####################
##### BOTH NNEs #####
#####################
### READING FROM FILE ###
bbb_both2 <- read.csv("NNEm1000_PROBIT_7_248.csv", header = T)

bbb_both2$NNE <- ifelse(bbb_both2$NNE > 20, NA, bbb_both2$NNE)
bbb_both2$NNE <- ifelse(bbb_both2$NNE < 1, NA, bbb_both2$NNE)

ggplot(data = bbb_both2, 
       aes(x = factor(n), y = NNE, fill = factor(TYPE))) + 
  geom_boxplot() +
  geom_hline(yintercept = 7.248, linetype = "dashed",
             col        = "red", size = 1)     +
  xlab("Sample size")        + 
  ylab("NNE")           +
  ggtitle("IV-based & Unadjusted Estimators of the NNE 
          as a Function of the Sample Size")  +
  scale_y_continuous(breaks = 1:10)          +
  theme_minimal() + 
  scale_fill_discrete(name = "NNE Type") + 
  theme(plot.title  = element_text(hjust = 0.5,     size = 25),
        axis.text.x = element_text(  size = 20),
        axis.text.y = element_text(  size = 20), 
        axis.title  = element_text(size  = 20), 
        legend.key.size = unit(1, 'cm'), 
        legend.text  = element_text(size = 20), 
        legend.title = element_text(size = 20))

ggsave("NNE_7_248_PROBITm1000R1.jpeg")
