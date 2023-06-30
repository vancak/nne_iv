###########################
##### BOXPLOTS LOGIT ######
###########################

#####################
##### BOTH NNTs #####
#####################
library(reshape2)
library(ggplot2)

m <- 1000

NNT_mat2 <- as.data.frame(NNT_mat)

colnames(NNT_mat2) <- c("NNT_UN", "NNT_IV")

NNT_mat2$id <- 1:m

NNT_psi$id <- 1:m

bbb2 <- melt(NNT_psi,
             id.var        = "id", 
             variable.name = "n", 
             value.name    = "NNT") 

NNT_undj$id <- 1:m
bbb2$TYPE <- "IV"

bbb_un <- melt(NNT_undj,
               id.var        = "id", 
               variable.name = "n", 
               value.name    = "NNT") 

# apply(NNT_psi, 2, mean)

bbb_un$NNT <- ifelse(bbb_un$NNT > 20, NA, bbb_un$NNT)
bbb_un$TYPE <- "Unadjusted"

bbb_both <- merge(x = bbb2, y = bbb_un, by = c("id", "n"))

bbb_both2 <- rbind(bbb2, bbb_un)

bbb_both2$NNT <- ifelse(bbb_both2$NNT > 15, NA, bbb_both2$NNT)
bbb_both2$NNT <- ifelse(bbb_both2$NNT < -25, NA, bbb_both2$NNT)

# write.csv(bbb_both2, "NNTm1000_LOGIT_4_65.csv", row.names = F)

### READING FROM FILE ###
bbb_both2 <- read.csv("NNTm1000_LOGIT_4_65.csv", header = T)

ggplot(data = bbb_both2, 
       aes(x = factor(n), y = NNT, fill = factor(TYPE))) + 
  geom_boxplot() +
  geom_hline(yintercept = 4.65, linetype = "dashed",
             col        = "red", size = 1)     +
  xlab("Sample size")        + 
  ylab("NNT")           +
  ggtitle("IV-based & Unadjusted Estimators of the NNT 
          as a Function of the Sample Size")  +
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

ggsave("NNT_4_65_LOGITm1000V1.jpeg")

###############################
### EIN BOXPLOTS - ONLY ADJ ###
###############################
EIN_mat2 <- as.data.frame(EIN_mat)

colnames(EIN_mat2) <- c("EIN_UN", "EIN_IV")

EIN_mat2$id <- 1:m

EIN_psi$id  <- 1:m

bbb2 <- melt(EIN_psi,
             id.var        = "id",
             variable.name = "n",
             value.name    = "EIN")

EIN_undj$id <- 1:m
bbb2$TYPE <- "IV"

bbb_un <- melt(EIN_undj,
               id.var        = "id",
               variable.name = "n",
               value.name    = "EIN")

# apply(NNT_psi, 2, mean)

bbb_un$EIN <- ifelse(bbb_un$EIN > 20, NA, bbb_un$EIN)
bbb_un$TYPE <- "Unadjusted"

bbb_both <- merge(x = bbb2, y = bbb_un, by = c("id", "n"))

bbb_both2 <- rbind(bbb2, bbb_un)

bbb_both2$EIN <- ifelse(bbb_both2$EIN > 10, NA, bbb_both2$EIN)
bbb_both2$EIN <- ifelse(bbb_both2$EIN < -6, NA, bbb_both2$EIN)

# write.csv(bbb_both2, "EINm1000_LOGIT_4_18.csv", row.names = F)

EIN_mat2 <- as.data.frame(EIN_psi)

EIN_mat2$id <- 1:m

library(reshape2)
library(ggplot2)

bbb2 <- melt(EIN_mat2,
             id.var        = "id", 
             variable.name = "n", 
             value.name    = "EIN") 

# apply(NNT_psi, 2, mean)

### READING FROM FILE ###
bbb2 <- read.csv("EINm1000_LOGIT_4_18.csv")

bbb2$EIN <- ifelse(bbb2$EIN > 10, NA, bbb2$EIN)
bbb2 <- bbb2[bbb2$TYPE == "IV",]

ggplot(data = bbb2, 
       aes(x = factor(n), y = EIN, group = n)) +
  geom_boxplot(fill = c("#F8766D")) +
  geom_hline(yintercept = 4.18, linetype = "dashed", col = "red", size = 1)     +
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

ggsave("EIN_4_18_LOGITm1000V1.jpeg")

#####################
##### BOTH NNEs #####
#####################

NNE_mat2 <- as.data.frame(NNE_mat)

colnames(NNE_mat2) <- c("NNE_UN", "NNE_IV")

NNE_mat2$id <- 1:m

NNE_psi$id  <- 1:m

bbb2 <- melt(NNE_psi,
             id.var        = "id", 
             variable.name = "n", 
             value.name    = "NNE") 

NNE_undj$id <- 1:m
bbb2$TYPE <- "IV"

bbb_un <- melt(NNE_undj,
               id.var        = "id", 
               variable.name = "n", 
               value.name    = "NNE") 

# apply(NNT_psi, 2, mean)

bbb_un$NNE <- ifelse(bbb_un$NNE > 20, NA, bbb_un$NNE)
bbb_un$TYPE <- "Unadjusted"

bbb_both <- merge(x = bbb2, y = bbb_un, by = c("id", "n"))

bbb_both2 <- rbind(bbb2, bbb_un)

#write.csv(bbb_both2, "NNE_NNTbet4.csv", row.names = F)


bbb_both2$NNE <- ifelse(bbb_both2$NNE > 10, NA, bbb_both2$NNE)
bbb_both2$NNE <- ifelse(bbb_both2$NNE < -7, NA, bbb_both2$NNE)

# write.csv(bbb_both2, "NNEm1000_LOGIT_5_60.csv", row.names = F)

### READING FROM FILE ###
bbb_both2 <- read.csv("NNEm1000_LOGIT_5_60.csv", header = T)

ggplot(data = bbb_both2, 
       aes(x = factor(n), y = NNE, fill = factor(TYPE))) + 
  geom_boxplot() +
  geom_hline(yintercept = 5.60, linetype = "dashed",
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

ggsave("NNE_5_60_LOGITm1000V1.jpeg")
