###code for reproducing results of Parallel copies of AA for Ames housing data set 
###with quarterly random forest experts built for the first year of the data

###Ames House data set was compiled by Dean De Cock and could be found at
###https://ww2.amstat.org/publications/jse/v19n3/decock.pdf

#install packages
library(plyr)
library(randomForest)
library(data.table)

#read df
dt <-read.csv("AmesHousing.csv",header=TRUE)

#total square footage
dt$Total.SQ.Footage <- dt$Total.Bsmt.SF + dt$Gr.Liv.Area

#remove outliers
dt <- subset(dt, dt$Gr.Liv.Area < 4000)
dt <- dt[!is.na(dt$Total.SQ.Footage),]

#number of quarterly experts
num_experts <- 4
#number of random shuffling
num_shuffles <- 500

#features for RF models
features <- colnames(dt)[colSums(is.na(dt)) == 0]
features <- features[!features %in% c("Order", "PID", "Mo.Sold", "Yr.Sold", "SalePrice")]

#divide data by packs
Quarter <- unique(dt[,c("Yr.Sold", "Mo.Sold")])
Quarter <- Quarter[order(Quarter$Yr.Sold, Quarter$Mo.Sold),]
Quarter$period <- ceiling(Quarter$Mo.Sold / 3) + 4*(Quarter$Yr.Sold - min(Quarter$Yr.Sold))
dt <- merge(dt, Quarter, by = c("Yr.Sold", "Mo.Sold"))

#create RF quarterly experts built on the first year of the data
for (i in 1: num_experts) {
  train <- subset(dt, dt$period == i) 
  set.seed(i)
  rf.fit <- randomForest(SalePrice~., data=train[, c("SalePrice", features)], importance=TRUE, ntree=100)
  dt[[paste("pred", i, sep="")]] <- predict(rf.fit, dt) / 10^6
  dt[[paste("loss", i, sep="")]] <- (dt$SalePrice / 10^6 - dt[[paste("pred", i, sep="")]])^2
}

#create shorter dataset
df <- dt[, c("PID", "Order", "SalePrice", "Total.SQ.Footage", "Yr.Sold", "Mo.Sold", "period",
             "pred1", "pred2", "pred3", "pred4", "loss1", "loss2", "loss3", "loss4")]

df$Sales <- df$SalePrice / 10^6

#calculate size of each pack
pack_size <- ddply(subset(df, df$Yr.Sold > min(dt$Yr.Sold)), .(Yr.Sold, Mo.Sold), summarise, K_current=length(PID))
for (j in 1:dim(pack_size)[1]) {
  pack_size[j, "K_inc"] <- max(pack_size[1:j, "K_current"])
}
pack_size$K_max <- max(pack_size$K_current)
pack_size$time <- seq(1, dim(pack_size)[1]) 

df <- merge(df, pack_size, by = c("Yr.Sold", "Mo.Sold"))

#set parameters for AAP
A <- 0
B <- max(dt$SalePrice) / 10^6 
N <- num_experts
eta <- 2/(B-A)^2

losses.shuffle <- matrix(0, ncol = num_shuffles+1, nrow = 1)
for (s in 1:(num_shuffles+1)) {
  df.shuffle <- df
  set.seed(s)
  #first cycle with the Order from the data set
  if (s == 1) {
    df.shuffle$Order_rand <- df.shuffle$Order
  } else {
    #random shuffling
    df.shuffle$Order_rand <- sample(seq(1, dim(df.shuffle)[1]), replace = FALSE)
  }
  df.shuffle <- df.shuffle[order(df.shuffle$Yr.Sold, df.shuffle$Mo.Sold, df.shuffle$Order_rand),]
  df.shuffle <- data.table(df.shuffle)
  df.shuffle[, C := sequence(.N), by = c("Yr.Sold", "Mo.Sold")]
  df.shuffle <- data.frame(df.shuffle)
  for (i in 1:max(pack_size$K_max)) {
    dt.pred <- subset(df.shuffle, (df.shuffle$C == i))
    T <- dim(dt.pred)[1]
    weights <- matrix(1, nrow = T+1, ncol = N)
    experts <- dt.pred[, grepl("pred" , names(dt.pred))] 
    outcomes <- dt.pred$Sales
    loss <- matrix(0, nrow = dim(dt.pred)[1], ncol = N)
    loss_AA <- matrix(0, nrow = dim(dt.pred)[1], ncol = 1)
    for (j in 1:T) {
      weights.norm <- weights[j, ] / sum(weights[j, ])
      g_A = -1 / eta * log(sum(weights.norm * exp(-eta * (experts[j,]-A)^2)))
      g_B = -1 / eta * log(sum(weights.norm * exp(-eta * (B-experts[j,])^2)))
      dt.pred$gamma[j] <- (A + B) / 2 - (g_B - g_A) / (2 * (B - A)) 
      loss[j, ] <- (experts[j,]- outcomes[j])^2
      weights[j+1, ] <- weights[j, ] * exp(-eta * loss[j, ])
    }
    if (i == 1) {
      dt.pred_tot <- dt.pred
    } else {
      dt.pred_tot <- rbind(dt.pred_tot, dt.pred)
    }
  }
  losses.shuffle[s] <- sum((dt.pred_tot$Sales - dt.pred_tot$gamma)^2)
}

#export to csv
write.csv(losses.shuffle, "losses_shuffle2.csv", row.names = FALSE)

