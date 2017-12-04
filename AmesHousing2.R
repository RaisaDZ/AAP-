###code for reproducing results of AAP for Ames housing data set 
###with quarterly random forest experts built for the first year of the data

###Ames House data set was compiled by Dean De Cock and could be found at
###https://ww2.amstat.org/publications/jse/v19n3/decock.pdf

#install packages
library(plyr)
library(randomForest)

#read df
dt <-read.csv("AmesHousing.csv",header=TRUE)

#total square footage
dt$Total.SQ.Footage <- dt$Total.Bsmt.SF + dt$Gr.Liv.Area

#remove outliers
dt <- subset(dt, dt$Gr.Liv.Area < 4000)
dt <- dt[!is.na(dt$Total.SQ.Footage),]

#number of quarterly experts
num_experts <- 4

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
T <- max(pack_size$time)
N <- num_experts
eta <- 2/(B-A)^2

#AAP
#learning_rate = 1: Aggregating Algorithm for Packs with an Unknown Maximum (AAP-incremental)
#learning_rate = 2: Aggregating Algorithm for Packs Averages (AAP-current)
#learning_rate = 3: Aggregating Algorithm for Packs with an Known Maximum (AAP-max)

for (learning_rate in 1:3) {
  weights <- matrix(1, nrow = T+1, ncol = N)
  for (i in 1:T) {
    dt.pred <- subset(df, df$time == i)  #subset the current pack
    experts <- dt.pred[, grepl("pred" , names(dt.pred))]  #select monthly experts
    #predictions within specified interval
    experts[experts < A] <- A
    experts[experts > B] <- B
    outcomes <- dt.pred$Sales
    gamma <- matrix(0, nrow = dim(dt.pred)[1], ncol = 1)  #matrix of predictions
    loss <- matrix(0, nrow = dim(dt.pred)[1], ncol = N)  #losses of experts
    loss_AA <- matrix(0, nrow = dim(dt.pred)[1], ncol = 1)  #losses of AAP
    for (j in 1:dim(dt.pred)[1]) {
      weights.norm <- weights[i, ] / sum(weights[i, ])  #normalise weights
      #calculate generalised predictions
      g_A = -1 / eta * log(sum(weights.norm * exp(-eta * (experts[j,]-A)^2)))
      g_B = -1 / eta * log(sum(weights.norm * exp(-eta * (B-experts[j,])^2)))
      #calculate prediction by AAP
      dt.pred[[paste("gamma", learning_rate, sep = "")]][j] <- (A + B) / 2 - (g_B - g_A) / (2 * (B - A)) 
      loss[j, ] <- (experts[j,] - outcomes[j])^2
    } 
    #weights update
    if (learning_rate == 1) {
      weights[i+1, ] <- weights[i, ] * exp(-(eta / subset(pack_size, pack_size$time == i)$K_inc) * apply(loss, 2, sum)) 
    }
    if (learning_rate == 2) {
      weights[i+1, ] <- weights[i, ] * exp(-(eta / subset(pack_size, pack_size$time == i)$K_current) * apply(loss, 2, sum)) 
    }
    if (learning_rate == 3) {
      weights[i+1, ] <- weights[i, ] * exp(-(eta / subset(pack_size, pack_size$time == i)$K_max) * apply(loss, 2, sum)) 
    }
    if (i == 1) {
      dt.pred_tot <- dt.pred
    } else {
      dt.pred_tot <- rbind(dt.pred_tot, dt.pred)
    }
  }
  if (learning_rate == 1) {
    dt.total <- dt.pred_tot[, c("PID", "gamma1")] 
  }
  if (learning_rate == 2) {
    dt.total <- merge(dt.total, dt.pred_tot[, c("PID", "gamma2")], by = c("PID"))
  }
  if (learning_rate == 3) {
    dt.total <- merge(dt.total, dt.pred_tot[, c("PID", "gamma3")], by = c("PID")) 
  }
}

df <- merge(df, dt.total, by = c("PID"))

#losses of AA
for (i in 1:3) {
  df[[paste("loss_AA", i, sep = "")]] <- (df$Sales - df[[paste("gamma", i, sep = "")]])^2
}

#aggregate losses by packs
losses <-aggregate(df[, grepl(c("loss") , names(df))], by=list(df$Yr.Sold, df$Mo.Sold), FUN=sum)
names(losses)[1:2] <- c("Yr.Sold", "Mo.Sold")
losses <- losses[order(losses$Yr.Sold, losses$Mo.Sold), ]

#calculate cumulative Losses
for (j in 1:num_experts) {
  losses[[paste("Loss", j, sep = "")]] <- cumsum(losses[[paste("loss", j, sep = "")]])
}
for (j in 1:3) {
  losses[[paste("Loss_AA", j, sep = "")]] <- cumsum(losses[[paste("loss_AA", j, sep = "")]])
}
losses <- merge(losses, pack_size, by = c("Yr.Sold", "Mo.Sold"))
losses <- losses[order(losses$time), ]

#export to csv
write.csv(losses, "losses_ames2.csv", row.names = FALSE)

