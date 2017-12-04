###code for reproducing results of AAP for Ames housing data set 
###with monthly linear regression experts built for the first year of the data

#install packages
library(plyr)

#read df
dt <-read.csv("AmesHousing.csv",header=TRUE)

#total square footage
dt$Total.SQ.Footage <- dt$Total.Bsmt.SF + dt$Gr.Liv.Area
dt$Sales <- dt$SalePrice / 10^6

#remove outliers
dt <- subset(dt, dt$Gr.Liv.Area < 4000)
dt <- dt[!is.na(dt$Total.SQ.Footage),]

#number of monthly experts
num_experts <- 12

dt <- dt[order(dt$Yr.Sold, dt$Mo.Sold), ]
dt <-  within(dt, period <- as.numeric(interaction(Yr.Sold, Mo.Sold, drop=TRUE, lex.order=TRUE)))
               
#create shorter dataset
df <- dt[, c("PID", "Order", "SalePrice", "Sales", "Total.SQ.Footage", "Neighborhood", "Yr.Sold", "Mo.Sold", "period")]
for (i in 1:num_experts) {
  df[[paste("flag", i, sep="")]] <- 0
  df[[paste("pred", i, sep="")]] <- 0
  df[[paste("loss", i, sep="")]] <- 0
}

neighborhood_num <- rep(0, num_experts)  #number of neighbourhoods in the model

#formulas for linear model
formula <- "Sales ~ Total.SQ.Footage + Neighborhood"
formula0 <- "Sales ~  Total.SQ.Footage"

#create monthly linear regression experts
#if the current neighbourhood is not presented in the model then use formula0, otherwise formula
for (i in 1:num_experts)  {
  lmfit <- lm(formula, data = subset(df, df$period == i))
  lmfit0 <- lm(formula0, data = subset(df, df$period == i))
  neighborhood_num[i] <- length(lmfit$coefficients) - 2
  Neighbor_list <- unique(subset(df, df$period == i)$Neighborhood)
  for (j in 1:dim(df)[1])  {
    df[[j, paste("flag", i, sep="")]] <- sum(grepl(paste("^",df$Neighborhood[j],"$", sep=""), Neighbor_list, fixed = FALSE))
    if (df[[j, paste("flag", i, sep="")]] == 1) { 
      df[[j, paste("pred", i, sep="")]] <- predict(lmfit, df[j,])
    } else {
      df[[j, paste("pred", i, sep="")]] <- predict(lmfit0, df[j,])
    }
  }
  df[[paste("loss", i, sep="")]] <- (df$Sales - df[[paste("pred", i, sep="")]])^2
}

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

#create seasonal and first year batch linear regressions for comparison
lm_year <- lm(formula, data = subset(dt, dt$Yr.Sold == min(dt$Yr.Sold)))
lm_year.0 <- lm(formula0, data = subset(dt, dt$Yr.Sold == min(dt$Yr.Sold)))
Neighbor_year <- unique(subset(dt, dt$Yr.Sold == min(dt$Yr.Sold))$Neighborhood)
for (j in 1:dim(df)[1]) {
  df$flag_year[j] <- sum(grepl(paste("^",df$Neighborhood[j],"$", sep=""), Neighbor_year, fixed = FALSE))
  if (df$flag_year[j] == 1) {
    df$pred_year[j] <- predict(lm_year, df[j, ])
  } else {
    df$pred_year[j] <- predict(lm_year.0, df[j, ])
  }
  df$pred_month[j] <- df[[paste("pred", df$Mo.Sold[j], sep = "")]][j]
} 

df$loss_year <- (df$Sales - df$pred_year)^2
df$loss_month <- (df$Sales - df$pred_month)^2

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
losses$Loss_year <- cumsum(losses$loss_year)
losses$Loss_month <- cumsum(losses$loss_month)

losses <- merge(losses, pack_size, by = c("Yr.Sold", "Mo.Sold"))
losses <- losses[order(losses$time), ]
#export to csv
write.csv(losses, "losses_ames1.csv", row.names = FALSE)

