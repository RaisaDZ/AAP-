###code for reproducing results of AAP for London housing data set 
###with quarterly random forest experts built for the first year of the data

#install packages
library(plyr)
library(FNN)

#read df
dt <-read.csv("Houseprice_2009_100km_London.csv",header=TRUE)
lsoa <- read.csv("NPSL_London100km.csv",header=TRUE)
metro <- read.csv("metro_and_railway_stations.csv",header=TRUE)
indexes <- read.csv("1871702.csv",header=TRUE)
names(indexes)[1] <- "Lsoa11"

#transform features
dt$Price <- gsub("_", "", dt$Price)
dt$Price <- as.numeric((dt$Price))
dt$Oseast1M <- gsub("_", "", dt$Oseast1M)
dt$Oseast1M <- as.numeric((dt$Oseast1M))
dt$Osnrth1M <- gsub("_", "", dt$Osnrth1M)
dt$Osnrth1M <- as.numeric((dt$Osnrth1M))
dt$name <- substr(dt$Postcode, 1, 2)
dt$name <- gsub('[0-9]+', '', dt$name)
dt$name2 <- gsub( " .*$", "", dt$Postcode)
dt$Price2 <- ifelse(dt$Price > 10^6, 1, dt$Price / 10^6) 
lsoa_unique <- unique(lsoa[,c("Oa11", "Lsoa11")])

#add deprivation data
df <- merge(dt, lsoa_unique[,c("Oa11", "Lsoa11")], by = "Oa11")
df <- merge(df, indexes, by = "Lsoa11")

#add metro data
nn = get.knnx(metro[, c("Easting", "Northing")], df[, c("Oseast1M", "Osnrth1M")], k=2)
#indexes of 2 nearest neibours
df$nn1 <- nn$nn.index[,1]
df$nn2 <- nn$nn.index[,2]
metro$nn1 <- seq(1, dim(metro)[1])
metro$nn2 <- seq(1, dim(metro)[1])
#add 1st nn
df <- merge(df, metro[, c("nn1", "StopType", "Easting", "Northing")], by = "nn1")
names(df)[names(df)=="StopType"] <- "StopType_1"
names(df)[names(df)=="Easting"] <- "Easting_1"
names(df)[names(df)=="Northing"] <- "Northing_1"
#add 2nd nn
df <- merge(df, metro[, c("nn2", "StopType", "Easting", "Northing")], by = "nn2")
names(df)[names(df)=="StopType"] <- "StopType_2"
names(df)[names(df)=="Easting"] <- "Easting_2"
names(df)[names(df)=="Northing"] <- "Northing_2"

#distance from stations
df$dist1 <- log(sqrt((df$Oseast1M - df$Easting_1)^2 + (df$Osnrth1M - df$Northing_1)^2))
df$dist21 <- log(sqrt((df$Oseast1M - df$Easting_2)^2 + (df$Osnrth1M - df$Northing_2)^2)
                 -sqrt((df$Oseast1M - df$Easting_1)^2 + (df$Osnrth1M - df$Northing_1)^2)+1)
#Stations
df$Station <- ""
df$Station <- ifelse((df$StopType_1 == "MET")&(df$StopType_2 == "MET"), "Tube", df$Station)
df$Station <- ifelse((df$StopType_1 == "MET")&(df$StopType_2 == "RLY"), "Both", df$Station)
df$Station <- ifelse((df$StopType_1 == "RLY")&(df$StopType_2 == "MET"), "Both", df$Station)
df$Station <- ifelse((df$StopType_1 == "RLY")&(df$StopType_2 == "RLY"), "RLY", df$Station)

#Create dummy variables for Property type
df$Property1 <- ifelse(df$Property_Type == "F", 1, 0)
df$Property2 <- ifelse(df$Property_Type == "S", 1, 0)
df$Property3 <- ifelse(df$Property_Type == "T", 1, 0)
#for station
df$Station1 <- ifelse(df$Station == "Tube", 1, 0)
df$Station2 <- ifelse(df$Station == "RLY", 1, 0)
#Leasehold
df$flag_leasehold <- ifelse(df$Freeorlease == "L", 1, 0)
#New building
df$flag_newbuilding <- ifelse(df$Newbuild == "N", 1, 0)

names(df)[names(df)=="INCOME.SCORE"] <- "Income_score"
names(df)[names(df)=="EMPLOYMENT.SCORE"] <- "Employment_score"
names(df)[names(df)=="HEALTH.DEPRIVATION.AND.DISABILITY.SCORE"] <- "HD_score"
names(df)[names(df)=="CRIME.AND.DISORDER.SCORE"] <- "Crime_score"
names(df)[names(df)=="LIVING.ENVIRONMENT.SCORE"] <- "LIV_score"
names(df)[names(df)=="Indoors.Sub.domain.Score"] <- "Indoors_score"
names(df)[names(df)=="Outdoors.Sub.domain.Score"] <- "Outdoors_score"
names(df)[names(df)=="Geographical.Barriers.Sub.domain.Score"] <- "GB_score"
names(df)[names(df)=="Wider.Barriers.Sub.domain.Score"] <- "WB_score"
names(df)[names(df)=="Children.Young.People.Sub.domain.Score"] <- "Child_score"
names(df)[names(df)=="Skills.Sub.domain.Score"] <- "Skills_score"
names(df)[names(df)=="IDACI.score"] <- "IDACI_score"
names(df)[names(df)=="IDAOPI.score"] <- "IDAOPI_score"

#create monthly linear regression models on the first year of the data
num_experts <- 12
df$ID <- seq(1, dim(df)[1])
df <- df[order(df$Month), ]
df <-  within(df, period <- as.numeric(interaction(Month, drop=TRUE, lex.order=TRUE)))

#features for RF
features <- c("Property1", "Property2", "Property3", "flag_leasehold", "flag_newbuilding",
              "dist1", "dist21", "Station1", "Station2", "Income_score",
              "Employment_score", "HD_score", "Crime_score", "LIV_score", "Indoors_score",
              "Outdoors_score", "GB_score", "WB_score", "Child_score", "Skills_score", 
              "IDACI_score", "IDAOPI_score")

for (i in 1:num_experts) {
  set.seed(i)
  rf.fit <- randomForest(Price2~., data=subset(df, df$period == i)[, c("Price2", features)], importance=TRUE, ntree=100)
  df[[paste("pred", i, sep="")]] <- predict(rf.fit, df)
  df[[paste("loss", i, sep="")]] <- (df["Price2"] - df[[paste("pred", i, sep="")]])^2 
}

#calculate size of each pack
pack_size <- ddply(subset(df, df$Year > min(dt$Year)), ~Month, summarise, K_current=length(ID))
for (j in 1:dim(pack_size)[1]) {
  pack_size[j, "K_inc"] <- max(pack_size[1:j, "K_current"])
}
pack_size$K_max <- max(pack_size$K_current)
pack_size$time <- seq(1, dim(pack_size)[1]) 

df <- merge(df, pack_size, by = c("Month"))

#set parameters for AAP
A <- 0
B <- 1 
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
    outcomes <- dt.pred$Price2
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
    dt.total <- dt.pred_tot[, c("ID", "gamma1")] 
  }
  if (learning_rate == 2) {
    dt.total <- merge(dt.total, dt.pred_tot[, c("ID", "gamma2")], by = c("ID"))
  }
  if (learning_rate == 3) {
    dt.total <- merge(dt.total, dt.pred_tot[, c("ID", "gamma3")], by = c("ID")) 
  }
}

df <- merge(df, dt.total, by = c("ID"))

#losses of AA
for (i in 1:3) {
  df[[paste("loss_AA", i, sep = "")]] <- (df$Price2 - df[[paste("gamma", i, sep = "")]])^2
}

#aggregate losses by packs
losses <-aggregate(df[, grepl(c("loss") , names(df))], by=list(df$Month), FUN=sum)
names(losses)[1] <- c("Month")
losses <- losses[order(losses$Month), ]

#calculate cumulative Losses
for (j in 1:num_experts) {
  losses[[paste("Loss", j, sep = "")]] <- cumsum(losses[[paste("loss", j, sep = "")]])
}
for (j in 1:3) {
  losses[[paste("Loss_AA", j, sep = "")]] <- cumsum(losses[[paste("loss_AA", j, sep = "")]])
}

losses <- merge(losses, pack_size, by = c("Month"))
losses <- losses[order(losses$time), ]
#export to csv
write.csv(losses, "losses_london2.csv", row.names = FALSE)
 