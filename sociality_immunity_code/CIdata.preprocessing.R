# This script covers the steps that input raw check-in data and produce the
# counts of check-ins in a more usable format.
# This requires two files: a mouse metadata files, for the mice released
# in the experiment, and a file containing all the check-ins.
# The experiment happened in two blocks, so we actually follow this process
# twice.
library(dplyr)

# The first step is to take the metadata and clean it up for use with the
# check-in data.
rawMice2021block1 <- read.csv("~/Data/mice.GxE.metadata.2021.block1.csv")
mice <- rawMice2021block1
mice <- filter(mice,Location=="SF") # This data frame includes metadata for mice that
# were lab-housed for a different element of the experiment.
# These mice are only analyzed in Oyesola et al.
# Each mouse had multiple forms of ID; we take the mouse ID from its left ear
# tag.
for (i in 1:nrow(mice)){
  if (is.na(mice$Left.ear.tag[i])==TRUE) {
    mice$Left.ear.tag[i] <- mice$Right.ear.tag[i]
  }
}
mice <- rename(mice,Mouse=Left.ear.tag,RFID.tag=RFID.Tag)%>%
  select(RFID.tag,Mouse,Wedge,Genotype,Sex,Cage,Parentage,
         Infected,Lost.RFID.1,Lost.RFID.2,Missing.1,Missing.2)
for(i in 1:nrow(mice)){
  mice$Wedge[i] <- paste("Wedge",mice$Wedge[i])
} # Note – Wedge is the same thing as enclosure in the text.
mice$Mouse <- as.character(mice$Mouse)
mice$time.control.1 <- ifelse(mice$Lost.RFID.1==TRUE|mice$Missing.1==TRUE,
                              mice$time.control.1 <- TRUE,mice$time.control.1 <- FALSE)
mice$time.control.2 <- ifelse(mice$Lost.RFID.2==TRUE|mice$Missing.2==TRUE,
                              mice$time.control.2 <- TRUE,mice$time.control.2 <- FALSE)
mice2021block1 <- mice
# The metadata is ready for use with the raw check-in data.

# Next we import the raw check-in data.
# It has four columns: Wedge (here meaning RFID reader; don't ask), timestamp,
# RFID, and mouse (which is translated by our system from the RFID).
rawCIdata2021block1 <- 
  read.csv("~/Data/2021.GxE.block1.clean.csv")
CIdata <- rawCIdata2021block1
# Some check-in records are either from test RFIDs or are nonsense produced by
# malfunctions, mostly water damage.
CIdata <- filter(CIdata,Mouse!="Unknown"&Mouse!="0"&is.na(Mouse)==FALSE)
CIdata$Mouse <- as.character(CIdata$Mouse)

# We want to convert the timestamps into more usable data
x=strptime(CIdata$Timestamp,format="%m/%d/%y %H:%M:%S")
x <- x-14400 # This adjusts the time by four hours (14400 seconds)
# Our system stores the time as UTC, but summer in NJ is UTC-04.
Date <- format(x,"%Y-%m-%d")
Hour <- format(x,"%H")
Minute <- format(x,"%M")
Second <- format(x,"%S")
Day <- format(x,"%m-%d")
CIdata <- cbind(CIdata,Date)%>%
  cbind(Hour)%>%
  cbind(Minute)%>%
  cbind(Second)%>%
  cbind(Day)
CIdata$Timestamp <- format(x,"%Y-%m-%d %H:%M:%S") # This is necessary for the association
# function to work properly - it needs the timestamp to be formatted in a certain way
rm(list=c("x","Date","Hour","Minute","Second","Day"))

Week1 <- as.factor(c("05-25","05-26","05-27","05-28","05-29","05-30","05-31"))
Week2 <- as.factor(c("06-01","06-02","06-03","06-04","06-05","06-06","06-07"))
Week3 <- as.factor(c("06-08","06-09","06-10","06-11","06-12","06-13","06-14"))
Week4 <- as.factor(c("06-15","06-16","06-17","06-18","06-19","06-20","06-21"))
Week5 <- as.factor(c("06-22","06-23","06-24","06-25","06-26","06-27","06-28"))
Week6 <- as.factor(c("06-29","06-30","07-01","07-02","07-03","07-04","07-05"))
CIdata$Week=ifelse(CIdata$Day %in% Week1, "Week 1",
                   ifelse(CIdata$Day %in% Week2, "Week 2",
                          ifelse(CIdata$Day %in% Week3, "Week 3",
                                 ifelse(CIdata$Day %in% Week4, "Week 4",
                                        ifelse(CIdata$Day %in% Week5, "Week 5",
                                               ifelse(CIdata$Day %in% Week6, "Week 6",
                                                      NA))))))
CIdata=rename(CIdata,Reader=Wedge) # Fixed it.

# Next we need to assign the locations and enclosures for each check-in.
# That requires listing which readers are at which locations and which
# readers are in which enclosures.
Wedge2 <- c(1,2,4,7,30)
Wedge3 <- c(3,5,6,13,71)
Wedge4 <- c(10,11,12,14,70)
Feeder <- c(30,70,71)
Tower <- c(1,13,14)
Left <- c(4,6,12)
Right <- c(2,5,10)
Base <- c(3,7,11)
CIdata$Wedge <- ifelse(CIdata$Reader %in% Wedge2,"Wedge 2",
                       ifelse(CIdata$Reader %in% Wedge3,"Wedge 3",
                              ifelse(CIdata$Reader %in% Wedge4,"Wedge 4",NA)))
CIdata$Location <- ifelse(CIdata$Reader %in% Feeder,"Feeder",
                           ifelse(CIdata$Reader %in% Tower, "Tower",
                                  ifelse(CIdata$Reader %in% Left, "Left side",
                                         ifelse(CIdata$Reader %in% Right,
                                                "Right side",
                                                ifelse(CIdata$Reader %in% Base,
                                                       "Base",NA)))))
rm(list=c("Wedge2","Wedge3","Wedge4",
          "Feeder","Tower","Left","Right","Base"))
# We also want to identify if check-ins are in the assigned wedge (for easier
# filtering later) and to attach the strain.
CIdata$AssignedWedge <- ""
CIdata$Genotype <- "" # Note that genotype = strain.
for(i in 1:nrow(CIdata)){
  ID <- CIdata$Mouse[i]
  CIdata$AssignedWedge[i] <- as.character(mice$Wedge[which(mice$Mouse==ID)])
  CIdata$Genotype[i] <- mice$Genotype[which(mice$Mouse==ID)]
} # This takes a little while to execute.

CIdata <- select(CIdata,Timestamp,RFID,Wedge,Reader,Location,
                 Mouse,Genotype,AssignedWedge,
                 Date,Week,Day,Hour,Minute,Second)
CIdata2021block1 <- CIdata

# Now we can repeat the process.  We will largely forego the annotations here.
rawMice2021block2 <- read.csv("~/Data/mice.GxE.metadata.2021.block2.csv")
mice <- rawMice2021block2
mice <- filter(mice,Location=="SF")
for (i in 1:nrow(mice)){
  if (is.na(mice$Left.ear.tag[i])==TRUE) {
    mice$Left.ear.tag[i] <- mice$Right.ear.tag[i]
  }
}
mice <- rename(mice,Mouse=Left.ear.tag,RFID.tag=RFID.Tag)%>%
  select(RFID.tag,Mouse,Wedge,Genotype,Sex,Cage,Parentage,
         Infected,Lost.RFID.1,Lost.RFID.2,Missing.1,Missing.2)
for(i in 1:nrow(mice)){
  mice$Wedge[i] <- paste("Wedge",mice$Wedge[i])
}
mice$Mouse <- as.character(mice$Mouse)
mice$time.control.1 <- ifelse(mice$Lost.RFID.1==TRUE|mice$Missing.1==TRUE,
                              mice$time.control.1 <- TRUE,mice$time.control.1 <- FALSE)
mice$time.control.2 <- ifelse(mice$Lost.RFID.2==TRUE|mice$Missing.2==TRUE,
                              mice$time.control.2 <- TRUE,mice$time.control.2 <- FALSE)
mice2021block2 <- mice

rawCIdata2021block2 <- 
  read.csv("~/Documents/Princeton EEB/Stony Ford 2021/SF 2021 Activity/2021.GxE.block2.clean.csv")
CIdata <- select(rawCIdata2021block2,Wedge:Mouse)
CIdata$Mouse[which(CIdata$Mouse=="304b")] <- "304"
CIdata <- filter(CIdata,Mouse!="107")
# Oddly 107 is recorded, but it shouldn't be there.  I'll have to look into this.

x=strptime(CIdata$Timestamp,format="%Y-%m-%d %H:%M:%S")
# x <- x-14400 # This adjusts the time by four hours (14400 seconds)
# It is, however, only necessary for some years (vagaries in how the messages were stored,
# I suppose?)
Date <- format(x,"%Y-%m-%d")
Hour <- format(x,"%H")
Minute <- format(x,"%M")
Second <- format(x,"%S")
Day <- format(x,"%m-%d")
CIdata <- cbind(CIdata,Date)%>%
  cbind(Hour)%>%
  cbind(Minute)%>%
  cbind(Second)%>%
  cbind(Day)
CIdata$Timestamp <- format(x,"%Y-%m-%d %H:%M:%S") # This is necessary for the association
# function to work properly - it needs the timestamp to be formatted in a certain way
rm(list=c("x","Date","Hour","Minute","Second","Day"))

# Note that the days for each week are different!
Week1 <- as.factor(c("07-20","07-21","07-22","07-23","07-24","07-25","07-26"))
Week2 <- as.factor(c("07-27","07-28","07-29","07-30","07-31","08-01","08-02"))
Week3 <- as.factor(c("08-03","08-04","08-05","08-06","08-07","08-08","08-09"))
Week4 <- as.factor(c("08-10","08-11","08-12","08-13","08-14","08-15","08-16"))
Week5 <- as.factor(c("08-17","08-18","08-19","08-20","08-21","08-22","08-23"))
CIdata$Week=ifelse(CIdata$Day %in% Week1, "Week 1",
                   ifelse(CIdata$Day %in% Week2, "Week 2",
                          ifelse(CIdata$Day %in% Week3, "Week 3",
                                 ifelse(CIdata$Day %in% Week4, "Week 4",
                                        ifelse(CIdata$Day %in% Week5, "Week 5",
                                                      NA)))))
CIdata=rename(CIdata,Reader=Wedge)
Wedge2 <- c(1,2,4,7,30)
Wedge3 <- c(3,5,6,13,71)
Wedge4 <- c(10,11,12,14,70)
Feeder <- c(30,70,71)
Tower <- c(1,13,14)
Left <- c(4,6,12)
Right <- c(2,5,10)
Base <- c(3,7,11)
CIdata$Wedge <- ifelse(CIdata$Reader %in% Wedge2,"Wedge 2",
                       ifelse(CIdata$Reader %in% Wedge3,"Wedge 3",
                              ifelse(CIdata$Reader %in% Wedge4,"Wedge 4",NA)))
CIdata$Location <- ifelse(CIdata$Reader %in% Feeder,"Feeder",
                          ifelse(CIdata$Reader %in% Tower, "Tower",
                                 ifelse(CIdata$Reader %in% Left, "Left side",
                                        ifelse(CIdata$Reader %in% Right,
                                               "Right side",
                                               ifelse(CIdata$Reader %in% Base,
                                                      "Base",NA)))))
rm(list=c("Wedge2","Wedge3","Wedge4",
          "Feeder","Tower","Left","Right","Base"))
CIdata$AssignedWedge <- ""
CIdata$Genotype <- ""
for(i in 1:nrow(CIdata)){
  ID <- CIdata$Mouse[i]
  CIdata$AssignedWedge[i] <- as.character(mice$Wedge[which(mice$Mouse==ID)])
  CIdata$Genotype[i] <- mice$Genotype[which(mice$Mouse==ID)]
} # This takes a little while to execute.

CIdata <- select(CIdata,Timestamp,RFID,Wedge,Reader,Location,
                 Mouse,Genotype,AssignedWedge,
                 Date,Week,Day,Hour,Minute,Second)
CIdata2021block2 <- CIdata

# We will also want the Julian night variable, because the mice are nocturnal.
##'   This function gives numeric values for each day and night of the dataset.  If
##'   certain dates are missing, they will not be given a number, so be aware that the
##'   number is not the number of days since the first.  Julian_Night is a variable that
##'   defines the night period as a unified period of time, rather than the day.  This
##'   fits better for nocturnal animals.
##'   
##'   @param ciData a data frame of check-in data with a "Date" variable
##'   
##'   @return The original data frame with two new columns for Julian Day and Night

julian.day.night <- function(ciData){
  ciData$Julian_Day <- as.numeric(as.factor(ciData$Date))
  ciData$Julian_Night<-ciData$Julian_Day
  for(i in 1:(nrow(ciData))){
    if(ciData[i,]$Hour %in% c("00","01","02","03","04","05","06",
                              "07","08","09","10","11")){
      ciData[i,]$Julian_Night <- ciData[i,]$Julian_Night-1
    }
  }
  # This takes a little while to do.
  return(ciData)
}
# And we want to filter out check-ins from 1) nights during which trapping took
# place and 2) nights with RFID reader malfunctions of some form (see Methods).
# First we do this for Block 1:
data <- julian.day.night(CIdata2021block1)
data <- filter(data,Wedge==AssignedWedge)
data.1 <- filter(data,Julian_Night<17)%>%
  filter(!((Wedge=="Wedge 3"|Wedge=="Wedge 4")&(Julian_Night==15|Julian_Night==16)))
data.2 <- filter(data,Julian_Night>15)%>%
  filter(Julian_Night<37)%>%
  filter(!((Wedge=="Wedge 2")&(Julian_Night==16|Julian_Night==17)))%>%
  filter(!(Wedge=="Wedge 4"&Julian_Night==36))

W2.CI.2021b1.1 <- filter(data.1,Wedge=="Wedge 2")
W3.CI.2021b1.1 <- filter(data.1,Wedge=="Wedge 3")
W4.CI.2021b1.1 <- filter(data.1,Wedge=="Wedge 4")

W2.CI.2021b1.2 <- filter(data.2,Wedge=="Wedge 2")
W3.CI.2021b1.2 <- filter(data.2,Wedge=="Wedge 3")
W4.CI.2021b1.2 <- filter(data.2,Wedge=="Wedge 4")
W4.CI.2021b1.2 <- filter(W4.CI.2021b1.2,!(Julian_Night%in%c(25:31)))

# And now for the Block 2 data:
data <- julian.day.night(CIdata2021block2)
data <- filter(data,Wedge==AssignedWedge)
data.1 <- filter(data,Julian_Night<9)%>%
  filter(!((Wedge=="Wedge 2"|Wedge=="Wedge 3")&(Julian_Night==8)))
data.2 <- filter(data,Julian_Night>8)%>%
  filter(Julian_Night<31)%>%
  filter(!((Wedge=="Wedge 4")&(Julian_Night==9)))%>% # Trap night for Wedge 4
  filter(!((Wedge=="Wedge 2"|Wedge=="Wedge 3")&(Julian_Night%in%c(29,30))))

# Set the break between 1 and 2 at 7/29, which is Julian_Night 9
# Trapout Julian_Night for Wedges 6 and 7: 27
# Trapout Julian_Night for Wedge 4: 31
# Trapout Julian_Night for Wedges 2 and 3: 29

W2.CI.2021b2.1 <- filter(data.1,Wedge=="Wedge 2")
W3.CI.2021b2.1 <- filter(data.1,Wedge=="Wedge 3")
W4.CI.2021b2.1 <- filter(data.1,Wedge=="Wedge 4")

W2.CI.2021b2.2 <- filter(data.2,Wedge=="Wedge 2")
W2.CI.2021b2.2 <- filter(W2.CI.2021b2.2,Julian_Night!=10)
W3.CI.2021b2.2 <- filter(data.2,Wedge=="Wedge 3")
W3.CI.2021b2.2 <- filter(W3.CI.2021b2.2,Julian_Night!=10)
W4.CI.2021b2.2 <- filter(data.2,Wedge=="Wedge 4")
W4.CI.2021b2.2 <- filter(W4.CI.2021b2.2,!(Julian_Night%in%c(10,11)))

# We've broken the data out into chunks for each enclosure during each part of
# the experiment.  This will be important for the social networks.
# But we can also make it one big file, if that's helpful.
temp <- rbind(W2.CI.2021b1.1,W2.CI.2021b1.2,
              W3.CI.2021b1.1,W3.CI.2021b1.2,
              W4.CI.2021b1.1,W4.CI.2021b1.2,
              W6.CI.2021b1.1,W6.CI.2021b1.2,
              W7.CI.2021b1.1,W7.CI.2021b1.2)
temp$Block <- "Block 1"
temp2 <- rbind(W2.CI.2021b2.1,W2.CI.2021b2.2,
              W3.CI.2021b2.1,W3.CI.2021b2.2,
              W4.CI.2021b2.1,W4.CI.2021b2.2,
              W6.CI.2021b2.1,W6.CI.2021b2.2,
              W7.CI.2021b2.1,W7.CI.2021b2.2)
temp2$Block <- "Block 2"
CIdata2021GxEfiltered <- filter(rbind(temp,temp2),
                                Wedge%in%c("Wedge 2","Wedge 3","Wedge 4"))

# At this point the data can be used for counting the number of check-ins of different
# types and for building social networks.

# The next script is ______


