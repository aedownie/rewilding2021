# This script takes the outputs of CIdata.preprocessing and generates all the
# measures of individual behavior.
# We don't import new files, just use the outputs from it (there are fully
# finished files included in the database though).

# Here I will count the number of check-ins by location, creating one column for each
# location, followed by aggregating them all to get the total number of check-ins.

##' A function to count the number of check-ins by location, creating one column for each
##' location, followed by aggregating them all to get the total number of check-ins.
##' Should accommodate different sets of locations.  Returns a data frame matching the
##' original "mice" data frame, but with new columns as above.
##' 
##' @param mice A data frame of the mice being examined
##' @param CIdata A data frame of the check-ins for those mice
##' @param assignedOnly A T/F parameter for whether only check-ins in assigned wedges count
##' 
##' @return The mice data frame, with new check-in columns for each location and the total.

location.counts <- function(mice,CIdata,assignedOnly){
  if(assignedOnly!=TRUE&assignedOnly!=FALSE){
    print("assignedOnly must be designated TRUE or FALSE.")
    return()
  }
  ifelse(assignedOnly==TRUE,
         data <- filter(CIdata,Wedge==AssignedWedge),
         data <- CIdata)
  CItable <- count(data,Mouse,Location)%>%
    pivot_wider(names_from=Location,values_from=n)
  mice <- left_join(mice,CItable,by="Mouse") # This keeps mice with no CIs
  rm(CItable)
  mice <- rename(mice,Left=`Left side`)%>%
    rename(Right=`Right side`)
  # If you want to adapt this code for some other set of locations, you may
  # need to expand things here.

  for (i in 1:nrow(mice)){
    wedge <- mice$Wedge[i]
    k <- ncol(mice)
    for (j in (k-4):k){
      ifelse(is.na(mice[i,j])==TRUE,mice[i,j] <- 0,NA)
    }
    next
  }
  rm(list=c("j","wedge")
  
  mice$Total <- 0
  k <- ncol(mice)
  for(i in 1:nrow(mice)){
    mice$Total[i] <- sum(mice[i,(k-5):(k-1)],na.rm=TRUE)
  }
  rm(list=c("i","k")
  return(mice)
}

# Now we can actually count all the check-ins.  For this we will need the
# metadata from the corresponding block.
temp <- rbind(W2.CI.2021b1.1,W2.CI.2021b1.2,
              W3.CI.2021b1.1,W3.CI.2021b1.2,
              W4.CI.2021b1.1,W4.CI.2021b1.2)
temp$Mouse <- as.character(temp$Mouse)
temp2 <- location.counts(select(mice2021block1,RFID.tag:time.control.2),
                         temp,
                         assignedOnly=TRUE)
temp2$Nights <- NA # Get the nights so that we can determine the number
# of check-ins per night for each individual.
for (i in 1:nrow(temp2)){
  if(temp2$Wedge[i]=="Wedge 2"){
    temp2$Nights[i] <- length(unique(rbind(W2.CI.2021b1.1,
                                           W2.CI.2021b1.2)$Julian_Night))
  }
  if(temp2$Wedge[i]=="Wedge 3"){
    temp2$Nights[i] <- length(unique(rbind(W3.CI.2021b1.1,
                                           W3.CI.2021b1.2)$Julian_Night))
  }
  if(temp2$Wedge[i]=="Wedge 4"){
    temp2$Nights[i] <- length(unique(rbind(W4.CI.2021b1.1,
                                           W4.CI.2021b1.2)$Julian_Night))
  }
}
mice2021block1 <- temp2

temp <- rbind(W2.CI.2021b2.1,W2.CI.2021b2.2,
              W3.CI.2021b2.1,W3.CI.2021b2.2,
              W4.CI.2021b2.1,W4.CI.2021b2.2)
temp$Mouse <- as.character(temp$Mouse)
temp2 <- location.counts(select(mice2021block2,RFID.tag:time.control.2),
                         temp,
                         assignedOnly=TRUE)
temp2$Nights <- NA
for (i in 1:nrow(temp2)){
  if(temp2$Wedge[i]=="Wedge 2"){
    temp2$Nights[i] <- length(unique(rbind(W2.CI.2021b2.1,
                                           W2.CI.2021b2.2)$Julian_Night))
  }
  if(temp2$Wedge[i]=="Wedge 3"){
    temp2$Nights[i] <- length(unique(rbind(W3.CI.2021b2.1,
                                           W3.CI.2021b2.2)$Julian_Night))
  }
  if(temp2$Wedge[i]=="Wedge 4"){
    temp2$Nights[i] <- length(unique(rbind(W4.CI.2021b2.1,
                                           W4.CI.2021b2.2)$Julian_Night))
  }
}
mice2021block2 <- temp2
mice2021total <- rbind(mice2021block1,mice2021block2)

# We need to calculate the number of check-ins per night
mice2021total$CIpN <- mice2021total$Total/mice2021total$Nights
# And we need to count the number of check-ins, nights, and check-ins per night
# for the second half of the experiment ("post-trap").
temp <- mice2021total
temp$PostTrapTotal <- NA
temp$PostTrapFeeder <- NA # For interest
temp$PostTrapNights <- NA
temp$PTCIpN <- NA
for (i in 1:nrow(temp)){
  if(temp$Block[i]=="Block 1"){
    if(temp$Wedge[i]=="Wedge 2"){
      temp$PostTrapTotal[i] <- nrow(filter(W2.CI.2021b1.2,Mouse==temp$Mouse[i]))
      temp$PostTrapFeeder[i] <- nrow(filter(W2.CI.2021b1.2,Mouse==temp$Mouse[i]&
                                              Location=="Feeder"))
      temp$PostTrapNights[i] <- length(unique(W2.CI.2021b1.2$Julian_Night))
    }
    if(temp$Wedge[i]=="Wedge 3"){
      temp$PostTrapTotal[i] <- nrow(filter(W3.CI.2021b1.2,Mouse==temp$Mouse[i]))
      temp$PostTrapFeeder[i] <- nrow(filter(W3.CI.2021b1.2,Mouse==temp$Mouse[i]&
                                              Location=="Feeder"))
      temp$PostTrapNights[i] <- length(unique(W3.CI.2021b1.2$Julian_Night))
    }
    if(temp$Wedge[i]=="Wedge 4"){
      temp$PostTrapTotal[i] <- nrow(filter(W4.CI.2021b1.2,Mouse==temp$Mouse[i]))
      temp$PostTrapFeeder[i] <- nrow(filter(W4.CI.2021b1.2,Mouse==temp$Mouse[i]&
                                              Location=="Feeder"))
      temp$PostTrapNights[i] <- length(unique(W4.CI.2021b1.2$Julian_Night))
    }
  }
  if(temp$Block[i]=="Block 2"){
    if(temp$Wedge[i]=="Wedge 2"){
      temp$PostTrapTotal[i] <- nrow(filter(W2.CI.2021b2.2,Mouse==temp$Mouse[i]))
      temp$PostTrapFeeder[i] <- nrow(filter(W2.CI.2021b2.2,Mouse==temp$Mouse[i]&
                                              Location=="Feeder"))
      temp$PostTrapNights[i] <- length(unique(W2.CI.2021b2.2$Julian_Night))
    }
    if(temp$Wedge[i]=="Wedge 3"){
      temp$PostTrapTotal[i] <- nrow(filter(W3.CI.2021b2.2,Mouse==temp$Mouse[i]))
      temp$PostTrapFeeder[i] <- nrow(filter(W3.CI.2021b2.2,Mouse==temp$Mouse[i]&
                                              Location=="Feeder"))
      temp$PostTrapNights[i] <- length(unique(W3.CI.2021b2.2$Julian_Night))
    }
    if(temp$Wedge[i]=="Wedge 4"){
      temp$PostTrapTotal[i] <- nrow(filter(W4.CI.2021b2.2,Mouse==temp$Mouse[i]))
      temp$PostTrapFeeder[i] <- nrow(filter(W4.CI.2021b2.2,Mouse==temp$Mouse[i]&
                                              Location=="Feeder"))
      temp$PostTrapNights[i] <- length(unique(W4.CI.2021b2.2$Julian_Night))
    }
  }
  temp$PTCIpN[i] <- temp$PostTrapTotal[i]/temp$PostTrapNights[i]
} # This takes a little bit, about 10 seconds
mice2021total <- temp

# In the data, we largely focus on whether the mice were challenged with T.
# muris.  However, you can also consider whether they developed an infection
# with worms present at the end of the experiment.
mice2021total$Wormy <- NA
d <- mice2021total
d$Wormy <- ifelse(d$Worm.count==0|is.na(d$Worm.count)==TRUE,"N","Y")
mice2021total <- d

# Now we want to calculate minimum distance traveled
# We have a function for calculating that quantity.

##' A function to determine the minimum distance traveled by the mice during
##' the experiment.  Note that the distances are very fiddly and vary from
##' enclosure to enclosure; the dimensions don't entirely match.
##'
##' @param ciData A table of check-in data with mouse, wedge, and location info
##' @param mice   A table of mouse metadata for the check-in data
##' @param wedge  The wedge/enclosure that you want to calculate minimum distance for
##'
##' @return a data frame containing a mouse column and a distance column

minimum.distance <- function(ciData,mice,wedge){
  travel <- rep(0,length.out=length(mice))
  if (!(wedge%in%c("Wedge 2","Wedge 3","Wedge 4"))){
    print("wedge argument must be 'Wedge 2', 'Wedge 3', or 'Wedge 4'")
    return(NULL)
  }
  for (mouse in mice){
    dist <- 0
    locations <- filter(ciData,Mouse==mouse)$Location
    for (i in 2:length(locations)){
      current <- locations[i]
      past <- locations[i-1]
      if (current==past){
        next
      }
      if (current!=past){
        if (current%in%c("Base","Feeder")&
             past%in%c("Base","Feeder")){
          if(wedge=="Wedge 2"){dist <- dist+3.0}
          if(wedge=="Wedge 3"|wedge=="Wedge 4"){dist <- dist+3.5}
          next
        }
        if (current%in%c("Tower","Feeder")&
            past%in%c("Tower","Feeder")){
          if(wedge=="Wedge 2"|wedge=="Wedge 4"){dist <- dist+11.5}
          if(wedge=="Wedge 3"){dist <- dist+10.0}
          next
        }
        if (current%in%c("Right side","Feeder")&
            past%in%c("Right side","Feeder")){
          if(wedge=="Wedge 2"|wedge=="Wedge 3"){dist <- dist+4.0}
          if(wedge=="Wedge 4"){dist <- dist+5.0}
          next
        }
        if (current%in%c("Left side","Feeder")&
            past%in%c("Left side","Feeder")){
          if(wedge=="Wedge 2"|wedge=="Wedge 3"){dist <- dist+5.0}
          if(wedge=="Wedge 4"){dist <- dist+4.5}
          next
        }
        if (current%in%c("Tower","Base")&
            past%in%c("Tower","Base")){
          if(wedge=="Wedge 2"|wedge=="Wedge 3"){dist <- dist+14.0}
          if(wedge=="Wedge 4"){dist <- dist+15.0}
          next
        }
        if (current%in%c("Right side","Base")&
            past%in%c("Right side","Base")){
          if(wedge=="Wedge 2"|wedge=="Wedge 4"){dist <- dist+7.5}
          if(wedge=="Wedge 3"){dist <- dist+8.0}
          next
        }
        if (current%in%c("Left side","Base")&
            past%in%c("Left side","Base")){
          if(wedge=="Wedge 2"){dist <- dist+6.5}
          if(wedge=="Wedge 3"|wedge=="Wedge 4"){dist <- dist+7.5}
          next
        }
        if (current%in%c("Right side","Tower")&
            past%in%c("Right side","Tower")){
          if(wedge=="Wedge 2"){dist <- dist+9.0}
          if(wedge=="Wedge 3"){dist <- dist+7.5}
          if(wedge=="Wedge 4"){dist <- dist+8.5}
          next
        }
        if (current%in%c("Left side","Tower")&
            past%in%c("Left side","Tower")){
          if(wedge=="Wedge 2"){dist <- dist+8.5}
          if(wedge=="Wedge 3"){dist <- dist+7.5}
          if(wedge=="Wedge 4"){dist <- dist+8.0}
          next
        }
        if (current%in%c("Left side","Right side")&
            past%in%c("Left side","Right side")){
          if(wedge=="Wedge 2"|wedge=="Wedge 3"){dist <- dist+7.0}
          if(wedge=="Wedge 4"){dist <- dist+6.5}
          next
        }
      }
    }
    travel[which(mice==mouse)] <- dist
  }
  return(data.frame(Mouse=mice,Distance=travel))
}

# We want to combine the two portions of CI data for each enclosure to make
# the process a little bit smoother.
W2.CI.2021b1 <- bind_rows(W2.CI.2021b1.1,W2.CI.2021b1.2)
W2.CI.2021b2 <- bind_rows(W2.CI.2021b2.1,W2.CI.2021b2.2)

W3.CI.2021b1 <- bind_rows(W3.CI.2021b1.1,W3.CI.2021b1.2)
W3.CI.2021b2 <- bind_rows(W3.CI.2021b2.1,W3.CI.2021b2.2)

W4.CI.2021b1 <- bind_rows(W4.CI.2021b1.1,W4.CI.2021b1.2)
W4.CI.2021b2 <- bind_rows(W4.CI.2021b2.1,W4.CI.2021b2.2)

# We run the minimum distance calculations 
W2.b1 <- minimum.distance(W2.CI.2021b1,unique(W2.CI.2021b1$Mouse),"Wedge 2")
W2.b2 <- minimum.distance(W2.CI.2021b2,unique(W2.CI.2021b2$Mouse),"Wedge 2")

W3.b1 <- minimum.distance(W3.CI.2021b1,unique(W3.CI.2021b1$Mouse),"Wedge 3")
W3.b2 <- minimum.distance(W3.CI.2021b2,unique(W3.CI.2021b2$Mouse),"Wedge 3")

W4.b1 <- minimum.distance(W4.CI.2021b1,unique(W4.CI.2021b1$Mouse),"Wedge 4")
W4.b2 <- minimum.distance(W4.CI.2021b2,unique(W4.CI.2021b2$Mouse),"Wedge 4")

minDist <- rbind(W2.b1,W3.b1,W4.b1,
                 W2.b2,W3.b2,W4.b2)
minDist$Mouse <- as.character(minDist$Mouse)
mice2021total <- left_join(mice2021total,minDist,
                  by="Mouse")
# Now we have calculated the minimum distance traveled for each mouse during
# the experiment.

# If you want, you can also calculate minimum distance for just the back
# portion of the experiment.
W2.b1.2 <- minimum.distance(W2.CI.2021b1.2,unique(W2.CI.2021b1.2$Mouse),"Wedge 2")
W3.b1.2 <- minimum.distance(W3.CI.2021b1.2,unique(W3.CI.2021b1.2$Mouse),"Wedge 3")
W4.b1.2 <- minimum.distance(W4.CI.2021b1.2,unique(W4.CI.2021b1.2$Mouse),"Wedge 4")
W2.b2.2 <- minimum.distance(W2.CI.2021b2.2,unique(W2.CI.2021b2.2$Mouse),"Wedge 2")
W3.b2.2 <- minimum.distance(W3.CI.2021b2.2,unique(W3.CI.2021b2.2$Mouse),"Wedge 3")
W4.b2.2 <- minimum.distance(W4.CI.2021b2.2,unique(W4.CI.2021b2.2$Mouse),"Wedge 4")

temp <- rbind(W2.b1.2,W3.b1.2,W4.b1.2,
              W2.b2.2,W3.b2.2,W4.b2.2)
temp$Mouse <- as.character(temp$Mouse)
temp <- rename(temp,Distance2=Distance)
mice2021total <- left_join(mice2021total,temp,
						   by="Mouse")

# The next individual statistic is roaming entropy, a measure of how spread
# through time and space are the mouse's check-ins.  The function for calculating
# this quantity takes a bit more time to run than minimum distance.  For example,
# calculating mean RE for 15 mice across ~22k CIs takes about ten minutes.

##' A function to determine the roaming entropy of mice during their release in
##' the enclosures.  It only calculates mean RE for the nighttime hours, but
##' the hours of interest could be adjusted.
##'
##' @param ciData A table of check-in data with mouse, wedge, and location info
##' @param mice   A table of mouse metadata for the check-in data
##'
##' @return a data frame containing a mouse column and a mean RE column

roaming.entropy <- function(ciData,mice){
  meanRE <- rep(0,length.out=length(mice))
  nightHours <- c("20","21","22","23","00","01","02","03","04","05","06","07")
  # These hours can be adjusted if one wants
  ciData <- filter(ciData,Hour%in%nightHours)
  readers <- c("Feeder","Right side","Left side","Tower","Base")
  hourList <- c(0,1,2,3,4,5,6,7,8,9,10,11)
  for(mouse in mice){
   CIs <- filter(ciData,Mouse==mouse)
   nights <- matrix(0,nrow=length(unique(CIs$Julian_Night)),ncol=5)
   # Any nights with no check-ins are not counted for the purposes of determining mean
   # roaming entropy.  Roaming entropy on such nights is not defined.
   if (nrow(CIs)==0){
     meanRE[which(mice==mouse)] <- NA
     next
   }
   for (i in 1:length(CIs$Hour)){
     CIs$Hour[i] <- hourList[which(nightHours==CIs$Hour[i])]
   }
   CIs$Hour <- as.numeric(CIs$Hour)
   CIs$bin <- ((CIs$Hour*60)+as.numeric(CIs$Minute))+1
   for (j in unique(CIs$Julian_Night)){
     times <- matrix(0,nrow=60*12,ncol=5)
     # We identify for each minute the location of the final check-in taking place during
     # that minute.  If no check-in takes place, then we identify the location of the most
     # recent check-in.
     successes=0
     for (i in 1:nrow(times)){
       set <- filter(CIs,bin==i&Julian_Night==j)
       if(nrow(set)==0){ # If there are no check-ins during that bin
         if (i==1){next} # Allow for 0s prior to check-ins during that night
         times[i,] <- times[i-1,] #Copy previous row if no new information
         next
       }
       if(nrow(set)>=1){
         successes=successes+1
         spot <- set$Location[nrow(set)] # The final check-in during that minute
         times[i,which(readers==spot)] <- 1
       }
       }
     first <- min(which(rowSums(times)==1)) # The first bin in which a check-in occurred
     # during a given night
     index <- which(unique(CIs$Julian_Night)==j)
     if(is.infinite(first)==TRUE){
       nights[index,] <- c(0,0,0,0,0)
       next
       }
     nights[index,] <- colSums(times)/(nrow(times)-(first-1)) # The number of times each
     # location is the sight of last check-in, divided by the number of bins with that
     # information.
     }
   # RE is -sum(log(p)*p)/log(k), where k is the number of readers (here k = 5)
   # When there are no visits to a given reader, the RE value is NaN, because x*log(x) is
   # not defined for x = 0.  However, the limit of x*log(x) approaches 0 as x -> 0.
   # So we will replace NaN values with 0 values to allow evaluation of RE.
   RE <- rep(0,length.out=nrow(nights))
   for (i in 1:length(RE)){
     components <- nights[i,]*log(nights[i,])
     components[which(is.nan(components)==TRUE)] <- 0
     RE[i] <- -1*sum(components)/log(5)
     } # The roaming entropy on each night for a given mouse
   RE_final <- mean(RE)
   meanRE[which(mice==mouse)] <- RE_final
   }
  return(data.frame(mice,meanRE))
}

# Now we can actually calculate mean RE for each enclosure during each block.
W2.b1.RE <- roaming.entropy(W2.CI.2021b1,unique(W2.CI.2021b1$Mouse))
W2.b2.RE <- roaming.entropy(W2.CI.2021b2,unique(W2.CI.2021b2$Mouse))
W3.b1.RE <- roaming.entropy(W3.CI.2021b1,unique(W3.CI.2021b1$Mouse))
W3.b2.RE <- roaming.entropy(W3.CI.2021b2,unique(W3.CI.2021b2$Mouse))
W4.b1.RE <- roaming.entropy(W4.CI.2021b1,unique(W4.CI.2021b1$Mouse))
W4.b2.RE <- roaming.entropy(W4.CI.2021b2,unique(W4.CI.2021b2$Mouse))
d <- rbind(W2.b1.RE,W3.b1.RE,W4.b1.RE,
           W2.b2.RE,W3.b2.RE,W4.b2.RE)
d <- rename(d,Mouse = mice)
d$Mouse <- as.character(d$Mouse)
d <- left_join(mice2021total,d,by="Mouse")
mice2021total <- d

