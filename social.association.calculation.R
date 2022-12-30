# This script does the calculation of social networks from the check-in data
# during the experiment.  This uses the products of CIdata.preprocessing.R;
# if you wish, you can also use the mouse metadata frames from
# individual.behavior.calculation.R, but it's not necessary to have done all
# of those calculations as well.  You can just use the ones after pre-
# processing.
# First we need to define a suite of functions that will be used.

##' A function that builds a GBI containing both individual associations and solitary
##' appearances.  It will also include the metadata of each row for the GBI (a gd), as
##' well as info about sampling periods and places.  The output will be formatted as a
##' list.  It can then easily be fed into the time.network.aura2.modified() function
##' that we updated from Raulo et al. (2021).
##' The code for this function is a slightly modified version of the code from Raulo et
##' al. (2021).  All credit to the original authors for developing the method and writing 
##' most of the code that is used herein.  We have retained their explanatory notes,
##' which start with a hash not followed by a space, and have added a few of our own,
##' which start with a hash that is followed by a space.
##' 
##' @param ciData The check-in data from which a GBI will be built
##' @param overlapTime The length of time for associations, defined in minutes
##' @param IDs A vector containing all the mouse IDs to be cataloged
##' 
##' @return a list with GBI, gd, sampling period-place combinations, and the night for
##' each sampling period-place combination (i.e. just the period)

association.info <- function(ciData,overlapTime,IDs){
  # This function uses a Julian_Night variable for defining activity periods.  This 
  # column is used because mice are nocturnal, meaning that it makes more sense to 
  # determine whether mice are associated during a given night (one activity period),
  # rather than a given day (part of two separate activity periods).
  # If you were to adapt this function for diurnal animals, you'd want to drop that 
  # variable and swap all the uses of it in the code for Julian_Day.
  # You can generate Julian_Night and Julian_Day from the julian.day.night function.

  ciData$Datetime <- ciData$Timestamp
  ciData$Datetime2 <- ciData$Datetime
  # Then we order the data by time, although this probably isn't strictly necessary.
  ciData <- ciData[order(ciData$Datetime),]
  
  # ids<-unique(ciData$Mouse) # This is the list of IDs for the mice in the network
  ids <- IDs
  
  #Now you want to make a blank 'follow' matrix in case you want to fill in who is following who i.e. a network:
  
  AM<-matrix(0,length(ids),length(ids),dimnames=list(ids,ids))
  
  #make a GBI  
  GBI<-matrix(0,1,length(ids),dimnames=list(0,ids))
  colnames(GBI)<-ids
  
  #and also make a GBI-info-dataframe (called gd) for storing info about these group events 
  gd<-data.frame(logger="",id1="",id2="",
                 id1.time="",id2.time="",
                 instance="",lognight="",
                 stringsAsFactors=F) 
  
  #Start a double loop, which is going to work its way through individual-indivdual dyads at a time
  for (i in 1:length(ids)){
    id1<-ids[i] #individual 1 is i
    for(j in 1:length(ids)){
      id2<-ids[j] #individual 2 is j
      
      #now, we only want to do the rest of the loop if individual i and j are different (could miss this out if wanted to fill in the diagonals)
      if(id1!=id2){
        
        #now cut the logger data so its only for these two individuals:
        ld.u<-ciData[as.character(ciData$Mouse) %in% c(as.character(id1),as.character(id2)),]
        
        #now cut the logger data so its only the loggers where both individuals where seen at
        shared.loggers<-unique(ld.u$Reader[ld.u$Mouse%in%id1][ld.u$Reader[ld.u$Mouse%in%id1] %in% ld.u$Reader[ld.u$Mouse%in%id2]])
        
        #now only carry on the process if they were at some point seen on a logger (if never seen at same logger then its obviously 0)
        if(length(shared.loggers)>0){
          
          #now cut the dyads logger data to the loggers they were both seen on
          ld.u<-ld.u[as.character(ld.u$Reader) %in% shared.loggers,]
          
          #now order the dataframe by logger and time and ID
          ld.u<-ld.u[order(ld.u$Reader,ld.u$Datetime,ld.u$Mouse),]
          
          #now find the instances that we see one individual, then the another, within x minutes at each logger
          
          first<-1:(nrow(ld.u)-1)
          second<-first+1
          #same logger:
          same.logger<-ld.u$Reader[second]==ld.u$Reader[first]
          #different id:
          diff.id<-ld.u$Mouse[second]!=ld.u$Mouse[first]
          #within x minutes:
          within.x<-abs(as.numeric(difftime(ld.u$Datetime[second],ld.u$Datetime2[first], units="mins"))) < overlapTime
          
          #an instance is when all of these are true
          instances<- same.logger & diff.id & within.x
          
          if(sum(instances)>0){
            
            # Add a column with logger and night as a single term
            ld.u$lognight_logger <- paste(ld.u$Julian_Night,ld.u$Reader,sep="-")
            
            #Get rid of duplicated instances within logger_lognight combination
            instances[instances][duplicated(ld.u$lognight_logger[second][instances])] <- F
            
            if(sum(instances)>0){
              AM[as.character(id1),as.character(id2)]<-sum(instances)
              
              #Now, from here all this part is only needed if we actually have some instances, so start another if command
              
              #and we can also store info in our GBI matrix in case that'll be helpful in future
              add.to<-nrow(GBI)+1
              GBI.u<-matrix(0,sum(instances),length(ids),dimnames=list(add.to:(nrow(GBI)+sum(instances)),ids))
              GBI.u[,c(as.character(id1),as.character(id2))]<-1
              
              #now store some info about these instances
              gd.u<-data.frame(logger=as.character(ld.u$Reader[second][instances]),
                               id1=id1,id2=id2,
                               id1.time=as.character(ld.u$Datetime2[second-1][instances]),
                               id2.time=as.character(ld.u$Datetime[second][instances]),
                               instance=nrow(gd):(nrow(gd)+(length(instances[instances])-1)),
                               lognight=as.character(ld.u$Julian_Night[second][instances]),
                               stringsAsFactors=F) 
              
              
              #now bind onto master objects:
              gd<-rbind(gd,gd.u)
              GBI<-rbind(GBI,GBI.u)
            } #ends the first 'if we some instances'
          } #ends the second 'if we have some instances'
        } #ends the 'if we have some shared loggers'
      } #ends the 'if we have 2 different individuals'
    } #ends the j loop
    #now before we end the i loop we might as well print some info that might be interesting to see how the loop is progressing
    print(paste0("individual_",i," total_instances_",nrow(gd)-1))
  } #ends the i loop
  
  #Now the whole loop is over, take off the 'ghost' first row of the gd, and then all done
  gd<-gd[-1,]
  GBI<-GBI[-1,]
  rownames(GBI)<-c(1:nrow(GBI))
  
  gd$lognight_logger<-paste(gd$lognight,gd$logger,sep="-")
  
  # C) CONSTRUCTING the GBI and GD matrices from loggerdata: Finding
  # all instances that each individuals was observed but not associated with anyone.
  #MAKE SOLITARY GD data frame
  ld.a<-ciData[,c("Mouse","Reader","Julian_Night")]
  ld.a<-unique(ld.a)
  ld.a$lognight_logger<-paste(ld.a$Julian_Night,ld.a$Reader,sep="-")
  ld.a$lognight_logger_ID<-paste(ld.a$lognight_logger,ld.a$Mouse)
  
  gd_all<-ld.a 
  length(ld.a$lognight_logger_ID)==length(unique(ld.a$lognight_logger_ID))
  gd_social1<-paste(gd$lognight_logger,gd$id1)
  gd_social2<-paste(gd$lognight_logger,gd$id2)
  all(gd_social2%in%gd_social1) #all true so can use just the other...
  
  gd_social<-paste(gd$lognight_logger,gd$id1)
  gd_social<-unique(gd_social)#
  
  gd_solitary<-gd_all[which(!gd_all$lognight_logger_ID%in%gd_social),]
  # It's possible that gd_solitary is empty.  Therefore we need a workaround.
  if(nrow(gd_solitary)>0){
    
    #double gd2 because pairwise GBI also has double observations per id
    library(splitstackshape)
    gd_solitary$freq<-2
    gd_solitary<-expandRows(gd_solitary, "freq")
    
    #Order gd_solitary to the same order as GBI
    colnames(gd_solitary)[1]<-"id1"
    iddf<-data.frame(id1=ids,idno=c(1:length(ids)))
    gd_solitary<-merge(gd_solitary,iddf,by="id1",all.x=T)
    gd_solitary<-gd_solitary[order(gd_solitary$idno),]
    gd_solitary<-gd_solitary[,1:4]
    
    gd_solitary$id2<-NA
    gd_solitary$id1.time<-NA
    gd_solitary$id2.time<-NA
    gd_solitary$instance<-c((length(gd$instance)+1):(length(gd$instance)+nrow(gd_solitary)))
    
    colnames(gd_solitary)<-c("id1","logger","lognight",
                             "lognight_logger","id2",
                             "id1.time","id2.time","instance")
    
    # MAKE SOLITARY GBI matrix
    GBI2<-matrix(0,1,length(ids),dimnames=list(0,ids))
    colnames(GBI2)<-ids
    suminstances<-as.data.frame(table(gd_solitary$id1))
    missing<-ids[which(!ids%in%suminstances$Var1)]
    suminstances$Var1 <- as.integer(as.character(suminstances$Var1)) # See note below
    if(length(missing)>0){
      dfm<-data.frame("Var1"=missing, "Freq"=0)
      suminstances<-rbind(suminstances,dfm)
    }
    # R was treating Var1 as a factor and not allowing the rbind() below to behave properly
    # due to invalid factor levels.  Replacing the factors led to IDs being associated
    # with the wrong frequencies, so I had to take the column out of factor status.
    
    for (i in 1:length(ids)){
      add.to<-nrow(GBI2)+1
      sumi<-(suminstances[which(suminstances$Var1==ids[i]),]$Freq) #doubling instances because pairvise GBI also has double observations per id
      if(length(sumi)==0){sumi <- 0}
      if(sumi>0){
        
        GBI2.u<-matrix(0,sumi,length(ids),dimnames=list(add.to:(nrow(GBI2)+sumi),ids))
        GBI2.u[,as.character(ids[i])]<-1
        GBI2<-rbind(GBI2,GBI2.u)
      }
    }
    if(nrow(GBI2)>1){ # Sometimes GBI2 ends up being empty, it seems, so I need to make
      # a workaround.
      GBI2<-GBI2[2:nrow(GBI2),]
      rownames(GBI2)<-gd_solitary$instance
      
      #Bind solitary GBI/gd to social GBI/gd
      
      GBI3<-rbind(GBI,GBI2)
      gd3<-rbind(gd,gd_solitary)
    } # remove braces if it doesn't work.
    if(nrow(GBI2)==1){ # New
      GBI3 <- GBI
      gd3 <- gd
    }# End new
  }
  if(nrow(gd_solitary)==0){
    GBI3 <- GBI
    gd3 <- gd
  } # Seems to work
  # To make networks from the GBI, the method requires two vectors.
  # The first gives the sampling period (i.e. which night) and place (i.e. which reader)
  # the overlap occurred in.
  # The second just gives the sampling period.
  # We make those here to allow easy movement from this function to the next.
  samp.period.place<-gd3$lognight_logger
  samp.period2<-as.numeric(as.character(gd3$lognight))
  
  return(list(GBI=GBI3,
              gd=gd3,
              samp.period.place=samp.period.place,
              samp.period2=samp.period2))
}

##' This function creates the actual social network based on GBI and gd data generated
##' with association.info().  It can use a few different edge weight metrics and can
##' handle individuals that lost RFID tags or went missing using the "time.controlled"
##' argument.  It can return the edge weights in a matrix, or it can return two separate
##' matrices, the numerator and denominator for calculating SRI.  This latter is useful
##' if adding together multiple different time frames for generating a network.
##' The code for this function is a slightly modified version of the code from Raulo et
##' al. (2021).  All credit to the original authors for developing the method and writing 
##' most of the code that is used herein.  We have retained their explanatory notes and 
##' added a few of our own.
##' 
##' @param gbi A group-by-identity matrix for the individuals of interest
##' @param index The index for calculating edge weights.  Values include "SRI",
##' "SRI_successes", "SRI_fails", "HWI", and "Fractions"
##' @param time.controlled A vector with TRUE values for all individuals that are
##' receiving time control because they lost an RFID or vanished without ever being
##' recaptured.
##' @param samp.period A vector for identifying sampling period of overlap for each row
##' in the GBI
##' @param samp.period.place A vector for identifying sampling period and place of overlap
##' for each row in the GBI
##' 
##' @return A matrix (or list with two matrices if index=="Fractions") with social
##' association strengths for all pairs.

time.network.aura2.modified<-function(gbi,index,time.controlled,samp.period,samp.period.place){
  am<-ya<-yb <- yab <-xab<-matrix(0,ncol(gbi),ncol(gbi))
  if(!is.numeric(samp.period)){
    print("samp.period should be numeric")
    samp.period<-as.numeric(samp.period)}
  
  for(i in 1:ncol(gbi)){
    sums<-gbi[which(gbi[,i]>0),]
    if(is.null(dim(sums))==T) {xab[,i]<-sums} else {xab[,i]<-(colSums(sums)/2)} #THIS IS "/2" because the GBI contains repeats of the grouping events, where ind 1 and ind 2 are replicated as ind 2 and ind 1
    i.samps.a<-unique(samp.period.place[which(gbi[,i]>0)]) # unique sampling periods and
    # places for individual i
    # if(time.controlled[i]==F){ya[,i]<-length(i.samps.a)}
    
    i.samps<-as.numeric(unique(samp.period[which(gbi[,i]>0)]))
    
    for(j in 1:ncol(gbi)){
      j.samps.a<-unique(samp.period.place[which(gbi[,j]>0)])
      # The following line is necessary for proper calculation of associations when both
      # mice appear but don't overlap
      yab[j,i] <- length(intersect(i.samps.a,j.samps.a))
      if(time.controlled[j]==FALSE&time.controlled[i]==FALSE){
        yb[j,i]<-length(j.samps.a)
        ya[j,i] <- length(i.samps.a)
      } # This code chunk is not in the original function
      j.samps<-as.numeric(unique(samp.period[which(gbi[,j]>0)]))
      
      
      if(time.controlled[i]==TRUE|time.controlled[j]==TRUE){
        # This condition was originally just "time.controlled==T"
        if(length(i.samps)==0)print(paste("No observations of individual",i))
        if(length(j.samps)==0)print(paste("No observations of individual",j))
        # I need to throw in some code here to smooth this over
        # Otherwise I end up with a bunch of warnings from the below code.
        shared.min<-max(c(min(i.samps),min(j.samps)))
        shared.max<-min(c(max(i.samps),max(j.samps)))
        
        if(!(sum(gbi[,i])==0 | sum(gbi[,j])==0)){
          if(shared.min>shared.max){ #if no overlap - marked as NA now, depends on
          	# biological question whether these zeros (never had a chance to interact)
          	# are inherently different from the "real zeros" (had a chance but didn't
          	# interact)
            ya[j,i]<-yb[j,i]<- yab[j,i] <- NA} else{
              shared.range<-shared.min:shared.max
              i.sum<-sum(i.samps.a%in%samp.period.place[samp.period%in%shared.range])
              j.sum<-sum(j.samps.a%in%samp.period.place[samp.period%in%shared.range])
              ya[j,i]<-i.sum
              yb[j,i]<-j.sum
              # This seems to work, but I'm a bit skittish about it.
            }
        }
      }
    }
  }
  yab <- yab-xab # To only get the number of locations at which individuals are double-counted
  # in ya and yb
  # if(time.controlled==F){yb <- t(yb)}
  # This ensures that ya and yb are not copies of each other.  Without this step, if
  # time.controlled==F, then the "SRI" output is asymmetrical, which makes no sense.
  if (index=="Fractions"){
    num <- xab
    denom <- (xab+(ya-xab)+(yb-xab)-yab)
    diag(num)<-0;num[is.na(num)]<-0;rownames(num)<-colnames(num)<-colnames(gbi)
    diag(denom)<-0;denom[is.na(denom)]<-0;rownames(denom)<-colnames(denom)<-colnames(gbi)
    # You may not wish to treat the NAs as 0s, depending on whether individuals did or
    # did not have some chance to overlap.
    return(list(Num=num,Denom=denom))
  }
  else{
    if(index=="SRI_successes"){
      am<-xab}
    if(index=="SRI_fails"){
      am<-(ya-xab)+((yb)-xab)-yab}
    if(index=="SRI"){
      am<-xab/(xab+(ya-xab)+((yb)-xab)-yab)}
    if(index=="HWI"){
      am<-xab/(xab+(0.5*((ya-xab)+((yb)-xab)-yab)))}
    diag(am)<-0;am[is.nan(am)]<-0;rownames(am)<-colnames(am)<-colnames(gbi)
    return(am)
  }
}

##' This function calculates an association matrix from data for the same individuals
##' across different time intervals.  Its principal purpose is to allow handling of
##' switches in tags and the variation in time control needed to handle such a scenario.
##' 
##' @param data A list of lists, each sub-list formatted with "Num" and "Denom" as in the
##' output from the "Fractions" option of time.network.aura2.modified()
##' 
##' @return The association matrix for the individuals across all time intervals, along
##' with the numerator and denominator for calculating said matrix.

SRI.intervals <- function(data){
  numerator <- denominator <- matrix(0,ncol=ncol(data[[1]]$Num),
                                     nrow=nrow(data[[1]]$Num))
  for(i in 1:length(data)){
    numerator <- numerator+data[[i]]$Num
    denominator <- denominator+data[[i]]$Denom
  }
  AM <- numerator/denominator
  diag(AM) <- 1 # Each individual perfectly associates with themselves
  return(list(AM=AM,Num=numerator,Denom=denominator))
}

##' The code here takes an association matrix and converts it into a data frame with
##' the IDs in columns and the association in a third column.  It should be useful
##' for various analyses.  It also has an option for use with an output from
##' time.network.aura2.modified() with the index "Fractions" – i.e. for us with an output
##' containing a numerator matrix and denominator matrix.
##' 
##' @param data An association matrix of any kind with row names and column names giving
##' IDs of individuals, or a list containing multiple matrices (AM, Num, Denom)
##' @param NumDenom A parameter for when you have a list with numerator and denominator
##' matrices
##' 
##' @return a data frame with two ID columns and one column denoting the association
##' 

collate.matrix <- function(data,NumDenom){
  if(!is.logical(NumDenom)){
    print("NumDenom must be logical.")
    return()
  }
  ifelse(NumDenom==FALSE,AM <- data,AM <- data$AM)
  id1=c()
  id2=c()
  SocAssoc=c()
  if(NumDenom==TRUE){
    Num=c()
    Denom=c()
  }
  for (i in 1:(ncol(AM)-1)){
    for (j in (i+1):nrow(AM)){
      id1 <- c(id1,colnames(AM)[i])
      id2 <- c(id2,rownames(AM)[j])
      SocAssoc <- c(SocAssoc,AM[j,i])
      if(NumDenom==TRUE){
        Num <- c(Num,data$Num[j,i])
        Denom <- c(Denom,data$Denom[j,i])
      }
    }
  }
  ifelse(NumDenom==TRUE,
         output <- data.frame(id1,id2,SocAssoc,Num,Denom),
         output <- data.frame(id1,id2,SocAssoc))
  return(output)
}

##' This function looks at whether the individuals were housed in the same cage and whether
##' they are of the same genotype.  This information can be tacked on to the data frame
##' for analysis of links.  It's quite specifically tailored to our data formatting.
##' 
##' @param assoc.data A data frame output of the sort from collate.matrix, although it can
##' feature other columns
##' @param metadata A data frame containing cage and genotype metadata for the IDs in the
##' assoc.data data frame; requires columns named "Genotype" and "Cage"
##' @param infectionCheck A logical TRUE/FALSE determining whether to look at infection
##' status for overlap in the dyad; requires a column named "Infected" in the metadata
##' 
##' @return a data frame with information about whether individuals share a genotype or
##' cage

overlap.check <- function(assoc.data,metadata,infectionCheck){
  rowNum <- nrow(assoc.data)
  genotypeOverlap=rep(TRUE,rowNum)
  cageOverlap=rep(TRUE,rowNum)
  if(infectionCheck==TRUE){infection=rep(NA,rowNum)}
  for (i in 1:rowNum){
    id1 <- assoc.data$id1[i]
    id2 <- assoc.data$id2[i]
    genotype1 <- metadata$Genotype[which(metadata$Mouse==id1)]
    genotype2 <- metadata$Genotype[which(metadata$Mouse==id2)]
    Cage1 <- metadata$Cage[which(metadata$Mouse==id1)]
    Cage2 <- metadata$Cage[which(metadata$Mouse==id2)]
    ifelse(genotype1==genotype2,genotypeOverlap[i] <- TRUE,genotypeOverlap[i] <- FALSE)
    ifelse(Cage1==Cage2,cageOverlap[i] <- TRUE,cageOverlap[i] <- FALSE)
    if(infectionCheck==TRUE){
      infected1 <- metadata$Infected[which(metadata$Mouse==id1)]
      infected2 <- metadata$Infected[which(metadata$Mouse==id2)]
      if(infected1=="Y"&infected2=="Y"){
        infection[i] <- "Both"
      }
      if(infected1!=infected2){
        infection[i] <- "One"
      }
      if(infected1=="N"&infected2=="N"){
        infection[i] <- "Neither"
      }
    }
  }
  ifelse(infectionCheck==TRUE,
         assoc.data <- cbind(assoc.data,genotypeOverlap,cageOverlap,infection),
         assoc.data <- cbind(assoc.data,genotypeOverlap,cageOverlap))
  return(assoc.data)

# Now that we have all these functions, we can actually calculate the networks.
# It's important to ensure that your check-in data has Julian night info,
# as well as info about when mice might have lost RFIDs (the Lost.RFID.1 and
# Lost.RFID.2 columns in the mice metadata data frames).

# Step 1: generate a GBI (here with a fifteen-minute overlap window)
W2.2021b1.1.info.15 <- association.info(W2.CI.2021b1.1,15,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.1.info.15 # simply for saving space
# Step 2: Get numerator and denominator for SRI
W2.2021b1.1.fracs.15 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
# Step 1: generate a GBI
W2.2021b1.2.info.15 <- association.info(W2.CI.2021b1.2,15,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.2.info.15
# Step 2: Get numerator and denominator for SRI
W2.2021b1.2.fracs.15 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
# Step 3: Calculate single association matrix from the two chunks of time.
W2.2021b1.AM.15 <- SRI.intervals(list(W2.2021b1.1.fracs.15,
                                   W2.2021b1.2.fracs.15))

# Step 4: Generate a data table from the association matrix
W2.2021b1.assoc.info <- collate.matrix(W2.20212b1.AM.15,NumDenom==TRUE)

# You can iterate Steps 1–4 out over and over to generate association matrices
# for each enclosure and each block, using different overlap windows and different
# sets of locations.  Sadly we did this all manually, for all the different
# overlap window lengths and enclosures.  It translates to a lot of code, all
# of which is shared below.  There is certainly an easier way to do it.

# The key is to generate data tables that join all this information together

# For the 2021 Block 1 data, a 60-minute overlap window
W2.2021b1.1.info.60 <- association.info(W2.CI.2021b1.1,60,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.1.info.60
W2.2021b1.1.fracs.60 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b1.2.info.60 <- association.info(W2.CI.2021b1.2,60,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.2.info.60
W2.2021b1.2.fracs.60 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b1.AM <- SRI.intervals(list(W2.2021b1.1.fracs.60,
                                   W2.2021b1.2.fracs.60))

W3.2021b1.1.info.60 <- association.info(W3.CI.2021b1.1,60,
                                        filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.1.info.60
W3.2021b1.1.fracs.60 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b1.2.info.60 <- association.info(W3.CI.2021b1.2,60,
                                        filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.2.info.60
W3.2021b1.2.fracs.60 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b1.AM <- SRI.intervals(list(W3.2021b1.1.fracs.60,
                                   W3.2021b1.2.fracs.60))

W4.2021b1.1.info.60 <- association.info(W4.CI.2021b1.1,60,
                                        filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.1.info.60
W4.2021b1.1.fracs.60 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b1.2.info.60 <- association.info(W4.CI.2021b1.2,60,
                                        filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.2.info.60
W4.2021b1.2.fracs.60 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b1.AM <- SRI.intervals(list(W4.2021b1.1.fracs.60,
                                   W4.2021b1.2.fracs.60))

# For the 2021 Block 1 data, but just feeders, a 60-minute overlap window
W2.2021b1.1.info.60.Feeder <- association.info(filter(W2.CI.2021b1.1,Location=="Feeder"),
                                               60,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.1.info.60.Feeder
W2.2021b1.1.fracs.60.Feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b1.2.info.60.Feeder <- association.info(filter(W2.CI.2021b1.2,Location=="Feeder"),
                                               60,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.2.info.60.Feeder
W2.2021b1.2.fracs.60.Feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b1.AM.Feeder <- SRI.intervals(list(W2.2021b1.1.fracs.60.Feeder,
                                          W2.2021b1.2.fracs.60.Feeder))

W3.2021b1.1.info.60.Feeder <- association.info(filter(W3.CI.2021b1.1,Location=="Feeder"),
                                               60,
                                        filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.1.info.60.Feeder
W3.2021b1.1.fracs.60.Feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b1.2.info.60.Feeder <- association.info(filter(W3.CI.2021b1.2,Location=="Feeder"),
                                                      60,
                                        filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.2.info.60.Feeder
W3.2021b1.2.fracs.60.Feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b1.AM.Feeder <- SRI.intervals(list(W3.2021b1.1.fracs.60.Feeder,
                                   W3.2021b1.2.fracs.60.Feeder))

W4.2021b1.1.info.60.Feeder <- association.info(filter(W4.CI.2021b1.1,Location=="Feeder"),
                                               60,
                                        filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.1.info.60.Feeder
W4.2021b1.1.fracs.60.Feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b1.2.info.60.Feeder <- association.info(filter(W4.CI.2021b1.2,Location=="Feeder"),
                                               60,
                                        filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.2.info.60.Feeder
W4.2021b1.2.fracs.60.Feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b1.AM.Feeder <- SRI.intervals(list(W4.2021b1.1.fracs.60.Feeder,
                                   W4.2021b1.2.fracs.60.Feeder))

# For the 2021, Block 2 data, a 60-minute overlap window
W2.2021b2.1.info.60 <- association.info(W2.CI.2021b2.1,60,
                                        filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.1.info.60
W2.2021b2.1.fracs.60 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b2.2.info.60 <- association.info(W2.CI.2021b2.2,60,
                                        filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.2.info.60
W2.2021b2.2.fracs.60 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b2.AM <- SRI.intervals(list(W2.2021b2.1.fracs.60,
                                   W2.2021b2.2.fracs.60))

W3.2021b2.1.info.60 <- association.info(W3.CI.2021b2.1,60,
                                        filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.1.info.60
W3.2021b2.1.fracs.60 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b2.2.info.60 <- association.info(W3.CI.2021b2.2,60,
                                        filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.2.info.60
W3.2021b2.2.fracs.60 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b2.AM <- SRI.intervals(list(W3.2021b2.1.fracs.60,
                                   W3.2021b2.2.fracs.60))

W4.2021b2.1.info.60 <- association.info(W4.CI.2021b2.1,60,
                                        filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.1.info.60
W4.2021b2.1.fracs.60 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b2.2.info.60 <- association.info(W4.CI.2021b2.2,60,
                                        filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.2.info.60
W4.2021b2.2.fracs.60 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b2.AM <- SRI.intervals(list(W4.2021b2.1.fracs.60,
                                   W4.2021b2.2.fracs.60))

# For the 2021 Block 2 data, but just feeders, 60-minute overlap window
W2.2021b2.1.info.60.Feeder <- association.info(filter(W2.CI.2021b2.1,Location=="Feeder"),
                                               60,
                                               filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.1.info.60.Feeder
W2.2021b2.1.fracs.60.Feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.1,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W2.2021b2.2.info.60.Feeder <- association.info(filter(W2.CI.2021b2.2,Location=="Feeder"),
                                               60,
                                               filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.2.info.60.Feeder
W2.2021b2.2.fracs.60.Feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.2,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W2.2021b2.AM.Feeder <- SRI.intervals(list(W2.2021b2.1.fracs.60.Feeder,
                                          W2.2021b2.2.fracs.60.Feeder))

W3.2021b2.1.info.60.Feeder <- association.info(filter(W3.CI.2021b2.1,Location=="Feeder"),
                                               60,
                                               filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.1.info.60.Feeder
W3.2021b2.1.fracs.60.Feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.1,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W3.2021b2.2.info.60.Feeder <- association.info(filter(W3.CI.2021b2.2,Location=="Feeder"),
                                               60,
                                               filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.2.info.60.Feeder
W3.2021b2.2.fracs.60.Feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.2,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W3.2021b2.AM.Feeder <- SRI.intervals(list(W3.2021b2.1.fracs.60.Feeder,
                                          W3.2021b2.2.fracs.60.Feeder))

W4.2021b2.1.info.60.Feeder <- association.info(filter(W4.CI.2021b2.1,Location=="Feeder"),
                                               60,
                                               filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.1.info.60.Feeder
W4.2021b2.1.fracs.60.Feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.1,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W4.2021b2.2.info.60.Feeder <- association.info(filter(W4.CI.2021b2.2,Location=="Feeder"),
                                               60,
                                               filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.2.info.60.Feeder
W4.2021b2.2.fracs.60.Feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.2,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W4.2021b2.AM.Feeder <- SRI.intervals(list(W4.2021b2.1.fracs.60.Feeder,
                                          W4.2021b2.2.fracs.60.Feeder))

# Here is the 2021, block 1, data, but without feeders, 60-minute overlap window
W2.2021b1.1.info.noF <- association.info(filter(W2.CI.2021b1.1,Location!="Feeder"),60,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.1.info.noF
W2.2021b1.1.fracs.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b1.2.info.noF <- association.info(filter(W2.CI.2021b1.2,Location!="Feeder"),60,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.2.info.noF
W2.2021b1.2.fracs.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b1.AM.noF <- SRI.intervals(list(W2.2021b1.1.fracs.noF,
                                   W2.2021b1.2.fracs.noF))

W3.2021b1.1.info.noF <- association.info(filter(W3.CI.2021b1.1,Location!="Feeder"),60,
                                        filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.1.info.noF
W3.2021b1.1.fracs.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b1.2.info.noF <- association.info(filter(W3.CI.2021b1.2,Location!="Feeder"),60,
                                        filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.2.info.noF
W3.2021b1.2.fracs.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b1.AM.noF <- SRI.intervals(list(W3.2021b1.1.fracs.noF,
                                   W3.2021b1.2.fracs.noF))

W4.2021b1.1.info.noF <- association.info(filter(W4.CI.2021b1.1,Location!="Feeder"),60,
                                        filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.1.info.noF
W4.2021b1.1.fracs.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b1.2.info.noF <- association.info(filter(W4.CI.2021b1.2,Location!="Feeder"),60,
                                        filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.2.info.noF
W4.2021b1.2.fracs.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b1.AM.noF <- SRI.intervals(list(W4.2021b1.1.fracs.noF,
                                   W4.2021b1.2.fracs.noF))

# Here are the 2021, Block 2, data without feeders, 60-minute overlap window
W2.2021b2.1.info.noF <- association.info(filter(W2.CI.2021b2.1,Location!="Feeder"),60,
                                         filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.1.info.noF
W2.2021b2.1.fracs.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                     filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.1,
                                                     d$samp.period2,
                                                     d$samp.period.place)
W2.2021b2.2.info.noF <- association.info(filter(W2.CI.2021b2.2,Location!="Feeder"),60,
                                         filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.2.info.noF
W2.2021b2.2.fracs.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                     filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.2,
                                                     d$samp.period2,
                                                     d$samp.period.place)
W2.2021b2.AM.noF <- SRI.intervals(list(W2.2021b2.1.fracs.noF,
                                   W2.2021b2.2.fracs.noF))

W3.2021b2.1.info.noF <- association.info(filter(W3.CI.2021b2.1,Location!="Feeder"),60,
                                         filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.1.info.noF
W3.2021b2.1.fracs.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                     filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.1,
                                                     d$samp.period2,
                                                     d$samp.period.place)
W3.2021b2.2.info.noF <- association.info(filter(W3.CI.2021b2.2,Location!="Feeder"),60,
                                         filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.2.info.noF
W3.2021b2.2.fracs.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                     filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.2,
                                                     d$samp.period2,
                                                     d$samp.period.place)
W3.2021b2.AM.noF <- SRI.intervals(list(W3.2021b2.1.fracs.noF,
                                   W3.2021b2.2.fracs.noF))

W4.2021b2.1.info.noF <- association.info(filter(W4.CI.2021b2.1,Location!="Feeder"),60,
                                         filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.1.info.noF
W4.2021b2.1.fracs.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                     filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.1,
                                                     d$samp.period2,
                                                     d$samp.period.place)
W4.2021b2.2.info.noF <- association.info(filter(W4.CI.2021b2.2,Location!="Feeder"),60,
                                         filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.2.info.noF
W4.2021b2.2.fracs.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                     filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.2,
                                                     d$samp.period2,
                                                     d$samp.period.place)
W4.2021b2.AM.noF <- SRI.intervals(list(W4.2021b2.1.fracs.noF,
                                   W4.2021b2.2.fracs.noF))

# For the 2021 Block 1 data, four-hour time window
W2.2021b1.1.info.240 <- association.info(W2.CI.2021b1.1,240,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.1.info.240
W2.2021b1.1.fracs.240 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b1.2.info.240 <- association.info(W2.CI.2021b1.2,240,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.2.info.240
W2.2021b1.2.fracs.240 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b1.AM.240 <- SRI.intervals(list(W2.2021b1.1.fracs.240,
                                   W2.2021b1.2.fracs.240))

W3.2021b1.1.info.240 <- association.info(W3.CI.2021b1.1,240,
                                        filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.1.info.240
W3.2021b1.1.fracs.240 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b1.2.info.240 <- association.info(W3.CI.2021b1.2,240,
                                        filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.2.info.240
W3.2021b1.2.fracs.240 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b1.AM.240 <- SRI.intervals(list(W3.2021b1.1.fracs.240,
                                   W3.2021b1.2.fracs.240))

W4.2021b1.1.info.240 <- association.info(W4.CI.2021b1.1,240,
                                        filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.1.info.240
W4.2021b1.1.fracs.240 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b1.2.info.240 <- association.info(W4.CI.2021b1.2,240,
                                        filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.2.info.240
W4.2021b1.2.fracs.240 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b1.AM.240 <- SRI.intervals(list(W4.2021b1.1.fracs.240,
                                   W4.2021b1.2.fracs.240))

# For the 2021, Block 2 data, with four-hour intervals
W2.2021b2.1.info.240 <- association.info(W2.CI.2021b2.1,240,
                                        filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.1.info.240
W2.2021b2.1.fracs.240 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b2.2.info.240 <- association.info(W2.CI.2021b2.2,240,
                                        filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.2.info.240
W2.2021b2.2.fracs.240 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b2.AM.240 <- SRI.intervals(list(W2.2021b2.1.fracs.240,
                                   W2.2021b2.2.fracs.240))

W3.2021b2.1.info.240 <- association.info(W3.CI.2021b2.1,240,
                                        filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.1.info.240
W3.2021b2.1.fracs.240 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b2.2.info.240 <- association.info(W3.CI.2021b2.2,240,
                                        filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.2.info.240
W3.2021b2.2.fracs.240 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b2.AM.240 <- SRI.intervals(list(W3.2021b2.1.fracs.240,
                                   W3.2021b2.2.fracs.240))

W4.2021b2.1.info.240 <- association.info(W4.CI.2021b2.1,240,
                                        filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.1.info.240
W4.2021b2.1.fracs.240 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b2.2.info.240 <- association.info(W4.CI.2021b2.2,240,
                                        filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.2.info.240
W4.2021b2.2.fracs.240 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b2.AM.240 <- SRI.intervals(list(W4.2021b2.1.fracs.240,
                                   W4.2021b2.2.fracs.240))

# For the 2021 Block 1 data, with 15-minute intervals
W2.2021b1.1.info.15 <- association.info(W2.CI.2021b1.1,15,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.1.info.15
W2.2021b1.1.fracs.15 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b1.2.info.15 <- association.info(W2.CI.2021b1.2,15,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.2.info.15
W2.2021b1.2.fracs.15 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b1.AM.15 <- SRI.intervals(list(W2.2021b1.1.fracs.15,
                                   W2.2021b1.2.fracs.15))

W3.2021b1.1.info.15 <- association.info(W3.CI.2021b1.1,15,
                                        filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.1.info.15
W3.2021b1.1.fracs.15 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b1.2.info.15 <- association.info(W3.CI.2021b1.2,15,
                                        filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.2.info.15
W3.2021b1.2.fracs.15 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b1.AM.15 <- SRI.intervals(list(W3.2021b1.1.fracs.15,
                                   W3.2021b1.2.fracs.15))

W4.2021b1.1.info.15 <- association.info(W4.CI.2021b1.1,15,
                                        filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.1.info.15
W4.2021b1.1.fracs.15 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b1.2.info.15 <- association.info(filter(W4.CI.2021b1.2,
                                               !(Julian_Night%in%c(25:31))),15,
                                        filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.2.info.15
W4.2021b1.2.fracs.15 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b1.AM.15 <- SRI.intervals(list(W4.2021b1.1.fracs.15,
                                   W4.2021b1.2.fracs.15))

# For the 2021, Block 2 data, with fifteen-minute intervals
W2.2021b2.1.info.15 <- association.info(W2.CI.2021b2.1,15,
                                        filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.1.info.15
W2.2021b2.1.fracs.15 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b2.2.info.15 <- association.info(W2.CI.2021b2.2,15,
                                        filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.2.info.15
W2.2021b2.2.fracs.15 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b2.AM.15 <- SRI.intervals(list(W2.2021b2.1.fracs.15,
                                   W2.2021b2.2.fracs.15))

W3.2021b2.1.info.15 <- association.info(W3.CI.2021b2.1,15,
                                        filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.1.info.15
W3.2021b2.1.fracs.15 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b2.2.info.15 <- association.info(W3.CI.2021b2.2,15,
                                        filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.2.info.15
W3.2021b2.2.fracs.15 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b2.AM.15 <- SRI.intervals(list(W3.2021b2.1.fracs.15,
                                   W3.2021b2.2.fracs.15))

W4.2021b2.1.info.15 <- association.info(W4.CI.2021b2.1,15,
                                        filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.1.info.15
W4.2021b2.1.fracs.15 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b2.2.info.15 <- association.info(W4.CI.2021b2.2,15,
                                        filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.2.info.15
W4.2021b2.2.fracs.15 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b2.AM.15 <- SRI.intervals(list(W4.2021b2.1.fracs.15,
                                   W4.2021b2.2.fracs.15))

# For the 2021 Block 1 data, with two-minute intervals
W2.2021b1.1.info.2 <- association.info(W2.CI.2021b1.1,2,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.1.info.2
W2.2021b1.1.fracs.2 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b1.2.info.2 <- association.info(W2.CI.2021b1.2,2,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.2.info.2
W2.2021b1.2.fracs.2 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b1.AM.2min <- SRI.intervals(list(W2.2021b1.1.fracs.2,
                                   W2.2021b1.2.fracs.2))

W3.2021b1.1.info.2 <- association.info(W3.CI.2021b1.1,2,
                                        filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.1.info.2
W3.2021b1.1.fracs.2 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b1.2.info.2 <- association.info(W3.CI.2021b1.2,2,
                                        filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.2.info.2
W3.2021b1.2.fracs.2 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b1.AM.2min <- SRI.intervals(list(W3.2021b1.1.fracs.2,
                                   W3.2021b1.2.fracs.2))

W4.2021b1.1.info.2 <- association.info(W4.CI.2021b1.1,2,
                                        filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.1.info.2
W4.2021b1.1.fracs.2 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b1.2.info.2 <- association.info(W4.CI.2021b1.2,2,
                                        filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.2.info.2
W4.2021b1.2.fracs.2 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b1.AM.2min <- SRI.intervals(list(W4.2021b1.1.fracs.2,
                                   W4.2021b1.2.fracs.2))

# For the 2021, Block 2 data, with two-minute intervals
W2.2021b2.1.info.2 <- association.info(W2.CI.2021b2.1,2,
                                        filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.1.info.2
W2.2021b2.1.fracs.2 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b2.2.info.2 <- association.info(W2.CI.2021b2.2,2,
                                        filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.2.info.2
W2.2021b2.2.fracs.2 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b2.AM.2min <- SRI.intervals(list(W2.2021b2.1.fracs.2,
                                   W2.2021b2.2.fracs.2))

W3.2021b2.1.info.2 <- association.info(W3.CI.2021b2.1,2,
                                        filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.1.info.2
W3.2021b2.1.fracs.2 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b2.2.info.2 <- association.info(W3.CI.2021b2.2,2,
                                        filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.2.info.2
W3.2021b2.2.fracs.2 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b2.AM.2min <- SRI.intervals(list(W3.2021b2.1.fracs.2,
                                   W3.2021b2.2.fracs.2))

W4.2021b2.1.info.2 <- association.info(W4.CI.2021b2.1,2,
                                        filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.1.info.2
W4.2021b2.1.fracs.2 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b2.2.info.2 <- association.info(W4.CI.2021b2.2,2,
                                        filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.2.info.2
W4.2021b2.2.fracs.2 <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b2.AM.2min <- SRI.intervals(list(W4.2021b2.1.fracs.2,
                                   W4.2021b2.2.fracs.2))

# For the 2021 Block 1 data, with 15-minute intervals and only feeders
W2.2021b1.1.info.15.feeder <- association.info(filter(W2.CI.2021b1.1,Location=="Feeder"),
                                               15,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.1.info.15.feeder
W2.2021b1.1.fracs.15.feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b1.2.info.15.feeder <- association.info(filter(W2.CI.2021b1.2,Location=="Feeder"),
                                               15,
                                        filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.2.info.15.feeder
W2.2021b1.2.fracs.15.feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b1.AM.15.feeder <- SRI.intervals(list(W2.2021b1.1.fracs.15.feeder,
                                      W2.2021b1.2.fracs.15.feeder))

W3.2021b1.1.info.15.feeder <- association.info(filter(W3.CI.2021b1.1,Location=="Feeder"),
                                               15,
                                        filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.1.info.15.feeder
W3.2021b1.1.fracs.15.feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b1.2.info.15.feeder <- association.info(filter(W3.CI.2021b1.2,Location=="Feeder"),
                                               15,
                                        filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.2.info.15.feeder
W3.2021b1.2.fracs.15.feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b1.AM.15.feeder <- SRI.intervals(list(W3.2021b1.1.fracs.15.feeder,
                                      W3.2021b1.2.fracs.15.feeder))

W4.2021b1.1.info.15.feeder <- association.info(filter(W4.CI.2021b1.1,Location=="Feeder"),
                                               15,
                                        filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.1.info.15.feeder
W4.2021b1.1.fracs.15.feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b1.2.info.15.feeder <- association.info(filter(W4.CI.2021b1.2,Location=="Feeder"),
                                               15,
                                        filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.2.info.15.feeder
W4.2021b1.2.fracs.15.feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b1.AM.15.feeder <- SRI.intervals(list(W4.2021b1.1.fracs.15.feeder,
                                      W4.2021b1.2.fracs.15.feeder))

# For the 2021, Block 2 data, with fifteen-minute intervals and feeders
W2.2021b2.1.info.15.feeder <- association.info(filter(W2.CI.2021b2.1,Location=="Feeder"),
                                               15,
                                        filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.1.info.15.feeder
W2.2021b2.1.fracs.15.feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b2.2.info.15.feeder <- association.info(filter(W2.CI.2021b2.2,Location=="Feeder"),
                                               15,
                                        filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.2.info.15.feeder
W2.2021b2.2.fracs.15.feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W2.2021b2.AM.15.feeder <- SRI.intervals(list(W2.2021b2.1.fracs.15.feeder,
                                      W2.2021b2.2.fracs.15.feeder))

W3.2021b2.1.info.15.feeder <- association.info(filter(W3.CI.2021b2.1,Location=="Feeder"),
                                               15,
                                        filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.1.info.15.feeder
W3.2021b2.1.fracs.15.feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b2.2.info.15.feeder <- association.info(filter(W3.CI.2021b2.2,Location=="Feeder"),
                                               15,
                                        filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.2.info.15.feeder
W3.2021b2.2.fracs.15.feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W3.2021b2.AM.15.feeder <- SRI.intervals(list(W3.2021b2.1.fracs.15.feeder,
                                      W3.2021b2.2.fracs.15.feeder))

W4.2021b2.1.info.15.feeder <- association.info(filter(W4.CI.2021b2.1,Location=="Feeder"),
                                               15,
                                        filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.1.info.15.feeder
W4.2021b2.1.fracs.15.feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.1,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b2.2.info.15.feeder <- association.info(filter(W4.CI.2021b2.2,Location=="Feeder"),
                                               15,
                                        filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.2.info.15.feeder
W4.2021b2.2.fracs.15.feeder <- time.network.aura2.modified(d$GBI,"Fractions",
                                                    filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.2,
                                                    d$samp.period2,
                                                    d$samp.period.place)
W4.2021b2.AM.15.feeder <- SRI.intervals(list(W4.2021b2.1.fracs.15.feeder,
                                      W4.2021b2.2.fracs.15.feeder))

# For the 2021 Block 1 data, with 15-minute intervals and only non-feeders
W2.2021b1.1.info.15.noF <- association.info(filter(W2.CI.2021b1.1,Location!="Feeder"),
                                               15,
                                               filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.1.info.15.noF
W2.2021b1.1.fracs.15.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.1,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W2.2021b1.2.info.15.noF <- association.info(filter(W2.CI.2021b1.2,Location!="Feeder"),
                                               15,
                                               filter(mice2021block1,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b1.2.info.15.noF
W2.2021b1.2.fracs.15.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block1,Wedge=="Wedge 2")$Lost.RFID.2,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W2.2021b1.AM.15.noF <- SRI.intervals(list(W2.2021b1.1.fracs.15.noF,
                                             W2.2021b1.2.fracs.15.noF))

W3.2021b1.1.info.15.noF <- association.info(filter(W3.CI.2021b1.1,Location!="Feeder"),
                                               15,
                                               filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.1.info.15.noF
W3.2021b1.1.fracs.15.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.1,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W3.2021b1.2.info.15.noF <- association.info(filter(W3.CI.2021b1.2,Location!="Feeder"),
                                               15,
                                               filter(mice2021block1,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b1.2.info.15.noF
W3.2021b1.2.fracs.15.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block1,Wedge=="Wedge 3")$Lost.RFID.2,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W3.2021b1.AM.15.noF <- SRI.intervals(list(W3.2021b1.1.fracs.15.noF,
                                             W3.2021b1.2.fracs.15.noF))

W4.2021b1.1.info.15.noF <- association.info(filter(W4.CI.2021b1.1,Location!="Feeder"),
                                               15,
                                               filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.1.info.15.noF
W4.2021b1.1.fracs.15.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.1,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W4.2021b1.2.info.15.noF <- association.info(filter(W4.CI.2021b1.2,Location!="Feeder"),
                                               15,
                                               filter(mice2021block1,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b1.2.info.15.noF
W4.2021b1.2.fracs.15.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block1,Wedge=="Wedge 4")$Lost.RFID.2,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W4.2021b1.AM.15.noF <- SRI.intervals(list(W4.2021b1.1.fracs.15.noF,
                                             W4.2021b1.2.fracs.15.noF))

# For the 2021, Block 2 data, with fifteen-minute intervals and only non-feeders
W2.2021b2.1.info.15.noF <- association.info(filter(W2.CI.2021b2.1,Location!="Feeder"),
                                               15,
                                               filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.1.info.15.noF
W2.2021b2.1.fracs.15.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.1,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W2.2021b2.2.info.15.noF <- association.info(filter(W2.CI.2021b2.2,Location!="Feeder"),
                                               15,
                                               filter(mice2021block2,Wedge=="Wedge 2")$Mouse)
d <- W2.2021b2.2.info.15.noF
W2.2021b2.2.fracs.15.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block2,Wedge=="Wedge 2")$Lost.RFID.2,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W2.2021b2.AM.15.noF <- SRI.intervals(list(W2.2021b2.1.fracs.15.noF,
                                             W2.2021b2.2.fracs.15.noF))

W3.2021b2.1.info.15.noF <- association.info(filter(W3.CI.2021b2.1,Location!="Feeder"),
                                               15,
                                               filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.1.info.15.noF
W3.2021b2.1.fracs.15.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.1,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W3.2021b2.2.info.15.noF <- association.info(filter(W3.CI.2021b2.2,Location!="Feeder"),
                                               15,
                                               filter(mice2021block2,Wedge=="Wedge 3")$Mouse)
d <- W3.2021b2.2.info.15.noF
W3.2021b2.2.fracs.15.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block2,Wedge=="Wedge 3")$Lost.RFID.2,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W3.2021b2.AM.15.noF <- SRI.intervals(list(W3.2021b2.1.fracs.15.noF,
                                             W3.2021b2.2.fracs.15.noF))

W4.2021b2.1.info.15.noF <- association.info(filter(W4.CI.2021b2.1,Location!="Feeder"),
                                               15,
                                               filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.1.info.15.noF
W4.2021b2.1.fracs.15.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.1,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W4.2021b2.2.info.15.noF <- association.info(filter(W4.CI.2021b2.2,Location!="Feeder"),
                                               15,
                                               filter(mice2021block2,Wedge=="Wedge 4")$Mouse)
d <- W4.2021b2.2.info.15.noF
W4.2021b2.2.fracs.15.noF <- time.network.aura2.modified(d$GBI,"Fractions",
                                                           filter(mice2021block2,Wedge=="Wedge 4")$Lost.RFID.2,
                                                           d$samp.period2,
                                                           d$samp.period.place)
W4.2021b2.AM.15.noF <- SRI.intervals(list(W4.2021b2.1.fracs.15.noF,
                                             W4.2021b2.2.fracs.15.noF))

# For the 2021 data in the first portion of the experiment with 15-minute windows
W2.2021b1.AM.1.15 <- SRI.intervals(list(W2.2021b1.1.fracs.15))
W3.2021b1.AM.1.15 <- SRI.intervals(list(W3.2021b1.1.fracs.15))
W4.2021b1.AM.1.15 <- SRI.intervals(list(W4.2021b1.1.fracs.15))
W6.2021b1.AM.1.15 <- SRI.intervals(list(W6.2021b1.1.fracs.15))
W7.2021b1.AM.1.15 <- SRI.intervals(list(W7.2021b1.1.fracs.15))

W2.2021b2.AM.1.15 <- SRI.intervals(list(W2.2021b2.1.fracs.15))
W3.2021b2.AM.1.15 <- SRI.intervals(list(W3.2021b2.1.fracs.15))
W4.2021b2.AM.1.15 <- SRI.intervals(list(W4.2021b2.1.fracs.15))
W6.2021b2.AM.1.15 <- SRI.intervals(list(W6.2021b2.1.fracs.15))
W7.2021b2.AM.1.15 <- SRI.intervals(list(W7.2021b2.1.fracs.15))

# For the 2021 data in the second portion of the experiment with 15-minute windows
W2.2021b1.AM.2.15 <- SRI.intervals(list(W2.2021b1.2.fracs.15))
W3.2021b1.AM.2.15 <- SRI.intervals(list(W3.2021b1.2.fracs.15))
W4.2021b1.AM.2.15 <- SRI.intervals(list(W4.2021b1.2.fracs.15))
W6.2021b1.AM.2.15 <- SRI.intervals(list(W6.2021b1.2.fracs.15))
W7.2021b1.AM.2.15 <- SRI.intervals(list(W7.2021b1.2.fracs.15))

W2.2021b2.AM.2.15 <- SRI.intervals(list(W2.2021b2.2.fracs.15))
W3.2021b2.AM.2.15 <- SRI.intervals(list(W3.2021b2.2.fracs.15))
W4.2021b2.AM.2.15 <- SRI.intervals(list(W4.2021b2.2.fracs.15))
W6.2021b2.AM.2.15 <- SRI.intervals(list(W6.2021b2.2.fracs.15))
W7.2021b2.AM.2.15 <- SRI.intervals(list(W7.2021b2.2.fracs.15))


# Now, having generated all of those networks, we can put all the network data
# together in table form.  When we do this, we create two columns for the two
# mice in the pair (id1 and id2) and place the social association level in a
# third.  In some cases we also want to put the numerator and denominator of
# the association metric in separate columns (if we want to model association
# as a binomially-distributed response variable, we need to do this).  In such
# cases, NumDenom=TRUE in collate.matrix().
# After doing this collation for each individual set of associations, we put
# everything together in data frames, adding new associations as new columns
# in the data frames with the pairs.  Ultimately we can combine everything into
# a single large data frame that will also get the immune similarity data.
# Again there is probably a faster way to do this, rather than manually.

# For the 2021, Block 1 data
W2.2021b1.assoc.info <- collate.matrix(W2.2021b1.AM,NumDenom=TRUE)
W2.2021b1.assoc.info <- overlap.check(W2.2021b1.assoc.info,mice2021block1,
                                      infectionCheck=TRUE)
W2.2021b1.assoc.info$Wedge <- "Wedge 2"
W2.2021b1.assoc.info$Sex <- "Female"

W3.2021b1.assoc.info <- collate.matrix(W3.2021b1.AM,NumDenom=TRUE)
W3.2021b1.assoc.info <- overlap.check(W3.2021b1.assoc.info,mice2021block1,
                                      infectionCheck=TRUE)
W3.2021b1.assoc.info$Wedge <- "Wedge 3"
W3.2021b1.assoc.info$Sex <- "Female"

W4.2021b1.assoc.info <- collate.matrix(W4.2021b1.AM,NumDenom=TRUE)
W4.2021b1.assoc.info <- overlap.check(W4.2021b1.assoc.info,mice2021block1,
                                      infectionCheck=TRUE)
W4.2021b1.assoc.info$Wedge <- "Wedge 4"
W4.2021b1.assoc.info$Sex <- "Female"

full.2021b1.assoc.info <- rbind(W2.2021b1.assoc.info,
                                W3.2021b1.assoc.info,
                                W4.2021b1.assoc.info)

# For the 2021, Block 2 data
W2.2021b2.assoc.info <- collate.matrix(W2.2021b2.AM,NumDenom=TRUE)
W2.2021b2.assoc.info <- overlap.check(W2.2021b2.assoc.info,mice2021block2,
                                      infectionCheck=TRUE)
W2.2021b2.assoc.info$Wedge <- "Wedge 2"
W2.2021b2.assoc.info$Sex <- "Female"

W3.2021b2.assoc.info <- collate.matrix(W3.2021b2.AM,NumDenom=TRUE)
W3.2021b2.assoc.info <- overlap.check(W3.2021b2.assoc.info,mice2021block2,
                                      infectionCheck=TRUE)
W3.2021b2.assoc.info$Wedge <- "Wedge 3"
W3.2021b2.assoc.info$Sex <- "Female"

W4.2021b2.assoc.info <- collate.matrix(W4.2021b2.AM,NumDenom=TRUE)
W4.2021b2.assoc.info <- overlap.check(W4.2021b2.assoc.info,mice2021block2,
                                      infectionCheck=TRUE)
W4.2021b2.assoc.info$Wedge <- "Wedge 4"
W4.2021b2.assoc.info$Sex <- "Female"

full.2021b2.assoc.info <- rbind(W2.2021b2.assoc.info,
                                W3.2021b2.assoc.info,
                                W4.2021b2.assoc.info)

# # For the 2021, Block 1, data in just the pre-trap period
# W2.2021b1.assoc.info.1 <- collate.matrix(W2.2021b1.AM.1,NumDenom=TRUE)
# W2.2021b1.assoc.info.1 <- overlap.check(W2.2021b1.assoc.info.1,mice2021block1,
#                                         infectionCheck=TRUE)
# W3.2021b1.assoc.info.1 <- collate.matrix(W3.2021b1.AM.1,NumDenom=TRUE)
# W3.2021b1.assoc.info.1 <- overlap.check(W3.2021b1.assoc.info.1,mice2021block1,
#                                         infectionCheck=TRUE)
# W4.2021b1.assoc.info.1 <- collate.matrix(W4.2021b1.AM.1,NumDenom=TRUE)
# W4.2021b1.assoc.info.1 <- overlap.check(W4.2021b1.assoc.info.1,mice2021block1,
#                                         infectionCheck=TRUE)
# full.2021b1.assoc.info.1 <- rbind(W2.2021b1.assoc.info.1,
#                                   W3.2021b1.assoc.info.1,
#                                   W4.2021b1.assoc.info.1)
# temp <- select(full.2021b1.assoc.info.1,id1:SocAssoc)
# d <- left_join(full.2021b1.assoc.info,temp,by=c("id1","id2"))
# d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc1=SocAssoc.y)
# full.2021b1.assoc.info <- d

# # For the 2021, Block 1, data in just the post-trap period
# W2.2021b1.assoc.info.2 <- collate.matrix(W2.2021b1.AM.2,NumDenom=TRUE)
# W2.2021b1.assoc.info.2 <- overlap.check(W2.2021b1.assoc.info.2,mice2021block1,
#                                         infectionCheck=TRUE)
# W3.2021b1.assoc.info.2 <- collate.matrix(W3.2021b1.AM.2,NumDenom=TRUE)
# W3.2021b1.assoc.info.2 <- overlap.check(W3.2021b1.assoc.info.2,mice2021block1,
#                                         infectionCheck=TRUE)
# W4.2021b1.assoc.info.2 <- collate.matrix(W4.2021b1.AM.2,NumDenom=TRUE)
# W4.2021b1.assoc.info.2 <- overlap.check(W4.2021b1.assoc.info.2,mice2021block1,
#                                         infectionCheck=TRUE)
# full.2021b1.assoc.info.2 <- rbind(W2.2021b1.assoc.info.2,
#                                   W3.2021b1.assoc.info.2,
#                                   W4.2021b1.assoc.info.2)
# temp <- select(full.2021b1.assoc.info.2,id1:SocAssoc)
# d <- left_join(full.2021b1.assoc.info,temp,by=c("id1","id2"))
# d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc2=SocAssoc.y)
# full.2021b1.assoc.info <- d

# # For the 2021, Block 2, data in just the pre-trap period (60-minute overlap)
# W2.2021b2.assoc.info.1 <- collate.matrix(W2.2021b2.AM.1,NumDenom=TRUE)
# W2.2021b2.assoc.info.1 <- overlap.check(W2.2021b2.assoc.info.1,mice2021block2,
#                                         infectionCheck=TRUE)
# W3.2021b2.assoc.info.1 <- collate.matrix(W3.2021b2.AM.1,NumDenom=TRUE)
# W3.2021b2.assoc.info.1 <- overlap.check(W3.2021b2.assoc.info.1,mice2021block2,
#                                         infectionCheck=TRUE)
# W4.2021b2.assoc.info.1 <- collate.matrix(W4.2021b2.AM.1,NumDenom=TRUE)
# W4.2021b2.assoc.info.1 <- overlap.check(W4.2021b2.assoc.info.1,mice2021block2,
#                                         infectionCheck=TRUE)
# full.2021b2.assoc.info.1 <- rbind(W2.2021b2.assoc.info.1,
#                                   W3.2021b2.assoc.info.1,
#                                   W4.2021b2.assoc.info.1)
# temp <- select(full.2021b2.assoc.info.1,id1:SocAssoc)
# d <- left_join(full.2021b2.assoc.info,temp,by=c("id1","id2"))
# d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc1=SocAssoc.y)
# full.2021b2.assoc.info <- d

# # For the 2021, Block 1, data in just the post-trap period (60-minute overlap)
# W2.2021b2.assoc.info.2 <- collate.matrix(W2.2021b2.AM.2,NumDenom=TRUE)
# W2.2021b2.assoc.info.2 <- overlap.check(W2.2021b2.assoc.info.2,mice2021block2,
#                                         infectionCheck=TRUE)
# W3.2021b2.assoc.info.2 <- collate.matrix(W3.2021b2.AM.2,NumDenom=TRUE)
# W3.2021b2.assoc.info.2 <- overlap.check(W3.2021b2.assoc.info.2,mice2021block2,
#                                         infectionCheck=TRUE)
# W4.2021b2.assoc.info.2 <- collate.matrix(W4.2021b2.AM.2,NumDenom=TRUE)
# W4.2021b2.assoc.info.2 <- overlap.check(W4.2021b2.assoc.info.2,mice2021block2,
#                                         infectionCheck=TRUE)
# full.2021b2.assoc.info.2 <- rbind(W2.2021b2.assoc.info.2,
#                                   W3.2021b2.assoc.info.2,
#                                   W4.2021b2.assoc.info.2)
# temp <- select(full.2021b2.assoc.info.2,id1:SocAssoc)
# d <- left_join(full.2021b2.assoc.info,temp,by=c("id1","id2"))
# d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc2=SocAssoc.y)
# full.2021b2.assoc.info <- d

# For the 2021, Block 1, data with four-hour overlap windows
W2.2021b1.assoc.info.240 <- collate.matrix(W2.2021b1.AM.240,NumDenom=TRUE)
W2.2021b1.assoc.info.240 <- overlap.check(W2.2021b1.assoc.info.240,mice2021block1,
                                        infectionCheck=TRUE)
W3.2021b1.assoc.info.240 <- collate.matrix(W3.2021b1.AM.240,NumDenom=TRUE)
W3.2021b1.assoc.info.240 <- overlap.check(W3.2021b1.assoc.info.240,mice2021block1,
                                        infectionCheck=TRUE)
W4.2021b1.assoc.info.240 <- collate.matrix(W4.2021b1.AM.240,NumDenom=TRUE)
W4.2021b1.assoc.info.240 <- overlap.check(W4.2021b1.assoc.info.240,mice2021block1,
                                        infectionCheck=TRUE)
full.2021b1.assoc.info.240 <- rbind(W2.2021b1.assoc.info.240,
                                  W3.2021b1.assoc.info.240,
                                  W4.2021b1.assoc.info.240)
temp <- select(full.2021b1.assoc.info.240,id1:SocAssoc)
d <- left_join(full.2021b1.assoc.info,temp,by=c("id1","id2"))
d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc240=SocAssoc.y)
full.2021b1.assoc.info <- d

# For the 2021, Block 2, data with four-hour overlap windows
W2.2021b2.assoc.info.240 <- collate.matrix(W2.2021b2.AM.240,NumDenom=TRUE)
W2.2021b2.assoc.info.240 <- overlap.check(W2.2021b2.assoc.info.240,mice2021block2,
                                        infectionCheck=TRUE)
W3.2021b2.assoc.info.240 <- collate.matrix(W3.2021b2.AM.240,NumDenom=TRUE)
W3.2021b2.assoc.info.240 <- overlap.check(W3.2021b2.assoc.info.240,mice2021block2,
                                        infectionCheck=TRUE)
W4.2021b2.assoc.info.240 <- collate.matrix(W4.2021b2.AM.240,NumDenom=TRUE)
W4.2021b2.assoc.info.240 <- overlap.check(W4.2021b2.assoc.info.240,mice2021block2,
                                        infectionCheck=TRUE)
full.2021b2.assoc.info.240 <- rbind(W2.2021b2.assoc.info.240,
                                  W3.2021b2.assoc.info.240,
                                  W4.2021b2.assoc.info.240)
temp <- select(full.2021b2.assoc.info.240,id1:SocAssoc)
d <- left_join(full.2021b2.assoc.info,temp,by=c("id1","id2"))
d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc240=SocAssoc.y)
full.2021b2.assoc.info <- d

# For the 2021, Block 1, data with fifteen-minute overlap windows
W2.2021b1.assoc.info.15 <- collate.matrix(W2.2021b1.AM.15,NumDenom=TRUE)
W2.2021b1.assoc.info.15 <- overlap.check(W2.2021b1.assoc.info.15,mice2021block1,
                                        infectionCheck=TRUE)
W3.2021b1.assoc.info.15 <- collate.matrix(W3.2021b1.AM.15,NumDenom=TRUE)
W3.2021b1.assoc.info.15 <- overlap.check(W3.2021b1.assoc.info.15,mice2021block1,
                                        infectionCheck=TRUE)
W4.2021b1.assoc.info.15 <- collate.matrix(W4.2021b1.AM.15,NumDenom=TRUE)
W4.2021b1.assoc.info.15 <- overlap.check(W4.2021b1.assoc.info.15,mice2021block1,
                                        infectionCheck=TRUE)
full.2021b1.assoc.info.15 <- rbind(W2.2021b1.assoc.info.15,
                                  W3.2021b1.assoc.info.15,
                                  W4.2021b1.assoc.info.15)
temp <- select(full.2021b1.assoc.info.15,id1:Denom)
d <- left_join(full.2021b1.assoc.info,temp,by=c("id1","id2"))
d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc15=SocAssoc.y,
            Num=Num.x,Num15=Num.y,Denom=Denom.x,Denom15=Denom.y)
full.2021b1.assoc.info <- d

# For the 2021, Block 2, data with fifteen-minute overlap windows
W2.2021b2.assoc.info.15 <- collate.matrix(W2.2021b2.AM.15,NumDenom=TRUE)
W2.2021b2.assoc.info.15 <- overlap.check(W2.2021b2.assoc.info.15,mice2021block2,
                                        infectionCheck=TRUE)
W3.2021b2.assoc.info.15 <- collate.matrix(W3.2021b2.AM.15,NumDenom=TRUE)
W3.2021b2.assoc.info.15 <- overlap.check(W3.2021b2.assoc.info.15,mice2021block2,
                                        infectionCheck=TRUE)
W4.2021b2.assoc.info.15 <- collate.matrix(W4.2021b2.AM.15,NumDenom=TRUE)
W4.2021b2.assoc.info.15 <- overlap.check(W4.2021b2.assoc.info.15,mice2021block2,
                                        infectionCheck=TRUE)
full.2021b2.assoc.info.15 <- rbind(W2.2021b2.assoc.info.15,
                                  W3.2021b2.assoc.info.15,
                                  W4.2021b2.assoc.info.15)
temp <- select(full.2021b2.assoc.info.15,id1:Denom)
d <- left_join(full.2021b2.assoc.info,temp,by=c("id1","id2"))
d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc15=SocAssoc.y,
            Num=Num.x,Num15=Num.y,Denom=Denom.x,Denom15=Denom.y)
full.2021b2.assoc.info <- d

# For the 2021, Block 1, data with two-minute overlap windows
W2.2021b1.assoc.info.2min <- collate.matrix(W2.2021b1.AM.2min,NumDenom=TRUE)
W2.2021b1.assoc.info.2min <- overlap.check(W2.2021b1.assoc.info.2min,mice2021block1,
                                        infectionCheck=TRUE)
W3.2021b1.assoc.info.2min <- collate.matrix(W3.2021b1.AM.2min,NumDenom=TRUE)
W3.2021b1.assoc.info.2min <- overlap.check(W3.2021b1.assoc.info.2min,mice2021block1,
                                        infectionCheck=TRUE)
W4.2021b1.assoc.info.2min <- collate.matrix(W4.2021b1.AM.2min,NumDenom=TRUE)
W4.2021b1.assoc.info.2min <- overlap.check(W4.2021b1.assoc.info.2min,mice2021block1,
                                        infectionCheck=TRUE)
full.2021b1.assoc.info.2min <- rbind(W2.2021b1.assoc.info.2min,
                                  W3.2021b1.assoc.info.2min,
                                  W4.2021b1.assoc.info.2min)
temp <- select(full.2021b1.assoc.info.2min,id1:SocAssoc)
d <- left_join(full.2021b1.assoc.info,temp,by=c("id1","id2"))
d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc2min=SocAssoc.y)
full.2021b1.assoc.info <- d

# For the 2021, Block 2, data with two-minute overlap windows
W2.2021b2.assoc.info.2min <- collate.matrix(W2.2021b2.AM.2min,NumDenom=TRUE)
W2.2021b2.assoc.info.2min <- overlap.check(W2.2021b2.assoc.info.2min,mice2021block2,
                                        infectionCheck=TRUE)
W3.2021b2.assoc.info.2min <- collate.matrix(W3.2021b2.AM.2min,NumDenom=TRUE)
W3.2021b2.assoc.info.2min <- overlap.check(W3.2021b2.assoc.info.2min,mice2021block2,
                                        infectionCheck=TRUE)
W4.2021b2.assoc.info.2min <- collate.matrix(W4.2021b2.AM.2min,NumDenom=TRUE)
W4.2021b2.assoc.info.2min <- overlap.check(W4.2021b2.assoc.info.2min,mice2021block2,
                                        infectionCheck=TRUE)
full.2021b2.assoc.info.2min <- rbind(W2.2021b2.assoc.info.2min,
                                  W3.2021b2.assoc.info.2min,
                                  W4.2021b2.assoc.info.2min)
temp <- select(full.2021b2.assoc.info.2min,id1:SocAssoc)
d <- left_join(full.2021b2.assoc.info,temp,by=c("id1","id2"))
d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc2min=SocAssoc.y)
full.2021b2.assoc.info <- d

# For the 2021, Block 1, data with fifteen-minute overlap windows and just feeders
W2.2021b1.assoc.info.15.feeder <- collate.matrix(W2.2021b1.AM.15.feeder,NumDenom=TRUE)
W2.2021b1.assoc.info.15.feeder <- overlap.check(W2.2021b1.assoc.info.15.feeder,
                                                mice2021block1,infectionCheck=TRUE)
W3.2021b1.assoc.info.15.feeder <- collate.matrix(W3.2021b1.AM.15.feeder,NumDenom=TRUE)
W3.2021b1.assoc.info.15.feeder <- overlap.check(W3.2021b1.assoc.info.15.feeder,
                                                mice2021block1,infectionCheck=TRUE)
W4.2021b1.assoc.info.15.feeder <- collate.matrix(W4.2021b1.AM.15.feeder,NumDenom=TRUE)
W4.2021b1.assoc.info.15.feeder <- overlap.check(W4.2021b1.assoc.info.15.feeder,
                                                mice2021block1,infectionCheck=TRUE)
full.2021b1.assoc.info.15.feeder <- rbind(W2.2021b1.assoc.info.15.feeder,
                                   W3.2021b1.assoc.info.15.feeder,
                                   W4.2021b1.assoc.info.15.feeder)
temp <- select(full.2021b1.assoc.info.15.feeder,id1:SocAssoc)
d <- left_join(full.2021b1.assoc.info,temp,by=c("id1","id2"))
d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc15feeder=SocAssoc.y)
full.2021b1.assoc.info <- d

# For the 2021, Block 2, data with fifteen-minute overlap windows and just feeders
W2.2021b2.assoc.info.15.feeder <- collate.matrix(W2.2021b2.AM.15.feeder,NumDenom=TRUE)
W2.2021b2.assoc.info.15.feeder <- overlap.check(W2.2021b2.assoc.info.15.feeder,
                                                mice2021block2,infectionCheck=TRUE)
W3.2021b2.assoc.info.15.feeder <- collate.matrix(W3.2021b2.AM.15.feeder,NumDenom=TRUE)
W3.2021b2.assoc.info.15.feeder <- overlap.check(W3.2021b2.assoc.info.15.feeder,
                                                mice2021block2,infectionCheck=TRUE)
W4.2021b2.assoc.info.15.feeder <- collate.matrix(W4.2021b2.AM.15.feeder,NumDenom=TRUE)
W4.2021b2.assoc.info.15.feeder <- overlap.check(W4.2021b2.assoc.info.15.feeder,
                                                mice2021block2,infectionCheck=TRUE)
full.2021b2.assoc.info.15.feeder <- rbind(W2.2021b2.assoc.info.15.feeder,
                                   W3.2021b2.assoc.info.15.feeder,
                                   W4.2021b2.assoc.info.15.feeder)
temp <- select(full.2021b2.assoc.info.15.feeder,id1:SocAssoc)
d <- left_join(full.2021b2.assoc.info,temp,by=c("id1","id2"))
d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc15feeder=SocAssoc.y)
full.2021b2.assoc.info <- d

# For the 2021, Block 1, data with fifteen-minute overlap windows and non-feeders
W2.2021b1.assoc.info.15.noF <- collate.matrix(W2.2021b1.AM.15.noF,NumDenom=TRUE)
W2.2021b1.assoc.info.15.noF <- overlap.check(W2.2021b1.assoc.info.15.noF,mice2021block1,
                                         infectionCheck=TRUE)
W3.2021b1.assoc.info.15.noF <- collate.matrix(W3.2021b1.AM.15.noF,NumDenom=TRUE)
W3.2021b1.assoc.info.15.noF <- overlap.check(W3.2021b1.assoc.info.15.noF,mice2021block1,
                                         infectionCheck=TRUE)
W4.2021b1.assoc.info.15.noF <- collate.matrix(W4.2021b1.AM.15.noF,NumDenom=TRUE)
W4.2021b1.assoc.info.15.noF <- overlap.check(W4.2021b1.assoc.info.15.noF,mice2021block1,
                                         infectionCheck=TRUE)
full.2021b1.assoc.info.15.noF <- rbind(W2.2021b1.assoc.info.15.noF,
                                   W3.2021b1.assoc.info.15.noF,
                                   W4.2021b1.assoc.info.15.noF)
temp <- select(full.2021b1.assoc.info.15.noF,id1:SocAssoc)
d <- left_join(full.2021b1.assoc.info,temp,by=c("id1","id2"))
d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc15noF=SocAssoc.y)
full.2021b1.assoc.info <- d

# For the 2021, Block 2, data with fifteen-minute overlap windows and non-feeders
W2.2021b2.assoc.info.15.noF <- collate.matrix(W2.2021b2.AM.15.noF,NumDenom=TRUE)
W2.2021b2.assoc.info.15.noF <- overlap.check(W2.2021b2.assoc.info.15.noF,mice2021block2,
                                         infectionCheck=TRUE)
W3.2021b2.assoc.info.15.noF <- collate.matrix(W3.2021b2.AM.15.noF,NumDenom=TRUE)
W3.2021b2.assoc.info.15.noF <- overlap.check(W3.2021b2.assoc.info.15.noF,mice2021block2,
                                         infectionCheck=TRUE)
W4.2021b2.assoc.info.15.noF <- collate.matrix(W4.2021b2.AM.15.noF,NumDenom=TRUE)
W4.2021b2.assoc.info.15.noF <- overlap.check(W4.2021b2.assoc.info.15.noF,mice2021block2,
                                         infectionCheck=TRUE)
full.2021b2.assoc.info.15.noF <- rbind(W2.2021b2.assoc.info.15.noF,
                                   W3.2021b2.assoc.info.15.noF,
                                   W4.2021b2.assoc.info.15.noF)
temp <- select(full.2021b2.assoc.info.15.noF,id1:SocAssoc)
d <- left_join(full.2021b2.assoc.info,temp,by=c("id1","id2"))
d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc15noF=SocAssoc.y)
full.2021b2.assoc.info <- d

# For the 2021, Block 1, data in just the pre-trap period (15-minute overlap)
W2.2021b1.assoc.info.1.15 <- collate.matrix(W2.2021b1.AM.1.15,NumDenom=TRUE)
W2.2021b1.assoc.info.1.15 <- overlap.check(W2.2021b1.assoc.info.1.15,mice2021block1,
                                        infectionCheck=TRUE)
W3.2021b1.assoc.info.1.15 <- collate.matrix(W3.2021b1.AM.1.15,NumDenom=TRUE)
W3.2021b1.assoc.info.1.15 <- overlap.check(W3.2021b1.assoc.info.1.15,mice2021block1,
                                        infectionCheck=TRUE)
W4.2021b1.assoc.info.1.15 <- collate.matrix(W4.2021b1.AM.1.15,NumDenom=TRUE)
W4.2021b1.assoc.info.1.15 <- overlap.check(W4.2021b1.assoc.info.1.15,mice2021block1,
                                        infectionCheck=TRUE)
full.2021b1.assoc.info.1.15 <- rbind(W2.2021b1.assoc.info.1.15,
                                  W3.2021b1.assoc.info.1.15,
                                  W4.2021b1.assoc.info.1.15)
temp <- select(full.2021b1.assoc.info.1.15,id1:SocAssoc)
d <- left_join(full.2021b1.assoc.info,temp,by=c("id1","id2"))
d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc1.15=SocAssoc.y)
full.2021b1.assoc.info <- d

# For the 2021, Block 1, data in just the post-trap period (15-minute overlap)
W2.2021b1.assoc.info.2.15 <- collate.matrix(W2.2021b1.AM.2.15,NumDenom=TRUE)
W2.2021b1.assoc.info.2.15 <- overlap.check(W2.2021b1.assoc.info.2.15,mice2021block1,
                                        infectionCheck=TRUE)
W3.2021b1.assoc.info.2.15 <- collate.matrix(W3.2021b1.AM.2.15,NumDenom=TRUE)
W3.2021b1.assoc.info.2.15 <- overlap.check(W3.2021b1.assoc.info.2.15,mice2021block1,
                                        infectionCheck=TRUE)
W4.2021b1.assoc.info.2.15 <- collate.matrix(W4.2021b1.AM.2.15,NumDenom=TRUE)
W4.2021b1.assoc.info.2.15 <- overlap.check(W4.2021b1.assoc.info.2.15,mice2021block1,
                                        infectionCheck=TRUE)
full.2021b1.assoc.info.2.15 <- rbind(W2.2021b1.assoc.info.2.15,
                                  W3.2021b1.assoc.info.2.15,
                                  W4.2021b1.assoc.info.2.15)
temp <- select(full.2021b1.assoc.info.2.15,id1:SocAssoc)
d <- left_join(full.2021b1.assoc.info,temp,by=c("id1","id2"))
d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc2.15=SocAssoc.y)
full.2021b1.assoc.info <- d

# For the 2021, Block 2, data in just the pre-trap period for 15-minute windows
W2.2021b2.assoc.info.1.15 <- collate.matrix(W2.2021b2.AM.1.15,NumDenom=TRUE)
W2.2021b2.assoc.info.1.15 <- overlap.check(W2.2021b2.assoc.info.1.15,mice2021block2,
                                        infectionCheck=TRUE)
W3.2021b2.assoc.info.1.15 <- collate.matrix(W3.2021b2.AM.1.15,NumDenom=TRUE)
W3.2021b2.assoc.info.1.15 <- overlap.check(W3.2021b2.assoc.info.1.15,mice2021block2,
                                        infectionCheck=TRUE)
W4.2021b2.assoc.info.1.15 <- collate.matrix(W4.2021b2.AM.1.15,NumDenom=TRUE)
W4.2021b2.assoc.info.1.15 <- overlap.check(W4.2021b2.assoc.info.1.15,mice2021block2,
                                        infectionCheck=TRUE)
full.2021b2.assoc.info.1.15 <- rbind(W2.2021b2.assoc.info.1.15,
                                  W3.2021b2.assoc.info.1.15,
                                  W4.2021b2.assoc.info.1.15)
temp <- select(full.2021b2.assoc.info.1.15,id1:SocAssoc)
d <- left_join(full.2021b2.assoc.info,temp,by=c("id1","id2"))
d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc1.15=SocAssoc.y)
full.2021b2.assoc.info <- d

# For the 2021, Block 2, data in just the post-trap period for 15-minute windows
W2.2021b2.assoc.info.2.15 <- collate.matrix(W2.2021b2.AM.2.15,NumDenom=TRUE)
W2.2021b2.assoc.info.2.15 <- overlap.check(W2.2021b2.assoc.info.2.15,mice2021block2,
                                        infectionCheck=TRUE)
W3.2021b2.assoc.info.2.15 <- collate.matrix(W3.2021b2.AM.2.15,NumDenom=TRUE)
W3.2021b2.assoc.info.2.15 <- overlap.check(W3.2021b2.assoc.info.2.15,mice2021block2,
                                        infectionCheck=TRUE)
W4.2021b2.assoc.info.2.15 <- collate.matrix(W4.2021b2.AM.2.15,NumDenom=TRUE)
W4.2021b2.assoc.info.2.15 <- overlap.check(W4.2021b2.assoc.info.2.15,mice2021block2,
                                        infectionCheck=TRUE)
full.2021b2.assoc.info.2.15 <- rbind(W2.2021b2.assoc.info.2.15,
                                  W3.2021b2.assoc.info.2.15,
                                  W4.2021b2.assoc.info.2.15)
temp <- select(full.2021b2.assoc.info.2.15,id1:SocAssoc)
d <- left_join(full.2021b2.assoc.info,temp,by=c("id1","id2"))
d <- rename(d,SocAssoc=SocAssoc.x,SocAssoc2.15=SocAssoc.y)
full.2021b2.assoc.info <- d

full.2021b1.assoc.info$Block <- "Block 1"
full.2021b2.assoc.info$Block <- "Block 2"
temp <- rbind(full.2021b1.assoc.info,
              full.2021b2.assoc.info)
full.2021.assoc.info <- temp
# We want to include both info about how many of the mice are infected and
# whether they have the same or different infection status.  Below we produce
# the latter.
d <- full.2021.assoc.info
d$InfDiff <- NA
d$InfDiff <- ifelse(d$infection=="One","Different","Same")
full.2021.assoc.info <- d
# This data frame here contains association info for every pair of mice during
# the experiment.  It also has some metadata, like their strain, whether they
# shared a cage, their infection statuses, etc.
# It can be used for studying social association behavior and for models of
# the relationships between said behavior and immune parameters.


