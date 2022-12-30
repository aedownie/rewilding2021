# This script collects the immune data into similarity measures.
# We do this for several aspects of the data, including CBC data, flow
# cytometry from the blood and MLNs, and for plasma cytokine levels and
# expression of cytokines in response to challenge by cells in the MLNs.
# We also calculate microbiome similarity.

# These metrics can be easily integrated with the table containing the social
# association metrics, from social.association.calculation.R

library(dplyr)
library(vegan)

##' A function that merges a pairwise association index of some sort into a data frame
##' for those individuals with two columns (id1 and id2) for the two members of the
##' dyad.  The pairwise association index must be in a data frame with the two ID
##' columns and one column entitled AI containing the association index of interest.
##' This function handles when the two names of a pair are id1 and id2 in one data
##' frame and id2 and id1 in the other.
##'
##' @param assocInfo A data frame with the association info between individuals
##'                  that you want to attach your new data to
##' @param data A data frame containing the new association info

merge.AI <- function(assocInfo,data){
  assocInfo$NewAI <- NA
  for (i in 1:nrow(assocInfo)){
    id1 <- assocInfo$id1[i]
    id2 <- assocInfo$id2[i]
    index <- which((data$id1==id1|data$id2==id1)&
                     (data$id1==id2|data$id2==id2))
    ifelse(length(index)==0,assocInfo$NewAI[i] <- NA,
           assocInfo$NewAI[i] <- data$AI[index])
  }
print("Don't forget to rename your NewAI column to something more specific.")
return(assocInfo)
}

# First we want to read in the flow cytometry data.  This is from the MLNs.
# We chose to work from relative abundances in different compartments.
# This is because the quantities in some compartments are nested within other
# compartments, which makes the full set of relationships complex and possibly
# a problem when using something like Euclidean distance.  Just focusing on
# relative abundance of different CD4 memory T cell phenotypes, for example,
# avoids confounds from different abundances of CD4 T cells overall.
fullMLN2021 <- read.csv("~/Block1_Block2.csv")
# A few IDs need to be corrected
fullMLN2021$Mouse.ID[which(fullMLN2021$Mouse.ID=="471/666")] <- "471"
fullMLN2021$Mouse.ID[which(fullMLN2021$Mouse.ID=="7010")] <- "436"
# Let's break it out by T cells and B cells
newTcelltemp <- fullMLN2021[,c(1,7:11,18:22)]
Bcelltemp <- fullMLN2021[,c(1,25:29)]

# We start with the B cell phenotypes.
Bcelltemp <- rename(Bcelltemp,Mouse=Mouse.ID,
                    B220percent=X..B220..cells.of.CD45.,
                    CD44hi=X..CD44..B.cells.of.CD45.,
                    CD62Lhi=X..CD62L..Bcells.of.CD45.,
                    CD62LhiCD44hi=X..CD62hiCD44..B.cells.of.CD45.,
                    CD62LhiCD44low=X..CD62hiCD44..B.cells.of.CD45..1)
# This doesn't have a double-negative column, or a CD62LlowCD44hi, so we get
# those manually.
Bcelltemp$CD62LlowCD44hi <- Bcelltemp$CD44hi-Bcelltemp$CD62LhiCD44hi
Bcelltemp$CD62LlowCD44low <- Bcelltemp$B220percent-(Bcelltemp$CD44hi+
                                                      Bcelltemp$CD62Lhi-
                                                      Bcelltemp$CD62LhiCD44hi)
# There's a little bit of rounding error, so for figuring out compartment percentages
# we want to use the sum.
BcellProp <- select(Bcelltemp,Mouse,CD62LhiCD44hi:CD62LlowCD44low)
for(i in 1:nrow(BcellProp)){
  total <- sum(BcellProp[i,2:5])
  for(j in 2:5){
    BcellProp[i,j] <- BcellProp[i,j]/total
  }
}

# We need to calculate the variation between each pair of mice with Jaccard
# index.
# NOTE: JI sometimes is rendered with 1 indicating perfect DISsimilarity and
# 0 indicating perfect similarity.  We use the reverse (see Raulo et al. 2021)
# for easier interpretation
Jdist <- as.matrix(vegdist(BcellProp[,2:5],method="jaccard"))
rownames(Jdist) <- colnames(Jdist) <- BcellProp$Mouse
BcelltempAI <- collate.matrix(Jdist,NumDenom=FALSE)
BcelltempAI <- rename(BcelltempAI,AI=SocAssoc)
# The matrix contains every single possible pair, but we only want pairs of
# mice that may have actually interacted through being in the same enclosure at
# the same time.  So we use merge.AI with the full.2021.assoc.info data frame
# to only select the right pairs.
# merge.AI can also handle variable location of ID numbers in columns. 
d <- merge.AI(full.2021.assoc.info,BcelltempAI)
d <- rename(d,BcellJI=NewAI)
d$BcellJI <- 1-d$BcellJI # Flipping the number as above
full.2021.assoc.info <- d
# And now we have immune similarity data, as we wanted.  From here we can
# conduct various analyses for factors influencing that similarity.
# We also have several other aspects of immune phenotype that we are interested
# in examining.

# Next we work with the MLN T cells
# In our first pass, we include both CD4 and CD8 T cells, scaled by their
# relative abundances.
colnames(newTcelltemp) <- c("Mouse","CD4","CD4effmem",
                            "CD4centmem","CD4naive","CD4doubleneg",
                            "CD8","CD8effmem","CD8centmem",
                            "CD8naive","CD8doubleneg")
for (i in 1:nrow(newTcelltemp)){
  totalCD4 <- sum(newTcelltemp[i,3:6])
  CD4prop <- newTcelltemp[i,2]/(newTcelltemp[i,2]+newTcelltemp[i,7])
  totalCD8 <- sum(newTcelltemp[i,8:11])
  CD8prop <- newTcelltemp[i,7]/(newTcelltemp[i,2]+newTcelltemp[i,7])
  for(j in 3:6){
    newTcelltemp[i,j] <- (newTcelltemp[i,j]/totalCD4)*CD4prop
  }
  for(j in 8:11){
    newTcelltemp[i,j] <- (newTcelltemp[i,j]/totalCD8)*CD8prop
  }
}
Jdist <- as.matrix(vegdist(newTcelltemp[,c(3:6,8:11)],method="jaccard"))
rownames(Jdist) <- colnames(Jdist) <- newTcelltemp$Mouse
TcelltempAI <- collate.matrix(Jdist,NumDenom=FALSE)
TcelltempAI <- rename(TcelltempAI,AI=SocAssoc)
d <- merge.AI(full.2021.assoc.info,TcelltempAI)
d <- rename(d,combTcellJI=NewAI)
d$combTcellJI <- 1-d$combTcellJI
full.2021.assoc.info <- d

# We can also do only the CD8 T cells
Jdist <- as.matrix(vegdist(newTcelltemp[,8:11],method="jaccard"))
rownames(Jdist) <- colnames(Jdist) <- newTcelltemp$Mouse
TcelltempAI <- collate.matrix(Jdist,NumDenom=FALSE)
TcelltempAI <- rename(TcelltempAI,AI=SocAssoc)
d <- merge.AI(full.2021.assoc.info,TcelltempAI)
d <- rename(d,CD8JI=NewAI)
d$CD8JI <- 1-d$CD8JI
full.2021.assoc.info <- d

# And only the CD4 T cells
Jdist <- as.matrix(vegdist(newTcelltemp[,3:6],method="jaccard"))
rownames(Jdist) <- colnames(Jdist) <- newTcelltemp$Mouse
TcelltempAI <- collate.matrix(Jdist,NumDenom=FALSE)
TcelltempAI <- rename(TcelltempAI,AI=SocAssoc)
d <- merge.AI(full.2021.assoc.info,TcelltempAI)
d <- rename(d,CD4JI=NewAI)
d$CD4JI <- 1-d$CD4JI
full.2021.assoc.info <- d

# We also want to assess combined CD4/CD8 T cell similarity in the whole blood
flowBlood2021 <- read.csv("~/2021.blood.flow.analysis.csv")
# First we need to ensure the metadata is clear
flowBlood2021$Genotype[which(flowBlood2021$Genotype=="C57/B6")] <- "C57"
flowBlood2021$Genotype[which(flowBlood2021$Genotype=="B6")] <- "C57"
flowBlood2021$Genotype[which(flowBlood2021$Genotype=="129SL")] <- "129"
flowBlood2021$Infection_Status[which(flowBlood2021$Infection_Status=="Tm")] <- "Y"
flowBlood2021$Infection_Status[which(flowBlood2021$Infection_Status=="No")] <- "N"
# Then we can select our columns of interest
Btemp <- filter(flowBlood2021,Location=="SF")%>%
  select(Sample.ID,Percent.CD4.Tcells,Percent.CD8.Tcells,
         Percent.CD4.effector.memory:Percent.CD8.double.negative)%>%
  rename(Mouse=Sample.ID)
Btemp$Mouse <- as.character(Btemp$Mouse)
# Here we want to calculate the right relative abundances, as with the MLNs.
for (i in 1:nrow(Btemp)){
  totalCD4 <- sum(Btemp[i,4:7])
  CD4prop <- Btemp[i,2]/(Btemp[i,2]+Btemp[i,3])
  totalCD8 <- sum(Btemp[i,8:11])
  CD8prop <- Btemp[i,3]/(Btemp[i,2]+Btemp[i,3])
  for(j in 4:7){
    Btemp[i,j] <- (Btemp[i,j]/totalCD4)*CD4prop
  }
  for(j in 8:11){
    Btemp[i,j] <- (Btemp[i,j]/totalCD8)*CD8prop
  }
}
# And then we can go through the standard process of merging the data.
tempDist <- as.matrix(vegdist(Btemp[,4:11],method="jaccard"))
rownames(tempDist) <- colnames(tempDist) <- Btemp$Mouse
BtempAI <- collate.matrix(tempDist,NumDenom=FALSE)
BtempAI <- rename(BtempAI,AI=SocAssoc)
d <- merge.AI(full.2021.assoc.info,BtempAI)
d <- rename(d,combTbloodJI=NewAI)
d$combTbloodJI <- 1-d$combTbloodJI
full.2021.assoc.info <- d

# The next slice of data for analysis is the CBC data.  This is in a different
# file, which contains CBC results for the three timepoints.  The results are
# both the absolute cell abundances and the relative abundances.
wbc2021 <- read.csv("~/2021.SF.CBC.data.relabeled.csv")
wbc2021 <- select(wbc2021,Unique_Sample_ID:RBC_count)
wbc2021 <- rename(wbc2021,Mouse=Princeton.NIH_Tag,Timepoint=Time_place_of_CBC_analysis)
wbc2021 <- filter(wbc2021,Location!="")
wbc2021 <- filter(wbc2021,Mouse%in%unique(mice2021total$Mouse))

# All cell count data is in 10^3 cells per µL, except for RBCs, which are 10^6/µL
# Accordingly we can calculate the logarithm and use that as our quantity.
# One issue that crops up is that basophil quantities are quite low, and some have 0.00
# for their value.  This then returns "-Inf" when I try to apply log10() to it, and
# that leads to NaNs after applying scale() to the column.  To address this, we 
# substitute a chosen value for 0.00, like 0.005, that can be evaluated with log10().
# This value is just below the detection threshold for our machine.  This seems
# unsatisfactory, but it retains information about low quantities.  This makes
# it our best choice.
for (i in 12:17){
  wbc2021[,i] <- log10(wbc2021[,i])+3
}
for(i in 1:nrow(wbc2021)){
  if(is.infinite(wbc2021[i,17])){wbc2021[i,17] <- log10(0.005)+3}
}
wbc2021[,23] <- log10(wbc2021[,23])+6
# We have three timepoints, and we don't want to scale distributions using data
# across the three timepoints.  We therefore need to break it down.
wbc2021NIH <- filter(wbc2021,Timepoint=="NIH")
for (i in c(12:17,23)){
  wbc2021NIH[,i] <- scale(wbc2021NIH[,i])
}
wbc2021Mid <- filter(wbc2021,Timepoint=="2 weeks")
for (i in c(12:17,23)){
  wbc2021Mid[,i] <- scale(wbc2021Mid[,i])
}
wbc2021Trap <- filter(wbc2021,Timepoint=="5 weeks")
for (i in c(12:17,23)){
  wbc2021Trap[,i] <- scale(wbc2021Trap[,i])
}
wbc2021scaled <- rbind(wbc2021NIH,wbc2021Mid,wbc2021Trap)
# This is useful for analysis of the relationships between individual behavior
# and cell concentrations.  But for similarity we again work with the
# proportion data.
wbc2021AItrapProp <- as.matrix(vegdist(filter(wbc2021,Timepoint=="5 weeks")[,18:22],
                                       method="jaccard"))
wbc2021AIstartProp <- as.matrix(vegdist(filter(wbc2021scaled,
                                               Timepoint=="NIH")[,18:22],
                                        method="jaccard"))
wbc2021AImidProp <- as.matrix(vegdist(filter(wbc2021scaled,
                                                    Timepoint=="2 weeks")[,18:22],
                                             method="jaccard"))

# Now we can place the similarity scores in the big data frame.
wbc2021AI <- collate.matrix(wbc2021AItrapProp,NumDenom=FALSE)
wbc2021AI <- rename(wbc2021AI,AI=SocAssoc)
d <- merge.AI(full.2021.assoc.info,wbc2021AI)
d <- rename(d,wbcPropJI=NewAI)
d$wbcPropJI <- 1-d$wbcPropJI

wbc2021AI <- collate.matrix(wbc2021AIstartProp,NumDenom=FALSE)
wbc2021AI <- rename(wbc2021AI,AI=SocAssoc)
d <- merge.AI(d,wbc2021AI)
d <- rename(d,wbcStartJ=NewAI)
d$wbcStartJ <- 1-d$wbcStartJ

wbc2021AI <- collate.matrix(wbc2021AImidProp,NumDenom=FALSE)
wbc2021AI <- rename(wbc2021AI,AI=SocAssoc)
d <- merge.AI(d,wbc2021AI)
d <- rename(d,wbcMidJ=NewAI)
d$wbcMidJ <- 1-d$wbcMidJ
full.2021.assoc.info <- d

# The last info we need to gather is the cytokine information.  This is in two
# forms: the plasma cytokine levels and the expression of cytokines by the MLNs
# in response to stimulus with various antigens.
# These two quantities are best analyzed with Manhattan distance

# First we do plasma cytokine levels
plasmaCyt <- read.csv("~/FINAL.csv")
plasmaCyt <- rename(plasmaCyt,Mouse=PrincetonID,
                    IFN=IFN....A4.,IL5=IL.5..A5.,TNF=TNF....A6.,
                    IL6=IL.6..A8.,IL17A=IL.17A..B4.,IL22=IL.22..B7.)
# We need to transform the plasma cytokine levels, because we expect them to
# follow a lognormal distribution.
for(i in 6:11){
  plasmaCyt[,i] <- scale(log10(plasmaCyt[,i]))
}
temp <- filter(plasmaCyt,Mouse%in%mice2021GxE$Mouse)
# After this the process is more or less the same as it is for cells.
plasmaCytDist <- as.matrix(dist(temp[,6:11],method="manhattan"))
rownames(plasmaCytDist) <- colnames(plasmaCytDist) <- temp$Mouse
fPCytAI <- collate.matrix(plasmaCytDist,NumDenom=FALSE)
fPCytAI <- rename(fPCytAI,AI=SocAssoc)
d <- merge.AI(full.2021.assoc.info,fPCytAI)
d <- rename(d,PCytDist=NewAI)
full.2021.assoc.info <- d

# The MLN cytokine expression is more complex.  Here we want to look at the
# expression relative to a control treatment (here exposure to PBS).  We
# normalize expression as log2((treatment/control)+1).  This is the same method
# as employed in Lin et al. (2020) and Yeung et al. (2020).
cytProd <- read.csv("Final final cytokines.csv")
temp <- cytProd
# Regrettably the normalization has to be done for each cytokine manually.
temp[,6:11] <- temp[,6:11]/temp[,5] # IFNg
temp[,13:18] <- temp[,13:18]/temp[,12] # IL5
temp[,20:25] <- temp[,20:25]/temp[,19] # TNFa
temp[,27:32] <- temp[,27:32]/temp[,26] # IL2
temp[,34:39] <- temp[,34:39]/temp[,33] # IL6
temp[,41:46] <- temp[,41:46]/temp[,40] # IL4
temp[,48:53] <- temp[,48:53]/temp[,47] # IL10
temp[,55:60] <- temp[,55:60]/temp[,54] # IL9
temp[,62:67] <- temp[,62:67]/temp[,61] # IL17A
temp[,69:74] <- temp[,69:74]/temp[,68] # IL22
temp[,76:81] <- temp[,76:81]/temp[,75] # IL13
# We then want to remove the control for each cytokine.
temp <- temp[,c(1,6:11,13:18,20:25,27:32,34:39,41:46,48:53,55:60,62:67,69:74,76:81)]
temp <- rename(temp,Mouse=X)
temp[,2:67] <- log2(temp[,2:67]+1)
# Normalization is now complete.
# Some of the mouse IDs have to be corrected.
temp$Mouse[which(temp$Mouse=="277/270")] <- "277"
temp$Mouse[which(temp$Mouse=="112/125/117")] <- "117"
temp$Mouse[which(temp$Mouse=="122/168")] <- "122"
temp$Mouse[which(temp$Mouse=="103/124")] <- "103"
temp$Mouse[which(temp$Mouse=="325/363")] <- "325"
temp$Mouse[which(temp$Mouse=="304/374")] <- "304"
temp$Mouse[which(temp$Mouse=="302/362")] <- "302"
temp$Mouse[which(temp$Mouse=="317/375")] <- "317"
# From here the process can proceed as usual.
cytProdDist <- as.matrix(dist(temp[,2:67],method="manhattan"))
colnames(cytProdDist) <- rownames(cytProdDist) <- temp$Mouse
cytProdAI <- collate.matrix(cytProdDist,NumDenom=FALSE)
cytProdAI <- rename(cytProdAI,AI=SocAssoc)
d <- merge.AI(full.2021.assoc.info,cytProdAI)
d <- rename(d,cytProdDist=NewAI)
full.2021.assoc.info <- d

# The last thing that needs to be added to the dataset is the microbiome data.
# This can be calculated in a relatively straightforward manner.
microbiome2021 <- read.delim("ASV_table_transposed.txt")
# First, some filtering to remove irrelevant mice (not SF) and the rows that
# have the taxonomic information.
microbiome2021 <- rbind(microbiome2021[1:7,],
                        filter(microbiome2021,is.na(Mouse)==FALSE))
taxonInfo <- microbiome2021[1:7,4:ncol(microbiome2021)]
mbData <- microbiome2021[8:nrow(microbiome2021),]
for(i in 4:ncol(mbData)){
  mbData[,i] <- as.numeric(mbData[,i])
} # Convert to numeric to make for easier work.
# Because we did our sequencing with both rewilded and laboratory mice, there
# are some taxa that don't show up in any of our rewilded mice.  We therefore
# remove those columns from the data set as below.
mbData <- mbData[,c(1:3,which(colSums(mbData[4:ncol(mbData)])!=0)+3)]
taxonInfo <- taxonInfo[,which(colnames(taxonInfo)%in%colnames(mbData))]

taxa <- rep(NA,ncol(taxonInfo))
greatestDepth <- rep(0,ncol(taxonInfo))
# We want to convert each column to a sensible-ish label
for(i in 1:ncol(taxonInfo)){
  depth <- 0
  for(j in 1:7){
    if(is.na(taxonInfo[j,i])==TRUE){
      depth <- j-1
      break
    }
  }
  taxon <- paste(c(colnames(taxonInfo)[i],
                            taxonInfo[1:depth,i]),collapse="_")
  taxa[i] <- taxon
  greatestDepth[i] <- depth
}

d <- filter(mbData,Timepoint=="5 weeks")
d <- filter(d,Mouse%in%mice2021GxE$Mouse) # 68 mice (fewer than I might have thought?)
colnames(d)[4:ncol(d)] <- taxa
d <- select(d,-Row.names,-Timepoint)
d <- d[,which(colSums(d[2:ncol(d)])!=0)] # More filtering
# We need to ensure that we are considering relative abundances, rather than
# absolute counts, as the latter will cause problems for our calculation of the
# similarities of pairs with Jaccard index.
for(i in 1:nrow(d)){
  d[i,2:ncol(d)] <- d[i,2:ncol(d)]/sum(d[i,2:ncol(d)])
} # Takes about 30 seconds to scale for relative abundance
RWmicrobiome2021 <- as.data.frame(d)
for(i in 2:ncol(RWmicrobiome2021)){
  RWmicrobiome2021[,i] <- as.numeric(RWmicrobiome2021[,i])
}
# From here it proceeds fairly smoothly.
temp <- filter(RWmicrobiome2021,Mouse!="446")
Jdist <- as.matrix(vegan::vegdist(temp[,2:ncol(temp)]),
                                  method="jaccard")
rownames(Jdist) <- colnames(Jdist) <- temp$Mouse
mbtempAI <- collate.matrix(Jdist,NumDenom=FALSE)
mbtempAI <- rename(mbtempAI,AI=SocAssoc)
d <- merge.AI(full.2021.assoc.info,mbtempAI) # Got some warnings, need to chase down
# May occur where we have multiple samples for a mouse?
# That seems to be correct.  Need to double-check the 446 results.
d <- rename(d,mbJI=NewAI)
d$mbJI <- 1-d$mbJI
full.2021.assoc.info <- d


