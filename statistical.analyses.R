# This script produces the statistical analyses in the manuscript.
# All data is prepared from the scripts "individual.behavior.calculation.R,"
# "social.association.calculation.R," and "immune.similarity.calculation.R"
# In general, in this file we will not include all the model selection
# done to identify final models.
library(dplyr)
library(tidyr)
library(brms)

# A key issue in assessment of some behavioral traits is that some mice either
# lost RFID tags during the experiment or escaped their home enclosure during
# the experiment.  These mice are not useful for assessing predictors of mouse
# behavior.  But they can be used for assessing the relationship between actual
# behavior and immune similarity (see Methods for more).
# In addition, several mice were released and either never recaptured or were
# recaptured at the mindpoint but not at the end.
# We must create a vector with these mice so that they can be filtered out when
# appropriate.
nonCaptures <- c("401","418","430","457","458","462","463","464","108",
                 "109","315","322","278","473")
escapees <- c("411","421","423","457","464","469","471","474","302","305")
RFIDloss <- c("407","408","429","436","471","104","304","288")
behaviorSelection <- union(nonCaptures,union(escapees,RFIDloss))
# behaviorSelection is our general-purpose filtering vector.

# Our first models explore predictors of individual behavior.
# Check-ins per night, by mouse strain, block, enclosure, infection status, cage of
# origin, parentage.
# Possible fixed effects: strain, infection status
# Possible random effects: block, wedge, cage of origin, parentage
# Results in Table S???
d <- mice2021GxE
# Some cage numbers are repeated across blocks, so a unique identifier including
# both is needed to avoid issues in using cage of origin as a random effect.
d$CageUnique <- paste(d$Cage,d$Block,sep="-")
d <- filter(d,!(Mouse%in%behaviorSelection))
d <- d[which(is.infinite(log10(d$CIpN))==FALSE),] # Because it has 0 check-ins
CIpN.model <- brm(log10(CIpN) ~ Genotype + Infected,
             data=d,family="gaussian",iter=3000,warmup=1000,
             chains=1,cores=1,inits="random")

# Post-trap check-ins per night, by mouse strain, block, enclosure, infection status,
# cage of origin, parentage.
# Possible fixed effects: strain, infection status
# Possible random effects: block, wedge, cage of origin, parentage
# Results in Table S3
d <- mice2021GxE
d$CageUnique <- paste(d$Cage,d$Block,sep="-")
d <- filter(d,!(Mouse%in%behaviorSelection))
d <- d[which(is.infinite(log10(d$PTCIpN))==FALSE),] # Because it has 0 check-ins
PTCIpN.model <- brm(log10(PTCIpN) ~ Genotype + Infected,
             data=d,family="gaussian",iter=3000,warmup=1000,
             chains=1,cores=1,inits="random")

# Proportion of feeder check-ins, by mouse strain, block, enclosure, infection status,
# cage of origin, parentage.
# Possible fixed effects: strain, infection status
# Possible random effects: block, wedge, cage of origin, parentage
# Test: linear mixed model, logistic
# Results in Table S1
d <- mice2021GxE
d$CageUnique <- paste(d$Cage,d$Block,sep="-")
d <- d[which(d$Total!=0&!(d$Mouse%in%behaviorSelection)),]
d <- d[complete.cases(d$Parentage),]
feederProp.model <- brm(Feeder|trials(Total) ~ Infected +
                                (1|Parentage) + (1|CageUnique),
                   data=d,family="binomial",warmup=1000,iter=3000,
                   chains=1,cores=1,inits="random")

# Minimum distance traveled as a function of strain, infection status, block, wedge,
# cage of origin, parentage
# Possible fixed effects: strain, infection status
# Possible random effects: block, wedge, cage of origin, parentage
# Test: linear mixed model, normal distribution
# Results in Table S???
d <- mice2021GxE
d$CageUnique <- paste(d$Cage,d$Block,sep="-")
d$Distance <- d$Distance/d$Nights
d <- d[which(d$Distance!=0),]
d <- d[which(!(d$Mouse%in%behaviorSelection)),]
d <- d[complete.cases(d$Parentage),]
dist.model <- brm(log10(Distance) ~ Genotype +
                     (1|Block) + (1|Wedge) + (1|CageUnique) + (1|Parentage),
                   family="gaussian",data=d,warmup=1000,iter=4000,
                   cores=1,chains=1,inits="random")

# Frequency of check-ins per night at each location as a function of location,
# strain, infection status, block, wedge, cage of origin, mouse, parentage
# Possible fixed effects: location, strain, infection status, block
# Possible random effects: wedge, cage of origin, parentage, mouse
# Test: linear mixed model, normal distribution
# Results in Table S???
d <- select(mice2021GxE,Mouse,Wedge,Genotype,Infected,Block,
            Feeder,Tower,Base,Left,Right,Nights)
d <- filter(d,!(Mouse%in%behaviorSelection))
d[,6:10] <- d[,6:10]/d[,11]
d2 <- pivot_longer(d,cols=c("Feeder","Tower","Base","Left","Right"),
                   names_to="Location",values_to="CIpN")
d2 <- filter(d2,CIpN!=0)
Location.model <- brm(log10(CIpN)~Location+Genotype+Infected+Block+
                        (1|Mouse),
                      data=d2,family="gaussian",warmup=1000,iter=3000,
                      cores=1,chains=1,inits="random")

# Next we examine predictors of social interaction
# Dyadic association strength as a function of strain, infection status, sibship, wedge,
# block, cage-sharing
# Possible fixed effects: strain, infection status, sibship, cage-sharing
# Possible random effects: wedge, block, individual ID
# Test: linear mixed model, logistic
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$SocAssoc15),]
d2 <- d2[complete.cases(d2$fullSib),]
d2 <- filter(d2,!(id1%in%behaviorSelection|id2%in%behaviorSelection))
assoc.model <- brm(Num15|trials(Denom15) ~ Genotype + cageOverlap + infection + 
                    (1|mm(id1,id2)),
                  data=d2,family="binomial",warmup=1000,iter=5000,
                  chains=1,cores=1,inits="random")
# We can do the same with other overlap windows.
assoc.model.1hr <- brm(Num1hr|trials(Denom1hr) ~ Genotype + infection + cageOverlap + fullSib + 
                      (1|mm(id1,id2)),
                    data=d2,family="binomial",warmup=1000,iter=5000,
                    chains=1,cores=1,inits="random")
assoc.model.4hr <- brm(Num4hr|trials(Denom4hr) ~ Genotype + infection + cageOverlap + fullSib + 
                         (1|mm(id1,id2)),
                       data=d2,family="binomial",warmup=1000,iter=5000,
                       chains=1,cores=1,inits="random")
assoc.model.2min <- brm(Num2min|trials(Denom2min) ~ Genotype + infection + cageOverlap + fullSib + 
                         (1|mm(id1,id2)),
                       data=d2,family="binomial",warmup=1000,iter=5000,
                       chains=1,cores=1,inits="random")
# Pearson's r for correlations between various overlap window metrics
cor(d2$SocAssoc15,d2$SocAssoc,method="pearson")
cor(d2$SocAssoc15,d2$SocAssoc240,method="pearson")
cor(d2$SocAssoc15,d2$SocAssoc2min,method="pearson")
cor(d2$SocAssoc,d2$SocAssoc240,method="pearson")
cor(d2$SocAssoc,d2$SocAssoc2min,method="pearson")
cor(d2$SocAssoc240,d2$SocAssoc2min,method="pearson")
# And 15-minute overlaps in the post-challenge period
PTassoc.model <- brm(Num15.2|trials(Denom15.2) ~ Genotype + infection + 
                       + cageOverlap + fullSib +
                       (1|Block) + (1|Wedge) + (1|mm(id1,id2)),
                     data=d2,family="binomial",warmup=1000,iter=3000,
                     chains=1,cores=1,inits="random")


# Next we examine relationships between individual behaviors and mouse immune
# parameters.
# In general, we want to use absolute quantities of immune cells, rather than
# the relative proportions, since those proportions are a little more complex
# for us to model.  The absolute quantities have all been log-transformed.
d <- left_join(mice2021GxE,
               select(wbc2021Trap,Mouse,WBC_count:Bas_percent),
               by="Mouse")
d <- filter(d,!(Mouse%in%behaviorSelection))
CIs.bas.model <- brm(Bas_count~log10(CIpN)+Genotype+Block,
                     data=d,silent=2,refresh=0,chains=1)
CIs.eos.model <- brm(Eos_count~log10(CIpN)+Genotype+Block,
                     data=d,silent=2,refresh=0,chains=1)
CIs.mon.model <- brm(Mon_count~log10(CIpN)+Genotype+Block,
                     data=d,silent=2,refresh=0,chains=1)
CIs.neu.model <- brm(Neu_count~log10(CIpN)+Genotype+Block,
                     data=d,silent=2,refresh=0,chains=1)
CIs.lym.model <- brm(Lym_count~log10(CIpN)+Genotype+Block,
                     data=d,silent=2,refresh=0,chains=1)
# Here we examine distance per night as a predictor of the immune cell
# concentrations.
d$DistpN <- d$Distance/d$Nights
Dist.bas.model <- brms::brm(Bas_count~log10(DistpN)+Genotype+Block,
                           data=d,silent=2,refresh=0,chains=1)
Dist.eos.model <- brms::brm(Eos_count~log10(DistpN)+Genotype+Block,
                           data=d,silent=2,refresh=0,chains=1)
Dist.mon.model <- brms::brm(Mon_count~log10(DistpN)+Genotype+Block,
                           data=d,silent=2,refresh=0,chains=1)
Dist.neu.model <- brms::brm(Neu_count~log10(DistpN)+Genotype+Block,
                           data=d,silent=2,refresh=0,chains=1)
Dist.lym.model <- brms::brm(Lym_count~log10(DistpN)+Genotype+Block,
                           data=d,silent=2,refresh=0,chains=1)

# Next we examine relationships between social associations and immune
# similarities.  For these analyses we prefer to consider whether the mice have
# the same or different infection status, rather than what their statuses are.
# A key feature of these models is a multiple-membership mixed effect for
# individual mouse ID.  This is achieved through the brms mm() function, which
# produces the right random effects structure.
# We also want to do beta regressions, since the distribution of similarities
# is going to be best approximated by such a distribution, as opposed to using
# a normal distribution.

# First we analyze the relationship between T cell similarities and social
# associations.
# Both CD4 and CD8 T cells in the MLNs (Fig. 3A)
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$combTcellJI,d2$fullSib,d2$SocAssoc15),]
Tmem.MLN.model.2 <-  brm(combTcellJI ~ SocAssoc15 + Genotype + InfDiff + fullSib +
                           (1|mm(id1,id2)),
                         data=d2,family="Beta",warmup=1000,iter=3000,
                         cores=1,chains=1,inits="random")
# CD4 T cells in the MLNs (Fig. 3A, 3B)
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$SocAssoc15,d2$CD4JI,d2$fullSib),]
CD4.model.2 <- brm(CD4JI ~ SocAssoc15 + Genotype + InfDiff + fullSib +
                     (1|mm(id1,id2)),
                   data=d2,family="Beta",warmup=1000,iter=3000,
                   cores=1,chains=1,inits="random")
# CD4 T cells in the MLNs without feeders and only at feeders (Fig. 4B)
CD4.noF.model <- brm(CD4JI ~ SocAssoc15noF + Genotype + InfDiff + fullSib +
                        (1|mm(id1,id2)),
                      data=d2,family="Beta",warmup=1000,iter=3000,
                      cores=1,chains=1,inits="random")
CD4.feeder.model <- brm(CD4JI ~ SocAssoc15feeder + Genotype + InfDiff + fullSib +
                           (1|mm(id1,id2)),
                         data=d2,family="Beta",warmup=1000,iter=5000,
                         cores=1,chains=1,inits="random")
# CD4 T cells in the MLNs with various overlap windows (Fig. 4C)
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$CD4JI,d2$fullSib,d2$SocAssoc),]
CD4.1hr.model <- brm(CD4JI ~ SocAssoc + Genotype + InfDiff + fullSib +
                        (1|mm(id1,id2)),
                      data=d2,family="Beta",warmup=1000,iter=3000,
                      cores=1,chains=1,inits="random")
CD4.4hr.model <- brm(CD4JI ~ SocAssoc240 + Genotype + InfDiff + fullSib +
                        (1|mm(id1,id2)),
                      data=d2,family="Beta",warmup=1000,iter=3000,
                      cores=1,chains=1,inits="random")
CD4.2min.model <- brm(CD4JI ~ SocAssoc2min + Genotype + InfDiff + fullSib +
                         (1|mm(id1,id2)),
                       data=d2,family="Beta",warmup=1000,iter=3000,
                       cores=1,chains=1,inits="random")
# CD8 T cells in the MLNs (Fig. 3A, 3B)
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$SocAssoc15,d2$CD4JI,d2$fullSib),]
CD8.model <- brm(CD4JI ~ SocAssoc15 + Genotype + InfDiff + fullSib +
                     (1|mm(id1,id2)),
                   data=d2,family="Beta",warmup=1000,iter=3000,
                   cores=1,chains=1,inits="random")
# CD4 and CD8 T cells in the peripheral blood (Fig. 3A)
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$combTbloodJI,d2$fullSib,d2$SocAssoc15),]
BloodT.model.2 <-  brm(combTbloodJI ~ SocAssoc15 + Genotype + InfDiff + fullSib +
                           (1|mm(id1,id2)),
                         data=d2,family="Beta",warmup=1000,iter=3000,
                         cores=1,chains=1,inits="random")
# B cells in the MLNs (Fig. 3A)
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$SocAssoc15,d2$BcellJI,d2$fullSib),]
Bcell.model.2 <- brm(BcellJI ~ SocAssoc15 + Genotype + InfDiff + fullSib +
               (1|mm(id1,id2)),
             data=d2,family="beta",warmup=1000,iter=3000,
             cores=1,chains=1,inits="random")

# CBC results from peripheral blood (Fig. 3A)
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$wbcJI,d2$SocAssoc15,d2$fullSib),]
CBC.JI.model.2 <- brm(wbcJI ~ SocAssoc15 + Genotype + InfDiff + fullSib +
                        (1|mm(id1,id2)),
                      data=d2,family="Beta",warmup=1000,iter=3000,
                      cores=1,chains=1,inits="random")
# CBC results for association at only feeders and only non-feeders (Fig. S4A)
CBC.noF.model <- brm(wbcJI ~ SocAssoc15noF + Genotype + InfDiff + fullSib +
                       (1|mm(id1,id2)),
                     data=d2,family="Beta",warmup=1000,iter=3000,
                     cores=1,chains=1,inits="random")
CBC.feeder.model <- brm(wbcJI ~ SocAssoc15feeder + Genotype + InfDiff + fullSib +
                          (1|mm(id1,id2)),
                        data=d2,family="Beta",warmup=1000,iter=3000,
                        cores=1,chains=1,inits="random")
# CBC results for association across different overlap windows (Fig. S4B)
CBC.1hr <- brm(wbcJI ~ SocAssoc + Genotype + InfDiff + fullSib +
                 (1|mm(id1,id2)),
               data=d2,family="Beta",warmup=1000,iter=3000,
               cores=1,chains=1,inits="random")
CBC.4hr <- brm(wbcJI ~ SocAssoc240 + Genotype + InfDiff + fullSib +
                 (1|mm(id1,id2)),
               data=d2,family="Beta",warmup=1000,iter=3000,
               cores=1,chains=1,inits="random")
CBC.2min <- brm(wbcJI ~ SocAssoc2min + Genotype + InfDiff + fullSib +
                  (1|mm(id1,id2)),
                data=d2,family="Beta",warmup=1000,iter=3000,
                cores=1,chains=1,inits="random")
# CBC results for different sample timings (Fig. 4A)
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$SocAssoc15,d2$fullSib),]
d3 <- d2[complete.cases(d2$wbcStartJ),]
CBC.start.model <- brm(wbcStartJ ~ SocAssoc15 + Genotype + InfDiff + fullSib +
                         (1|mm(id1,id2)),
                       data=d3,family="beta",iter=3000,warmup=1000,
                       cores=1,chains=1,inits="random")
d3 <- d2[complete.cases(d2$wbcMidJ),]
CBC.mid.model <- brm(wbcMidJ ~ SocAssoc15 + Genotype + InfDiff + fullSib +
                       (1|mm(id1,id2)),
                     data=d3,family="beta",iter=3000,warmup=1000,
                     cores=1,chains=1,inits="random")

# We also want to assess the influence of the microbiome for each of these
# different aspects of immune similarity, as well as how microbiome similarity
# correlates with social association.
# First we check the relationship between microbiome and association.
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$SocAssoc15,d2$mbJI,d2$fullSib),]
Micro.model <- brm(mbJI ~ SocAssoc15 + Genotype + InfDiff + fullSib +
                     (1|mm(id1,id2)),
                   data=d2,family="beta",warmup=1000,iter=3000,
                   cores=1,chains=1,inits="random")
# Next we consider how microbiome fits with a few aspects of T cell similarity
# in the MLNs (Fig. 4D for the CD4 results).
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$SocAssoc15,d2$combTcellJI,d2$mbJI,d2$fullSib),]
Tmem.micro.model <- brm(combTcellJI ~ SocAssoc15 + Genotype + InfDiff + fullSib + mbJI +
                          (1|mm(id1,id2)),
                        data=d2,family="beta",warmup=1000,iter=3000,
                        cores=1,chains=1,inits="random")
CD4.micro.model <- brm(CD4JI ~ SocAssoc15 + Genotype + InfDiff + fullSib + mbJI +
                         (1|mm(id1,id2)),
                       data=d2,family="beta",warmup=1000,iter=3000,
                       cores=1,chains=1,inits="random")
# Then we consider how microbiome fits with CBC similarity (Fig. S4C).
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$SocAssoc15,d2$wbcJI,d2$mbJI,d2$fullSib),]
CBC.micro.model <- brm(wbcJI ~ SocAssoc15 + Genotype + InfDiff + fullSib + mbJI +
                         (1|mm(id1,id2)),
                       data=d2,family="beta",warmup=1000,iter=3000,
                       cores=1,chains=1,inits="random")
# And lastly we consider how microbiome fits with B cell similarity.
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$SocAssoc15,d2$BcellJI,d2$mbJI,d2$fullSib),]
Bcell.micro.model <- brm(BcellJI ~ SocAssoc15 + Genotype + InfDiff + fullSib + mbJI +
                           (1|mm(id1,id2)),
                         data=d2,family="beta",warmup=1000,iter=3000,
                         cores=1,chains=1,inits="random")



