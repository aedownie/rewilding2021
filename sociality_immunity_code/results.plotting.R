# This script contains methods for plotting the model results and other figures
# in the manuscript.  To do this, we need the outputs from previous scripts:
# individual.behavior.calculation.R, social.association.calculation.R, and
# statistical.analyses.R
library(dplyr)
library(tidyr)
library(ggplot2)
library(MetBrewer)
library(ggridges)

# We start with the plots showing behavior

# The plots for Fig. 1 (except 1A)

# Fig. 1B
d <- select(mice2021total,Mouse,Total,Nights,Genotype,Infected)
ggplot(filter(d,!(Mouse%in%behaviorSelection)),
       aes(Total/Nights,Distance/Nights,color=Genotype))+
  geom_point(size=1.8)+
  scale_color_manual(values=c("#FFA533","#FF3366","#006699"),
                     labels=c("129S1","C57BL/6","PWK/PhJ"),
                     name="Mouse strain")+
  theme_classic()+theme(axis.line=element_line(lineend="round",linewidth=0.75),
                        legend.position=c(0.75,0.75),
                        legend.background=element_rect(color="grey10"),
                        legend.text=element_text(size=10),
                        legend.title=element_text(hjust=0.5),
                        axis.text=element_text(size=11),
                        axis.title=element_text(size=14),
                        axis.ticks=element_line(linewidth = 0.75))+
  labs(x="Check-ins per night during experiment",
       y="Minimum distance traveled per night (m)")

# Fig. 1C
ggplot(filter(mice2021total,!(Mouse%in%behaviorSelection)),
       aes(Genotype,Feeder/Total,color=Genotype))+
  geom_boxplot(outlier.shape=NA)+geom_jitter(size=1.8)+
  theme_classic()+theme(axis.line=element_line(lineend="round",linewidth = 0.75),
                        axis.text=element_text(size=11),
                        axis.title=element_text(size=14),
                        axis.ticks=element_line(linewidth = 0.75))+
  labs(x="Mouse strain",
       y="Proportion of check-ins at feeder")+
  scale_color_manual(values=c("#FFA533","#FF3366","#006699"))+
  scale_x_discrete(labels=c("129S1","C57BL/6","PWK/PhJ"))+
  guides(color="none")

# Fig. 1D
ggplot(filter(mice2021total,!(Mouse%in%behaviorSelection)),
       aes(Genotype,meanRE,color=Genotype))+
  geom_boxplot(outlier.shape=NA)+geom_jitter(size=2.0)+
  theme_classic()+
  theme(axis.line=element_line(lineend="round",linewidth = 0.75),
        axis.text=element_text(size=11),
        axis.title=element_text(size=14),
        axis.ticks=element_line(linewidth = 0.75))+
  labs(x="Mouse strain",
       y="Mean nightly roaming entropy")+
  scale_color_manual(values=c("#FFA533","#FF3366","#006699"))+
  scale_x_discrete(labels=c("129S1","C57BL/6","PWK/PhJ"))+
  guides(color="none")

# Now the plots for Fig. 2.  First, some reconfiguring of the data.
d <- full.2021.assoc.info
d$pair <- NA
for(i in 1:nrow(d)){
  d$pair[i] <- paste(d$id1[i],d$id2[i],sep="-")
}
d <- select(filter(d,Wedge!="Wedge 6"&Wedge!="Wedge 7"&
                     (!(id1%in%behaviorSelection)&!(id2%in%behaviorSelection))),
            SocAssoc15feeder,SocAssoc15,SocAssoc15noF,
            pair,Wedge,Block,Genotype,infection)
d <- as.data.frame(pivot_longer(d,1:3,names_to="Location"))
# For Fig. 2A
ggplot(filter(d,Location=="SocAssoc15"),aes(value,fill=Genotype))+
  geom_histogram()+facet_wrap(~Genotype,
                              labeller=as_labeller(c('129'="129S1",
                                                     'C57'="C57BL/6",
                                                     'PWK'="PWK/PhJ")))+
  labs(x="Spatiotemporal association of pair",y="Number of pairs")+
  scale_fill_manual(values=c("#FFA533","#FF3366","#006699"))+
  theme_classic()+
  theme(axis.line=element_line(lineend="round",linewidth = 0.75),
        axis.text=element_text(size=11),
        axis.title=element_text(size=14),
        axis.ticks=element_line(linewidth = 0.75),
        strip.background=element_rect(linewidth=1.25),
        strip.text = element_text(size=11))+
  guides(fill="none")

# For Fig. 2B
ggplot(filter(d,Location!="SocAssoc15"),aes(Location,value,color=Genotype,group=pair))+
  geom_point()+
  geom_line(alpha=0.4)+facet_wrap(~Genotype,
                                  labeller=as_labeller(c('129'="129S1",
                                                         'C57'="C57BL/6",
                                                         'PWK'="PWK/PhJ")))+
  labs(x="Set of RFID reader locations",
       y="Pairwise association between individuals")+
  theme_classic()+
  scale_x_discrete(limits=c("SocAssoc15feeder","SocAssoc15noF"),
                   labels=c("Feeder\nonly","Non-\nfeeder"))+
  scale_color_manual(values=c("#FFA533","#FF3366","#006699"))+
  theme(axis.text.x=element_text(angle=45,hjust=1.0),
        axis.text=element_text(size=11),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75),
        strip.background=element_rect(linewidth=1.25),
        strip.text = element_text(size=11))+
  guides(color="none")

# Different time frames, shown in Fig. S1.
d <- full.2021.assoc.info
d$pair <- NA
for(i in 1:nrow(d)){
  d$pair[i] <- paste(d$id1[i],d$id2[i],sep="-")
}
d <- select(filter(d,Wedge!="Wedge 6"&Wedge!="Wedge 7"&
                     (!(id1%in%behaviorSelection)&!(id2%in%behaviorSelection))),
            SocAssoc240,SocAssoc,SocAssoc15,SocAssoc2min,
            pair,Wedge,Block,Genotype,infection)
d <- as.data.frame(pivot_longer(d,1:4,names_to="Overlap"))

ggplot(d,aes(Overlap,value,color=Genotype,group=pair))+
  geom_point()+
  geom_line(alpha=0.4)+facet_wrap(~Genotype,
                                  labeller=as_labeller(c('129'="129S1",
                                                         'C57'="C57BL/6",
                                                         'PWK'="PWK/PhJ")))+
  labs(x="Length of overlap window for associations",
       y="Spatiotemporal association between individuals")+
  theme_classic()+
  scale_x_discrete(limits=c("SocAssoc240","SocAssoc","SocAssoc15","SocAssoc2min"),
                   labels=c("4 hours","1 hour","15 minutes","2 minutes"))+
  scale_color_manual(values=c("#FFA533","#FF3366","#006699"))+
  theme(axis.text.x=element_text(angle=45,hjust=1.0),
        axis.text=element_text(size=11),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75),
        strip.background=element_rect(linewidth=1.25),
        strip.text = element_text(size=11))+
  guides(color="none")

# Effects of sibling relationships and cage-sharing on associations in Fig. S2
# Fig. S2A: sibling relationships
d <- filter(full.2021.assoc.info,
              !(id1%in%behaviorSelection|id2%in%behaviorSelection)&
              is.na(fullSib)==FALSE)
ggplot(d,aes(fullSib,SocAssoc15,fill=Genotype))+
  geom_boxplot(outlier.shape=NA)+
  theme_classic()+
  theme(axis.line=element_line(lineend="round",linewidth=0.75),
        axis.text=element_text(size=11),
        axis.title=element_text(size=14),
        axis.ticks=element_line(linewidth = 0.75),
        legend.position=c(0.85,0.2),
        legend.background=element_rect(color="grey10"),
        legend.text=element_text(size=10),
        legend.title=element_text(hjust=0.5))+
  labs(x="Sibling relationship",
       y="Strength of pairwise association")+
  scale_fill_manual(values=c("#FFA533","#FF3366","#006699"),
                    labels=c("129S1","C57BL/6","PWK/PhJ"),
                    name="Strain")+
  scale_x_discrete(labels=c("Not siblings\n(n = 262)","Siblings\n(n = 28)"))

# Fig. S2B: cage-sharing
ggplot(d,aes(cageOverlap,SocAssoc15,fill=Genotype))+
  geom_boxplot(outlier.shape=NA)+
  theme_classic()+
  theme(axis.line=element_line(lineend="round",linewidth = 0.75),
        axis.text=element_text(size=11),
        axis.title=element_text(size=14),
        axis.ticks=element_line(linewidth = 0.75))+
  labs(x="Cage-sharing?",
       y="Strength of pairwise association")+
  scale_fill_manual(values=c("#FFA533","#FF3366","#006699"))+
  scale_x_discrete(labels=c("Different cages\n(n = 251)","Same cage\n(n = 39)"))+
  guides(fill="none")

# Average strength of an individual's associations (Fig. S3)
ggplot(filter(mice2021GxE,is.na(StrengthAv15adj)==FALSE),
       aes(Genotype,StrengthAv15adj,color=Genotype,shape=Block))+
  #geom_boxplot(outlier.shape=NA)+
  geom_jitter(size=2.0)+
  theme_classic()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12,hjust=0.5),
        legend.position=c(0.2,0.85),
        legend.background=element_rect(color="grey10"))+
  labs(x="Mouse strain",
       y="Average strength of an individual's\nspatiotemporal associations")+
  scale_color_manual(values=c("#FFA533","#FF3366","#006699"))+
  scale_x_discrete(labels=c("129S1","C57BL/6","PWK/PhJ"))+
  guides(color="none")

# We want to show some of the basic similarity data plotted against pairwise
# associations.  This we do in Fig. S5.
# CD4 similarity vs. social association (Fig. S5A)
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$CD4JI,d2$fullSib,d2$SocAssoc15),]
ggplot(d2,aes(SocAssoc15,CD4JI,color=Genotype))+
  geom_point(size=1.8,alpha=0.5)+
  theme_classic()+
  theme(axis.line=element_line(lineend="round",linewidth = 0.75),
        axis.text=element_text(size=11),
        axis.title=element_text(size=14),
        axis.ticks=element_line(linewidth = 0.75),
        legend.position=c(0.8,0.25),
        legend.background=element_rect(color="grey10"),
        legend.text=element_text(size=10),
        legend.title=element_text(hjust=0.5))+
  labs(x="Pairwise spatiotemporal association with\n15-minute overlap",
       y="Pairwise similarity of MLN CD4 T cell\nmemory phenotypes")+
  scale_color_manual(values=c("#FFA533","#FF3366","#006699"),
                     name="Mouse strain",labels=c("129S1","C57BL/6","PWK/PhJ"))

# CBC similarity vs. social association (Fig. S5B)
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$CD4JI,d2$fullSib,d2$SocAssoc15),]
ggplot(d2,aes(SocAssoc15,wbcJI,color=Genotype))+
  geom_point(size=1.8,alpha=0.5)+
  theme_classic()+
  theme(axis.line=element_line(lineend="round",linewidth = 0.75),
        axis.text=element_text(size=11),
        axis.title=element_text(size=14),
        axis.ticks=element_line(linewidth = 0.75),
        legend.position=c(0.8,0.2),
        legend.background=element_rect(color="grey10"),
        legend.text=element_text(size=10),
        legend.title=element_text(hjust=0.5))+
  labs(x="Pairwise spatiotemporal association with\n15-minute overlap",
       y="Pairwise similarity of CBC phenotypes in\nperipheral blood")+
  scale_color_manual(values=c("#FFA533","#FF3366","#006699"),
                     name="Mouse strain",labels=c("129S1","C57BL/6","PWK/PhJ"))

# We want to show the variation in pairwise CBC similarity across the duration
# of the experiment, from pre-release to endpoint.  We do this in Fig. S6.
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$wbcJI,d2$wbcStartJ,d2$fullSib),]
ggplot(d2,aes(wbcStartJ,wbcJI,color=Genotype))+
  geom_point(alpha=0.5)+
  geom_rug(alpha=0.2)+
  theme_classic()+
  labs(x="CBC pairwise similarity (Jaccard index)\npre-release",
       y="CBC pairwise similarity (Jaccard index)\nat endpoint")+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75),
        legend.position=c(0.25,0.75),
        legend.background=element_rect(color="grey10"),
        legend.text=element_text(size=10),
        legend.title=element_text(hjust=0.5))+
  scale_color_manual(values=c("#FFA533","#FF3366","#006699"),
                     labels=c("129S1","C57BL/6","PWK/PhJ"),
                     name="Strain")

# Next we want to plot model results.
# First we show the association coefficients for the various aspects of
# cellular immune similarity (Fig. 3A).
post1 <- as.data.frame(CBC.JI.model.2)$b_SocAssoc15[1:1000]
post2 <- as.data.frame(Tmem.MLN.model.2)$b_SocAssoc15[1:1000]
post3 <- as.data.frame(CD4.model.2)$b_SocAssoc15[1:1000]
post4 <- as.data.frame(CD8.model.2)$b_SocAssoc15[1:1000]
post5 <- as.data.frame(BloodT.model.2)$b_SocAssoc15[1:1000]
post6 <- as.data.frame(Bcell.model.2)$b_SocAssoc15[1:1000] # B cells in MLNs
posts <- as.data.frame(cbind(post1,post2,post3,post4,post5,post6))
posts <- pivot_longer(posts,1:6,names_to="Model")
ggplot(posts,
       aes(value, Model))+
  geom_density_ridges(aes(fill=Model),scale=1.0,rel_min_height=0.0001)+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=11),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75))+
  labs(y="Immune cell dataset analyzed",
       x="Regression coefficient for relationship between\nspatiotemporal association and immune similarity")+
  scale_y_discrete(limits=c("post6","post5","post4",
                            "post3","post2","post1"),
                   labels=c("B cells\nin MLNs","CD4/CD8 T\nin blood",
                            "CD8 T in\nMLNs","CD4 T in\nMLNs",
                            "CD4/CD8 T\nin MLNs","CBC in\nblood"))+
  scale_fill_manual(values=rev(met.brewer("VanGogh3",6,type="discrete")))+
  guides(fill="none")

# Full results for analysis of CD4 T cells in the MLNs (Fig. 3B).
temp <- as.data.frame(CD4.model.2)
post1 <- temp$b_SocAssoc15[1:1000]
post2 <- temp$b_GenotypeC57[1:1000]
post3 <- temp$b_GenotypePWK[1:1000]
post4 <- temp$b_InfDiffSame[1:1000]
post5 <- temp$b_fullSibTRUE[1:1000]
posts <- as.data.frame(cbind(post1,post2,post3,post4,post5))
posts <- pivot_longer(posts,1:5,names_to="Parameter")
ggplot(posts,aes(value,Parameter))+
  geom_density_ridges(aes(fill=Parameter),scale=1.0)+
  stat_summary(fun="mean",size=0.5)+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=11),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75))+
  labs(y="Parameter in model",
       x="Regression-estimated coefficient value for\nrelationship with immune similarity")+
  scale_y_discrete(labels=c("Spatiotemporal\nassociation",
                            "C57 vs. 129\nstrain","PWK vs. 129\nstrain",
                            "Shared\ninfection\nstatus","Full\nsiblings"))+
  scale_fill_manual(values=met.brewer("VanGogh3",5,type="discrete"))+
  guides(fill="none")

# Social association coefficients for cytokine similarity analyses (Fig. 3C).
post1 <- as.data.frame(PCyt.mnhttn)$b_SocAssoc15[1:1000]
post2 <- as.data.frame(cytProd.model)$b_SocAssoc15[1:1000]
posts <- as.data.frame(cbind(post1,post2))
posts <- pivot_longer(posts,1:2,names_to="Model")
ggplot(posts,
       aes(-value, Model))+
  geom_density_ridges(aes(fill=Model),scale=1.0,rel_min_height=0.0001)+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=11),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75))+
  labs(y="Cytokine dataset analyzed",
       x="Regression coefficient for relationship between\nspatiotemporal association and immune similarity")+
  scale_y_discrete(labels=c("Plasma\ncytokines","MLN\ncytokine\nproduction"))+
  scale_fill_manual(values=met.brewer("VanGogh3",2,type="discrete"))+
  guides(fill="none")

# Full model results for CBC similarity analysis
temp <- as.data.frame(CBC.JI.model.2)
post1 <- temp$b_SocAssoc15[1:1000]
post2 <- temp$b_GenotypeC57[1:1000]
post3 <- temp$b_GenotypePWK[1:1000]
post4 <- temp$b_InfDiffSame[1:1000]
post5 <- temp$b_fullSibTRUE[1:1000]
posts <- as.data.frame(cbind(post1,post2,post3,post4,post5))
posts <- pivot_longer(posts,1:5,names_to="Parameter")
ggplot(posts,aes(value,Parameter))+
  geom_density_ridges(aes(fill=Parameter),scale=1.0)+
  stat_summary(fun="mean",size=0.5)+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=11),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75))+
  labs(y="Parameter in model",
       x="Regression-estimated coefficient value for\nrelationship with immune similarity")+
  scale_y_discrete(labels=c("Spatiotemporal\nassociation",
                            "C57 vs. 129\nstrain","PWK vs. 129\nstrain",
                            "Shared\ninfection\nstatus","Full\nsiblings"))+
  scale_fill_manual(values=met.brewer("VanGogh3",5,type="discrete"))+
  guides(fill="none")

# Social association coefficients from analyses of CBC similarity at different
# times during the experiment (Fig. 4A).
post1 <- as.data.frame(CBC.start.model)$b_SocAssoc15[1:1000]
post2 <- as.data.frame(CBC.mid.model)$b_SocAssoc15[1:1000]
post3 <- as.data.frame(CBC.JI.model.2)$b_SocAssoc15[1:1000]
posts <- as.data.frame(cbind(post1,post2,post3))
posts <- pivot_longer(posts,1:3,names_to="Model")
ggplot(posts,aes(value,Model))+
  geom_density_ridges(aes(fill=Model),scale=1.0,rel_min_height=0.0001)+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=11),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title.y=element_text(size=14),
        axis.title.x=element_text(size=12),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75))+
  labs(y="White blood cell dataset analyzed",
       x="Regression coefficient for relationship between\nspatiotemporal association and immune similarity")+
  scale_y_discrete(labels=c("CBC pre-\nrelease","CBC at\nmidpoint",
                            "CBC at\nendpoint"))+
  scale_fill_manual(values=met.brewer("VanGogh3",3,type="discrete"))+
  guides(fill="none")

# Social association coefficients from different spatial datasets for CD4 T
# cell similarities (Fig. 4B)
post1 <- as.data.frame(CD4.model.2)$b_SocAssoc15[1:1000]
post2 <- as.data.frame(CD4.feeder.model)$b_SocAssoc15feeder[1:1000]
post3 <- as.data.frame(CD4.noF.model)$b_SocAssoc15noF[1:1000]
posts <- as.data.frame(cbind(post1,post2,post3))
posts <- pivot_longer(posts,1:3,names_to="Model")
ggplot(posts,
       aes(value,Model))+
  geom_density_ridges(aes(fill=Model),scale=1.0,rel_min_height=0.0001)+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=11),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title.y=element_text(size=14),
        axis.title.x=element_text(size=12),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75))+
  labs(y="RFID stations for calculating\nassociation in model",
       x="Regression coefficient for relationship between\nspatiotemporal association and immune similarity")+
  scale_y_discrete(labels=c("All\nstations","Feeders\nonly","Non-feeders"))+
  scale_fill_manual(values=met.brewer("VanGogh3",3,type="discrete"))+
  guides(fill="none")

# Social association coefficients from different spatial datasets for CBC cell
# similarities (Fig. S6A)
post1 <- as.data.frame(CBC.JI.model.2)$b_SocAssoc15[1:1000]
post2 <- as.data.frame(CBC.feeder.model)$b_SocAssoc15feeder[1:1000]
post3 <- as.data.frame(CBC.noF.model)$b_SocAssoc15noF[1:1000]
posts <- as.data.frame(cbind(post1,post2,post3))
ggplot(posts,
       aes(value,Model))+
  geom_density_ridges(aes(fill=Model),scale=1.0,rel_min_height=0.0001)+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=11),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title.y=element_text(size=14),
        axis.title.x=element_text(size=12),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75))+
  labs(y="RFID stations for calculating\nassociation in model",
       x="Regression coefficient for relationship between\nspatiotemporal association and immune similarity")+
  scale_y_discrete(labels=c("All\nstations","Feeders\nonly","Non-feeders"))+
  scale_fill_manual(values=met.brewer("VanGogh3",3,type="discrete"))+
  guides(fill="none")

# Social association coefficients from different overlap windows for CD4 T cell
# similarities (Fig. 4D)
post1 <- as.data.frame(CD4.space.model)$b_SpaceUse[1:1000]
post2 <- as.data.frame(CD4.4hr.model)$b_SocAssoc240[1:1000]
post3 <- as.data.frame(CD4.1hr.model)$b_SocAssoc[1:1000]
post4 <- as.data.frame(CD4.model.2)$b_SocAssoc15[1:1000]
post5 <- as.data.frame(CD4.2min.model)$b_SocAssoc2min[1:1000]
posts <- as.data.frame(cbind(post1,post2,post3,post4,post5))
posts <- pivot_longer(posts,1:5,names_to="Model")
ggplot(posts,
       aes(value, Model))+
  geom_density_ridges(aes(fill=Model),scale=1.0,rel_min_height=0.0001)+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=11),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title.y=element_text(size=14),
        axis.title.x=element_text(size=12),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75))+
  labs(y="Method for calculating social network",
       x="Regression coefficient for relationship between\nspatiotemporal association and immune similarity")+
  scale_y_discrete(labels=c("Shared\nspace use","4 hour\noverlap","1 hour\noverlap",
                            "15 minute\noverlap",
                            "2 min\noverlap"))+
  scale_fill_manual(values=met.brewer("VanGogh3",5,type="discrete"))+
  guides(fill="none")

# Social association coefficients from different overlap windows for CBC
# similarities (Fig. S6B)
post1 <- as.data.frame(CBC.space.model)$b_SpaceUse[1:1000]
post2 <- as.data.frame(CBC.4hr)$b_SocAssoc240[1:1000]
post3 <- as.data.frame(CBC.1hr)$b_SocAssoc[1:1000]
post4 <- as.data.frame(CBC.JI.model.2)$b_SocAssoc15[1:1000]
post5 <- as.data.frame(CBC.2min)$b_SocAssoc2min[1:1000]
posts <- as.data.frame(cbind(post1,post2,post3,post4,post5))
posts <- pivot_longer(posts,1:4,names_to="Model")
ggplot(posts,
ggplot(posts,
       aes(value, Model))+
  geom_density_ridges(aes(fill=Model),scale=1.0,rel_min_height=0.0001)+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=11),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title.y=element_text(size=14),
        axis.title.x=element_text(size=12),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75))+
  labs(y="Method for calculating social network",
       x="Regression coefficient for relationship between\nspatiotemporal association and immune similarity")+
  scale_y_discrete(labels=c("Shared\nspace use","4 hour\noverlap","1 hour\noverlap",
                            "15 minute\noverlap",
                            "2 min\noverlap"))+
  scale_fill_manual(values=met.brewer("VanGogh3",5,type="discrete"))+
  guides(fill="none")

# Models including microbiome similarity
# First CD4 T cell memory similarity with microbiome
post1 <- as.data.frame(CD4.micro.model)$b_SocAssoc15[1:1000] # Social association
post2 <- as.data.frame(CD4.micro.model)$b_GenotypeC57[1:1000] # C57 vs. 129
post3 <- as.data.frame(CD4.micro.model)$b_GenotypePWK[1:1000] # PWK vs. 129
post4 <- as.data.frame(CD4.micro.model)$b_InfDiffSame[1:1000] # Shared infection status
post5 <- as.data.frame(CD4.micro.model)$b_fullSibTRUE[1:1000] # Full siblings
post6 <- as.data.frame(CD4.micro.model)$b_mbJI[1:1000] # Microbiome
posts <- as.data.frame(cbind(post1,post2,post3,post4,
                             post5,post6))
posts <- pivot_longer(posts,1:6,names_to="Parameter")
ggplot(posts,aes(value,Parameter))+
  geom_density_ridges(aes(fill=Parameter),scale=1.0,rel_min_height=0.0001)+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=11),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title.y=element_text(size=14),
        axis.title.x=element_text(size=12),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75))+
  labs(y="Parameter in model",
       x="Regression coefficient value for\nrelationship with immune similarity")+
  scale_y_discrete(labels=c("Social\nassociation","C57 vs.\n129","PWK vs.\n129",
                            "Same\ninfection\nstatus","Full\nsiblings","Microbiome"))+
  scale_fill_manual(values=met.brewer("VanGogh3",6,type="discrete"))+
  guides(fill="none")

# Then CBC similarity with microbiome
post1 <- as.data.frame(CBC.micro.model)$b_SocAssoc15[1:1000] # Social association
post2 <- as.data.frame(CBC.micro.model)$b_GenotypeC57[1:1000] # C57 vs. 129
post3 <- as.data.frame(CBC.micro.model)$b_GenotypePWK[1:1000] # PWK vs. 129
post4 <- as.data.frame(CBC.micro.model)$b_InfDiffSame[1:1000] # Shared infection status
post5 <- as.data.frame(CBC.micro.model)$b_fullSibTRUE[1:1000] # Full siblings
post6 <- as.data.frame(CBC.micro.model)$b_mbJI[1:1000] # Microbiome
posts <- as.data.frame(cbind(post1,post2,post3,post4,
                             post5,post6))
posts <- pivot_longer(posts,1:6,names_to="Parameter")
ggplot(posts,aes(value,Parameter))+
  geom_density_ridges(aes(fill=Parameter),scale=1.0)+
  stat_summary(fun="mean",size=0.5)+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=11),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75))+
  labs(y="Parameter in model",
       x="Regression-estimated coefficient value for\nrelationship with immune similarity")+
  scale_y_discrete(labels=c("Spatiotemporal\nassociation","C57 vs.\n129","PWK vs.\n129",
                            "Same\ninfection\nstatus","Full\nsiblings","Microbiome"))+
  scale_fill_manual(values=met.brewer("VanGogh3",6,type="discrete"))+
  guides(fill="none")

# Here is the comparison of association during the first portion of the
# experiment with association during the second portion, after the midpoint
# trapping session (Fig. S9).
d2 <- filter(full.2021.assoc.info,Wedge%in%c("Wedge 2","Wedge 3","Wedge 4"))
d2 <- d2[complete.cases(d2$SocAssoc15),]
ggplot(filter(d2,!(id1%in%behaviorSelection|id2%in%behaviorSelection)),
       aes(SocAssoc1.15,SocAssoc2.15,color=Genotype))+
  geom_point(alpha=0.5)+
  labs(x="Spatiotemporal association between release and\nmidpoint trapping session",
       y="Spatiotemporal association between midpoint and\nendpoint trapping session")+
  geom_abline(slope=1,intercept=0,linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(linewidth=0.75),
        legend.position=c(0.75,0.25),
        legend.background=element_rect(color="grey10"),
        legend.text=element_text(size=10),
        legend.title=element_text(hjust=0.5))+
  scale_color_manual(values=c("#FFA533","#FF3366","#006699"),
                     labels=c("129S1","C57BL/6","PWK/PhJ"),
                     name="Mouse strain")

# This figure shows the comparison between space use and 15-minute
# spatiotemporal association. (Fig. S8)
d2 <- full.2021.assoc.info
d2 <- d2[complete.cases(d2$SpaceUse,d2$SocAssoc15),]
ggplot(d4,aes(SpaceUse,SocAssoc15,color=Genotype))+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(axis.line=element_line(lineend="round",linewidth=0.75),
        legend.position=c(0.2,0.75),
        legend.background=element_rect(color="grey10"),
        legend.text=element_text(size=10),
        legend.title=element_text(hjust=0.5),
        axis.text=element_text(size=11),
        axis.title=element_text(size=14),
        axis.ticks=element_line(linewidth = 0.75))+
  scale_color_manual(values=c("#FFA533","#FF3366","#006699"),
                     labels=c("129S1","C57BL/6","PWK/PhJ"),
                     name="Mouse strain")+
  coord_cartesian(ylim=c(0,0.75))+
  labs(x="Pairwise space use similarity",
       y="Spatiotemporal association from\n15-minute overlap window")

# Finally, we have some figures showing trends of aggregate individual check-in
# patterns across the experiment.
# First we show variation in check-ins recorded by night (Fig. S10A and B).
d <- CIdata2021GxEfiltered
d <- filter(d,!Mouse%in%behaviorSelection)
ggplot(filter(d,Block=="Block 1"),
       aes(Julian_Night,fill=Genotype))+
  geom_bar()+
  facet_wrap(~Genotype,scales="free",
             labeller=as_labeller(c('129'="129S1",
                                    'C57'="C57BL/6",
                                    'PWK'="PWK/PhJ")))+
  labs(x="Night of the Experiment",y="Number of total check-ins on night",
       title="Block 1")+
  theme_classic()+
  theme(axis.text.x=element_blank(),plot.title=element_text(hjust=0.5),
        axis.line=element_line(lineend="round",linewidth = 0.75),
        axis.text.y=element_text(size=11),
        axis.title=element_text(size=14),
        axis.ticks=element_line(linewidth = 0.75))+
  scale_fill_manual(values=c("#FFA533","#FF3366","#006699"))+
  guides(fill="none")

ggplot(filter(d,Block=="Block 2"),
       aes(Julian_Night,fill=Genotype))+
  geom_bar()+
  facet_wrap(~Genotype,scales="free",
             labeller=as_labeller(c('129'="129S1",
                                    'C57'="C57BL/6",
                                    'PWK'="PWK/PhJ")))+
  labs(x="Night of the Experiment",y="Number of total check-ins on night",
       title="Block 2")+
  theme_classic()+
  theme(axis.text.x=element_blank(),plot.title=element_text(hjust=0.5),
        axis.line=element_line(lineend="round"))+
  scale_fill_manual(values=c("#FFA533","#FF3366","#006699"))+
  guides(fill="none")

# We also show how the aggregate total of check-ins varies by hour (Fig. S11).
ggplot(d,aes(Hour,fill=Genotype))+
  geom_bar()+
  facet_wrap(~Genotype,scales="free",
             labeller=as_labeller(c('129'="129S1",
                                    'C57'="C57BL/6",
                                    'PWK'="PWK/PhJ")))+
  labs(x="Hour of the day",y="Number of total check-ins during hour\nacross experiment")+
  theme_classic()+
  theme(axis.text.x=element_blank(),
        axis.line=element_line(lineend="round"))+
  scale_fill_manual(values=c("#FFA533","#FF3366","#006699"))+
  guides(fill="none")

ggplot(filter(d,Genotype=="PWK"),
       aes(Hour))+
  geom_bar(fill="#006699")+
  facet_wrap(~Block,scales="free")+
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.line=element_line(lineend="round"))+
  labs(x="Hour of the day",y="Number of total check-ins during hour\nacross experiment")

# And the breakdown of check-ins by location (Fig. S12).
ggplot(d,aes(Location,fill=Genotype))+
  geom_bar()+facet_grid(rows=vars(Block),cols=vars(Genotype),scales="free",
                        labeller=as_labeller(c('129'="129S1",
                                               'C57'="C57BL/6",
                                               'PWK'="PWK/PhJ",
                                               'Block 1'="Block 1",
                                               'Block 2'="Block 2")))+
  labs(x="Location of reader in enclosure",
       y="Number of check-ins across experiment")+
  theme_classic()+
  theme(axis.line=element_line(lineend="round"),
        axis.text.x=element_text(angle=45,hjust=1.0))+
  scale_fill_manual(values=c("#FFA533","#FF3366","#006699"))+
  guides(fill="none")

