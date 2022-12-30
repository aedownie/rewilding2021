# This script contains methods for plotting the model results and other figures
# in the manuscript.  To do this, we need the outputs from previous scripts:
# individual.behavior.calculation.R, social.association.calculation.R, and
# statistical.analyses.R
library(dplyr)
library(tidyr)
library(ggplot2)

# We start with the plots showing behavior

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
posts <- tidyr::pivot_longer(posts,1:6,names_to="Model")
ggplot(posts,
       aes(value, Model))+
  geom_violin(aes(fill=Model))+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(size=0.75))+
  labs(y="Immune cell dataset analyzed",
       x="Estimated social association coefficient")+
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
posts <- tidyr::pivot_longer(posts,1:5,names_to="Parameter")
ggplot(posts,aes(value,Parameter))+
  geom_violin(aes(fill=Parameter),scale="width",linewidth=0.6)+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(size=0.75))+
  labs(y="Parameter in model",
       x="Estimated coefficient value in model")+
  scale_y_discrete(labels=c("Social\nassociation",
                            "C57 vs. 129\nstrain","PWK vs. 129\nstrain",
                            "Shared\ninfection\nstatus","Full\nsiblings"))+
  scale_fill_manual(values=met.brewer("VanGogh3",5,type="discrete"))+
  guides(fill="none")

# Social association coefficients for cytokine similarity analyses (Fig. 3C).
post1 <- as.data.frame(PCyt.mnhttn)$b_SocAssoc15[1:1000]
post2 <- as.data.frame(cytProd.model)$b_SocAssoc15[1:1000]
posts <- as.data.frame(cbind(post1,post2))
posts <- tidyr::pivot_longer(posts,1:2,names_to="Model")
ggplot(posts,
       aes(-value, Model))+
  geom_violin(aes(fill=Model))+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(size=0.75))+
  labs(y="Cytokine dataset analyzed",
       x="Estimated social association coefficient")+
  scale_y_discrete(labels=c("Plasma\ncytokines","MLN cytokine\nproduction"))+
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
posts <- tidyr::pivot_longer(posts,1:5,names_to="Parameter")
ggplot(posts,aes(value,Parameter))+
  geom_violin(aes(fill=Parameter))+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(size=0.75))+
  labs(y="Parameter in model",
       x="Estimated coefficient value in model")+
  scale_y_discrete(labels=c("Social\nassociation",
                            "C57 vs. 129\nstrain","PWK vs. 129\nstrain",
                            "Shared\ninfection\nstatus","Full siblings"))+
  scale_fill_manual(values=met.brewer("VanGogh3",5,type="discrete"))+
  guides(fill="none")

# Social association coefficients from analyses of CBC similarity at different
# times during the experiment (Fig. 4A).
post1 <- as.data.frame(CBC.start.model)$b_SocAssoc15[1:1000]
post2 <- as.data.frame(CBC.mid.model)$b_SocAssoc15[1:1000]
post3 <- as.data.frame(CBC.JI.model.2)$b_SocAssoc15[1:1000]
posts <- as.data.frame(cbind(post1,post2,post3))
posts <- tidyr::pivot_longer(posts,1:3,names_to="Model")
ggplot(posts,aes(value,Model))+
  geom_violin(aes(fill=Model))+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(size=0.75))+
  labs(y="White blood cell dataset analyzed",
       x="Estimated social association coefficient")+
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
posts <- tidyr::pivot_longer(posts,1:3,names_to="Model")
ggplot(posts,
       aes(value,Model))+
  geom_violin(aes(fill=Model))+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(size=0.75))+
  labs(y="RFID stations for association in model",
       x="Estimated social association coefficient")+
  scale_y_discrete(labels=c("All\nstations","Feeders\nonly","Non-feeders"))+
  scale_fill_manual(values=met.brewer("VanGogh3",3,type="discrete"))+
  guides(fill="none")

# Social association coefficients from different spatial datasets for CBC cell
# similarities (Fig. S6A)
post1 <- as.data.frame(CBC.JI.model.2)$b_SocAssoc15[1:1000]
post2 <- as.data.frame(CBC.feeder.model)$b_SocAssoc15feeder[1:1000]
post3 <- as.data.frame(CBC.noF.model)$b_SocAssoc15noF[1:1000]
posts <- as.data.frame(cbind(post1,post2,post3))
posts <- tidyr::pivot_longer(posts,1:3,names_to="Model")
ggplot(posts,
       aes(value,Model))+
  geom_violin(aes(fill=Model))+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(size=0.75))+
  labs(y="RFID stations for association in model",
       x="Estimated social association coefficient")+
  scale_y_discrete(labels=c("All\nstations","Feeders\nonly","Non-feeders"))+
  scale_fill_manual(values=met.brewer("VanGogh3",3,type="discrete"))+
  guides(fill="none")

# Social association coefficients from different overlap windows for CD4 T cell
# similarities (Fig. 4D)
post1 <- as.data.frame(CD4.4hr.model)$b_SocAssoc240[1:1000]
post2 <- as.data.frame(CD4.1hr.model)$b_SocAssoc[1:1000]
post3 <- as.data.frame(CD4.model.2)$b_SocAssoc15[1:1000]
post4 <- as.data.frame(CD4.2min.model)$b_SocAssoc2min[1:1000]
posts <- as.data.frame(cbind(post1,post2,post3,post4))
posts <- tidyr::pivot_longer(posts,1:4,names_to="Model")
ggplot(posts,
       aes(value, Model))+
  geom_violin(aes(fill=Model))+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(size=0.75))+
  labs(y="Length of overlap window in model",
       x="Estimated social association coefficient")+
  scale_y_discrete(labels=c("4 hours","1 hour","Base\n(15 min)",
                            "2 min"))+
  scale_fill_manual(values=met.brewer("VanGogh3",4,type="discrete"))+
  guides(fill="none")

# Social association coefficients from different overlap windows for CBC
# similarities (Fig. S6B)
post1 <- as.data.frame(CBC.4hr)$b_SocAssoc240[1:1000]
post2 <- as.data.frame(CBC.1hr)$b_SocAssoc[1:1000]
post3 <- as.data.frame(CBC.JI.model.2)$b_SocAssoc15[1:1000]
post4 <- as.data.frame(CBC.2min)$b_SocAssoc2min[1:1000]
posts <- as.data.frame(cbind(post1,post2,post3,post4))
posts <- tidyr::pivot_longer(posts,1:4,names_to="Model")
ggplot(posts,
       aes(value, Model))+
  geom_violin(aes(fill=Model))+
  geom_vline(aes(xintercept=0),linetype=2)+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.text.y=element_text(angle=45,hjust=0.5),
        axis.title=element_text(size=14),
        axis.line=element_line(lineend="round",linewidth=0.75),
        axis.ticks=element_line(size=0.75))+
  labs(y="Length of overlap window in model",
       x="Estimated social association coefficient")+
  scale_y_discrete(labels=c("4 hours","1 hour","Base\n(15 min)",
                            "2 min"))+
  scale_fill_manual(values=met.brewer("VanGogh3",4,type="discrete"))+
  guides(fill="none")


