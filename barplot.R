

# define the working directory
setwd("C:\\Users\\acer\\Desktop\\writing\\5. result\\home\\homozygous line")


# read.table allows to specify the options of the file format
#Matr <- read.table("GermanFlora_forR.csv", header = TRUE, sep = ";", dec=",")

ROS1 <- read.csv("A1.csv")


attach(ROS1)
str(ROS1)

#separate the data and set up a data frame
library(tidyr)
library(dplyr)
library(plyr)
PreROS <- ROS1[,2:32]
TRLU <- rowSums(PreROS)
TRLUe3 <- TRLU/1000
SumROS <- data.frame(ROS1$Cycle.Nr.,TRLU,TRLUe3)
seROS <- separate(SumROS,ROS1.Cycle.Nr., into = c("group", "sample"),sep=1)
str(seROS)

treatments <- as.factor(seROS$group)
treatments <- revalue(treatments, c("A"="Col-0+Water","B"="Col-0+Flagellin22","C"="SPE2-GFP+Water","D"="SPE2-GFP+Flagellin22","E"="Col-01+Water","F"="Col-01+Flagellin22","G"="SPE7-GFP+Water","H"="SPE7-GFP+Flagellin22"))
FinalTROS <- data.frame(seROS, treatments)

write.table(FinalTROS,file="A1 TROS.csv",sep = ",") # export a table

#description statistics
library(Rmisc)
StaROS <- summarySE(FinalTROS, measurevar="TRLUe3", groupvars=c("treatments"))
str(StaROS)

#set up a data frame for legand
NaROS <- separate(StaROS,treatments, into = c("plant", "elicitor"),extra= "merge",remove = FALSE)
plant <- as.factor(NaROS$plant)
Plant <- revalue(plant, c("Col"="Col-0","SPE2"="SPE2-GFP","SPE7"="SPE7-GFP"))
elicitor <- as.factor(NaROS$elicitor)
Elicitor <- revalue(elicitor, c("0+Water"="Water","GFP+Water"="Water","01+Water"="Water","0+Flagellin22"="Flagellin22","GFP+Flagellin22"="Flagellin22","01+Flagellin22"="Flagellin22"))

NameROS <- data.frame(select(NaROS, -plant,-elicitor),Plant,Elicitor)
str(NameROS)
NameROS$Plant <- ordered(NameROS$Plant,
                           levels = c("Col-0","SPE2-GFP","SPE7-GFP"))
NameROS$Elicitor <- ordered(NameROS$Elicitor,
                            levels = c("Water", "Flagellin22"))


#anova
# Compute the analysis of variance
res.aov <- aov(TRLUe3 ~ treatments, data = FinalTROS)
# Summary of the analysis
summary(res.aov)

TukeyHSD(res.aov)
siROS <- TukeyHSD(res.aov)
sigROS <- as.data.frame(siROS[1:1])
str(sigROS)

write.table(sigROS,file="A1 TukeyTROS.csv",sep = ",") 

# TukeyHSD 
library(agricolae)
lev <- HSD.test(res.aov, "treatments", group=TRUE)
print(lev)
#example:LevROS <- as.data.frame(lev$groups)
#write.table(Lev,file="LevTROS.csv",sep = ",")
## open file and ename by hand "treatment"
#LevTROS<- read.csv("LevTROS.csv")
#str(LevTROS)

#MergeROS <- merge(NameROS, LevTROS, by.x = "treatments", by.y = "treatment") 
sig <-c("d","bc","d","bcd","d","b","cd","a") #manually add
plotROS <- data.frame(NameROS,sig)
plotROS$siglevel <- as.character(factor(plotROS$sig))
str(plotROS)

write.table(plotROS,file="plotROS.csv",sep = ",")
plotROS1 <- read.csv("plotROS.csv")
str(plotROS1)
plotROS1$siglevel1 <- as.character(factor(plotROS1$siglevel))
str(plotROS1)

#PLOT
library(ggplot2)
library(ggpubr)
library(ggsci)

pd <- position_dodge(0.9)

Plot <-
  ggplot(plotROS1, aes(x=Elicitor,y=TRLUe3, fill = Plant))+
  geom_bar(position=pd,stat="identity", width = 0.6,colour="black",size= 0.7)+      # Thinner lines
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(x=1, y=500, label="Stretch it"), vjust=-1)+
  geom_errorbar(aes(ymin=TRLUe3-se, ymax=TRLUe3+se), colour="black", width=0.3,size=0.7,position=pd) +
  ylab("Total Luminesence (RLU*1000)") +
  xlab ("") +
  ggtitle("ROS Brust") + 
  fill_palette(palette ="grey")+
  theme_bw() +
  
  theme(panel.border=element_rect(size=1.5),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line (colour = "black"),
        title= element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold",size = 12),
        axis.text.x = element_text(face = "bold" ,size = 12, colour = "black"),
        axis.text.y = element_text(face = "bold",size = 12, colour = "black"),
        legend.title = element_blank(),
        legend.background=element_rect(fill = NA, colour = NA),
        legend.text=element_text(face = 4),
        legend.justification=c(1,1),
        legend.position=c(1,1))
Plot

FinalPlot<-
  Plot + geom_text(aes(label = plotROS1$siglevel1,y=TRLUe3+se), size=5,stat = "identity",position = pd,vjust = -0.2)
  

FinalPlot 
#"SPE2-1 barplot with tukey 450x427"
