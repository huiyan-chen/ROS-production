#ROS-Huiyan

detach()

# define the working directory
setwd("C:\\Users\\acer\\Desktop\\writing\\5. result\\home\\homozygous line")


# read.table allows to specify the options of the file format
ROS1 <- read.csv("A1.csv")
attach(ROS1)
str(ROS1)

#separate the data  
library(tidyr)
TROS <- data.frame(ROS1$Cycle.Nr.,gather(ROS1, "Cycle.Nr.", "RLU", 2:32))
seROS <- separate(TROS,ROS1.Cycle.Nr., into = c("group", "sample"),sep=1)
sepROS <- separate(seROS,Cycle.Nr., into = c("useless", "cycles"),sep=1)
str(sepROS)

#set up new data frame
library(dplyr)
library(plyr)
cycless <- as.numeric(sepROS$cycles)
time <- cycless*2-2
summary(time)
Treatments <- as.factor(sepROS$group)
Treatments <- revalue(Treatments, c("A"="Col-0+Water","B"="Col-0+Flagellin22","C"="SPE2-GFP+Water","D"="SPE2-GFP+Flagellin22","E"="Col-0+Water1","F"="Col-0+Flagellin221","G"="SPE7-GFP+Water","H"="SPE7-GFP+Flagellin22"))
FinalROS <- data.frame(select(sepROS, -useless,-group),time,Treatments)
str(FinalROS)

#description statistics
library(Rmisc)
StaROS <- summarySE(FinalROS, measurevar="RLU", groupvars=c("Treatments","time"))
str(StaROS)

# draw a liner plot
library(ggplot2)
library(ggpubr)
library(ggsci)
#PROS <- ggplot(StaROS, aes(x=time, y=RLU,colour=treatments, group=treatments)) + 
#  geom_errorbar(aes(ymin=RLU-se, ymax=RLU+se), width=.1) +
#  geom_line() +
#  geom_point()
#ggpar(PROS, palette = "npg")

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

  ggplot(StaROS, aes(x=time, y=RLU, colour=Treatments, group=Treatments)) + 
  geom_errorbar(aes(ymin=RLU-se, ymax=RLU+se), colour="black", width=4, position=pd) +
  geom_line(position=pd,size=1.2) +
  geom_point(position=pd, size=1, shape=20) + # 21 is filled circle
  xlab("Time (min)") +
  ylab("Luminesence (RLU)") +
  color_palette("npg") +
  ggtitle("ROS Brust") +
  expand_limits(x = c(0, 60), y = c(0, 40000)) +                        # Expand y range
  theme_bw() +
  
  theme(panel.border=element_rect(size=2),
        title= element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold",size = 12),
        axis.text = element_text(face = "bold",size = 12),
        legend.background=element_rect(fill = NA, colour = NA),
        legend.text=element_text(face = 3),
        legend.justification=c(1,1),
        legend.position=c(1,1))  # Position legend in bottom right

write.table(FinalROS,file="SPE2-1 RLU.csv",sep = ",")
##580*415
