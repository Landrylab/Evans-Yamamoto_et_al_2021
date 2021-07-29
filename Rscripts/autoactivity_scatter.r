#!/usr/bin/env Rscript
if(require("ggplot2")){
  #print("ggplot2 is loaded correctly")
} else {
  print("trying to install ggplot2...")
  install.packages("ggplot2")
  if(require("ggplot2")){
    print("ggplot2 installed and loaded")
  } else {
    stop("could not install ggplot2")
  }
}
if(require("scales")){
  #print("scales is loaded correctly")
} else {
  print("trying to install scales...")
  install.packages("scales")
  if(require("scales")){
    print("scales installed and loaded")
  } else {
    stop("could not install scales")
  }
}
if(require("cowplot")){
  #print("cowplot is loaded correctly")
} else {
  print("trying to install cowplot...")
  install.packages("cowplot")
  if(require("cowplot")){
    print("cowplot installed and loaded")
  } else {
    stop("could not install cowplot")
  }
}
if(require("grid")){
  #print("grid is loaded correctly")
} else {
  print("trying to install grid...")
  install.packages("grid")
  if(require("grid")){
    print("grid installed and loaded")
  } else {
    stop("could not install grid")
  }
}
if(require("gridExtra")){
  #print("gridExtra is loaded correctly")
} else {
  print("trying to install gridExtra...")
  install.packages("gridExtra")
  if(require("gridExtra")){
    print("gridExtra installed and loaded")
  } else {
    stop("could not install gridExtra")
  }
}

if(require("GGally")){
  #print("gridExtra is loaded correctly")
} else {
  print("trying to install gridExtra...")
  install.packages("GGally")
  if(require("GGally")){
    print("GGally installed and loaded")
  } else {
    stop("could not install GGally")
  }
}




library(ggplot2)
library(scales) # for muted function
library(gridExtra)
library(grid)
library(reshape2)
library(egg)

args = commandArgs(trailingOnly=TRUE)


data = read.csv(file=args[1],header=TRUE)
out       = args[2]
width_cm  = as.numeric(args[3])
height_cm = as.numeric(args[4])
p_width   = as.numeric(args[5])
p_height  = as.numeric(args[6])

data$cond = paste(data$Strain_type,data$Condition)
data$group = paste(data$cond,data$Replicate)
scatter  = 
  ggplot(data=data,
         aes(x=AA_median_rank,y=log10(AA_median+1)),color="#000000")+ 
  geom_point(alpha=1,size=0.5) + 
  scale_x_continuous(expand=c(0,0)) +  
  scale_y_continuous(expand=c(0,0),limits=c(-0.1,ceiling(log10(max(data$AA_median)))*1.1 ))+ 
  theme(
    axis.title=element_text(size=9.5,color="#000000"), 
    legend.position='none',
    legend.key = element_blank(), 
    strip.background = element_rect(color="#FFFFFF",fill="#FFFFFF"),
    panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
    panel.spacing = unit(0.75, "lines"),
    panel.background = element_blank(),
    aspect.ratio=0.5,
    axis.ticks =element_line(color = "#000000", size = 0.5),
    axis.text.x = element_text(color="#000000",size=9.5),
    axis.text.y = element_text(color="#000000",size=9.5) ,
    panel.border = element_rect(color="#000000", fill = NA,size=1) 
    )+
  xlab("Rank")+
  ylab(bquote(~"Median"~ '(Raw interaction signal)'))+
  facet_wrap(~group)

p <-set_panel_size(scatter,width  = unit(width_cm, "cm"),height = unit(height_cm, "cm"))
grid.newpage()
grid.draw(p)
ggsave(plot=p,out,height=p_height,width=p_width,units = "cm")
