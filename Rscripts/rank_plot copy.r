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


data150 = data[data$Rank < 101,]
data150$th = as.factor(data150$TH)

plot = ggplot(data=data150,aes(Rank,log10(score))) + geom_jitter(aes(colour = th),size=0.25,alpha=0.5)+  scale_color_manual(breaks = c("0", "1"),values=c("#ababab", "#FF0000"))+theme(legend.position='none',legend.key = element_blank(), strip.background = element_rect(color="#FFFFFF",fill="#FFFFFF"),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.spacing = unit(0.75, "lines"),panel.background = element_blank(),aspect.ratio=0.02,axis.ticks =element_line(color = "#000000", size = 0.5),axis.text.x = element_blank(),axis.text.y = element_text(color="#000000",size=9.5) ,axis.title=element_blank(),axis.line.x=element_line(size=0.5),axis.line.y =element_line(size=0.5) ) +scale_x_continuous(expand=c(0,0),limits=c(0.7,101),breaks=c())+xlab("")+scale_y_continuous(expand=c(0,0),limits=c(-0.1,4),breaks=c(0,1,2,3,4,5,6,7))

P <-set_panel_size(plot,width  = unit(13.5, "cm"),height = unit(3, "cm"))
grid.newpage()
grid.draw(P)
name = paste("./plot/",gsub(".csv", "points_pca.pdf", args[1]))
nm   = gsub(" ", "", name)
ggsave(plot=P,nm ,height=2)



data150$Method  = NULL
data150$Selection  = NULL
data150$n              = NULL
data150$Norm_meth  = NULL
data150$Th_method      = NULL
data150$s  = NULL
data150$Precision  = NULL
data150$Recall         = NULL
data150$MCC  = NULL
data150$score       = NULL
data150$TH       = NULL
data150$th       = NULL

data150hm = melt(data150,id.var=c("Rank","ORF_pair"))

data150hm2 = data150hm[data150hm$variable!= "X4_ALL",]

plot =ggplot(data150hm2)+ geom_tile(aes(Rank,variable,fill=value))+ theme_bw()   + theme(panel.border= element_rect(size=0.5), legend.position='none',legend.key = element_blank(), strip.background = element_rect(colour="#FFFFFF", color="#FFFFFF",fill="#FFFFFF"),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_blank(),aspect.ratio=0.02,axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank()  ,axis.title = element_blank()     )+ scale_x_discrete(expand=c(0,0),limits=c(0,101)) + scale_y_discrete(expand=c(0,0),limits=rev) +scale_fill_gradientn(colors=c("#FFFFFF","#ababab","#000000","#000000"),na.value = "#FFFFFF")
P <-set_panel_size(plot,width  = unit(13.5, "cm"),height = unit(1, "cm"))
grid.newpage()
grid.draw(P)
name = paste("./plot/",gsub(".csv", "hmap.pdf", args[1]))
nm   = gsub(" ", "", name)
ggsave(plot=P,nm ,height=1)
