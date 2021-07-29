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


data$reported = as.factor(data$BinUNION)
christien = data
christien$Tarassov = NULL
christien$BFG.Y2H = NULL
CR_2 = na.omit(christien)


scatter  = ggplot(data=CR_2,aes(x=log10(BFG.PCA),y=Christien_2_2,color=reported))+scale_color_manual(breaks = c("0", "1"),values=c("#000000", "#FF0000")) + geom_point(alpha=0.3,size=1) + scale_x_continuous(expand=c(0,0),limits=c(0.45,2.75) ) +  scale_y_continuous(expand=c(0,0),limits=c(-13.5,13.5) )+ theme(legend.position='none',legend.key = element_blank(), strip.background = element_rect(color="#FFFFFF",fill="#FFFFFF"),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.spacing = unit(0.75, "lines"),panel.background = element_blank(),aspect.ratio=1.0,axis.ticks =element_line(color = "#000000", size = 0.5),axis.text.x = element_text(color="#000000",size=9.5),axis.text.y = element_text(color="#000000",size=9.5) , panel.border = element_rect(color="#000000", fill = NA,size=1),axis.title=element_text(size=9.5) )  +geom_vline(xintercept=log10(9),linetype="dashed",color="#ababab")+geom_hline(yintercept=2.5,linetype="dashed",color="#ababab")
p <-set_panel_size(scatter,width  = unit(3.5, "cm"),height = unit(3.5, "cm"))
grid.newpage()
grid.draw(p)
ggsave(plot=p,"comparison_Christien_2_2.pdf" ,width=2.5,height=2.5)
scatter  = ggplot(data=CR_2,aes(x=log10(BFG.PCA),y=Christien_4_2,color=reported))+scale_color_manual(breaks = c("0", "1"),values=c("#000000", "#FF0000")) + geom_point(alpha=0.3,size=1) + scale_x_continuous(expand=c(0,0),limits=c(0.45,2.75) ) +  scale_y_continuous(expand=c(0,0),limits=c(-13.5,13.5) )+ theme(legend.position='none',legend.key = element_blank(), strip.background = element_rect(color="#FFFFFF",fill="#FFFFFF"),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.spacing = unit(0.75, "lines"),panel.background = element_blank(),aspect.ratio=1.0,axis.ticks =element_line(color = "#000000", size = 0.5),axis.text.x = element_text(color="#000000",size=9.5),axis.text.y = element_text(color="#000000",size=9.5) , panel.border = element_rect(color="#000000", fill = NA,size=1),axis.title=element_text(size=9.5) )  +geom_vline(xintercept=log10(9),linetype="dashed",color="#ababab")+geom_hline(yintercept=2.5,linetype="dashed",color="#ababab")
p <-set_panel_size(scatter,width  = unit(3.5, "cm"),height = unit(3.5, "cm"))
grid.newpage()
grid.draw(p)
ggsave(plot=p,"comparison_Christien_4_2.pdf" ,width=2.5,height=2.5)

scatter  = ggplot(data=CR_2,aes(x=log10(BFG.PCA),y=Christien_2_4,color=reported))+scale_color_manual(breaks = c("0", "1"),values=c("#000000", "#FF0000")) + geom_point(alpha=0.3,size=1) + scale_x_continuous(expand=c(0,0),limits=c(0.45,2.75) ) +  scale_y_continuous(expand=c(0,0),limits=c(-13.5,13.5) )+ theme(legend.position='none',legend.key = element_blank(), strip.background = element_rect(color="#FFFFFF",fill="#FFFFFF"),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.spacing = unit(0.75, "lines"),panel.background = element_blank(),aspect.ratio=1.0,axis.ticks =element_line(color = "#000000", size = 0.5),axis.text.x = element_text(color="#000000",size=9.5),axis.text.y = element_text(color="#000000",size=9.5) , panel.border = element_rect(color="#000000", fill = NA,size=1),axis.title=element_text(size=9.5) )  +geom_vline(xintercept=log10(9),linetype="dashed",color="#ababab")+geom_hline(yintercept=2.5,linetype="dashed",color="#ababab")
p <-set_panel_size(scatter,width  = unit(3.5, "cm"),height = unit(3.5, "cm"))
grid.newpage()
grid.draw(p)
ggsave(plot=p,"comparison_Christien_2_4.pdf" ,width=2.5,height=2.5)

scatter  = ggplot(data=CR_2,aes(x=log10(BFG.PCA),y=Christien_4_4,color=reported))+scale_color_manual(breaks = c("0", "1"),values=c("#000000", "#FF0000")) + geom_point(alpha=0.3,size=1) + scale_x_continuous(expand=c(0,0),limits=c(0.45,2.75) ) +  scale_y_continuous(expand=c(0,0),limits=c(-13.5,13.5) )+ theme(legend.position='none',legend.key = element_blank(), strip.background = element_rect(color="#FFFFFF",fill="#FFFFFF"),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.spacing = unit(0.75, "lines"),panel.background = element_blank(),aspect.ratio=1.0,axis.ticks =element_line(color = "#000000", size = 0.5),axis.text.x = element_text(color="#000000",size=9.5),axis.text.y = element_text(color="#000000",size=9.5) , panel.border = element_rect(color="#000000", fill = NA,size=1),axis.title=element_text(size=9.5) )  +geom_vline(xintercept=log10(9),linetype="dashed",color="#ababab")+geom_hline(yintercept=2.5,linetype="dashed",color="#ababab")
p <-set_panel_size(scatter,width  = unit(3.5, "cm"),height = unit(3.5, "cm"))
grid.newpage()
grid.draw(p)
ggsave(plot=p,"comparison_Christien_4_4.pdf" ,width=2.5,height=2.5)





tarasov = data
tarasov$Christien_2_2 = NULL
tarasov$Christien_2_4 = NULL
tarasov$Christien_4_2 = NULL
tarasov$Christien_4_4 = NULL
tarasov$BFG.Y2H = NULL
TR_2 = na.omit(tarasov)

TR_3 = TR_2[TR_2$Tarassov!=0,]
scatter  = ggplot(data=TR_3,aes(x=log10(BFG.PCA),y=log10(Tarassov),color=reported))+scale_color_manual(breaks = c("0", "1"),values=c("#000000", "#FF0000")) + geom_point(alpha=0.3,size=1) + scale_x_continuous(expand=c(0,0),limits=c(0.45,2.75) ) +  scale_y_continuous(expand=c(0,0),limits=c(2.45,5.05) )+ theme(legend.position='none',legend.key = element_blank(), strip.background = element_rect(color="#FFFFFF",fill="#FFFFFF"),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.spacing = unit(0.75, "lines"),panel.background = element_blank(),aspect.ratio=1.0,axis.ticks =element_line(color = "#000000", size = 0.5),axis.text.x = element_text(color="#000000",size=9.5),axis.text.y = element_text(color="#000000",size=9.5) , panel.border = element_rect(color="#000000", fill = NA,size=1),axis.title=element_text(size=9.5) )  +geom_vline(xintercept=log10(9),linetype="dashed",color="#ababab")+geom_hline(yintercept=log10(23000),linetype="dashed",color="#ababab")
p <-set_panel_size(scatter,width  = unit(3.5, "cm"),height = unit(3.5, "cm"))
grid.newpage()
grid.draw(p)
ggsave(plot=p,"comparison_tarassov.pdf" ,width=2.5,height=2.5)
