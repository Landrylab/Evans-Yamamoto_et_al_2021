#!/usr/bin/env Rscript
if(require("ggplot2")){
    print("ggplot2 is loaded correctly")
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
    print("scales is loaded correctly")
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
    print("cowplot is loaded correctly")
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
    print("grid is loaded correctly")
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
    print("gridExtra is loaded correctly")
} else {
    print("trying to install gridExtra...")
    install.packages("gridExtra")
    if(require("gridExtra")){
        print("gridExtra installed and loaded")
    } else {
        stop("could not install gridExtra")
    }
}

if(require("stringr")){
    print("stringr is loaded correctly")
} else {
    print("trying to install stringr...")
    install.packages("stringr")
    if(require("stringr")){
        print("stringr installed and loaded")
    } else {
        stop("could not install stringr")
    }
}

library(ggplot2)
library(scales) # for muted function
library(cowplot)
library(gridExtra)
library(grid)
library(stringr)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

args = commandArgs(trailingOnly=TRUE)

hap = read.csv(file=args[1],header=TRUE)
hap_freq  = ggplot(hap) +   geom_histogram(aes(x=log10(Cells)),bins=30,fill="#2f4893",color="#999999",size=0.1) + theme_bw()  +scale_x_continuous(expand=c(0,0),breaks=seq(0,max(Cells),1)) +  scale_y_continuous(expand=c(0,0))  + theme(panel.spacing = unit(2, "lines"),legend.position='bottom',legend.key = element_blank(), strip.background = element_rect(colour="#FFFFFF", color="#FFFFFF",fill="#FFFFFF"),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_blank(),aspect.ratio=0.8,axis.ticks = element_line(colour = "#000000", size = 0.25),axis.text.x = element_text(color="#000000",size=9.5),axis.text.y = element_text(color="#000000",size=9.5),strip.text.y = element_text(size = 9.5, colour = "#000000", angle = 270),strip.text.x = element_text(size = 9.5, colour = "#000000"))+xlab(bquote(~Log[10]~ '('  ~italic(Cells)[Haploid]~ ')'))+ylab("Number of strains")+facet_wrap(~Pool,nrow=2)

P <-set_panel_size(hap,width  = unit(2.5, "cm"),height = unit(5, "cm"))
grid.newpage()
grid.draw(P)

hap_freq_png= paste(str_sub(args[1], start = 1, end = -5),".png")

ggsave(hap_freq_png,plot=P,height=4,limitsize=FALSE)

hap_cov = read.csv(file=args[2],header=TRUE)
plot = ggplot(data=data,aes(Threshold,Percentage),group=Pool) + geom_line(aes( color=Pool))+scale_color_manual(breaks = c("x", "y"),values=c("#ababab", "#2f4893")) + theme_bw()+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0),limits=c(0,105))+theme(legend.position="bottom")+theme(legend.key = element_blank(), strip.background = element_rect(colour="#FFFFFF", color="#FFFFFF",fill="#FFFFFF"),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_blank(),aspect.ratio=1.0,axis.ticks = element_line(colour = "#000000", size = 0.25),axis.text.x = element_text(color="#000000",size=9.5),axis.text.y = element_text(color="#000000",size=10),strip.text.y = element_text(size = 9.5, colour = "#000000", angle = 270),strip.text.x = element_text(size = 9.5, colour = "#000000"))

P <-set_panel_size(hap_cov,width  = unit(2.5, "cm"),height = unit(2.5, "cm"))
grid.newpage()
grid.draw(P)
hap_cov_png = paste(str_sub(args[2], start = 1, end = -5),".png")
ggsave(hap_cov_png,plot=P,height=4,limitsize=FALSE)
