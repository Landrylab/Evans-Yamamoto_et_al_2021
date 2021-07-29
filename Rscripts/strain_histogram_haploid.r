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

if(require("egg")){
  print("egg is loaded correctly")
} else {
  print("trying to install egg...")
  install.packages("egg")
  if(require("egg")){
    print("egg")
  } else {
    stop("egg")
  }
}

library(ggplot2)
library(scales) # for muted function
library(cowplot)
library(gridExtra)
library(grid)
library(stringr)
library(egg)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


args = commandArgs(trailingOnly=TRUE)


data      = read.csv(file=args[1],header=TRUE)
out       = args[2]
width_cm  = args[3]
height_cm = args[4]
p_height  = args[5]
p_width   = args[6]

P  = ggplot(data) +  
  geom_histogram(aes(x=log10(F)),bins=25,fill="#2f4893",color="#999999",size=0.1) +
  theme_bw()  +
  scale_x_continuous( expand=c(0,0),
                      limits=c(round(log10(min(data$F))-2,digit=0)-0.1,round(log10(max(data$F))+2,digits=0)+0.1),
                      breaks=seq(round(log10(min(data$F))-2,digits=0),round(log10(max(data$F))+2,digits=0),1) 
                      )+  
  theme( panel.spacing = unit(2, "lines"),
         legend.key = element_blank(),
         strip.background = element_rect( color="#FFFFFF",fill="#FFFFFF"),
         panel.border = element_rect(size=1.0),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         panel.background = element_blank(),
         aspect.ratio=0.8,
         axis.ticks = element_line(color = "#000000", size = 0.5),
         axis.text.x = element_text(color="#000000",size=9.5),
         axis.text.y = element_text(color="#000000",size=9.5),
         strip.text.y = element_text(size = 9.5, color = "#000000", angle = 270), 
         strip.text.x = element_text(size = 9.5, color = "#000000")
         )+
  xlab(bquote(~Log["10"]~ '('  ~italic("F")["Marginal barcode"]~ ')'))+
  ylab("Number of barcodes")+
  facet_grid(Screening~Ori)
P2 = P +scale_y_continuous(expand=c(0,0),limits=c(0-(ceiling(max(ggplot_build(P)$data[[1]]$count))*0.01),ceiling(max(ggplot_build(P)$data[[1]]$count)*1.1))) 
ggsave(out,plot=P,height=4,limitsize=FALSE)

p <-set_panel_size(P2,width  = unit(width_cm, "cm"),height = unit(height_cm, "cm"))
grid.newpage()
grid.draw(p)
print(out)
ggsave(plot=p,out,height=as.numeric(p_height),width=as.numeric(p_width),units = "cm")
