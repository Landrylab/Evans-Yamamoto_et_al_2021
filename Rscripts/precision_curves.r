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



group.colors = c("Precision" = "#FF0000","Recall" = "#0000FF","MCC" = "#000000")

#Delete except Rank, Precision, Recall , MCC

#Melt data for plot
#dat_m = melt(dat,id.vars=c("Rank"))


plot = ggplot(data=data) + geom_line(aes(log10(Rank),Precision),color="#FF0000")+ geom_line(aes(log10(Rank),Recall),color="#0000FF")+ geom_line(aes(log10(Rank),MCC),color="#000000")+ theme_bw()+scale_x_continuous(expand=c(0,0),limits=c(-0.05,max(log10(data$Rank)))) + scale_y_continuous(expand=c(0,0),limits=c(-0.01,1.05))+theme(legend.position="none",panel.border = element_rect(size=1), legend.key = element_blank(), strip.background = element_rect(colour="#FFFFFF", color="#FFFFFF",fill="#FFFFFF"),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),panel.background = element_blank(),aspect.ratio=1.0,axis.ticks = element_line(colour = "#000000", size = 0.5),axis.text.x = element_text(color="#000000",size=9.5),axis.text.y = element_text(color="#000000",size=9.5),strip.text.y = element_text(size = 9.5, colour = "#000000", angle = 270),strip.text.x = element_text(size = 9.5, colour = "#000000"))


p <-set_panel_size(plot,width  = unit(3.5, "cm"),height = unit(3.5, "cm"))
grid.newpage()
grid.draw(p)
name = paste("./",gsub(".csv", ".pdf", args[1]))
nm   = gsub(" ", "", name)
ggsave(plot=p,nm ,width=2.5,height=2.5)
