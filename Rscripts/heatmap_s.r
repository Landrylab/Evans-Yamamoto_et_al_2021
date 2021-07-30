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

heatmap = ggplot(data )+   
    geom_tile(aes(prey,bait,fill = log10(s)))  +
    theme_bw()   + 
    theme( panel.spacing = unit(1, "lines"),
           legend.position='bottom',
           legend.key = element_blank(), 
           strip.background = element_rect(colour="#FFFFFF", color="#FFFFFF",fill="#FFFFFF"),
           panel.border = element_rect(fill=NA, colour = "black", size=1),
           panel.grid.minor = element_blank(),
           panel.grid.major = element_blank(),
           panel.background = element_blank(),
           aspect.ratio=1.0,
           axis.ticks = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_blank()    ,
           axis.title =  element_text(size=9.5)       )+
    scale_x_discrete(expand=c(0,0),position = "top") + 
    scale_y_discrete(expand=c(0,0),limits=rev) +
    scale_fill_gradientn(colors=c("#FFFFFF","#FFFFFF","#FFFFFF","#ADD8E6","#001144","#000000"),limits=c(ceiling(log10(min(data$s))-2),ceiling(log10(max(data$s))+0.1)),na.value = "#000000")+
    xlab("DHFR F[1,2] barcodes")+
    ylab("DHFR F[3] barcodes")+facet_grid(Replicate~Condition)


p <-set_panel_size(heatmap,width  = unit(width_cm, "cm"),height = unit(height_cm, "cm"))
grid.newpage()
grid.draw(p)
ggsave(plot=p,out,height=p_height,width=p_width,units = "cm")


th = theme(aspect.ratio=0.5,
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(fill=NA, colour = "black", size=1),
    axis.title=element_text(size=9.5,color="#000000"), 
    legend.key=element_blank(),
    strip.background = element_rect(colour="#FFFFFF", color="#FFFFFF",fill="#FFFFFF"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(colour = "#000000", size = 0.5),
    axis.text = element_text(size=9.5,color="#000000"),
    strip.text.y = element_text(size=9.5,color="#000000")
    )

histogram  = ggplot(data) +   
    geom_histogram(aes(x=log10(s),y=log10(..count..+1)),bins=30,fill="#2f4893",color="#999999",size=0.1) + 
    theme_bw()  +
    xlab(bquote(~Log["10"]~ '(Raw interaction signal)'))+
    ylab(bquote(~Log["10"]~ '(Barcode counts + 1)'))+
    facet_grid(Replicate~Condition,scales = "free")

histogram2 = histogram+  
    scale_y_continuous(expand=c(0,0),
                       limits=c( 0-(ceiling(log10(max(ggplot_build(histogram)$data[[1]]$count))))*0.01,
                                 ceiling(log10(max(ggplot_build(histogram)$data[[1]]$count)))*1.1) )+
    scale_x_continuous(expand=c(0,0),limits=c(round(log10(min(data$s))-1,digit=0)-0.1,round(log10(max(data$s))+1,digits=0)+0.1),breaks=seq(round(log10(min(data$s))-1,digit=0),round(log10(max(data$s))+1,digits=0)+0.05,1))  +th
    


hs   = gsub("heatmap", "histogram", out)

pl = histogram2
p <-set_panel_size(pl,width  = unit(width_cm, "cm"),height  = unit(height_cm/2.5, "cm"))
grid.newpage()
grid.draw(p)
ggsave(plot=p,hs,width=p_width,height=p_height,units = "cm")
