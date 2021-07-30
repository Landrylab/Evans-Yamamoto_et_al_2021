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


curve_data = data[data$TH==1,]
plot = ggplot(data=curve_data) +
  geom_line(aes(log10(Rank),Precision),color="#FF0000")+ 
  geom_line(aes(log10(Rank),Recall),color="#0000FF")+ 
  geom_line(aes(log10(Rank),MCC),color="#000000")+ 
  theme_bw()+
  scale_x_continuous(expand=c(0,0),limits=c(-0.05,max(log10(data$Rank)))) + 
  scale_y_continuous(expand=c(0,0),limits=c(-0.01,1.05))+
  theme(legend.position="none",
        panel.border = element_rect(size=1), 
        legend.key = element_blank(), 
        strip.background = element_rect(colour="#FFFFFF", color="#FFFFFF",fill="#FFFFFF"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        aspect.ratio=1.0,
        axis.ticks = element_line(colour = "#000000", size = 0.5),
        axis.text.x = element_text(color="#000000",size=9.5),
        axis.text.y = element_text(color="#000000",size=9.5),
        strip.text.y = element_text(size = 9.5, colour = "#000000", angle = 270),
        strip.text.x = element_text(size = 9.5, colour = "#000000"))+
  ylab("Value")+
  xlab(bquote(~Log["10"]~ 'Rank'))
texts = ggplot()+annotate("text", x = 0,y =2, label = "Precision", angle=0, size=3, color='#FF0000', face="bold")+
  annotate("text", x = 0,y =1, label = "Recall", angle=0, size=3, color='#0000FF', face="bold")+
  annotate("text", x = 0,y =0, label = "MCC", angle=0, size=3, color='#000000', face="bold")+
  theme(legend.position="none",
        panel.border = element_blank(), 
        legend.key = element_blank(), 
        strip.background = element_rect(colour="#FFFFFF", color="#FFFFFF",fill="#FFFFFF"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+scale_x_continuous(limits=c(-0.5,1))

p <-set_panel_size(plot,width  = unit(3.5, "cm"),height = unit(3.5, "cm"))
p2 <-set_panel_size(texts+theme(plot.margin = unit(c(0,0,0,0), "cm")),width  = unit(2, "cm"),height = unit(3.5, "cm"))
P = grid.draw(plot_grid(p, p2, align = "h", rel_widths= c(1, 0.25),nrow = 1))

grid.newpage()
grid.draw(P)
nm   = gsub(".pdf","Precision_Recall_curve.pdf", out)
ggsave(plot=P,nm ,width=10,height=5,units = "cm")
print(paste("Exported: ",nm))

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

histogram  = ggplot(curve_data) +   
  geom_histogram(aes(x=log2(n),y=log2(..count..+1)),bins=30,fill="#2f4893",color="#999999",size=0.1) + 
  theme_bw()  +
  xlab("Number of replicates \nper protein pair")+
  ylab("Protein pairs")

histogram2 = histogram+  
  scale_y_continuous(expand=c(0,0),
                     limits=c( 0-(ceiling(log2(max(ggplot_build(histogram)$data[[1]]$count))))*0.01,
                               ceiling(log2(max(ggplot_build(histogram)$data[[1]]$count)))*1.1) )+
  scale_x_continuous(expand=c(0,0),limits=c(round(log2(min(data$n))-1,digit=0)-0.1,round(log2(max(data$n))+1,digits=0)+0.1),breaks=seq(round(log2(min(data$n))-1,digit=0),round(log2(max(data$n))+1,digits=0)+0.05,1))  +th

p <-set_panel_size(histogram2,width  = unit(3.5, "cm"),height = unit(3.5, "cm"))

grid.newpage()
grid.draw(p)
nm   = gsub(".pdf","replicates.pdf", out)
ggsave(plot=p,nm ,width=5,height=5,units = "cm")

print(paste("Exported: ",nm))

best_rank = data[data$MCC == max(data$MCC),]$Rank[1]
print(paste("Optimal rank threashold for calling PPIs are:",best_rank))
plot_area = min(ceiling(data[data$MCC == max(data$MCC),]$Rank[1]*1.2),max(data$Rank))
print(paste("Plotting scores for the top ",plot_area," PPIs."))

data_ex = data[data$Rank < plot_area,]
data_ex$th = as.factor(data_ex$TH)

plot = ggplot(data=data_ex,aes(Rank-0.5,log10(score))) + 
  geom_jitter(aes(colour = th),size=0.25,alpha=0.5)+  
  scale_color_manual(breaks = c("0", "1"),values=c("#ababab", "#FF0000"))+
  theme(
    legend.position='none',
    legend.key = element_blank(), 
    strip.background = element_rect(color="#FFFFFF",fill="#FFFFFF"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(0.75, "lines"),
    panel.background = element_blank(),
    aspect.ratio= 1.5/plot_area,
    axis.ticks =element_line(color = "#000000", size = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color="#000000",size=9.5),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9.5,color="#000000"),
    axis.line.x=element_line(size=0.5),
    axis.line.y =element_line(size=0.5)
    ) +
  scale_x_continuous(expand=c(0,0),limits=c(0,plot_area),breaks=c())+
  xlab("")+
  ylab(bquote(~Log["10"]~ '(PPI score)'))+
  scale_y_continuous(expand=c(0,0),limits=c(-0.5,ceiling(log10(max(data$score)))*1.1),breaks=c(0,1,2,3,4,5,6,7))+
  annotate("segment", x = best_rank+0.025-0.25, xend = best_rank-1.015-0.25, y = ceiling(log10(max(data$score)))-1, yend =ceiling(log10(max(data$score)))-1, colour = "#000000", size=0.25, alpha=1)+
  annotate("text", x = best_rank-1,y = ceiling(log10(max(data$score))), label = "Optimal MCC", angle=0, size=3, color='#000000', face="bold")+
  annotate("segment", x = best_rank-0.25, xend = best_rank-0.25, y = ceiling(log10(max(data$score)))-2, yend =  ceiling(log10(max(data$score)))-1, colour = "#000000", size=0.25, alpha=1)+    
  annotate("segment", x = best_rank-1-0.25, xend = best_rank-0.5-0.25, y = ceiling(log10(max(data$score)))-1, yend =  ceiling(log10(max(data$score)))-1.2, colour = "#000000", size=0.25, alpha=1)+
  annotate("segment", x = best_rank-1-0.25, xend = best_rank-0.5-0.25, y = ceiling(log10(max(data$score)))-1, yend =  ceiling(log10(max(data$score)))-0.8, colour = "#000000", size=0.25, alpha=1)+theme(plot.margin = unit(c(0,0,0,0), "cm"))


data_ex$Method  = NULL
data_ex$Selection  = NULL
data_ex$n              = NULL
data_ex$Norm_meth  = NULL
data_ex$Th_method      = NULL
data_ex$s  = NULL
data_ex$Precision  = NULL
data_ex$Recall         = NULL
data_ex$MCC  = NULL
data_ex$score       = NULL
data_ex$TH       = NULL
data_ex$th       = NULL
data_ex$score_TH       = NULL
data_ex$X1_Two.hybrid = NULL
data_ex$X2_PCA = NULL
data_ex$X4_ALL = NULL

data_ex_hm = melt(data_ex,id.var=c("Rank","ORF_pair"))

data_ex_hm2 = data_ex_hm[data_ex_hm$variable!= "X4_ALL",]

data_ex_hm2$Reported_Binary_PPIs = data_ex_hm2$value
 
hm =ggplot(data_ex_hm2)+ 
  geom_tile(aes(Rank,variable,fill=Reported_Binary_PPIs))+ 
  theme_bw()   + 
  theme(panel.border= element_rect(size=1), 
        legend.position='none',
        legend.key = element_blank(), 
        strip.background = element_rect( color="#FFFFFF",fill="#FFFFFF"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        aspect.ratio=0.05,
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()  ,
        axis.title.x = element_text(size=9.5,color="#000000"),
        axis.title.y = element_blank()     )+ 
  scale_x_discrete(expand=c(0,0),limits=c(0.5,plot_area+0.5)) + 
  scale_y_discrete(expand=c(0,0),limits=rev) +
  scale_fill_gradientn(colors=c("#FFFFFF","#ababab","#000000","#000000"),limits=c(0,max(data_ex_hm2$Reported_Binary_PPIs)),na.value = "#FFFFFF")+
  xlab("Protein pairs ordered by rank")+
  ylab("")
Citations = seq(0,max(data_ex_hm2$Reported_Binary_PPIs)+2,1)
legend = as.data.frame(Citations)
legend$Y = 1
leg =ggplot(legend)+ 
  geom_tile(aes(Citations,Y,fill=Citations))+ 
  theme_bw()   + 
  theme(panel.border= element_rect(size=1), 
        legend.position='none',
        legend.key = element_blank(), 
        strip.background = element_rect( color="#FFFFFF",fill="#FFFFFF"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        aspect.ratio=0.05,
        axis.ticks.y = element_blank(),
        axis.ticks.length.x.bottom = unit(.15, "cm"),
        axis.text.x = element_text(color="#000000",size=9.5),
        axis.text.y = element_blank()  ,
        axis.title.x = element_text(size=9.5,color="#000000"),
        axis.title.y = element_blank())+ 
  scale_x_continuous(expand=c(0,0),breaks = seq(0,max(data_ex_hm2$Reported_Binary_PPIs)+2,2)) + 
  scale_y_discrete(expand=c(0,0)) +
  scale_fill_gradientn(colors=c("#FFFFFF","#ababab","#000000","#000000"),limits=c(0,max(data_ex_hm2$Reported_Binary_PPIs)+2),na.value = "#FFFFFF")+xlab("Binary PPIs reported in BioGRID")

P1 <-set_panel_size(plot+theme(plot.margin = unit(c(1,1,0.5,1), "cm")),width  = unit(width_cm*plot_area/10, "cm"),height = unit(3, "cm"))
P2 <-set_panel_size(hm+theme(plot.margin = unit(c(0,0,1,0), "cm")),width  = unit(width_cm*plot_area/10, "cm"),height = unit(0.3, "cm"))
P3 <-set_panel_size(leg+theme(plot.margin = unit(c(0,1,1,1), "cm")),width  = unit(3, "cm"),height = unit(0.3, "cm"))
grid.newpage()
P = grid.draw(plot_grid(P1, P2, align = "v", rel_heights  = c(1, 0.1),ncol = 1))
nm   = gsub(".pdf","rankplot.pdf", out)
ggsave(plot=P,nm ,width=(p_width*plot_area/100)*4,height=p_height,units = "cm")
P = grid.draw(plot_grid(P3))
nm   = gsub(".pdf","rankplot_legend.pdf", out)
ggsave(plot=P,nm ,width=10,height=2.5,units = "cm")
print(paste("Exported: ",nm))