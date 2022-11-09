#library(ggplot2)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(MicrobiotaProcess)
source("../src/plotpoint.R")

#######################################################################
##################### different mean ##################################
#######################################################################
data <- read.table("../result/kwfprfnr_1024_norm_mean.60.txt", header=T, sep="\t")
#data2 <- read.table("../result/kwfprfnr_512_norm_mean.txt", header=T, sep="\t")
#data <- rbind(data, data2)
data$methods <- "KW"
data <- data %>% select(c("t","FPR","FNR","methods"))
dat <- read.table("../result/dafprfnr_1024_norm_mean.60.txt", header=T, sep="\t")
#dat2 <- read.table("../result/dafprfnr_512_norm_mean.txt", header=T, sep="\t")
#dat <- rbind(dat, dat2)
dat$methods <- "mp_diff_analysis"
dat <- dat %>% select(c("t","FPR","FNR","methods"))

dt <- read.table("../result/ancombcfprfnr_1024_norm_mean.60.txt", header=T, sep="\t")
#dt2 <- read.table("../result/daglmfprfnr_512_norm_mean.txt", header=T, sep="\t")
#dt <- rbind(dt, dt2)
dt$methods <- "ANCOMBC"
dt <- dt %>% select(c("t", "FPR", "FNR", "methods"))

#udata <- read.table("../result/dafprfnr_unstrict_1024_norm_mean.txt", header=T, sep="\t")
#udata2 <- read.table("../result/dafprfnr_unstrict_512_norm_mean.txt", header=T, sep="\t")
#udata <- rbind(udata, udata2)
#udata$methods <- "diff_analysis(KW,WX)[UN]"
#udata <- udata %>% select(c("t","FPR","FNR","methods"))

#udt <- read.table("../result/daglmfprfnr_unstrict_1024_norm_mean.txt", header=T, sep="\t")
#udt2 <- read.table("../result/daglmfprfnr_unstrict_512_norm_mean.txt", header=T, sep="\t")
#udt <- rbind(udt, udt2)
#udt$methods <- "diff_analysis(glm)[UN]"
#udt <- udt %>% select(c("t","FPR","FNR","methods"))

meseq <- read.table("../result/metafprfnr_1024_norm_mean.60.txt", header=T, sep="\t")
#meseq2 <- read.table("../result/metafprfnr_512_norm_mean.txt", header=T, sep="\t")
#meseq <- rbind(meseq, meseq2)
meseq$methods <- "metagenomeSeq"
meseq <- meseq %>% select(c("t","FPR","FNR","methods"))

anmean <- read.table("../result/ancfprfnr_1024_norm_mean.60.txt", header=T, sep="\t")
anmean$methods <- "ANCOM"
anmean <- anmean %>% select(c("t","FPR","FNR","methods"))

lefse <- read.table("../result/lefsefprfnr_1024_norm_mean.60.txt", header=T, sep="\t")
lefse$methods <- "LEfSe"
lefse <- lefse %>% select(c("t", "FPR", "FNR", "methods"))

lindamean <- read.table('../result/lindafprfnr_1024_norm_mean.60.txt', header=T, sep='\t')
lindamean$methods <- 'LinDA'
lindamean <- lindamean %>% select(c('t', 'FPR', 'FNR', 'methods'))

ziconseqmean <- read.table('../result/zicoseqfprfnr_1024_norm_mean.60.txt', header=T, sep='\t')
ziconseqmean$methods <- 'ZicoSeq'
ziconseqmean <- ziconseqmean %>% select(t, FPR, FNR, methods)

data <- rbind(data, dat, dt, meseq, anmean, lefse, lindamean, ziconseqmean)
data$t <- data$t * 2

dat <- data[,c("FPR", "FNR")] 
sampletmp <- data[,c("methods"),drop=FALSE]

data <- melt(data,id=c("t","methods"), variable.name="type")
head(data)
p1 <- plotpoint(data=data, x="t",y="value",
                group="methods",type="type",
                pointsize=0.8) +
      ylab(NULL) +
      xlab("Different Mean") +
      scale_linetype_manual(values=c("solid", "dashed"))+
      scale_shape_manual(values=c(21, 4))+
      scale_color_manual(values = colors) +
      #scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3")) +
      guides(color = guide_legend(keywidth = 0.8, keyheight = 0.8),
             shape = guide_legend(keywidth = 0.8, keyheight = 0.8))+ 
      theme(legend.position="bottom",
            legend.box="horizontal",
            legend.spacing.x = unit(0.02, "cm"),
            legend.text = element_text(size=8),
            legend.title=element_text(size=9,face="bold"),
            legend.background=element_rect(fill=NA),
            axis.text=element_text(size=7))

p1_box <- ggbox(obj=dat, sampleda=sampletmp, p_textsize = 1.4,
	        indexNames=c("FPR", "FNR"), factorNames="methods", testmethod=NULL) +
           ggtitle("different mean")+
           scale_fill_manual(values = colors) +
           guides(fill=guide_legend(keywidth = 0.8, keyheight = 0.8)) +
           theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                 strip.background=element_rect(colour = NULL,fill = "grey"),
		 legend.position="bottom",legend.box="horizontal",legend.spacing.x = unit(0.02, "cm"),
		 legend.background=element_rect(fill=NA), plot.title=element_text(size=8, face="bold", hjust = 0.5)) 
######################################################################

######################################################################
##################### different sd  ##################################
######################################################################
dasd <- read.table("../result/kwfprfnr_1024_norm_sd.60.txt", header=T, sep="\t")
#dasd2 <- read.table("../result/kwfprfnr_1024_norm_sd.txt", header=T, sep="\t")
#dasd <- rbind(dasd, dasd2)
dasd$methods <- "KW"
dasd <- dasd %>% select(c("t", "FPR", "FNR", "methods"))

datsd <- read.table("../result/dafprfnr_1024_norm_sd.60.txt", header=T, sep="\t")
#datsd2 <- read.table("../result/dafprfnr_1024_norm_sd.txt", header=T, sep="\t")
#datsd <- rbind(datsd, datsd2)
datsd$methods <- "mp_diff_analysis"
datsd <- datsd %>% select(c("t", "FPR", "FNR", "methods"))

dtsd <- read.table("../result/ancombcfprfnr_1024_norm_sd.60.txt", header=T, sep="\t")
#dtsd2 <- read.table("../result/daglmfprfnr_1024_norm_sd.txt", header=T, sep="\t")
#dtsd <- rbind(dtsd, dtsd2)
dtsd$methods <- "ANCOMBC"
dtsd <- dtsd %>% select(c("t", "FPR", "FNR", "methods"))

#udatsd <- read.table("../result/dafprfnr_unstrict_1024_norm_sd.txt", header=T, sep="\t")
#udatsd2 <- read.table("../result/dafprfnr_unstrict_512_norm_sd.txt", header=T, sep="\t")
#udatsd <- rbind(udatsd, udatsd2)
#udatsd$methods <- "diff_analysis(KW,WX)[UN]"
#udatsd <- udatsd %>% select(c("sd", "FPR", "FNR", "methods"))

#udtsd <- read.table("../result/daglmfprfnr_unstrict_1024_norm_sd.txt", header=T, sep="\t")
#udtsd2 <- read.table("../result/daglmfprfnr_unstrict_512_norm_sd.txt", header=T, sep="\t")
#udtsd <- rbind(udtsd, udtsd2)
#udtsd$methods <- "diff_analysis(glm)[UN]"
#udtsd <- udtsd %>% select(c("sd", "FPR", "FNR", "methods"))

meseqsd <- read.table("../result/metafprfnr_1024_norm_sd.60.txt", header=T, sep="\t")
#meseqsd2 <- read.table("../result/metafprfnr_512_norm_sd.txt", header=T, sep="\t")
#meseqsd <- rbind(meseqsd, meseqsd2)
meseqsd$methods <- "metagenomeSeq"
meseqsd <- meseqsd %>% select(c("t", "FPR", "FNR", "methods")) 

ansd <- read.table("../result/ancfprfnr_1024_norm_sd.60.txt", header=T, sep="\t")
ansd$methods <- "ANCOM"
ansd <- ansd %>% select(c("t", "FPR", "FNR", "methods"))

lefsesd <- read.table("../result/lefsefprfnr_1024_norm_sd.60.txt", header=T, sep="\t")
lefsesd$methods <- "LEfSe"
lefsesd <- lefsesd %>% select(c("t", "FPR", "FNR", "methods"))

lindasd <- read.table('../result/lindafprfnr_1024_norm_sd.60.txt', header=T, sep='\t')
lindasd$methods <- 'LinDA'
lindasd <- lindasd %>% select(c('t', 'FPR', 'FNR', 'methods'))

ziconseqsd <- read.table('../result/zicoseqfprfnr_1024_norm_sd.60.txt', header=T, sep='\t')
ziconseqsd$methods <- 'ZicoSeq'
ziconseqsd <- ziconseqsd %>% select(t, FPR, FNR, methods)

dasd <- rbind(dasd, datsd, dtsd, meseqsd, ansd, lefsesd, lindasd, ziconseqsd)

datsd <- dasd[,c("FPR", "FNR")]
sampletmpsd <- dasd[,c("methods"), drop=FALSE]

dasd <- melt(dasd, id=c("t", "methods"), variable.name="type")

p2 <- plotpoint(data=dasd,x="t",y="value",
                group="methods",type="type",
                pointsize=0.8) +
       scale_linetype_manual(values=c("solid", "dashed"))+
       scale_shape_manual(values=c(21, 4))+
       scale_color_manual(values = colors) +
       ylab(NULL)+
       xlab("Standard Deviation")+
       #scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"))+
       theme(legend.position="none",
             axis.text=element_text(size=7))

p2_box <- ggbox(obj=datsd, sampleda=sampletmpsd, p_textsize = 1.4,
                indexNames=c("FPR", "FNR"), factorNames="methods", testmethod=NULL) +
           ggtitle("different standard deviation in class")+
           scale_fill_manual(values = colors) +
           guides(fill=guide_legend(keywidth = 0.8, keyheight = 0.8)) +
           theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                strip.background=element_rect(colour = NULL,fill = "grey"),
		legend.position="none", plot.title=element_text(size=8, face="bold", hjust = 0.5))
#####################################################################

#####################################################################
##################### mean sd       #################################
#####################################################################

damsd <- read.table("../result/kwfprfnr_1024_norm_meansd.60.txt", header=T, sep="\t")
#damsd2 <- read.table("../result/kwfprfnr_1024_norm_meansd.txt", header=T, sep="\t")
#damsd <- rbind(damsd, damsd2)
damsd$methods <- "KW"
damsd <- damsd %>% select(c("t", "FPR", "FNR", "methods"))

datmsd <- read.table("../result/dafprfnr_1024_norm_meansd.60.txt", header=T, sep="\t")
#datmsd2 <- read.table("../result/dafprfnr_1024_norm_meansd.txt", header=T, sep="\t")
#datmsd <- rbind(datmsd, datmsd2)
datmsd$methods <- "mp_diff_analysis"
datmsd <- datmsd %>% select(c("t", "FPR", "FNR", "methods"))

dtmsd <- read.table("../result/ancombcfprfnr_1024_norm_meansd.60.txt", header=T, sep="\t")
#dtmsd2 <- read.table("../result/daglmfprfnr_1024_norm_meansd.txt", header=T, sep="\t")
#dtmsd <- rbind(dtmsd, dtmsd2)
dtmsd$methods <- "ANCOMBC"
dtmsd <- dtmsd %>% select(c("t", "FPR", "FNR", "methods"))

#udatmsd <- read.table("../result/dafprfnr_unstrict_1024_norm_meansd.txt", header=T, sep="\t")
#udatmsd2 <- read.table("../result/dafprfnr_unstrict_512_norm_meansd.txt", header=T, sep="\t")
#udatmsd <- rbind(udatmsd, udatmsd2)
#udatmsd$methods <- "diff_analysis(KW,WX)[UN]"
#udatmsd <- udatmsd %>% select(c("sd", "FPR", "FNR", "methods"))

#udtmsd <- read.table("../result/daglmfprfnr_unstrict_1024_norm_meansd.txt", header=T, sep="\t")
#udtmsd2 <- read.table("../result/daglmfprfnr_unstrict_512_norm_meansd.txt", header=T, sep="\t")
#udtmsd <- rbind(udtmsd, udtmsd2)
#udtmsd$methods <- "diff_analysis(glm)[UN]"
#udtmsd <- udtmsd %>% select(c("sd", "FPR", "FNR", "methods"))

meseqmsd <- read.table("../result/metafprfnr_1024_norm_meansd.60.txt", header=T, sep="\t")
#meseqmsd2 <- read.table("../result/metafprfnr_512_norm_meansd.txt", header=T, sep="\t")
meseqmsd$methods <- "metagenomeSeq"
meseqmsd <- meseqmsd %>% select(c("t", "FPR", "FNR", "methods"))

anmsd <- read.table("../result/ancfprfnr_1024_norm_meansd.60.txt", header=T, sep="\t")
anmsd$methods <- "ANCOM"
anmsd <- anmsd %>% select(c("t", "FPR", "FNR", "methods"))

lefsemsd <- read.table("../result/lefsefprfnr_1024_norm_meansd.60.txt", header=T, sep="\t")
lefsemsd$methods <- "LEfSe"
lefsemsd <- lefsemsd %>% select(c("t", "FPR", "FNR", "methods"))

lindamsd <- read.table('../result/lindafprfnr_1024_norm_meansd.60.txt', header=T, sep='\t')
lindamsd$methods <- 'LinDA'
lindamsd <- lindamsd %>% select(c('t', 'FPR', 'FNR', 'methods'))

ziconseqmsd <- read.table('../result/zicoseqfprfnr_1024_norm_meansd.60.txt', header=T, sep='\t')
ziconseqmsd$methods <- 'ZicoSeq'
ziconseqmsd <- ziconseqmsd %>% select(t, FPR, FNR, methods)

damsd <- rbind(damsd, datmsd, dtmsd, meseqmsd, anmsd, lefsemsd, lindamsd, ziconseqmsd)

datmsd <- damsd[,c("FPR", "FNR")]
sampletmpmsd <- damsd[,c("methods"), drop=FALSE]

damsd <- melt(damsd, id=c("t", "methods"), variable.name="type")
p3 <- plotpoint(data=damsd,x="t",y="value",
                group="methods",type="type",
                pointsize=0.8) +
      scale_linetype_manual(values=c("solid", "dashed"))+
      scale_shape_manual(values=c(21, 4)) +
      scale_color_manual(values = colors) +
      ylab(NULL) + 
      xlab("Standard Deviation") + 
      #scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3"))+
      theme(legend.position="none",
            axis.text=element_text(size=7))

p3_box <- ggbox(obj=datmsd, sampleda=sampletmpmsd, p_textsize=1.5,
	        indexNames=c("FPR", "FNR"), factorNames="methods", testmethod=NULL) +
           ggtitle("different standard deviation in subclass")+
           scale_fill_manual(values = colors) +
           guides(fill=guide_legend(keywidth = 0.8, keyheight = 0.8)) +
           theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                  strip.background=element_rect(colour = NULL,fill = "grey"),
		  legend.position="none", 
                  plot.title=element_text(size=8, face="bold", hjust = 0.5))

####################################################################
############### figure assemble #############
####################################################################
legend <- get_legend(p1)
p1 <- p1 + theme(legend.position="none")
p <- ggdraw(plot_grid(plot_grid(p1, p2, p3, labels = c("A", "B", "C"), align = "h", nrow=1),
     plot_grid(NULL,legend,NULL, ncol=3, align="h"), ncol=1, rel_heights=c(1,  0.2)))

svg("../figures/norm_compare_point.60.svg", width=7.5, height=3)
p
dev.off()


#########################################################
################ BOXPLOT FIGURE ##############
#########################################################

legendbox <- get_legend(p1_box)
p1_box <- p1_box + theme(legend.position="none")
p_box <- ggdraw(plot_grid(plot_grid(p1_box, p2_box, p3_box, labels = c("A", "B", "C"), align = "h", nrow=1),
			  plot_grid(NULL, legendbox, NULL, ncol=3, align="h"), ncol=1, rel_heights=c(1,  0.2)))

svg("../figures/norm_compare_boxplot.60.svg", width=8, height=4)
p_box
dev.off()
