library(MicrobiotaProcess)
#suppressPackageStartupMessages(library(tidyverse))
#suppressPackageStartupMessages(library(metagenomeSeq))
#library(ANCOMBC)
library(tibble)
library(dplyr)
library(tidyr)
seednum <- 1024
set.seed(seednum)
sampleda <- data.frame(group=c(rep("Case", 50), rep("Control", 50)),
					   subgroup=c(rep("Case1",25),rep("Case2",25), 
								  rep("Control1",25), rep("Control2", 25)))
rownames(sampleda) <- paste0('S', seq_len(nrow(sampleda)))
source('./src/getFPR_FNR.R')
source("./src/rel_ancom.R")
source('./src/daa_funs.R')
kwfprfnr <- list()
dafprfnr <- list()
#dafprfnr2 <- list()
#daglmfprfnr2 <- list() 
metafprfnr <- list()
ancfprfnr <- list()
ancombcfprfnr <- list()
lefsefprfnr <- list()
lindafprfnr <- list()
zicoseqfprfnr <- list()
for (i in seq(0.01, 3, by=0.01)){
    print(i)
    negda <- data.frame(matrix(rlnorm(500*100, meanlog=10,sdlog=1), nrow=100, byrow=TRUE))
    colnames(negda) <- paste0("N", seq_len(ncol(negda)))
    poida1 <- data.frame(matrix(rlnorm(500*50, meanlog=10-i, sdlog=1), 
								nrow=50, byrow=TRUE))
    poida2 <- data.frame(matrix(rlnorm(500*50, meanlog=10+i, sdlog=1),
								nrow=50, byrow=TRUE))
    poida <- rbind(poida1,poida2)
    colnames(poida) <- paste0("P", seq_len(ncol(poida)))
    totalda <- cbind(negda,poida)
    rownames(totalda) <- paste0('S', seq_len(nrow(totalda)))
    ref <- data.frame(f=colnames(totalda),ref=c(rep("N",500),rep("P",500)))
    ref$ref <- factor(ref$ref, levels=c("P","N")) 
#    ######################### kruskal.test ###############################
#    kwres <- multi_compare(fun="kruskal.test", data=merge(totalda, sampleda, by=0),
#    					   feature=colnames(totalda),factorNames="group") 
#    resp <- unlist(lapply(kwres,function(x)x$p.value)) 
#    kwres <- data.frame(f=colnames(totalda), pvalue=resp) 
#    kwres <- as.vector(kwres[kwres$pvalue<=0.05, 1])
#    kwres <- confusiontab(kwres, ref)
#    kwfprfnr[[as.character(i)]] <-  kwres
#    ######################### ANCOM-II ###################################
#    anres <- tryCatch(rel_ancom(OTUdat=cbind(Sample.ID=rownames(totalda), totalda), 
#                                Vardat=cbind(Sample.ID=rownames(sampleda), sampleda), 
#                                main.var="group"),
#                      error=function(e)NULL)
#    anres <- tryCatch(confusiontab(anres$detected_microbes, ref),error=function(e)confusiontab(NULL, ref))
#    ancfprfnr[[as.character(i)]] <- anres
#    ######################### diff_analysis(KW,WX) ########################
    mpse <- MPSE(assays=list(Abundance=t(totalda)), colData=sampleda)
#    res <- tryCatch(mp_diff_analysis(mpse,
#                                 .group = group,
#                                 .sec.group = subgroup,
#                                 .abundance = Abundance,
#                                 force = TRUE,
#                                 relative = FALSE,
#                                 ldascore = 0,
#                                 action = 'only'
#                                 ),error=function(e)NULL)
#    res <- tryCatch(confusiontab(res[!is.na(res$group),'f'], ref),
#    				error=function(e)confusiontab(NULL,ref))
#    dafprfnr[[as.character(i)]] <- res
#    #########################  ANCOMBC            ########################
#    ps <- as.phyloseq(mpse)
#    res.bc <- tryCatch(ancombc(ps, formula = 'group', group='group', global=TRUE),error=function(e)NULL)
#    res.bc <- tryCatch(confusiontab(rownames(res.bc$res$diff_abn[res.bc$res$diff_abn[[1]],,drop=F]), ref),error=function(e)confusiontab(NULL,ref))
#    ancombcfprfnr[[as.character(i)]] <- res.bc
#    #########################  metagenomeSeq      ########################
#    saa <- AnnotatedDataFrame(sampleda)
#    nMR <- newMRexperiment(t(totalda), phenoData=saa)
#    nMR <- cumNorm(nMR,p=0.5)
#    nMRpd <- pData(nMR)
#    mod <- model.matrix(~1 + group, data = nMRpd)
#    #res5 <- fitFeatureModel(nMR, mod)
#    res5 <- fitZig(nMR, mod)
#    res5 <- MRcoefs(res5)
#    res5 <- tryCatch(confusiontab(rownames(res5[res5$pvalues<=0.05,]),ref),
#                                  error=function(e)confusiontab(NULL, ref))
#    metafprfnr[[as.character(i)]] <- res5
#    ##################### LEfSe ###############################
#    lefseda <- merge(sampleda, totalda, by=0)
#    lefseda$Row.names <- NULL
#    lefseda <- data.frame(t(lefseda), check.names=FALSE) %>% rownames_to_column(var="feature")
#    tmpfile1 <- paste0(tempfile(), i)
#    tmpfile2 <- paste0(tmpfile1, ".format")
#    outfile <- paste0(tmpfile1, ".out")
#    write.table(lefseda, tmpfile1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
#    CMD1 <- paste("format_input.py", tmpfile1, tmpfile2, "-c 1 -s 2 -o -1", sep=" ")
#    CMD2 <- paste("run_lefse.py", tmpfile2, outfile, "--min_c 3 -f 1 -b 1 -l 0 -y 1", sep=" ")
#    system(CMD1)
#    system(CMD2)
#    lefseout <- read.table(outfile, sep="\t", header=F, row.names=1)
#    flags <- suppressWarnings(ifelse(is.na(as.numeric(lefseout$V5) <= 0.05),FALSE,TRUE))
#    lefseres <- tryCatch(confusiontab(rownames(lefseout[flags,]), ref),error=function(e)confusiontab(NULL, ref))
#    lefsefprfnr[[as.character(i)]] <- lefseres
    ###################### LinDA ##############################
    linda.res <- obtain_linda_out(mpse)
    linda.res <- run_linda(linda.res)
    lindafprfnr[[as.character(i)]] <- tryCatch(confusiontab(linda.res, ref),error=function(e)confusiontab(NULL, ref))

    ##################### ZiCoSeq ############################
    zicoseq.res <- obtain_zicoseq_out(mpse)
    zicoseq.res <- run_zicoseq(zicoseq.res)
    zicoseqfprfnr[[as.character(i)]] <- tryCatch(confusiontab(zicoseq.res, ref),error=function(e)confusiontab(NULL, ref))
}
#
#
#kwfprfnr <- do.call("rbind", kwfprfnr) %>% rownames_to_column(var="t")
#write.table(kwfprfnr,paste0("./result/kwfprfnr_", seednum,"_rlnorm_mean.txt"), sep="\t",
#			quote=FALSE, row.names=FALSE, col.names=TRUE)
#
#
#ancfprfnr <- do.call("rbind", ancfprfnr) %>% rownames_to_column(var="t")
#write.table(ancfprfnr, paste0("./result/ancfprfnr_", seednum, "_rlnorm_mean.txt"), sep="\t",
#			quote=FALSE, row.names=FALSE, col.names=TRUE)
#
#ancombcfprfnr <- do.call("rbind", ancombcfprfnr) %>% rownames_to_column(var="t")
#write.table(ancombcfprfnr, paste0("./result/ancombcfprfnr_", seednum, "_rlnorm_mean.txt"), sep="\t",
#            quote=FALSE, row.names=FALSE, col.names=TRUE)    
#
#dafprfnr <- do.call("rbind", dafprfnr) %>% rownames_to_column(var="t")
#write.table(dafprfnr,paste0("./result/dafprfnr_", seednum,"_rlnorm_mean.txt"), sep="\t",
#			quote=FALSE, row.names=FALSE, col.names=TRUE)
#
#metafprfnr <- do.call("rbind", metafprfnr) %>% rownames_to_column(var="t")
#write.table(metafprfnr, paste0("./result/metafprfnr_", seednum, "_rlnorm_mean.txt"), sep="\t",
#			quote=FALSE, row.names=FALSE, col.names=TRUE)
#
#lefsefprfnr <- do.call("rbind", lefsefprfnr) %>% rownames_to_column(var="t")
#write.table(lefsefprfnr, paste0("./result/lefsefprfnr_", seednum, "_rlnorm_mean.txt"), sep="\t",
#	                            quote=FALSE, row.names=FALSE, col.names=TRUE)

lindafprfnr <- do.call('rbind', lindafprfnr) %>% rownames_to_column(var='t')
write.table(lindafprfnr, paste0("./result/lindafprfnr_", seednum, "_rlnorm_mean.txt"), sep="\t",
            quote=FALSE, row.names=FALSE, col.names=TRUE)


zicoseqfprfnr <- do.call('rbind', zicoseqfprfnr) %>% rownames_to_column(var='t')
write.table(zicoseqfprfnr, paste0("./result/zicoseqfprfnr_", seednum, "_rlnorm_mean.txt"), sep="\t",
            quote=FALSE, row.names=FALSE, col.names=TRUE)
