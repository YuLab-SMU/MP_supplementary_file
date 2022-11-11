library(MicrobiotaProcess)
library(tibble)
library(MicrobiomeStat)
library(GUniFrac)
library(ANCOMBC)
library(metagenomeSeq)
library(matrixStats)


## MicrobiotaProcess
run_mpse_daa <- function(mpse, p.adj = TRUE, method = 'bonferroni'){
    if (p.adj){
        filter.p <- 'fdr'
    }else{
        filter.p <- 'pvalue'
    }
    sample.da <- mpse |> mp_extract_sample() |> tibble::column_to_rownames(var='Sample')
    if (ncol(sample.da) > 1){
        group <- colnames(sample.da)[1]
        adj.group <- colnames(sample.da)[2]
    }else{
        group <- colnames(sample.da)
        adj.group <- NULL
    }
    mpse.diff <- tryCatch(mpse |> mp_filter_taxa(min.abun=1, min.prop=.1, include.lowest=T) |>
        mp_diff_analysis(
          .abundance = "total", 
          force = T,
          relative = F,
          .group = !!as.symbol(group),
          .sec.group = {{adj.group}},
          normalization = 1e8,
          p.adjust = method,
          filter.p = filter.p
        ),
        error=function(e)NULL
    )
    if (is.null(mpse.diff)){
        return(NA)
    }else{
       mpse.diff |>
        mp_extract_feature() |>
        dplyr::filter(!is.na(!!as.symbol(paste0("Sign_", group)))) |>
        dplyr::pull("OTU") -> mp.diff
    }
    return(mp.diff)
}


## linda
obtain_linda_out <- function(mpse, method = 'bonferroni'){
    ps <- mpse |> as.phyloseq(.abundance = 'Abundance')
    sample.da <- mpse |> mp_extract_sample() |> tibble::column_to_rownames(var='Sample')
    if (ncol(sample.da) > 1){
        group <- colnames(sample.da)[1]
        adj.group <- colnames(sample.da)[2]
    }else{
        group <- colnames(sample.da)
        adj.group <- NULL
    }    
    linda.diff <- tryCatch(linda(phyloseq.obj=ps, 
                        formula=paste0("~", paste0(c(group, adj.group), collapse='+')),
                        prev.filter=0.1,
                        p.adj.method = method),
                        error = function(e)NULL)
    if (is.null(linda.diff)){
        return(NULL)
    }
    ind1 <- which(grepl(group, names(linda.diff$output)))
    return(linda.diff$output[[ind1]])
}

run_linda <- function(linda.diff, p.adj = TRUE){
    #linda.diff <- obtain_linda_out(mpse)
    if (p.adj){
        p.adj <- 'padj'
    }else{
        p.adj <- 'pvalue'
    } 
    linda.diff[linda.diff[[p.adj]] <= 0.05,] |>
        rownames() -> linda.diff
    return(linda.diff)
    
}

## ZicoSeq
obtain_zicoseq_out <- function(mpse, method='bonferroni'){
    otu.da <- mpse |> mp_extract_assays(.abundance = "total") %>% as.matrix()
    sample.da <- mpse |> mp_extract_sample() |> tibble::column_to_rownames(var='Sample')
    if (ncol(sample.da) > 1){
        group <- colnames(sample.da)[1]
        adj.group <- colnames(sample.da)[2]
    }else{
        group <- colnames(sample.da)
        adj.group <- NULL
    }    
    feature.dat.type = 'proportion'
    otu.da <- otu.da[rowSds(otu.da) != 0, ]
    zicoseq.diff <- GUniFrac::ZicoSeq(feature.dat = otu.da, 
                                      meta.dat = sample.da, 
                                      grp.name = group, 
                                      adj.name = adj.group,
                                      prev.filter = .1, 
                                      perm.no = 999,
                                      feature.dat.type = feature.dat.type
                    )
    zicoseq.diff$p.adj.df <- p.adjust(zicoseq.diff$p.raw, method)
    return(zicoseq.diff)
}

run_zicoseq <- function(zicoseq.diff, p.adj = TRUE){
    if (p.adj){
        p.adj <- 'p.adj.fdr'
    }else{
        p.adj <- 'p.raw'
    }
    names(zicoseq.diff[[p.adj]][zicoseq.diff[[p.adj]] < .05]) -> zicoseq.diff
    return(zicoseq.diff)
}

## ANCOMBC
obtain_ancombc_out <- function(mpse, method='bonferroni'){
    ps <- mpse %>% as.phyloseq(.abundance='total')
    sample.da <- mpse |> mp_extract_sample() |> tibble::column_to_rownames(var='Sample')
    if (ncol(sample.da) > 1){
        group <- colnames(sample.da)[1]
        adj.group <- colnames(sample.da)[2]
    }else{
        group <- colnames(sample.da)
        adj.group <- NULL
    }    
    ancombc.diff <- ANCOMBC::ancombc(ps,
                                 formula = paste0(c(group, adj.group), collapse="+"),
                                 p_adj_method= method,
                                 prv_cut = 0.1)
    return(list(xx=ancombc.diff$res, group = group))
}

run_ancombc <- function(res, p.adj = TRUE){
    #res <- obtain_ancombc_out(mpse)
    ind1 <- which(grepl(res$group, names(res$xx$diff_abn)))
    if (p.adj){
        p.adj <- 'q_val'
    }else{
        p.adj <- 'p_val'
    }
    ancombc.diff <- rownames(res$xx[[p.adj]][res$xx[[p.adj]][[ind1]] <= 0.05,,drop=FALSE]) 
    
    return(ancombc.diff)
}             

## metagenomeSeq
obtain_metaseq_out <- function(mpse, method='bonferroni'){
    otu.da <- mpse |> mp_extract_assays(.abundance = "Abundance")
    sample.da <- mpse |> mp_extract_sample() |> tibble::column_to_rownames(var='Sample')
    if (ncol(sample.da) > 1){
        group <- colnames(sample.da)[1]
        adj.group <- colnames(sample.da)[2]
    }else{
        group <- colnames(sample.da)
        adj.group <- NULL
    }    
    design <- model.matrix(as.formula(paste('~', group)), sample.da)
    data <- newMRexperiment(otu.da, phenoData = AnnotatedDataFrame(sample.da))
    data <- cumNorm(data, p = cumNormStatFast(data))
    res <- fitFeatureModel(data, design, coef = 2)
    pval <- res@pvalues
    res <- list(pval=pval, fdr=p.adjust(pval, method = method))
    return(res)
}

run_metaseq <- function(pval, p.adj = TRUE){
    #pval <- obtain_metaseq_out(mpse)
    if (p.adj){
        pval <- pval$fdr
    }else{
        pval <-pval$pval 
    }
    res <- names(pval[pval <= 0.05])
    res <- res[!is.na(res)]
    return(res)
}


## LEfSe
obtain_lefse_out <- function(mpse){
    otu.da <- mpse |> mp_extract_assays(.abundance = 'total', byRow = F)
    sample.da <- mpse |> mp_extract_sample() %>% column_to_rownames(var='Sample')
    if (ncol(sample.da) > 1){
        sample.da <- sample.da[,c(1, 2)]
        s <- 2
    }else{
        s <- -1
    }    
    lefseda <- merge(sample.da, otu.da, by=0)
    lefseda$Row.names <- NULL
    lefseda <- data.frame(t(lefseda), check.names=FALSE) %>% rownames_to_column(var="feature")
    tmpfile1 <- tempfile()
    tmpfile2 <- paste0(tmpfile1, ".format")
    outfile <- paste0(tmpfile1, ".out")
    write.table(lefseda, tmpfile1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    CMD1 <- paste("format_input.py", tmpfile1, tmpfile2, "-c 1 -s ", s," -o 1000000", sep=" ")
    CMD2 <- paste("run_lefse.py", tmpfile2, outfile, "--min_c 3 -f 1 -b 1 -l 2 -y 1", sep=" ")
    system(CMD1)
    system(CMD2)
    lefse.out <- read.table(outfile, sep="\t", header=F, row.names=1)
    lefse.out$V5 <- as.numeric(lefse.out$V5)
    lefse.out <- lefse.out[!is.na(lefse.out$V5),]
    return(lefse.out)
}

run_lefse <- function(lefse.out, p.adj = TRUE){
    #lefse.out <- obtain_lefse_out(mpse)
    if (p.adj){
        res <- rownames(lefse.out[lefse.out$V5 * nrow(mpse) <= 0.05,,drop=FALSE])
    }else{
        res <- rownames(lefse.out[lefse.out$V5 <= 0.05,,drop=FALSE])
    }
    res <- gsub("^f_", "", res)
    return(res)
}

cal_sparsity <- function(mpse){
    xx <- mpse %>% mp_extract_assays(.abundance = 'Abundance')
    xx <- sum(xx==0)/prod(dim(xx))
    return(xx)
}

for (file in list.files(path='./dataset/', pattern='*_mpse.rds')){
    prefix <- gsub("_mpse.rds", "", file)
    print (prefix)
    mpse <- readRDS(paste0("./dataset/",file))
    mpse %<>% mp_filter_taxa(.abundance=Abundance, min.abun=1, min.prop=.1, include.lowest=T)
    mpse %<>% mp_decostand(.abundance=Abundance, method='total')
    prefix <- gsub("_mpse.rds", "", file)
    print ("Running LEfSe (default, No p.adjust) ...............................")
    lefse <- obtain_lefse_out(mpse)
    lefse1 <- run_lefse(lefse, p.adj = FALSE)
    print ("Running LEfSe (p.adjust) ...............................")
    lefse2 <- run_lefse(lefse, p.adj = TRUE)
    print ("Running mp_diff_analysis (No p.adjust) ....................")
    mp.diff1 <- run_mpse_daa(mpse, p.adj = FALSE)    
    print ("Running mp_diff_analysis (default, p.adjust) ....................")
    mp.diff2 <- run_mpse_daa(mpse, p.adj = TRUE)
    print ("Running LinDA (No p.adjust) ..............................")
    linda.out <- obtain_linda_out(mpse)
    linda.diff1 <- run_linda(linda.out, p.adj = FALSE)
    print ("Running LinDA (default, p.adjust) ..............................")
    linda.diff2 <- run_linda(linda.out, p.adj = TRUE)    
    print ('Running ANCOMBC (No p.adjust) .............................')
    ancombc.out <- obtain_ancombc_out(mpse)
    ancombc.diff1 <- run_ancombc(ancombc.out, p.adj = FALSE)    
    print ('Running ANCOMBC (default, p.adjust) .............................')
    ancombc.diff2 <- run_ancombc(ancombc.out, p.adj = TRUE)
    print ('Running metagenomSeq (default, No p.adjust) ........................')
    metaseq.out <- obtain_metaseq_out(mpse)
    metaseq.diff1 <- run_metaseq(metaseq.out, p.adj = FALSE)
    print ('Running metagenomSeq (p.adjust) ........................')
    metaseq.diff2 <- run_metaseq(metaseq.out, p.adj = TRUE)
    print ('Running ZicoSeq (No p.adjust) .............................')
    zicoseq.out <- obtain_zicoseq_out(mpse)
    zicoseq.diff1 <- run_zicoseq(zicoseq.out, p.adj = FALSE)
    print ('Running ZicoSeq (default, p.adjust) .............................')
    zicoseq.diff2 <- run_zicoseq(zicoseq.out, p.adj = TRUE)
    res <- list(LEfSe = lefse1,
                LEfSe.adj = lefse2, 
                mp_diff_analysis = mp.diff1,
                mp_diff_analysis.adj = mp.diff2,
                LinDA = linda.diff1, 
                LinDA.adj = linda.diff2,
                ANCOMBC = ancombc.diff1,
                ANCOMBC.adj = ancombc.diff2,
                metagenomSeq = metaseq.diff1,
                metagenomeSeq.adj = metaseq.diff2,
                ZicoSeq = zicoseq.diff1,
                ZicoSeq.adj = zicoseq.diff2
    ) 
    res <- list(Diff=res, sample=ncol(mpse), features=nrow(mpse), sparsity=cal_sparsity(mpse))
    saveRDS(res, paste0("./result/", prefix, "_common.diff.res.prop.adj.noadj.rds"))
}
