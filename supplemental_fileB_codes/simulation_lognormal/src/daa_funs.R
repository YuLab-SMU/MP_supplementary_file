library(MicrobiomeStat)
library(GUniFrac)
library(matrixStats)
obtain_linda_out <- function(mpse){
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
                        prev.filter=0.1),
                        error = function(e)NULL)
    if (is.null(linda.diff)){
        return(NULL)
    }
    ind1 <- which(grepl(paste0("^",group), names(linda.diff$output)))
    return(linda.diff$output[[ind1]])
}

run_linda <- function(linda.diff, p.adj = TRUE){
    if (p.adj){
        p.adj <- 'padj'
    }else{
        p.adj <- 'pvalue'
    }
    linda.diff[linda.diff[[p.adj]] <= 0.05,] |>
        rownames() -> linda.diff
    return(linda.diff)

}

obtain_zicoseq_out <- function(mpse){
    otu.da <- mpse |> mp_extract_assays(.abundance = "Abundance") %>% as.matrix()
    sample.da <- mpse |> mp_extract_sample() |> tibble::column_to_rownames(var='Sample')
    adj.group <- NULL
    group <- colnames(sample.da)[1]
    otu.da <- otu.da[rowSds(otu.da) != 0, ]
    otu.da <- otu.da[!apply(otu.da,1,function(x)sum(is.na(x))>1),]
    zicoseq.diff <- GUniFrac::ZicoSeq(feature.dat = otu.da,
                                      meta.dat = sample.da,
                                      grp.name = group,
                                      adj.name = adj.group,
                                      prev.filter = .1,
                                      perm.no = 999,
                                      excl.pct = 0.5,
				      feature.dat.type = "other"
                    )
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
