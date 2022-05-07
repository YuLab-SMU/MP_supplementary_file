# genelist and reference gene P (Positive) or N (Negative)
confusiontab <- function(lists, ref){
    totaln <- nrow(ref)/2
    if (is.null(lists)||length(lists)==0){
        pred <- c(rep("N", nrow(ref)))
        pred <- factor(pred, levels=c("P", "N"))
        ref <- ref[,2]
        res <- table(pred, ref)
        res <- data.frame(res)
    }else{
        others <- setdiff(as.vector(ref[,1]), lists)
        pred <- data.frame(f=c(lists, others),
                           pred=c(rep("P", length(lists)),
                           rep("N", length(others))))
        pred$pred <- factor(pred$pred, levels=c("P","N"))
        ref <- ref[match(pred$f,as.vector(ref[,1])),2]
        pred <- pred$pred
        res <- table(pred, ref)
        res <- data.frame(res)
    }
    res <- data.frame(FPR=100*res[3,3]/totaln, FNR=100*res[2,3]/totaln, 
                      TRP=100*res[1,3]/totaln, TNR=100*res[4,3]/totaln)
    return(res)
}

