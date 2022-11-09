library(ComplexUpset)
library(magrittr)
library(ggplot2)
#set_size = function(w, h, factor=1.5) {
#    s = 1 * factor
#    options(
#        repr.plot.width=w * s,
#        repr.plot.height=h * s,
#        repr.plot.res=100 / factor,
#        jupyter.plot_mimetypes='image/png',
#        jupyter.plot_scale=1
#    )
#}
#
#for (file in list.files(path="./result/",  pattern="*diff.res.prop.adj.noadj.rds")){
#    prefix <- gsub('_common.diff.res.prop.adj.noadj.rds', '', file)
#    print (prefix)
#    res <- readRDS(paste0("./result/",file))
#    res <- res$Diff 
#    upset.data <- data.frame(Features=Reduce(union, res))
#    upset.data <- upset.data[!is.na(upset.data$Features),,drop=FALSE]
#    print (dim(upset.data))
#    lapply(res, function(i) ifelse(upset.data$Features %in% i, 1, 0)) %>% 
#        do.call('cbind', .) -> yy
#    print (dim(yy))
#    res <- data.frame(cbind(upset.data, yy))
#    #p <- ggplotify::as.ggplot(upset(res, nsets = 20, order.by='freq', point.size = 1.2, line.size = 0.4, text.scale=.8))
#    p <- upset(res, intersect=colnames(res)[-1], set_sizes=upset_set_size(position='right'), 
#               width_ratio=0.15, 
#               base_annotations=list('Intersection size'=intersection_size(
#                   text_colors = c(bar_number_threshold=2, on_background = "black", on_bar = "white"),
#                   text = list(angle=90, vjust=.5, hjust=0.5, position=position_nudge(y=5), size=2),
#                   width=4/5)
#               ), 
#               matrix=intersection_matrix(geom=geom_point(size=1), segment = geom_segment(size=.1)))
#    
#    ggplot2::ggsave(paste0('./result/', prefix, '_upset.prop.adj.noadj.svg'), p, device='svg', width=w, height=h)
#    #pdf(paste0('./result/', prefix, '_upset.svg'), width=6, height=6, onefile=F)
#    #upset(res, nsets = 6)
#    #dev.off()
#}


res <- readRDS('./result/ArcticFireSoils_common.diff.res.prop.adj.noadj.rds')
res <- res$Diff
upset.data <- data.frame(Features=Reduce(union, res))
lapply(res, function(i) ifelse(upset.data$Features %in% i, 1, 0)) %>% do.call('cbind', .) -> yy
res <- data.frame(cbind(upset.data, yy))
upset(res, intersect=colnames(res)[-1], set_sizes=upset_set_size(position='right'), width_ratio=0.1, base_annotations=list('Intersection size'=intersection_size(text_colors=c(bar_number_threshold=2, on_background = "black", on_bar = "white"), text = list(size=2, angle=90, vjust=.5, hjust=0.5, position=position_nudge(y=5)), width=3/5)), matrix=intersection_matrix(geom=geom_point(size=1), segment = geom_segment(size=.1)), name = "ArcticFireSoils") -> p1
p1$patches$plots[[1]] <- p1$patches$plots[[1]] + scale_y_continuous(expand=c(0,0,0,1)) + aplot:::theme_no_margin()
p1$patches$plots[[2]] <- p1$patches$plots[[2]] + aplot:::theme_no_margin()
svg('./result/ArcticFireSoils_upset.prop.adj.noadj.svg', width=13, height=5)
p1
dev.off()


res <- readRDS('./result/BISCUIT_common.diff.res.prop.adj.noadj.rds')
res <- res$Diff
upset.data <- data.frame(Features=Reduce(union, res))
lapply(res, function(i) ifelse(upset.data$Features %in% i, 1, 0)) %>% do.call('cbind', .) -> yy
res <- data.frame(cbind(upset.data, yy))
upset(res, intersect=colnames(res)[-1], set_sizes=upset_set_size(position='right'), width_ratio=0.1, base_annotations=list('Intersection size'=intersection_size(text_colors=c(bar_number_threshold=2, on_background = "black", on_bar = "white"), text = list(size=2, angle=90, vjust=.5, hjust=0.5, position=position_nudge(y=.21)), width=3/5)), matrix=intersection_matrix(geom=geom_point(size=1), segment = geom_segment(size=.1)), name = 'BISCUIT') -> p2
p2$patches$plots[[1]] <- p2$patches$plots[[1]] + scale_y_continuous(expand=c(0,0,0,1)) + aplot:::theme_no_margin()
p2$patches$plots[[2]] <- p2$patches$plots[[2]] + aplot:::theme_no_margin()
svg('./result/BISCUIT_upset.prop.adj.noadj.svg', width=7, height=4.5)
p2
dev.off()


res <- readRDS('./result/Blueberry_common.diff.res.prop.adj.noadj.rds')
res <- res$Diff
upset.data <- data.frame(Features=Reduce(union, res))
lapply(res, function(i) ifelse(upset.data$Features %in% i, 1, 0)) %>% do.call('cbind', .) -> yy
res <- data.frame(cbind(upset.data, yy))
upset(res, intersect=colnames(res)[-1], set_sizes=upset_set_size(position='right'), width_ratio=0.15, base_annotations=list('Intersection size'=intersection_size(text_colors=c(bar_number_threshold=2, on_background = "black", on_bar = "white"), text = list(size=2, angle=90, vjust=.5, hjust=0.5, position=position_nudge(y=3)), width=3/5)), matrix=intersection_matrix(geom=geom_point(size=1), segment = geom_segment(size=.1)), name = 'Blueberry') -> p3
p3$patches$plots[[1]] <- p3$patches$plots[[1]] + scale_y_continuous(expand=c(0,0,0,1)) + aplot:::theme_no_margin()
p3$patches$plots[[2]] <- p3$patches$plots[[2]] + aplot:::theme_no_margin()
svg('./result/Blueberry_upset.prop.adj.noadj.svg', width=8, height=4.5)
p3
dev.off()

res <- readRDS('./result/crc_baxter_common.diff.res.prop.adj.noadj.rds')
res <- res$Diff
upset.data <- data.frame(Features=Reduce(union, res))
lapply(res, function(i) ifelse(upset.data$Features %in% i, 1, 0)) %>% do.call('cbind', .) -> yy
res <- data.frame(cbind(upset.data, yy))
upset(res, intersect=colnames(res)[-1], set_sizes=upset_set_size(position='right'), width_ratio=0.15, base_annotations=list('Intersection size'=intersection_size(text_colors=c(bar_number_threshold=2, on_background = "black", on_bar = "white"), text = list(size=2, angle=90, vjust=.5, hjust=0.5, position=position_nudge(y=.5)), width=3/5)), matrix=intersection_matrix(geom=geom_point(size=1), segment = geom_segment(size=.1)), name = 'crc_baxter') -> p4
p4$patches$plots[[2]] <- p4$patches$plots[[2]] + aplot:::theme_no_margin()
p4$patches$plots[[1]] <- p4$patches$plots[[1]] + scale_y_continuous(expand=c(0,0,0,1)) + aplot:::theme_no_margin()
svg('./result/crc_baxter_upset.prop.adj.noadj.svg', width=8, height=4.5)
p4
dev.off()

res <- readRDS('./result/Crohn_common.diff.res.prop.adj.noadj.rds')
res <- res$Diff
upset.data <- data.frame(Features=Reduce(union, res))
lapply(res, function(i) ifelse(upset.data$Features %in% i, 1, 0)) %>% do.call('cbind', .) -> yy
res <- data.frame(cbind(upset.data, yy))
upset(res, intersect=colnames(res)[-1], set_sizes=upset_set_size(position='right'), width_ratio=0.15, base_annotations=list('Intersection size'=intersection_size(text = list(size=2.3), width=4/5)), matrix=intersection_matrix(geom=geom_point(size=1), segment = geom_segment(size=.1)), name = 'Crohn') -> p5
p5$patches$plots[[2]] <- p5$patches$plots[[2]] + aplot:::theme_no_margin()
p5$patches$plots[[1]] <- p5$patches$plots[[1]] + scale_y_continuous(expand=c(0,0,0,1)) + aplot:::theme_no_margin()
svg('./result/Crohn_upset.prop.adj.noadj.svg', width=8, height=5)
p5
dev.off()

res <- readRDS('./result/Diarrheal_common.diff.res.prop.adj.noadj.rds')
res <- res$Diff
upset.data <- data.frame(Features=Reduce(union, res))
lapply(res, function(i) ifelse(upset.data$Features %in% i, 1, 0)) %>% do.call('cbind', .) -> yy
res <- data.frame(cbind(upset.data, yy))
upset(res, intersect=colnames(res)[-1], set_sizes=upset_set_size(position='right'), width_ratio=0.15, base_annotations=list('Intersection size'=intersection_size(text_colors=c(bar_number_threshold=2, on_background = "black", on_bar = "white"), text = list(size=2, angle=90, vjust=.5, hjust=0.5, position=position_nudge(y=.5)), width=3/5)), matrix=intersection_matrix(geom=geom_point(size=1), segment = geom_segment(size=.1)), name = 'Diarrheal') -> p6
p6$patches$plots[[2]] <- p6$patches$plots[[2]] + aplot:::theme_no_margin()
p6$patches$plots[[1]] <- p6$patches$plots[[1]] + scale_y_continuous(expand=c(0,0,0,1)) + aplot:::theme_no_margin()
svg('./result/Diarrheal_upset.prop.adj.noadj.svg', width=8, height=4.5)
p6
dev.off()


res <- readRDS('./result/HLT_common.diff.res.prop.adj.noadj.rds')
res <- res$Diff
upset.data <- data.frame(Features=Reduce(union, res))
lapply(res, function(i) ifelse(upset.data$Features %in% i, 1, 0)) %>% do.call('cbind', .) -> yy
res <- data.frame(cbind(upset.data, yy))
upset(res, intersect=colnames(res)[-1], set_sizes=upset_set_size(position='right'), width_ratio=0.15, base_annotations=list('Intersection size'=intersection_size(text_colors=c(bar_number_threshold=2, on_background = "black", on_bar = "white"), text = list(size=2, angle=90, vjust=.5, hjust=0.5, position=position_nudge(y=.5)), width=3/5)), matrix=intersection_matrix(geom=geom_point(size=1), segment = geom_segment(size=.1)), name = 'HLT') -> p7
p7$patches$plots[[1]] <- p7$patches$plots[[1]] + scale_y_continuous(expand=c(0,0,0,1)) + aplot:::theme_no_margin()
p7$patches$plots[[2]] <- p7$patches$plots[[2]] + aplot:::theme_no_margin()
svg('./result/HLT_upset.prop.adj.noadj.svg', width=8, height=4.5)
p7
dev.off()


res <- readRDS('./result/ob_turnbaugh_common.diff.res.prop.adj.noadj.rds')
res <- res$Diff
upset.data <- data.frame(Features=Reduce(union, res))
lapply(res, function(i) ifelse(upset.data$Features %in% i, 1, 0)) %>% do.call('cbind', .) -> yy
res <- data.frame(cbind(upset.data, yy))
upset(res, intersect=colnames(res)[-1], set_sizes=upset_set_size(position='right'), width_ratio=0.15, base_annotations=list('Intersection size'=intersection_size(text_colors=c(bar_number_threshold=2, on_background = "black", on_bar = "white"), text = list(size=2, angle=90, vjust=.5, hjust=0.5, position=position_nudge(y=.5)), width=3/5)), matrix=intersection_matrix(geom=geom_point(size=1), segment = geom_segment(size=.1)), name = 'ob_turnbaugh') -> p8
p8$patches$plots[[2]] <- p8$patches$plots[[2]] + aplot:::theme_no_margin()
p8$patches$plots[[1]] <- p8$patches$plots[[1]] + scale_y_continuous(expand=c(0,0,0,1)) + aplot:::theme_no_margin()
svg('./result/ob_turnbaugh_upset.prop.adj.noadj.svg', width=8, height=4.5)
p8
dev.off()


res <- readRDS('./result/seston_plastic_mccormick_common.diff.res.prop.adj.noadj.rds')
res <- res$Diff
upset.data <- data.frame(Features=Reduce(union, res))
lapply(res, function(i) ifelse(upset.data$Features %in% i, 1, 0)) %>% do.call('cbind', .) -> yy
res <- data.frame(cbind(upset.data, yy))
upset(res, intersect=colnames(res)[-1], set_sizes=upset_set_size(position='right'), width_ratio=0.15, base_annotations=list('Intersection size'=intersection_size(text_colors=c(bar_number_threshold=2, on_background = "black", on_bar = "white"), text = list(size=2, angle=90, vjust=.5, hjust=0.5, position=position_nudge(y=.5)), width=3/5)), matrix=intersection_matrix(geom=geom_point(size=1), segment = geom_segment(size=.1)), name = 'seston_plastic_mccormick') -> p9
p9$patches$plots[[1]] <- p9$patches$plots[[1]] + scale_y_continuous(expand=c(0,0,0,1)) + aplot:::theme_no_margin()
p9$patches$plots[[2]] <- p9$patches$plots[[2]] + aplot:::theme_no_margin()
svg('./result/seston_plastic_mccormick_upset.prop.adj.noadj.svg', width=10, height=5)
p9
dev.off()


res <- readRDS('./result/Smoke_common.diff.res.prop.adj.noadj.rds')
res <- res$Diff
upset.data <- data.frame(Features=Reduce(union, res))
lapply(res, function(i) ifelse(upset.data$Features %in% i, 1, 0)) %>% do.call('cbind', .) -> yy
res <- data.frame(cbind(upset.data, yy))
upset(res, intersect=colnames(res)[-1], set_sizes=upset_set_size(position='right'), width_ratio=0.15, base_annotations=list('Intersection size'=intersection_size(text_colors=c(bar_number_threshold=2, on_background = "black", on_bar = "white"), text = list(size=2, angle=90, vjust=.5, hjust=0.5, position=position_nudge(y=.5)), width=3/5)), matrix=intersection_matrix(geom=geom_point(size=1), segment = geom_segment(size=.1)), name = 'Smoke') -> p10
p10$patches$plots[[2]] <- p10$patches$plots[[2]] + aplot:::theme_no_margin()
p10$patches$plots[[1]] <- p10$patches$plots[[1]] + scale_y_continuous(expand=c(0,0,0,1)) + aplot:::theme_no_margin()
svg('./result/Smoke_upset.prop.adj.noadj.svg', width=8, height=4.5)
p10
dev.off()

design <- "
  11
  23
"

f <- aplot::plot_list(p1, p3, p9, design = design, widths = c(1, 1.6), heights=c(1, .7))

svg('./result/soil_root_env_upset.svg', width=16, height=10)
f
dev.off()


design2 <- '
  12
  34
  56
'

f2 <- aplot::plot_list(p2, p4, p5, p6, p7, p8, p10, ncol=2)

svg('./result/human_gut.svg', width = 14, height = 18)
f2
dev.off()
