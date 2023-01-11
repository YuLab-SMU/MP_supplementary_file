library(dplyr)
library(ggtree)
library(ggplot2)
library(aplot)
library(shadowtext)
num.list <- list()
sample.list <- list()
n <- 0
for (file in list.files(path="./result/",  pattern="*diff.res.prop.adj.noadj.rds")){
    n <- n + 1
    prefix <- gsub('_common.diff.res.prop.adj.noadj.rds', '', file)
    res <- readRDS(paste0("./result/",file))
    num.list[[n]] <- purrr::map(res$Diff, function(x)length(x[!is.na(x)])) %>%
           tibble::as_tibble() %>% dplyr::mutate(study=prefix)
    sample.list[[n]] <- data.frame(study=prefix, sample=res$sample, features=res$features, sparsity=res$sparsity) 
}


num.list <- num.list %>% dplyr::bind_rows() 

dat <- num.list %>% tibble::column_to_rownames(var='study') #%>% scale() %>% data.frame()
dat <- apply(dat, 1, function(x){x - mean(x) / sd(x)}) %>% t() %>% data.frame()

dat %>% dist() %>% hclust(method='ave') %>% ape::as.phylo() -> tr1
dat %>% t() %>% dist() %>% hclust(method='ave') %>% ape::as.phylo() -> tr2

df <- data.frame(method=tr2$tip.label, adj.method=ifelse(grepl(".adj$", tr2$tip.label), 'bonferroni', "raw pvalue"))
#df[df$method %in% c('LEfSe', 'metagenomSeq'), 'default'] <- 'Yes'
#df[df$method %in% c('LEfSe.adj', 'metagenomeSeq.adj'), 'default'] <- 'No'

p1 <- ggtree(tr1)

p2 <- ggtree(tr2, layout='den')


num.dat <- num.list %>% tidyr::pivot_longer(
               cols=!'study', names_to = 'method', values_to = 'num'
            ) %>% 
           left_join(
               dat %>% 
               tibble::rownames_to_column(var = 'study') %>%
               tidyr::pivot_longer(cols=!'study', names_to = 'method', values_to = 'scale.num'))


sample.df <- sample.list %>% dplyr::bind_rows()
sample.df$type <- ifelse(sample.df$study %in% c('ArcticFireSoils', 'Blueberry', 'seston_plastic_mccormick'), 'Soil', 'Human Stool')
sample.df[sample.df$study == 'seston_plastic_mccormick', 'type'] <- 'Marine plastic' 
p3 <- ggplot(data=num.dat, aes(x=method, y=study)) +
      geom_tile(aes(fill=scale.num), color='grey90') +
      geom_shadowtext(aes(label = num), 
                      size = 2.2, 
                      color='black', 
                      bg.colour='white'
      ) +
      labs(x = NULL, y = NULL) +
      scale_fill_viridis_c(option = 'magma', direction = -1, guide = 'none') + 
      theme_bw() +
      coord_cartesian(expand = F) +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = -45, hjust = 0),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.key.width = unit(.3, 'cm'),
        legend.key.height = unit(.3, 'cm'),
        legend.box.margin=ggplot2::margin()       
      )

p4 <- ggplot(df, aes(y='adj.method', x=method, fill=adj.method)) +
      geom_tile() +
      theme_bw() +
      scale_fill_manual(values=c("#CC79A7", "#0072B2")) +
      xlab(NULL) +
      ylab(NULL) +
      scale_y_discrete(position = 'right') +
      coord_cartesian(expand = F) +
      #coord_flip(expand = F) + 
      theme_bw() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.key.width = unit(.3, 'cm'),
        legend.key.height = unit(.3, 'cm'),
        legend.box.margin=ggplot2::margin()
      ) 

p5 <- ggplot(sample.df, aes(x='samples', y=study, fill = sample)) + 
      geom_tile() +
      xlab(NULL) +
      ylab(NULL) +
      theme_bw() +
      coord_cartesian(expand = F) +
      #scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
      scale_fill_gradient2() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = .5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.key.width = unit(.3, 'cm'),
        legend.key.height = unit(.3, 'cm'),
        legend.box.margin=ggplot2::margin()
      )

p6 <- ggplot(sample.df, aes(x="features", y=study, fill = features)) +
      geom_tile() +
      xlab(NULL) +
      ylab(NULL) +
      coord_cartesian(expand = F) +
      theme_bw() +
      scale_fill_viridis_c() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = .5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.key.width = unit(.3, 'cm'),
        legend.key.height = unit(.3, 'cm'),
        legend.box.margin=ggplot2::margin()
      )

p7 <- ggplot(sample.df, aes(x="type", y=study, fill = type)) +
      geom_tile() +
      xlab(NULL) +
      ylab(NULL) +
      scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
      coord_cartesian(expand = F) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = .5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.key.width = unit(.3, 'cm'),
        legend.key.height = unit(.3, 'cm'),
        legend.box.margin=ggplot2::margin()
      )

p8 <- ggplot(sample.df, aes(x="sparsity", y=study, fill = sparsity)) +
      geom_tile() +
      xlab(NULL) +
      ylab(NULL) +
      scale_fill_viridis_c(option='E') +
      coord_cartesian(expand = F) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = .5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.key.width = unit(.3, 'cm'),
        legend.key.height = unit(.3, 'cm'),
        legend.box.margin = ggplot2::margin()
      )

pp <- p3 %>% 
      insert_left(p7, width = .025) %>%
      insert_left(p8, width = .0251) %>%
      insert_left(p5, width = .0252) %>%
      insert_left(p6, width = .0253) %>%
      insert_left(p1, width = .1) %>% 
      insert_top(p4, height = .022) %>% 
      insert_top(p2, height = .08)

svg('./result/num_plot.svg', width=6, height = 5.5)
pp
dev.off()
