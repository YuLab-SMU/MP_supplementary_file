library(ggplot2)
plotpoint <- function(data,
                      x,
                      y,
                      group,
                      type,pointsize=1,
		      pointstroke=0.2,
                      ...){
    p <- ggplot(data=data, 
                aes_string(x=x,y=y, color=group, shape=type))+
	 geom_point(size=pointsize, stroke=pointstroke, alpha=0.7,...) +
         theme_bw()+
         #scale_y_continuous(#expand=c(0,0),
         #                   limits=c(0,110))+
         scale_y_continuous(breaks=seq(0, 100, by=10),
                            limits=c(0,110)) +
	 theme(panel.grid=element_blank(),
               legend.box.spacing=unit(0.02,"cm"))
    p
}


plotline <- function(data, 
                     x,
                     y,
                     group,
                     type,...){

    p <- ggplot(data=data,
                aes_string(x=x,
                           y=y,
                           color=group,
                           linetype=type))+
         #geom_path()+
         #geom_point(size=0.5,...) +
         geom_line()+
         #geom_smooth(method="lm") +
         theme_bw() +
         scale_y_continuous(#expand=c(0,0.1),
                            limits=c(0,105))+
         theme(panel.grid=element_blank(),
               legend.box.spacing=unit(0.02,"cm"))
    p

}
