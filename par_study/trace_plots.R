library(ggplot2)
library(reshape2)
library(dplyr)
library(ggh4x)
library(gridExtra)
library(RColorBrewer)
library(grid)

cols <- c('stat','env_rat1', 'm', 'T', 'K1', 'K2', 'Alpha1', 'Alpha2', 'M1', 'M2')
data <- read.csv('data_mT_col2_w07.csv', header=FALSE, col.names = cols)
data$rat <- data$M1/(data$M1 + data$M2)

myColors <- brewer.pal(8,"Spectral")
names(myColors) <- levels(data$m)
colScale <- scale_colour_manual(name = "m",values = myColors)


data_rat <- select(filter(data, stat=='mean'), c('T','m','rat'))
dat.mrat <- melt(data_rat, id.var=c('m','T'))
prat <- ggplot(dat.mrat) + geom_line(aes(x=log10(T), y=value, group=m, colour=as.factor(m)), lwd=1.25) + 
  geom_point(aes(x=log10(T), y=value, group=m, colour=as.factor(m)),lwd=2) +
  theme_bw() + colScale + ylim(0,1) + ggtitle(expression(frac(M[1],M[1]+M[2])))

data_I <- select(filter(data, stat=='mean'), c('T','m','K1','K2','Alpha1', 'Alpha2'))
dat.mI <- melt(data_I, id.var=c('m','T'))
pKI <- ggplot(dat.mI) + geom_line(aes(x=log10(T), y=value, group=m, colour=as.factor(m)), lwd=1.25) + 
  geom_point(aes(x=log10(T), y=value, group=m, colour=as.factor(m)),lwd=2) +
  facet_wrap(~variable, scales='free') +
  theme_bw() + colScale + facetted_pos_scales(
    y = rep(list(
      scale_y_continuous(limits = c(0, 10000)),
      scale_y_continuous(limits = c(0, 1))), each = 2))

grid.arrange(prat, pKI, nrow=1)
