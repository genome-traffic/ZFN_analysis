#!/usr/bin/env Rscript
library("devtools")
#install_github("omarwagih/ggseqlogo")
library (ggplot2)
#library (data.table)
require(ggseqlogo)

setwd("/media/nikolai/DATADRIVE1/BioProjects/TAL_generator/logo")
logodata <- read.table("input.logo")

cs1 = make_col_scheme(chars=c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X'),
cols=c('red','red4','blue','blue4','green','green4','yellow','yellow4','magenta','magenta4','brown','brown4','gray','gray4','pink','pink4','cyan','cyan4','gold','gold4','aquamarine','aquamarine4','tomato','tomato4'))


cs2 = make_col_scheme(chars=c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X'),  values=5:28)

val_col <- colorRampPalette(c('firebrick','cornflowerblue'))(24)
cs3 = make_col_scheme(chars=c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X'),
                      cols=(val_col))


ggplot() + geom_logo(as.vector(logodata$V1), col_scheme=cs3, method='p',seq_type='other', namespace=c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X')) + theme_logo()
