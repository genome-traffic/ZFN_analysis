#!/usr/bin/env Rscript
library("devtools")
library (ggplot2)
require(gtools)
library(plyr)
library (data.table)
#require(ggseqlogo)

setwd("/media/nikolai/DATADRIVE1/BioProjects/TAL_generator/halfsams/")

#odd<-fread("odd.filtered.sam", select = c(1,3), header=F, sep="\t", stringsAsFactors=FALSE)
#even<-fread("even.filtered.sam", select = c(1,3), header=F, sep="\t", stringsAsFactors=FALSE)

odd<-fread("odd.f1.sam", select = c(1,3), header=F, sep="\t", stringsAsFactors=FALSE)
even<-fread("even.f1.sam", select = c(1,3), header=F, sep="\t", stringsAsFactors=FALSE)


merged <- merge(odd,even,by="V1")
merged[] <- lapply(merged, gsub, pattern='-', replacement='')
merged$pattern <- paste(
substr(merged$V3.x,1,1),substr(merged$V3.y,1,1),
substr(merged$V3.x,2,2),substr(merged$V3.y,2,2),
substr(merged$V3.x,3,3),substr(merged$V3.y,3,3),
substr(merged$V3.x,4,4),substr(merged$V3.y,4,4),
substr(merged$V3.x,5,5),substr(merged$V3.y,5,5),
substr(merged$V3.x,6,6),substr(merged$V3.y,6,6),sep="")

patternsonly <- merged$pattern
write.table(merged, "result.txt", sep="\t")
write.table(patternsonly, "resultpattern.txt", sep="\t",col.names = F, row.names = F)

alpha <- LETTERS[seq( from = 1, to = 24 )]

t1 <- table(c(alpha,substr(merged$pattern, 1, 1)))-1
t2 <- table(c(alpha,substr(merged$pattern, 2, 2)))-1
t3 <- table(c(alpha,substr(merged$pattern, 3, 3)))-1
t4 <- table(c(alpha,substr(merged$pattern, 4, 4)))-1
t5 <- table(c(alpha,substr(merged$pattern, 5, 5)))-1
t6 <- table(c(alpha,substr(merged$pattern, 6, 6)))-1
t7 <- table(c(alpha,substr(merged$pattern, 7, 7)))-1
t8 <- table(c(alpha,substr(merged$pattern, 8, 8)))-1
t9 <- table(c(alpha,substr(merged$pattern, 9, 9)))-1
t10 <- table(c(alpha,substr(merged$pattern, 10, 10)))-1
t11 <- table(c(alpha,substr(merged$pattern, 11, 11)))-1
t12 <- table(c(alpha,substr(merged$pattern, 12, 12)))-1

combined <- cbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12)
colnames(combined) <- c(1:12)

library(devtools)
install_github('kkdey/Logolas')
library("Logolas")
library(RColorBrewer)

color_profile <- list("type" = "per_row","col" = RColorBrewer::brewer.pal(dim(combined)[1],name ="Spectral"))
color_profile$col[12] <- "#000000"
color_profile$col <- c(rbind(color_profile$col,color_profile$col))
color_profile$col <- c("#9E0142","#9E0142","#3288BD","#3288BD","#9E0142","#9E0142","#3288BD","#3288BD","#9E0142","#9E0142","#3288BD","#3288BD","#9E0142","#9E0142","#3288BD","#3288BD","#9E0142","#9E0142","#3288BD","#3288BD","#9E0142","#9E0142","#3288BD","#3288BD")

color_profile2 <- list("type" = "per_row","col" = colorRampPalette(brewer.pal(8,"Dark2"))(22))

logomaker(combined,hist = FALSE,
          frame_width = 1,
          color_profile = color_profile,
          ic.scale = FALSE,
          alpha = 2,
          yscale_change=TRUE,
          xlab="position",
          col_line_split = "grey80")