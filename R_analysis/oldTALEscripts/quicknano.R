#!/usr/bin/env Rscript
library (ggplot2)
library (data.table)

system.time(c1<-fread(paste("/home/nikolai/Desktop/s.txt", sep=""), header=T, sep="\t",stringsAsFactors=FALSE))

#plot(c1$sequence_length_template,c1$mean_qscore_template)

plot(c1$sequence_length_template, c1$sequence_length_complement, cex=0.5, pch=21, bg = rainbow(length(c1$sequence_length_2d))[rank(c1$sequence_length_2d)], log="xy")
abline(v=3600)
abline(v=13572)
abline(h=3600)
abline(h=13572)

#plot(c1$sequence_length_template,c1$mean_qscore_template, ylim=c(5, 15), xlim=c(1,20000))

plot(c1$sequence_length_2d,c1$mean_qscore_2d)

plot(c1$sequence_length_2d,c1$sequence_length_complement, cex=0.5, pch=21, bg = heat.colors(length(c1$mean_qscore_2d))[rank(c1$mean_qscore_2d)], log="xy")
plot(c1$sequence_length_2d,c1$sequence_length_template, cex=0.5, pch=21, bg = heat.colors(length(c1$mean_qscore_2d))[rank(c1$mean_qscore_2d)], log="xy")
plot(c1$sequence_length_template,c1$sequence_length_complement, cex=0.5, pch=21, bg = heat.colors(length(c1$mean_qscore_2d))[rank(c1$mean_qscore_2d)], log="xy")

# 
#plot(c1$sequence_length_template, c1$sequence_length_complement)
# 
plot(density(c1$sequence_length_template), xlim=c(0,20000))
abline(v=3600)
abline(v=13572)

#c1 <- c1[which( c1$sequence_length_template > 11000 & c1$sequence_length_template < 14000),]
#c1 <- c1[which( c1$sequence_length_2d > 0),]
#c1 <- c1[which( c1$sequence_length_2d > 12000),]
