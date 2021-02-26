#!/usr/bin/env Rscript
library("devtools")
library (ggplot2)
require(gtools)
library(plyr)
library(devtools)
#install_github('kkdey/Logolas')
#library("Logolas")
library(RColorBrewer)
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)


setwd("/media/nikiwind/DATADRIVE2/ZFN_analysis/R_analysis/")

#qdata <- fread("logs.txt")
qdata <- fread("logs2.txt")
setnames(qdata, c("read", "tool","score"))
weird <- table(qdata$read)
bad <- names(weird[weird < 17])
qdata <- qdata[! read %in% bad]
qdata$position <- seq(17)
qdata$set <- c("odd","even","odd","even","odd","even","odd","even","odd","even","odd","even","odd","even","odd","even","odd")

#pdata <- fread("p.txt", header=F)
pdata <- fread("p2.txt", header=F)
setnames(pdata, c("read", "1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17"))
pdata <-melt(pdata, id.vars = "read", variable.name = "position",value.name = "ZF")

data <- merge(pdata, qdata)
data <- na.omit(data)
datastats <- aggregate(score ~ position + tool, data=data, summary)

#ggplot(data = data, aes(x=position, y=score)) + geom_boxplot() + facet_wrap(tool ~ ZF, scales="free")
#ggplot(data = data, aes(x=position, y=score)) + geom_boxplot(aes(fill=ZF)) + facet_wrap( ~ tool, scales="free")
#ggplot(data = data, aes(x=position, y=ZF)) + geom_bar(aes(fill=ZF), stat="identity") 
#ggplot(data = data) +  geom_bar(aes(x=position, y=ZF,fill=ZF), stat="identity") + facet_wrap( ~ set, scales="free")

#cutoff per quartile
# data[which(tool=="mmap" & score=="0")]$ZF <- NA
# for (p in 1:17) {
# data[which(tool=="gmap" & position==p & score<=datastats$score[p,5])]$ZF <- NA
# }
# for (p in 1:17) {
#   data[which(tool=="mmap" & position==p & score>=datastats$score[p+17,2])]$ZF <- NA
# }


#cutoff per median
data[which(tool=="mmap" & score=="0")]$ZF <- NA
for (p in 1:17) {
  data[which(tool=="gmap" & position==p & score<=datastats$score[p,3])]$ZF <- NA
}
for (p in 1:17) {
  data[which(tool=="mmap" & position==p & score>=datastats$score[p+17,3])]$ZF <- NA
}

dataclean <- as.data.frame.matrix((table(data[,2:3])))
datacleanper <- sweep(dataclean, 1, rowSums(dataclean), FUN="/")
datacleanper <- datacleanper*100
datacleanper$position <- seq(17)
datacleanper <- melt(datacleanper, id.vars="position")
setnames(datacleanper, c("position", "module","percent"))

#ordervector <- tibble(module = c("o1f","o1r","e2f","e2r","o3f","o3r","e4f","e4r","o5f","o5r","e6f","e6r","o7f","o7r","e8f","e8r","o9f","o9r","e10f","e10r","o11f","o11r","e12f","e12r","o13f","o13r","e14f","e14r","o15f","o15r","e16f","e16r","o17fAD","o17rAD"))
ordervector <- tibble(module = c("o1f","o1r","e2f","e2r","o3f","o3r","e4f","e4r","o5f","o5r","e6f","e6r","o7f","o7r","e8f","e8r","o9f","o9r","e10f","e10r","o11f","o11r","e12f","e12r","o13f","o13r","e14f","e14r","o15f","o15r","e16f","e16r","o17fAD","o17rAD"))
#ordervector <- tibble(module = c("o1f","o1r","o3f","o3r","o5f","o5r","o7f","o7r","o9f","o9r","o11f","o11r","o13f","o13r","o15f","o15r","o17fAD","o17rAD","e2f","e2r","e4f","e4r","e6f","e6r","e8f","e8r","e10f","e10r","e12f","e12r","e14f","e14r","e16f","e16r"))

datacleanper <- left_join(data.frame(name=ordervector),datacleanper,by="module")
datacleanper$module.ordered <- factor(datacleanper$module, levels=ordervector$module)

coul = brewer.pal(11, "Spectral")
coul = rev(colorRampPalette(coul)(34))
ggplot(data = datacleanper, aes(x=position, y=percent)) + geom_bar(aes(fill=module.ordered), colour="black",stat="identity") +  scale_x_continuous(breaks = seq(1, 17)) +  scale_fill_manual(values = coul)





coul1 = brewer.pal(9, "Greens")
coul1 = colorRampPalette(coul1)(17)
coul2 = brewer.pal(9, "Oranges")
coul2 = colorRampPalette(coul2)(17)
#coul3 <- c(rbind(coul1,rev(coul2))) 
coul3 <- c(rbind(coul1,coul2)) 
datacleanper <- ddply(datacleanper, .(position), transform, pos = cumsum(percent) - (0.5 * percent))
ggplot(data = datacleanper, aes(x=position, y=percent)) + geom_bar(aes(fill=module.ordered), colour="black",stat="identity") +  scale_x_continuous(breaks = seq(1, 17)) +  scale_fill_manual(values = coul3)






ggplot(data = datacleanper, aes(x=position, y=percent)) + geom_bar(aes(fill=module.ordered), colour="black",stat="identity") +  scale_x_continuous(breaks = seq(1, 17)) 

activ <- datacleanper[(datacleanper$module=="o17fAD" | datacleanper$module=="o17rAD") ,]
ggplot(data = activ, aes(x=position, y=percent)) + geom_bar(aes(fill=module.ordered), colour="black",stat="identity") +  scale_x_continuous(breaks = seq(1, 17)) + scale_fill_hue(l=60, c=50)


first <- datacleanper[(datacleanper$module=="o1f" | datacleanper$module=="o1r") ,]
ggplot(data = first, aes(x=position, y=percent)) + geom_bar(aes(fill=module.ordered), colour="black",stat="identity") +  scale_x_continuous(breaks = seq(1, 17)) + scale_fill_hue(l=50, c=45)

middle <- datacleanper[(datacleanper$module=="e8f" | datacleanper$module=="e8r" | datacleanper$module=="o9f" | datacleanper$module=="o9r") ,]
ggplot(data = middle, aes(x=position, y=percent)) + geom_bar(aes(fill=module.ordered), colour="black",stat="identity") +  scale_x_continuous(breaks = seq(1, 17)) + scale_fill_hue(l=50, c=45)


# lastpos <- lastpos[tool=="gmap"]
# p <- ggplot(data = data, aes(x=tool, y=score)) + geom_boxplot(aes(fill=ZF)) + facet_wrap( ~ tool, scales="free")

# data1 <- data[ZF=="o1f" | ZF=="o1r"]
# ggplot(data = data1, aes(x=position, y=ZF)) + geom_bar(aes(fill=ZF), stat="identity") 
# data2 <- data[ZF=="e2f" | ZF=="e2r"]
# ggplot(data = data2, aes(x=position, y=ZF)) + geom_bar(aes(fill=ZF), stat="identity") 
# data3 <- data[ZF=="o3f" | ZF=="o3r"]
# ggplot(data = data3, aes(x=position, y=ZF)) + geom_bar(aes(fill=ZF), stat="identity") 
# data3 <- data[ZF=="e4f" | ZF=="e4r"]
# ggplot(data = data3, aes(x=position, y=ZF)) + geom_bar(aes(fill=ZF), stat="identity") 
# data17 <- data[ZF=="o17fAD" | ZF=="o17rAD"]
# ggplot(data = data17, aes(x=position, y=ZF)) + geom_bar(aes(fill=ZF), stat="identity") 
# 
# 
# datapos <- data[ZF=="o1f" | ZF=="o1r"]$position
# 
# data17 <- data[ZF=="o17fAD" | ZF=="o17rAD"]$position

#list.files(pattern=".txt$") # use the pattern argument to define a common pattern  for import files with regex. Here: .csv
#list.filenames<-list.files(pattern=".txt$")
#list.filenames <- mixedsort(list.filenames)

# 
# for (i in 1:length(list.filenames))
# {
#   alldata <- do.call(cbind,lapply(list.filenames,read.table, header = FALSE,sep="\t"))
# }
# 
# alldata[(seq(4, 35, by = 3))] <- list(NULL)
# colnames(alldata) <- c("readname",seq(1, 12,by=0.5))
# colnames(alldata)[c(3,5,7,9,11,13,15,17,19,21,23,25)] <- "MAPQ"
# 
# alldata$string <- gsub(" ","", paste(substr(alldata[,2], 0, 1),
#                                      substr(alldata[,4], 0, 1),
#                                      substr(alldata[,6], 0, 1),
#                                      substr(alldata[,8], 0, 1),
#                                      substr(alldata[,10], 0, 1),
#                                      substr(alldata[,12], 0, 1),
#                                      substr(alldata[,14], 0, 1),
#                                      substr(alldata[,16], 0, 1),
#                                      substr(alldata[,18], 0, 1),
#                                      substr(alldata[,20], 0, 1),
#                                      substr(alldata[,22], 0, 1),
#                                      substr(alldata[,24], 0, 1),
#                                      collapse = NULL))
# 
# alldata[alldata == '*'] <- NA
# alldata[alldata == 255] <- 0
# alldata$sum <- rowSums(alldata[,c(3,5,7,9,11,13,15,17,19,21,23,25)])
# 
# alldata$missing <- rowSums(is.na(alldata))
# #plot(hist(alldata$missing, breaks=16))
# 
# alpha <- LETTERS[seq( from = 1, to = 24 )]
# 
# chq1 <- substr(alldata[which(alldata[,3]=="40"),2],0,1)
#   chq1 <- append(chq1, alpha)
#   chq1RVD <- mapvalues(chq1, from=as.vector(RVD2$module), to=as.vector(RVD2$RVD))
# chq2 <- substr(alldata[which(alldata[,5]=="40"),4],0,1)
#   chq2 <- append(chq2, alpha)
#   chq2RVD <- mapvalues(chq2, from=as.vector(RVD2$module), to=as.vector(RVD2$RVD))
# chq3 <- substr(alldata[which(alldata[,7]=="40"),6],0,1)
#   chq3 <- append(chq3, alpha)
#   chq3RVD <- mapvalues(chq3, from=as.vector(RVD2$module), to=as.vector(RVD2$RVD))
# chq4 <- substr(alldata[which(alldata[,9]=="40"),8],0,1)
#   chq4 <- append(chq4, alpha)
#   chq4RVD <- mapvalues(chq4, from=as.vector(RVD2$module), to=as.vector(RVD2$RVD))
# chq5 <- substr(alldata[which(alldata[,11]=="40"),10],0,1)
#   chq5 <- append(chq5, alpha)
#   chq5RVD <- mapvalues(chq5, from=as.vector(RVD2$module), to=as.vector(RVD2$RVD))
# chq6 <- substr(alldata[which(alldata[,13]=="40"),12],0,1)
#   chq6 <- append(chq6, alpha)
#   chq6RVD <- mapvalues(chq6, from=as.vector(RVD2$module), to=as.vector(RVD2$RVD))
# chq7 <- substr(alldata[which(alldata[,15]=="40"),14],0,1)
#   chq7 <- append(chq7, alpha)
#   chq7RVD <- mapvalues(chq7, from=as.vector(RVD2$module), to=as.vector(RVD2$RVD))
# chq8 <- substr(alldata[which(alldata[,17]=="40"),16],0,1)
#   chq8 <- append(chq8, alpha)
#   chq8RVD <- mapvalues(chq8, from=as.vector(RVD2$module), to=as.vector(RVD2$RVD))
# chq9 <- substr(alldata[which(alldata[,19]=="40"),18],0,1)
#   chq9 <- append(chq9, alpha)
#   chq9RVD <- mapvalues(chq9, from=as.vector(RVD2$module), to=as.vector(RVD2$RVD))
# chq10 <- substr(alldata[which(alldata[,21]=="40"),20],0,1)
#   chq10 <- append(chq10, alpha)
#   chq10RVD <- mapvalues(chq10, from=as.vector(RVD2$module), to=as.vector(RVD2$RVD))
# chq11 <- substr(alldata[which(alldata[,23]=="40"),22],0,1)
#   chq11 <- append(chq11, alpha)
#   chq11RVD <- mapvalues(chq11, from=as.vector(RVD2$module), to=as.vector(RVD2$RVD))
# chq12 <- substr(alldata[which(alldata[,25]=="40"),24],0,1)
#   chq12 <- append(chq12, alpha)
#   chq12RVD <- mapvalues(chq12, from=as.vector(RVD2$module), to=as.vector(RVD2$RVD))
# 
# t1 <- table(chq1)-1
# t2 <- table(chq2)-1
# t3 <- table(chq3)-1
# t4 <- table(chq4)-1
# t5 <- table(chq5)-1
# t6 <- table(chq6)-1
# t7 <- table(chq7)-1
# t8 <- table(chq8)-1
# t9 <- table(chq9)-1
# t10 <- table(chq10)-1
# t11 <- table(chq11)-1
# t12 <- table(chq12)-1
# 
# rt1 <- table(chq1RVD)-1
#   frt1 <- rt1/sum(rt1)
# rt2 <- table(chq2RVD)-1
#   frt2 <- rt2/sum(rt2)
# rt3 <- table(chq3RVD)-1
#   frt3 <- rt3/sum(rt3)
# rt4 <- table(chq4RVD)-1
#   frt4 <- rt4/sum(rt4)
# rt5 <- table(chq5RVD)-1
#   frt5 <- rt5/sum(rt5)
# rt6 <- table(chq6RVD)-1
#   frt6 <- rt6/sum(rt6)
# rt7 <- table(chq7RVD)-1
#   frt7 <- rt7/sum(rt7)
# rt8 <- table(chq8RVD)-1
#   frt8 <- rt8/sum(rt8)
# rt9 <- table(chq9RVD)-1
#   frt9 <- rt9/sum(rt9)
# rt10 <- table(chq10RVD)-1
#   frt10 <- rt10/sum(rt10)
# rt11 <- table(chq11RVD)-1
#   frt11 <- rt11/sum(rt11)
# rt12 <- table(chq12RVD)-1
#   frt12 <- rt12/sum(rt12)
# 
# ncombined <- cbind(length(chq1),length(chq2),length(chq3),length(chq4),length(chq5),length(chq6),length(chq7),length(chq8),length(chq9),length(chq10),length(chq11),length(chq12))
# combined <- cbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12)
# combinedRVD <- cbind(rt1,rt2,rt3,rt4,rt5,rt6,rt7,rt8,rt9,rt10,rt11,rt12)
# colnames(combined) <- c(1:12)
# colnames(combinedRVD) <- c(1:12)
# 
# # Caluclate nucleotide frequencies
# 
# nucfreq <- matrix(0, ncol = 4, nrow = 12)
# colnames(nucfreq) <- c("A","T","G","C")
# 
# for (rvdtocheck in c("AA","AA","CI","CP","CT","DN","DT","EN","HA","HD","HG","HN","HS","HT","HT","KT","NN","NS","NT","RA","RD","RI","VA","YT"))
#   {
#   nucfreq[1,] <- nucfreq[1,] + as.numeric(RVD2[RVD2$RVD==rvdtocheck,4:7][1,] * frt1[frt1 = rvdtocheck])
#   nucfreq[2,] <- nucfreq[2,] + as.numeric(RVD2[RVD2$RVD==rvdtocheck,4:7][1,] * frt2[frt2 = rvdtocheck])
#   nucfreq[3,] <- nucfreq[3,] + as.numeric(RVD2[RVD2$RVD==rvdtocheck,4:7][1,] * frt3[frt3 = rvdtocheck])
#   nucfreq[4,] <- nucfreq[4,] + as.numeric(RVD2[RVD2$RVD==rvdtocheck,4:7][1,] * frt4[frt4 = rvdtocheck])
#   nucfreq[5,] <- nucfreq[5,] + as.numeric(RVD2[RVD2$RVD==rvdtocheck,4:7][1,] * frt5[frt5 = rvdtocheck])
#   nucfreq[6,] <- nucfreq[6,] + as.numeric(RVD2[RVD2$RVD==rvdtocheck,4:7][1,] * frt6[frt6 = rvdtocheck])
#   nucfreq[7,] <- nucfreq[7,] + as.numeric(RVD2[RVD2$RVD==rvdtocheck,4:7][1,] * frt7[frt7 = rvdtocheck])
#   nucfreq[8,] <- nucfreq[8,] + as.numeric(RVD2[RVD2$RVD==rvdtocheck,4:7][1,] * frt8[frt8 = rvdtocheck])
#   nucfreq[9,] <- nucfreq[9,] + as.numeric(RVD2[RVD2$RVD==rvdtocheck,4:7][1,] * frt9[frt9 = rvdtocheck])
#   nucfreq[10,] <- nucfreq[10,] + as.numeric(RVD2[RVD2$RVD==rvdtocheck,4:7][1,] * frt10[frt10 = rvdtocheck])
#   nucfreq[11,] <- nucfreq[11,] + as.numeric(RVD2[RVD2$RVD==rvdtocheck,4:7][1,] * frt11[frt11 = rvdtocheck])
#   nucfreq[12,] <- nucfreq[12,] + as.numeric(RVD2[RVD2$RVD==rvdtocheck,4:7][1,] * frt12[frt12 = rvdtocheck])
# }
# 
# nucfreq <- (nucfreq/rowSums(nucfreq))*100
# 
# # Plotting...
# 
# color_profile <- list("type" = "per_row","col" = RColorBrewer::brewer.pal(dim(combined)[1],name ="Spectral"))
# color_profile$col[12] <- "#000000"
# color_profile$col <- c(rbind(color_profile$col,color_profile$col))
# color_profile$col <- c("#9E0142","#9E0142","#3288BD","#3288BD","#9E0142","#9E0142","#3288BD","#3288BD","#9E0142","#9E0142","#3288BD","#3288BD","#9E0142","#9E0142","#3288BD","#3288BD","#9E0142","#9E0142","#3288BD","#3288BD","#9E0142","#9E0142","#3288BD","#3288BD")
# 
# color_profile2 <- list("type" = "per_row","col" = colorRampPalette(brewer.pal(8,"Dark2"))(22))
# 
# logomaker(combined,
#           frame_width = 1,
#           color_profile = color_profile,
#           ic.scale = FALSE,
#           yscale_change=TRUE,
#           xlab="position",
#           xlab_fontsize=16,
#           xaxis_fontsize=16,
#           pop_name = "module frequency",
#           col_line_split = "grey80")
# 
# combined0 <- combined 
# combined0[combined0 > 0] <- 0
# combined0 = combined0[c(seq(1,24,2)), ]
# diag(combined0) <- 100
# 
# color_profile3 <- color_profile
# color_profile3$col <- color_profile3$col[1:12]
# 
# logomaker(combined0,
#           frame_width = 1,
#           color_profile = color_profile3,
#           ic.scale = FALSE,
#           yscale_change=TRUE,
#           xlab="position",
#           xlab_fontsize=16,
#           xaxis_fontsize=16,
#           pop_name = "module frequency",
#           col_line_split = "grey80")
# 
# 
# logomaker(combinedRVD,
#           frame_width = 1,
#           color_profile = color_profile2,
#           ic.scale = FALSE,
#           yscale_change=TRUE,
#           xlab="position",
#           xlab_fontsize=16,
#           xaxis_fontsize=16,
#           pop_name = "RVD frequency",
#           col_line_split = "grey80")
# 
# 
# 
# 
# color_profile$col <- c("darkred","darkblue","darkgrey","darkgreen")
# 
# tnucfreq <- t(nucfreq)
# colnames(tnucfreq) <- c(1:12)
# logomaker(tnucfreq,hist = FALSE,
#           frame_width = 1,
#           color_profile = color_profile,
#           ic.scale = FALSE,
#           alpha = 2,
#           yscale_change=TRUE,
#           xlab_fontsize=16,
#           xaxis_fontsize=16,
#           xlab="position",
#           pop_name = "nucleotide specificity",
#           col_line_split = "grey80")
# 
# logomaker(tnucfreq,
#           frame_width = 1,
#           color_profile = color_profile,
#           ic.scale = FALSE,
#           yscale_change=TRUE,
#           xlab_fontsize=16,
#           xaxis_fontsize=16,
#           xlab="position",
#           pop_name = "nucleotide specificity",
#           col_line_split = "grey80")
# 
# 

