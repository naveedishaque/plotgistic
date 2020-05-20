#

# Author: Naveed Ishaque
# Date: 19th Feb 2019

# USAGE: R -f plot_gistic.R --args [output from gistic_2_dataframe_for_plotting.pl] [BED file for custom genes, 4 column]

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (output from gistic_2_dataframe_for_plotting.pl)", call.=FALSE)
}
if (length(args)>0) {
  # default output file
  file = args[1]
}
custom_bed=""
if (length(args)>1) {
  # default output file
  custom_bed = args[2]
}

file
custom_bed

gistic_results <- read.table(file, header=T, sep="")

amp<- gistic_results[gistic_results[,1]=="Amp",]
del<- gistic_results[gistic_results[,1]=="Del",]
chr<- gistic_results[gistic_results[,1]=="Chr",]
gene_amp<- gistic_results[gistic_results[,1]=="Gene_amp",]
gene_del<- gistic_results[gistic_results[,1]=="Gene_del",]

custom_genes<-NULL

if (nchar(custom_bed)>1){
  custom_genes<-read.table(custom_bed, sep="", header=F)
  for(i in 1:nrow(custom_genes)){
    chr_details = chr[chr[,2]==custom_genes[i,1],]
    custom_genes[i,2] = custom_genes[i,2] + chr_details$Position_offset - chr_details$Position
  }
}
custom_genes

pdf(paste(file,".pdf", sep=""), width=15, height=7.5)

par(mar=c(1.5,5,1,1))
par(fig=c(0,1,0.66,1))
plot(amp$Position_offset, amp$X.log10.q.value., col=rgb(191/255,95/255,0), lwd=2, ylab="", type="l", xaxt="null",xaxs="i")
mtext(side=2, line=3, "Amplification")
mtext(side=2,expression(paste("-log" ["10"], italic(" q-value"))),line=2)
abline(v=c(0,chr$Position_offset))   
abline(h=2, col="darkgrey", lty=2)
abline(h=1, col="darkgrey", lty=2)

abline(v=gene_amp$Position_offset, lty=2, col=rgb(191/255,95/255,0))
text(gene_amp$Position_offset-10000000,max(amp$X.log10.q.value.),labels=gene_amp$average_amplitude, adj=1, srt=90, cex=0.6)

if (nchar(custom_bed)>1){
  abline(v=custom_genes[,2], , lty=2, col="grey")
  text(custom_genes[,2]-10000000,max(amp$X.log10.q.value.),labels=custom_genes[,4], adj=1, srt=90, cex=0.6)
}

par(mar=c(1,5,1,1))   
par(fig=c(0,1,0.33,0.66), new=TRUE)
plot(amp$Position_offset, amp$frequency,ylim=c(-1,1),col=NA,xaxt="null",ylab="Frequency of CNAs", xlab="Chromosome",xaxs="i", yaxs="i")
mtext(side=2,c("Deletions                     \n\n","                         Amplifications\n\n"),col=c(rgb(0,102/255,204/255),rgb(191/255,95/255,0)))
polygon(c(0,amp$Position_offset,max(amp$Position_offset)), c(0,amp$frequency,0),col=rgb(191/255,95/255,0),border="NA")
polygon(c(0,del$Position_offset,max(del$Position_offset)), -c(0,del$frequency,0),col=rgb(0,102/255,204/255),border="NA") 
#polygon(c(0,amp$Position_offset,max(amp$Position_offset)), c(0,amp$average_amplitude,0),border="darkred",col="NA") 
#polygon(c(0,del$Position_offset,max(del$Position_offset)), -c(0,del$average_amplitude,0),border="darkblue",col="NA")
abline(v=c(0,chr$Position_offset))
abline(h=0.25, col="grey", lty=2)
abline(h=0.5, col="grey", lty=2)
abline(h=0.75, col="grey", lty=2)
abline(h=1, col="grey", lty=2)
abline(h=-0.25, col="grey", lty=2)
abline(h=-0.5, col="grey", lty=2)
abline(h=-0.75, col="grey", lty=2)
abline(h=-1, col="grey", lty=2)
box()
axis(1,chr$Chromosome, at=(chr$X.log10.q.value.),tick=F,pos=-0.97,las=2)
axis(3,chr$Chromosome, at=(chr$X.log10.q.value.),tick=F,pos=0.97,las=2)

par(mar=c(1,5,1.5,1))
par(fig=c(0,1,0,0.33), new=TRUE)
plot(del$Position_offset, -del$X.log10.q.value., col=rgb(0,102/255,204/255), lwd=2, ylab="", type="l", xaxt="null",xaxs="i")
mtext(side=2, line=3, "Deletion")
mtext(side=2,expression(paste("-log" ["10"], italic(" q-value"))),line=2)
abline(v=c(0,chr$Position_offset))   
abline(h=-2, col="darkgrey", lty=2)
abline(h=-1, col="darkgrey", lty=2)

abline(v=gene_del$Position_offset, lty=2, col=rgb(0,102/255,204/255))
text(gene_del$Position_offset-10000000,-max(del$X.log10.q.value.),labels=gene_del$average_amplitude, adj=0, srt=90, cex=0.6)

if (nchar(custom_bed)>1){
  abline(v=custom_genes[,2], , lty=2, col="grey")
  text(custom_genes[,2]-10000000,-max(del$X.log10.q.value.),labels=custom_genes[,4], adj=0, srt=90, cex=0.6)
}

dev.off()
