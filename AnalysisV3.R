# Robertson, Alvarez, et al. 2020
# Code for working with Reynoutria epiGBS data
# Feb 1, 2020

# Set working directory here
setwd("/Volumes/ANALYSIS3/2018Knotweed/")
library(tidyverse)
options(contrasts = c("contr.treatment","contr.poly"))  # The default


# This function formats the SNP table from the VCF file to create allele frequencies. It also implements a depth cutoff of 10x
PrepGenetic<-function(Depth=10){
  #bcftools query -f '%CHROM %POS  %REF  %ALT [ %SAMPLE=%RO] [ %SAMPLE=%DP]\n' snp.vcf.gz > AD.table
  
  ## Preliminary steps to get data organized and ready to convert
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  
  allele.depth2 <-fread("AD.table", header=FALSE,data.table = FALSE)
  
  temp<-data.frame(t(allele.depth2[1,-c(1:4)]))
  temp<-data.frame(str_split_fixed(temp$X1,"=",2))
  colnames(allele.depth2)<-c("Chr","Pos","Ref","Alt",as.character(temp$X1))
  sampleNames<-unique(as.character(temp$X1))
  
  sampleSize<-length(sampleNames)
  
  allele.count <- cbind(allele.depth2[1:4],allele.depth2[,c(5:(sampleSize+4))])
  allele.depth <- cbind(allele.depth2[1:4],allele.depth2[,-c(1:(sampleSize+4))])
  fwrite(allele.depth,"depthMatrix.txt",sep="\t",quote=FALSE,row.names = FALSE)
  
  for(i in 5:ncol(allele.count)){
    print(i)
    allele.count[,i]<-gsub(paste(colnames(allele.count)[i],"=",sep=""),"",allele.count[,i])
    allele.depth[,i]<-gsub(paste(colnames(allele.depth)[i],"=",sep=""),"",allele.depth[,i])
  }
  
  #Depth=10
  
  tempmat<-data.matrix(allele.depth[,-c(1:4)])
  tempmat[tempmat < Depth] <- NA
  
  AlleleFreq<-data.frame(allele.count[,c(1:4)],(data.matrix(allele.count[,-c(1:4)])/tempmat))
  
  fwrite(AlleleFreq,paste("Allele.frequency.",Depth,".txt",sep=""),quote=FALSE,row.names = FALSE,sep="\t")
  
}
# This function filters loci and samples with too much missing data, then imputes the rest.
FilterGenetic<-function(){
  library(data.table)
  library(tidyr)
  ######
  
  AlleleFreq<-fread("Allele.frequency.10.txt",data.table = FALSE)
  
  AlleleFreq$LocusMissing<-rowSums(is.na(AlleleFreq[,-c(1:4)]))/ncol(AlleleFreq[,-c(1:4)])
  
  AlleleFreq<-dplyr::filter(AlleleFreq,LocusMissing<1)
  AlleleFreq$LocusMissing<-NULL
  
  InvariantFilter<-function(){
    # Filter out fully methylated invariant sites
    print("Filtering 100% methylated....")
    AlleleFreq2<-AlleleFreq[,-c(1:4)]
    AlleleFreq2[is.na(AlleleFreq2)] <- 1
    AlleleFreq$Invariant<-rowSums(AlleleFreq2==1)/ncol(AlleleFreq2)
    rm(AlleleFreq2)
    AlleleFreq<-dplyr::filter(AlleleFreq,!Invariant==1)
    AlleleFreq$Invariant<-NULL
    
    # Filter out not methylated invariant sites
    print("Filtering 0% methylated...")
    AlleleFreq2<-AlleleFreq[,-c(1:4)]
    AlleleFreq2[is.na(AlleleFreq2)] <- 0
    AlleleFreq$Invariant<-rowSums(AlleleFreq2==0)/ncol(AlleleFreq2)
    rm(AlleleFreq2)
    AlleleFreq<-dplyr::filter(AlleleFreq,!Invariant==1)
    AlleleFreq$Invariant<-NULL
    return(AlleleFreq)
  }
  AlleleFreq<-InvariantFilter()
  
  SampleFilter<-function(MissingMax=0.50){
    SampMissing<-c()
    for(i in 5:ncol(AlleleFreq)){
      SampMissing<-c(SampMissing,sum(is.na(AlleleFreq[,i]))/nrow(is.na(AlleleFreq)))
    }
    SampMissing<-c(0,0,0,0,SampMissing)
    SampMissing
    sum(SampMissing>.80)
    AlleleFreq<-AlleleFreq[,!(SampMissing>MissingMax)]
    return(AlleleFreq)
  }
  #AlleleFreq<-SampleFilter(0.8)
  
  LocusFilter<-function(MissingMax=0.5){
    AlleleFreq$LocusMissing<-rowSums(is.na(AlleleFreq[,-c(1:4)]))/ncol(AlleleFreq[,-c(1:4)])
    summary(AlleleFreq$LocusMissing)
    AlleleFreq<-dplyr::filter(AlleleFreq,LocusMissing<=MissingMax)
    AlleleFreq$LocusMissing<-NULL
    return(AlleleFreq)
  }
  #AlleleFreq<-LocusFilter(0.5)
  
  SampleFilter<-function(MissingMax=0.80){
    SampMissing<-c()
    for(i in 5:ncol(AlleleFreq)){
      SampMissing<-c(SampMissing,sum(is.na(AlleleFreq[,i]))/nrow(is.na(AlleleFreq)))
    }
    SampMissing<-c(0,0,0,0,SampMissing)
    SampMissing
    sum(SampMissing>.80)
    AlleleFreq<-AlleleFreq[,!(SampMissing>MissingMax)]
    return(AlleleFreq)
  }
  AlleleFreq<-SampleFilter(0.8)
  
  LocusFilter<-function(MissingMax=0.5){
    AlleleFreq$LocusMissing<-rowSums(is.na(AlleleFreq[,-c(1:4)]))/ncol(AlleleFreq[,-c(1:4)])
    summary(AlleleFreq$LocusMissing)
    AlleleFreq<-dplyr::filter(AlleleFreq,LocusMissing<=MissingMax)
    AlleleFreq$LocusMissing<-NULL
    return(AlleleFreq)
  }
  AlleleFreq<-LocusFilter(0.4)
  
  AlleleFreq<-InvariantFilter()
  
  PlotMissing<-function(){
    PlotDat<-AlleleFreq
    PlotDat<-data.frame(ChrPos=paste(PlotDat$Chr,PlotDat$Pos,sep=":"),PlotDat[,-c(1:4)])
    library(reshape2)
    PlotDat2<-reshape2::melt(PlotDat,id.vars=colnames(PlotDat)[1])
    PlotDat2$variable<-gsub("SPALT_","",PlotDat2$variable)
    PlotDat2$value<-is.na(PlotDat2$value)
    PlotDat2$value<-gsub(FALSE,"Present",PlotDat2$value)
    PlotDat2$value<-gsub(TRUE,"Imputed",PlotDat2$value)
    PlotDat3<-data.frame(table(PlotDat2$variable,PlotDat2$value))
    PlotDat4 <- PlotDat3  %>%
      group_by(Var1, Var2) %>%
      summarise(n = sum(Freq)) %>%
      mutate(percentage = n / sum(n))
    
    order<-unique(PlotDat4$Var1)
    
    library(ggplot2)
    library(ggthemes)
    library(ggsci)
    library(tidyverse)
    
    #pdf("CrazyPlot.pdf",height=10,width=8)
    #ggplot(PlotDat2,aes(y=as.numeric(as.factor(ChrPos)),x=variable,color=value)) + 
    #  geom_point(shape=15,alpha=0.5) + #coord_polar() +
    #  theme_bw() + scale_color_npg(name="Value handling") +
    #  theme(axis.title.y=element_blank(),
    #        axis.text.y=element_blank(),
    #        axis.ticks.y=element_blank())
    #dev.off()
    
    # Samples
    #pdf("SampleCoverage.pdf",height=8,width=8)
    GenSamp<-ggplot(PlotDat4,aes(x=Var1,y=percentage,fill=Var2)) + 
      geom_bar(stat="identity") + labs(x="Sample",y="Percentage of loci") +
      geom_hline(yintercept = 0.2,linetype = "dashed",color="white")+
      theme_minimal() + scale_fill_hc(name="Data handling") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "bottom") +
      theme(axis.title.x=element_text(vjust=-2)) +
      theme(axis.title.y=element_text(angle=90, vjust=-0.5)) +
      theme(plot.title=element_text(size=15, vjust=3)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"))
    #dev.off()
    
    PlotDat3<-data.frame(table(PlotDat2$ChrPos,PlotDat2$value))
    PlotDat4 <- PlotDat3  %>%
      group_by(Var1, Var2) %>%
      summarise(n = sum(Freq)) %>%
      mutate(percentage = n / sum(n))
    
    #pdf("LocusCoverage.pdf",height=8,width=8)
    GenLoc<-ggplot(PlotDat4[1:10000,],aes(x=as.numeric(as.factor(Var1)),y=percentage,fill=Var2)) + 
      geom_area() + labs(x="Locus",y="Percentage of data at a locus") +
      geom_hline(yintercept = 0.5,linetype = "dashed",color="white") +
      theme_minimal() + scale_fill_hc(guide=FALSE) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(axis.title.x=element_text(vjust=-2)) +
      theme(axis.title.y=element_text(angle=90, vjust=-0.5)) +
      theme(plot.title=element_text(size=15, vjust=3)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"))
    #dev.off()
    
    library(cowplot)
    #pdf("FigS2.pdf",height=8,width=8)
    #print(plot_grid(Samp,Loc,rel_heights = c(1.2, 1),
    #                #rel_widths = c(1, 1.2),
    #                nrow = 2,labels="AUTO"))
    #dev.off()
    save(GenLoc,GenSamp,order,file="GenPlots.RData")
    
  }
  PlotMissing()
  load("GenPlots.RData")
  print(GenLoc)
  print(GenSamp)
  
  
  
  library(impute)
  tempmat<-data.matrix(AlleleFreq[,-c(1:4)])
  imputed<-impute.knn(tempmat,colmax=0.9)
  imputed<-data.frame(AlleleFreq[,c(1:4)],imputed$data)
  
  fwrite(imputed,"AlleleFreq.txt",quote=FALSE,sep="\t",row.names = FALSE)
  
}
# This function uses STAMPP to get the genetic relatedness matrix (based on Nei's d).
PopGenPrep_stampp<-function(){
  ###
  
  library(StAMPP)
  library(data.table)
  freq<-fread("AlleleFreq.txt")
  Pheno<-fread("pheno.csv",data.table = FALSE)
  colnames(Pheno)[1]<-"Sample"
  
  stamped<-data.frame(Sample=colnames(freq)[-c(1:4)],Pop=NA,Ploidy=rep(6,ncol(freq)-4),Format="freq",t(freq[,-c(1:4)]))
  tmp<-stamped[,c(1,3:4)]
  tmp<-merge(tmp,Pheno,all.x = TRUE,by="Sample")
  stamped<-merge(tmp[,c(1,4,2,3)],stamped[,-c(2:4)],by="Sample")
  #stamped$Pop<-as.numeric(as.factor(stamped$Pop))
  
  allele.freq <- stamppConvert(stamped, type='r')
  allele.GL <- stampp2genlight(allele.freq, TRUE)
  
  ClonalID<-function(){
    library(poppr)
    library(philentropy)
    distplot<-as.matrix(distance(t(freq[,-c(1:4)]),method = "euclidean",use.row.names = TRUE))
    rownames(distplot) <- gsub("FALLOPIA_","",rownames(distplot))
    
    tmp<-stamped[,-c(1:4)]
    rownames(tmp)<-gsub("FALLOPIA_","",stamped$Sample)
    hc <- hclust(dist(tmp), "ave")
    
    library(rafalib)
    myplclust(hc, labels=gsub("FALLOPIA_","",rownames(tmp)),
              lab.col=as.fumeric(stamped$pop), cex=1.5,
              hang=0.1)
    
    
    pdf("GenotypeDistances.pdf",height=8,width=10)
    #print(plot(hclust(as.dist(distplot))))
    print(myplclust(hc, labels=gsub("FALLOPIA_","",rownames(tmp)),
                    lab.col=as.fumeric(stamped$pop), cex=1.5,
                    hang=0.1))
    dev.off()
    af.dist<-data.matrix(distplot)
    #af.gc<-recode_polyploids(allele.GL,newploidy = TRUE)
    
    e<-t(stamped[,-c(1:4)])
    
    library(genefilter)
    rv <- rowVars(e)
    idx <- order(-rv)[1:500]
    
    
    library(RColorBrewer) 
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(nrow(e))
    
    library(gplots) ##Available from CRAN
    
    colourCount = length(unique(stamped$pop))
    getPalette = colorRampPalette(brewer.pal(20, "Dark2"))
    cols <- getPalette(colourCount)[as.fumeric(stamped$pop)]
    head(cbind(colnames(e),cols))
    
    pdf("Heatmap.truncated.pdf",height=10,width=10)
    print(
    heatmap.2(data.matrix(e)[idx,], labCol=stamped$pop,
              trace="none", 
              ColSideColors=cols, 
              col=hmcol)
    )
    dev.off()
    
    pdf("Heatmap.jpeg",height=1000,width=1000)
    print(
      heatmap.2(data.matrix(e), labCol=stamped$pop,
                trace="none", 
                ColSideColors=cols, 
                col=hmcol)
    )
    dev.off()
    
    
    
    
    a<-mlg.filter(allele.GL,distance=af.dist,
                  stats="ALL",
                  threshold = 1e+06 + .Machine$double.eps^0.5,
                  algorithm = "average_neighbor")
    f<-mlg.filter(allele.GL,distance=af.dist,
                  stats="ALL",
                  threshold = 1e+06 + .Machine$double.eps^0.5,
                  algorithm = "f")
    n<-mlg.filter(allele.GL,distance=af.dist,
                  stats="ALL",
                  threshold = 1e+06 + .Machine$double.eps^0.5,
                  algorithm = "n")
    
    fanlist <- list(farthest = f, average = a, nearest = n)
    
    plot_filter_stats(allele.GL, fanlist, af.dist, cols=NULL, NULL, 
                      "Scott")
    #print(farthest_thresh <- cutoff_predictor(fanlist$farthest$THRESHOLDS))
    print(farthest_thresh <- cutoff_predictor(fanlist$average$THRESHOLDS))
    #print(farthest_thresh <- cutoff_predictor(fanlist$nearest$THRESHOLDS))
    
    al.filt<-mlg.filter(allele.GL,distance=af.dist,
                        stats="ALL",
                        threshold = farthest_thresh,
                        algorithm = "a")
    al.filt$MLGS
    print(paste("Found",length(unique(al.filt$MLGS)),"unique genotypes"))
    
    
    # Updated cutoff
    al.filt<-mlg.filter(allele.GL,distance=af.dist,
                        stats="ALL",
                        threshold = 13.5,
                        algorithm = "a")
    al.filt$MLGS
    print(paste("Found",length(unique(al.filt$MLGS)),"unique genotypes"))
    
    tmp<-melt(af.dist)
    tmp$Var2<-gsub("FALLOPIA_","",tmp$Var2)
    tmp<-dplyr::filter(tmp,!value==0)
    
    pdf("Histogram.pdf",height=7,width = 7)
    ggplot(tmp,aes(x=value)) + geom_histogram(bins = 50) + theme_bw() +
      geom_vline(xintercept = 8.225) +
      labs(x="Euclidean distance",y="Pairwise similarity value")
    dev.off()
    
    return(af.dist)
  }
  af.dist<-ClonalID()
  fwrite(af.dist,"Distances.euclidean.csv",quote=FALSE,row.names = TRUE,col.names = TRUE)
  
  allele.D.ind <- stamppNeisD(allele.freq, pop=FALSE)
  allele.D.pop <- stamppNeisD(allele.freq, pop=TRUE)
  allele.relationship<-stamppGmatrix(allele.freq)
  fwrite(allele.relationship,"Gmatrix.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names = FALSE)
  
  allele.D.ind<-data.frame(allele.D.ind)
  colnames(allele.D.ind) <- c(row.names(allele.D.ind))
  #names(allele.D.pop)[1:45] <- c(row.names(allele.D.pop))
  fwrite(allele.D.ind, "genetic.distance.ind.txt", quote=FALSE, sep="\t", row.names = TRUE, col.names = FALSE)
  
  allele.D.pop<-data.frame(allele.D.pop)
  colnames(allele.D.pop)<-c("GIN1","GIN2","GIO1","GIO2","MSN","MSO")
  fwrite(allele.D.pop, "genetic.distance.pop.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
  
  #Amova <- stamppAmova(allele.D.ind, allele.freq, 10000)
  
  allele.fst <- stamppFst(allele.freq, nboots=999, percent=95, nclusters=6)
  save(list = c("allele.D.ind","allele.D.pop","allele.fst"), file = "GeneticDistanceAF.RData",compress = TRUE)
  load("GeneticDistanceAF.RData")
}
# This function formats the methylation table from the BED file to create methylation frequencies. It also implements a depth cutoff of 10x
PrepMeth<-function(Depth=10){
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  allele.depth2<-fread("methylation.bed",data.table = FALSE)
  
  allele.count <- cbind(allele.depth2[1:4],dplyr::select(allele.depth2,contains("methylated")))
  allele.depth <- cbind(allele.depth2[1:4],dplyr::select(allele.depth2,contains("total")))
  
  
  
  tempmat<-data.matrix(allele.depth[,-c(1:4)])
  tempmat[tempmat < Depth] <- NA
  
  AlleleFreq<-data.frame(allele.count[,c(1:4)],(data.matrix(allele.count[,-c(1:4)])/tempmat))
  
  fwrite(AlleleFreq,paste("Meth.frequency.",Depth,".txt",sep=""),quote=FALSE,row.names = FALSE,sep="\t")
  
  
  
}
# This function filters loci and samples with too much missing data, then imputes the rest.
FilterMeth<-function(){
  ######
  
  library(data.table)
  library(dplyr)
  library(stringr)
  
  AlleleFreq<-fread("Meth.frequency.10.txt",data.table = FALSE)
  Genetic<-fread("AlleleFreq.txt",data.table = FALSE)
  
  
  AlleleFreq$LocusMissing<-rowSums(is.na(AlleleFreq[,-c(1:4)]))/ncol(AlleleFreq[,-c(1:4)])
  
  AlleleFreq<-dplyr::filter(AlleleFreq,LocusMissing<1)
  AlleleFreq$LocusMissing<-NULL
  
  InvariantFilter<-function(){
    # Filter out fully methylated invariant sites
    print("Filtering 100% methylated....")
    AlleleFreq2<-AlleleFreq[,-c(1:4)]
    AlleleFreq2[is.na(AlleleFreq2)] <- 1
    AlleleFreq$Invariant<-rowSums(AlleleFreq2==1)/ncol(AlleleFreq2)
    rm(AlleleFreq2)
    AlleleFreq<-dplyr::filter(AlleleFreq,!Invariant==1)
    AlleleFreq$Invariant<-NULL
    
    # Filter out not methylated invariant sites
    print("Filtering 0% methylated...")
    AlleleFreq2<-AlleleFreq[,-c(1:4)]
    AlleleFreq2[is.na(AlleleFreq2)] <- 0
    AlleleFreq$Invariant<-rowSums(AlleleFreq2==0)/ncol(AlleleFreq2)
    rm(AlleleFreq2)
    AlleleFreq<-dplyr::filter(AlleleFreq,!Invariant==1)
    AlleleFreq$Invariant<-NULL
    return(AlleleFreq)
  }
  AlleleFreq<-InvariantFilter()
  
  colnames(AlleleFreq)<-gsub("_methylated","",colnames(AlleleFreq))
  
  AlleleFreq<-data.frame(AlleleFreq[,c(1:4)],dplyr::select(AlleleFreq, one_of(colnames(Genetic)[-c(1:4)])))
  
  SampMissing<-c()
  for(i in 5:ncol(AlleleFreq)){
    SampMissing<-c(SampMissing,sum(is.na(AlleleFreq[,i]))/nrow(is.na(AlleleFreq)))
  }
  SampMissing<-c(0,0,0,0,SampMissing)
  SampMissing
  sum(SampMissing>.80)
  #AlleleFreq<-AlleleFreq[,!(SampMissing>.80)]
  
  AlleleFreq$LocusMissing<-rowSums(is.na(AlleleFreq[,-c(1:4)]))/ncol(AlleleFreq[,-c(1:4)])
  summary(AlleleFreq$LocusMissing)
  AlleleFreq<-dplyr::filter(AlleleFreq,LocusMissing<=.4) ## ********
  AlleleFreq$LocusMissing<-NULL
  
  AlleleFreq<-InvariantFilter()
  
  AlleleFreq2<-data.frame(ID=paste(AlleleFreq$chr,AlleleFreq$pos,sep=":"),AlleleFreq)
  
  PlotMissing<-function(){
    load("GenPlots.RData")
    PlotDat<-AlleleFreq2[,-c(2:5)]
    #PlotDat<-data.frame(ChrPos=paste(PlotDat$Chr,PlotDat$Pos,sep=":"),PlotDat[,-c(1:4)])
    library(reshape2)
    PlotDat2<-reshape2::melt(PlotDat,id.vars=colnames(PlotDat)[1])
    #PlotDat2$variable<-gsub("SPALT_","",PlotDat2$variable)
    PlotDat2$value<-is.na(PlotDat2$value)
    PlotDat2$value<-gsub(FALSE,"Present",PlotDat2$value)
    PlotDat2$value<-gsub(TRUE,"Imputed",PlotDat2$value)
    PlotDat3<-data.frame(table(PlotDat2$variable,PlotDat2$value))
    PlotDat4 <- PlotDat3  %>%
      group_by(Var1, Var2) %>%
      summarise(n = sum(Freq)) %>%
      mutate(percentage = n / sum(n))
    
    PlotDat4$Var1<-factor(PlotDat4$Var1,levels=order)
    
    #PlotDat4<-PlotDat4[order(colnames(Genetic)[-c(1:4)]),] 
    
    library(ggplot2)
    library(ggthemes)
    library(ggsci)
    
    #pdf("CrazyPlot.pdf",height=10,width=8)
    #ggplot(PlotDat2,aes(y=as.numeric(as.factor(ChrPos)),x=variable,color=value)) + 
    #  geom_point(shape=15,alpha=0.5) + #coord_polar() +
    #  theme_bw() + scale_color_npg(name="Value handling") +
    #  theme(axis.title.y=element_blank(),
    #        axis.text.y=element_blank(),
    #        axis.ticks.y=element_blank())
    #dev.off()
    
    # Samples
    #pdf("SampleCoverage.pdf",height=8,width=8)
    MethSamp<-ggplot(PlotDat4,aes(x=Var1,y=percentage,fill=Var2)) + 
      geom_bar(stat="identity") + labs(x="Sample",y="Percentage of loci") +
      geom_hline(yintercept = 0.2,linetype = "dashed",color="white")+
      theme_minimal() + scale_fill_hc(name="Data handling") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "bottom") +
      theme(axis.title.x=element_text(vjust=-2)) +
      theme(axis.title.y=element_text(angle=90, vjust=-0.5)) +
      theme(plot.title=element_text(size=15, vjust=3)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"))
    #dev.off()
    
    PlotDat3<-data.frame(table(PlotDat2$ID,PlotDat2$value))
    PlotDat4 <- PlotDat3  %>%
      group_by(Var1, Var2) %>%
      summarise(n = sum(Freq)) %>%
      mutate(percentage = n / sum(n))
    PlotDat4<-data.frame(PlotDat4)
    
    #pdf("LocusCoverage.pdf",height=8,width=8)
    MethLoc<-ggplot(PlotDat4,aes(x=as.numeric(as.factor(Var1)),y=percentage,fill=Var2)) + 
      geom_area() + labs(x="Locus",y="Percentage of data at a locus") +
      geom_hline(yintercept = 0.5,linetype = "dashed",color="white") +
      theme_minimal() + scale_fill_hc(guide=FALSE) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(axis.title.x=element_text(vjust=-2)) +
      theme(axis.title.y=element_text(angle=90, vjust=-0.5)) +
      theme(plot.title=element_text(size=15, vjust=3)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"))
    #dev.off()
    
    library(cowplot)
    #pdf("FigS2.pdf",height=8,width=8)
    #print(plot_grid(Samp,Loc,rel_heights = c(1.2, 1),
    #                #rel_widths = c(1, 1.2),
    #                nrow = 2,labels="AUTO"))
    #dev.off()
    
    save(MethLoc,MethSamp,file="MethPlots.RData")
    
  }
  PlotMissing()
  #load("MethPlots.RData")
  
  out<-data.frame()
  library(impute)
  
  tempmat<-data.matrix(AlleleFreq[,-c(1:4)])
  imputed<-impute.knn(tempmat,colmax=0.9)
  imputed<-data.frame(AlleleFreq[,c(1:4)],imputed$data)
  
  fwrite(imputed,"MethFreq.txt",quote=FALSE,sep="\t",row.names = FALSE)
  #####
  
  # macau
  
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  allele.depth2<-fread("methylation.bed",data.table = FALSE)
  
  #freq<-fread("MethFreq.txt",data.table = FALSE)[,-c(1:2)]
  freq<-fread("MethFreq.txt",data.table = FALSE)
  allele.depth <- cbind(allele.depth2[1:4],dplyr::select(allele.depth2,contains("total")))
  #rm(allele.depth2)
  colnames(allele.depth)<-gsub("_total","",colnames(allele.depth))
  colnames(allele.depth)<-gsub("-",".",colnames(allele.depth),fixed = TRUE)
  #colnames(allele.depth)
  allele.depth<-dplyr::select(allele.depth, one_of(colnames(freq)))
  allele.depth$chrpos<-paste(allele.depth$chr,allele.depth$pos,allele.depth$context,sep=":")
  #allele.depth<-data.frame(chrpos=paste(allele.depth$chr,allele.depth$pos,allele.depth$context,sep=":"),allele.depth[,-c(1:4)])
  freq$chrpos<-paste(freq$chr,freq$pos,freq$context,sep=":")
  #freq<-data.frame(chrpos=paste(freq$chr,freq$pos,freq$context,sep=":"),freq[,-c(1:4)])
  allele.depth<-dplyr::filter(allele.depth, chrpos %in% freq$chrpos)
  #allele.depth$chrpos<-NULL
  
  tmp<-data.matrix(allele.depth[,-c(1:4,ncol(allele.depth))])
  allele.depth[,-c(1:4,ncol(allele.depth))]<-tmp
  rm(tmp)
  #for(i in 2:ncol(allele.depth)){
  #  print(i)
  #  allele.depth[,i]<-as.numeric(as.character(allele.depth[,i]))
  #}
  
  #allele.depth<-allele.depth %>% mutate_all(~replace(., is.na(.), 10))
  allele.depth[is.na(allele.depth)]<-10
  allele.depth<-allele.depth[order(match(allele.depth$chrpos,freq$chrpos)),]
  
  allele.depth<-allele.depth[,-c(1:4)]
  allele.depth<-allele.depth[,c(ncol(allele.depth),1:(ncol(allele.depth)-1))]
  freq<-freq[,-c(1:5)]
  freq<-freq[,c(ncol(freq),1:(ncol(freq)-1))]
  
  fwrite(allele.depth,"counts.macau.txt",quote=FALSE,row.names = FALSE,sep = "\t")
  
  tempmat<-data.matrix(freq[,-1])
  tempdep<-data.matrix(allele.depth[,-1])
  tempmeth<-data.matrix(tempmat*tempdep)
  tempmeth<-round(tempmeth)
  #sum(tempmeth>tempdep)
  freq[,-1]<-tempmeth
  fwrite(freq,"mcounts.macau.txt",quote=FALSE,row.names = FALSE,sep = "\t")
  
}
# This function makes Figure S2
MakeFilterPlots<-function(){
  load("GenPlots.RData")
  load("MethPlots.RData")
  
  library(cowplot)
  pdf("FigS2.pdf",height=8,width=12)
  print(plot_grid(GenSamp,MethSamp,GenLoc,MethLoc,nrow = 2,ncol=2,labels="AUTO",rel_heights = c(1.2,1.2,1,1)))
  dev.off()
}


# Load and format analysis files - SNP matrix, methylation frequencies, phenotypes

library(data.table)
allele<-fread("AlleleFreq.txt",data.table = FALSE)
pheno<-fread("fallopia_sample_inf.txt",data.table=FALSE)
loc<-fread("locations.csv",data.table = FALSE)
pheno<-merge(pheno,loc,by="site",all.x=TRUE)
pheno<-pheno[,-c(7:9)]
pheno<-pheno[,c(2,1,3:ncol(pheno))]

colnames(allele)<-gsub(".","-",colnames(allele),fixed = TRUE)
pheno<-dplyr::filter(pheno,sampleID %in% colnames(allele))
pheno<-pheno[match(colnames(allele)[-c(1:4)], pheno$sampleID),]
meth<-fread("MethFreq.txt",data.table = FALSE)#[,-1]

pheno$sampleID<-gsub("FALLOPIA_","",pheno$sampleID)
colnames(meth)<-gsub("FALLOPIA_","",colnames(meth))
colnames(allele)<-gsub("FALLOPIA_","",colnames(allele))

# for coding region subsets - defunct
#frags<-fread("frags.filtered.csv")
#allele<-dplyr::filter(allele,Chr %in% frags$queryID)
#meth<-dplyr::filter(meth,chr %in% frags$queryID)


# C/T removal - defunct

CTremoval<-function(){
  allele$ChrPos<-paste(allele$Chr,allele$Pos,sep=":")
  allele.t<-dplyr::filter(allele, grepl("A",Ref))
  allele.t<-dplyr::filter(allele, grepl("G",Alt))
  ToRemove<-allele.t$ChrPos
  allele.t<-dplyr::filter(allele, grepl("T",Ref))
  allele.t<-dplyr::filter(allele, grepl("C",Alt))
  ToRemove<-c(ToRemove,allele.t$ChrPos)
  allele<-dplyr::filter(allele,!ChrPos %in% ToRemove)
  allele$ChrPos<-NULL
  #remove(allele.t)
  return(allele)
}
#allele<-CTremoval()


# Identify latent factors

GetLFs<-function(df=allele,K=3){
  library(lfmm)
  Y <- df[,-c(1:4)]
  pc <- prcomp(t(Y))
  plot((pc$sdev[1:20])^2, xlab = 'PC', ylab = "Variance explained")
  points(4,(pc$sdev[4])^2, type = "h", lwd = 3, col = "blue")
  # two latent factors
  
  modmat<-model.matrix(~treatment,data=pheno)
  
  Y<-data.frame(t(df[,-c(1:4)]))
  mod.lfmm <- lfmm_ridge(Y = Y, 
                         X = modmat, 
                         K = K)
  return(mod.lfmm$U)
}
LF.allele<-GetLFs(df=allele,K=3) 
LF.meth<-GetLFs(df=meth,K=4)


# MLL plots
MLL<-function(){
  freq<-allele
  
  ClonalID<-function(){
    library(poppr)
    library(philentropy)
    distplot<-as.matrix(distance(t(freq[,-c(1:4)]),method = "euclidean",use.row.names = TRUE))
    rownames(distplot) <- gsub("FALLOPIA_","",rownames(distplot))
    
    tmp<-t(freq[,-c(1:4)])
    tmp<-data.frame(sampleID=rownames(tmp),tmp)
    tmp<-merge(pheno,tmp,by="sampleID",all.y = TRUE)
   
    
    
    e<-t(tmp[,-c(1:6)])
    
    library(genefilter)
    rv <- rowVars(e)
    idx <- order(-rv)[1:1000]
    tmp.filt<-t(data.matrix(e)[idx,])
    tmp.filt<-data.frame(tmp[,1:6],tmp.filt)
    
    library(philentropy)
    m<-data.matrix(dist(tmp[,-c(1:6)],method = "binary"))
    xy <- t(combn(colnames(m), 2))
    toplot<-data.frame(xy, dist=m[xy])
    pdf("Histogram.pdf",height=6,width=6)
    print(ggplot(toplot,aes(x=1-dist)) + 
            geom_histogram(binwidth = 0.0001,aes(color=)) + 
            theme_bw() +
      labs(x="Jaccard similarity",y="Frequency"))
    dev.off()
    
    hc <- hclust(dist(tmp[,-c(1:6)]), "ave")
    #hc <- hclust(dist(tmp.filt[,-c(1:6)]), "ave")
    
    #library(pvclust)
    #phc<-pvclust(tmp[,-c(1:6)], method.hclust = "average",
    #        method.dist = "correlation", nboot = 999)
    
    library(rafalib)
    myplclust(hc, labels=gsub("FALLOPIA_","",tmp$sampleID),
              lab.col=as.fumeric(tmp$species), cex=1.5,
              hang=0.1)
    
    
    pdf("GenotypeDistances.pdf",height=8,width=10)
    #print(plot(hclust(as.dist(distplot))))
    print(myplclust(hc, labels=gsub("FALLOPIA_","",tmp$sampleID),
                    lab.col=as.fumeric(tmp$species), cex=1.5,
                    hang=0.1))
    dev.off()
    af.dist<-data.matrix(distplot)
    #af.gc<-recode_polyploids(allele.GL,newploidy = TRUE)
    
    e<-t(stamped[,-c(1:4)])
    
    library(genefilter)
    rv <- rowVars(e)
    idx <- order(-rv)[1:500]
    
    
    library(RColorBrewer) 
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(nrow(e))
    
    library(gplots) ##Available from CRAN
    
    colourCount = length(unique(stamped$pop))
    getPalette = colorRampPalette(brewer.pal(20, "Dark2"))
    cols <- getPalette(colourCount)[as.fumeric(stamped$pop)]
    head(cbind(colnames(e),cols))
    
    pdf("Heatmap.truncated.pdf",height=10,width=10)
    print(
      heatmap.2(data.matrix(e)[idx,], labCol=stamped$pop,
                trace="none", 
                ColSideColors=cols, 
                col=hmcol)
    )
    dev.off()
    
    pdf("Heatmap.jpeg",height=1000,width=1000)
    print(
      heatmap.2(data.matrix(e), labCol=stamped$pop,
                trace="none", 
                ColSideColors=cols, 
                col=hmcol)
    )
    dev.off()
    
    
    
    
    a<-mlg.filter(allele.GL,distance=af.dist,
                  stats="ALL",
                  threshold = 1e+06 + .Machine$double.eps^0.5,
                  algorithm = "average_neighbor")
    f<-mlg.filter(allele.GL,distance=af.dist,
                  stats="ALL",
                  threshold = 1e+06 + .Machine$double.eps^0.5,
                  algorithm = "f")
    n<-mlg.filter(allele.GL,distance=af.dist,
                  stats="ALL",
                  threshold = 1e+06 + .Machine$double.eps^0.5,
                  algorithm = "n")
    
    fanlist <- list(farthest = f, average = a, nearest = n)
    
    plot_filter_stats(allele.GL, fanlist, af.dist, cols=NULL, NULL, 
                      "Scott")
    #print(farthest_thresh <- cutoff_predictor(fanlist$farthest$THRESHOLDS))
    print(farthest_thresh <- cutoff_predictor(fanlist$average$THRESHOLDS))
    #print(farthest_thresh <- cutoff_predictor(fanlist$nearest$THRESHOLDS))
    
    al.filt<-mlg.filter(allele.GL,distance=af.dist,
                        stats="ALL",
                        threshold = farthest_thresh,
                        algorithm = "a")
    al.filt$MLGS
    print(paste("Found",length(unique(al.filt$MLGS)),"unique genotypes"))
    
    
    # Updated cutoff
    al.filt<-mlg.filter(allele.GL,distance=af.dist,
                        stats="ALL",
                        threshold = 13.5,
                        algorithm = "a")
    al.filt$MLGS
    print(paste("Found",length(unique(al.filt$MLGS)),"unique genotypes"))
    
    tmp<-melt(af.dist)
    tmp$Var2<-gsub("FALLOPIA_","",tmp$Var2)
    tmp<-dplyr::filter(tmp,!value==0)
    
    pdf("Histogram.pdf",height=7,width = 7)
    ggplot(tmp,aes(x=value)) + geom_histogram(bins = 50) + theme_bw() +
      geom_vline(xintercept = 8.225) +
      labs(x="Euclidean distance",y="Pairwise similarity value")
    dev.off()
    
    return(af.dist)
  }
  af.dist<-ClonalID()
  fwrite(af.dist,"Distances.euclidean.csv",quote=FALSE,row.names = TRUE,col.names = TRUE)
  
}

# Depth PCA
DepthPCA<-function(){
  depth<-fread("depthMatrix.txt",data.table = FALSE)
  depth$chrpos<-paste(depth$Chr,depth$Pos,sep=":")
  allele$chrpos<-paste(allele$Chr,allele$Pos,sep=":")
  depth<-dplyr::filter(depth,chrpos %in% allele$chrpos)
  depth<-dplyr::select(depth, one_of(colnames(allele)))
  depth$chrpos<-NULL
  depth[,-c(1:4)]<-data.matrix(depth[,-c(1:4)])
  depth[is.na(depth)]<-10
  
  topc<-data.matrix(t(depth[,-c(1:4)]))
  pcdata<-prcomp(topc)
  sum<-summary(pcdata)
  sum<-data.frame((sum$importance))
  sum<-data.frame(var=rownames(sum),sum)[2,]
  sum<-reshape2::melt(sum,id.vars="var")
  pcdata<-data.frame(pheno,pcdata$x)
  pcdata$site<-gsub("MSO9","MSO",pcdata$site)
  
  library(RColorBrewer)
  colourCount = length(unique(sum$variable))
  getPalette = colorRampPalette(brewer.pal(8, "Set1"))
  
  library(ggplot2)
  library(ggthemes)
  PC<-ggplot(pcdata,aes(x=PC1,y=PC2,color=site,shape=treatment)) + 
    geom_point() + theme_bw() + scale_color_few() + guides(color=guide_legend(title="Site")) +
    guides(shape=guide_legend(title="Exposure"))
  PVE<-ggplot(sum,aes(x=variable,y=value,fill=variable)) + 
    geom_bar(stat="identity",aes(factor(variable)),fill=getPalette(colourCount)) +
    theme_bw() + labs(x="Principal component",y="Proportion of variance explained")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  library(cowplot)
  
  pdf("DepthPCA.pdf",height=8,width=8)
  print(plot_grid(PC,PVE,nrow=2,labels="AUTO"))
  dev.off()
}
#DepthPCA()

ClusteringAndHeatMap<-function(){
  library(WGCNA)
  
  
  # Clustering
  a.d<-dist(t(allele[,-c(1:4)]),method = "euclidean")
  
  dat<-t(allele[,-c(1:4)])
  
  library(shipunov)
  hc<-Bclust(dat,mc.cores = 6,bootstrap = TRUE)
  plot(hc)
  #library(pvclust)
  #hc<-pvclust(dat,parallel = 8)
  
  pdf("GenotypeDistances.bootstraps.pdf",width=15,height=15)
  plot(hc)
  dev.off()
  
  
  hc<-hclust(a.d,"ave")
  plot(hc)
  clusts<-cutree(hc,3)
  clusts[which(clusts==3)]<-2
  
  dhc <- as.dendrogram(hc)
  ddata <- dendro_data(dhc, type = "rectangle")
  p2 <- ggplot(segment(ddata)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
  
  labs <- label(ddata)
  ph2<-pheno
  colnames(ph2)[1]<-"label"
  labs <- merge(labs,ph2,by="label",all=TRUE)
  
  p2 + geom_text(data=label(ddata),
                 aes(label=label, x=x, y=0, colour=labs$treatment))
  
  
  library(rafalib)
  pdf("GenotypeDistances.pdf",height=8,width=10)
  print(myplclust(hc, labels=colnames(allele)[-c(1:4)],
                  lab.col=as.fumeric(as.character(paste("C",clusts,sep=""))), cex=1.5,
                  hang=0.1))
  dev.off()
  
  clusts<-data.frame(sampleID=names(clusts),Cluster=clusts)
  pheno2<-merge(pheno,clusts,by="sampleID")
  #fwrite(pheno2,"updated_pheno2.csv",quote=FALSE,row.names = FALSE)
  
  library(genefilter)
  e<-(allele[,-c(1:4)])
  rv <- rowVars(e)
  idx <- order(-rv)[1:500]
  tmp.filt<-(data.matrix(e)[idx,])
  #tmp.filt<-data.frame(tmp[,1:6],tmp.filt)
  
  library(gplots)
  hr <- hclust(dist(tmp.filt,method = "binary"), method="ave")
  #hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete") 
  ## Tree cutting
  mycl <- cutree(hr, h=0.99)
  mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
  ## Plot heatmap 
  mycol <- redgreen(75)#colorpanel(40, "darkblue", "yellow", "white") # or try redgreen(75)
  heatmap.2(data.matrix(allele[,-c(1:4)]), Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), 
            col=mycol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc) 
  
  
  a.d<-dist(t(allele[,-c(1:4)]),method = "binary")
  
  m<-data.matrix(a.d)
  xy <- t(combn(colnames(m), 2))
  toplot<-data.frame(xy, dist=m[xy])
  
  tmp<-dplyr::filter(pheno2,Cluster==1)
  jap<-dplyr::filter(toplot,X1 %in% tmp$sampleID)
  jap<-dplyr::filter(toplot,X2 %in% tmp$sampleID)
  colnames(jap)[1]<-"sampleID"
  jap<-merge(jap,pheno2,by="sampleID",all.x = TRUE)
  colnames(jap)[1]<-"X1"
  colnames(jap)[2]<-"sampleID"
  jap<-merge(jap,pheno2,by="sampleID",all.x=TRUE)
  colnames(jap)[2]<-"X2"
  jap$within<-jap$site.x==jap$site.y
  jap<-jap[,c(1:3,16)]
  
  
  tmp<-dplyr::filter(pheno2,Cluster==2)
  boh<-dplyr::filter(toplot,X1 %in% tmp$sampleID)
  boh<-dplyr::filter(toplot,X2 %in% tmp$sampleID)
  colnames(boh)[1]<-"sampleID"
  boh<-merge(boh,pheno2,by="sampleID",all.x = TRUE)
  colnames(boh)[1]<-"X1"
  colnames(boh)[2]<-"sampleID"
  boh<-merge(boh,pheno2,by="sampleID",all.x=TRUE)
  colnames(boh)[2]<-"X2"
  boh$within<-boh$site.x==boh$site.y
  boh<-boh[,c(1:3,16)]
  
  tmp<-toplot
  colnames(tmp)[1]<-"sampleID"
  tmp<-merge(tmp,pheno2,by="sampleID",all.x = TRUE)
  colnames(tmp)[1]<-"X1"
  colnames(tmp)[2]<-"sampleID"
  tmp<-merge(tmp,pheno2,by="sampleID",all.x=TRUE)
  colnames(tmp)[2]<-"X2"
  tmp$within<-"Between species"
  tmp$species<-"Between species"
  tmp$withinspecies<-tmp$Cluster.x==tmp$Cluster.y
  tmp<-tmp[,c(1:3,16:18)]
  tmp<-dplyr::filter(tmp,withinspecies==FALSE)
  
  
  jap$species<-"japonica"
  boh$species<-"bohemica"
  jap$withinspecies<-TRUE
  boh$withinspecies<-TRUE
  
  tp2<-rbind(jap,boh)
  
  tp2$within<-gsub("TRUE","Within site",tp2$within)
  tp2$within<-gsub("FALSE","Among sites",tp2$within)
  tp2$within<-factor(tp2$within)
  tp2$within<-relevel(tp2$within,"Within site")
  
  library(ggplot2)
  library(ggthemes)
  library(cowplot)
  library(ggsci)
  
  withinspecies<-ggplot(tp2,aes(x=1-dist)) + 
    geom_histogram(binwidth = 0.0005,aes(fill=species),position="dodge",alpha=0.75) + 
    facet_grid(within~.,scales = "free_y") +
    xlim(0.98,1) +
    theme_bw() + scale_fill_tableau(name="Species") +
    labs(x="Jaccard similarity",y="Frequency") +
    theme(legend.position = "bottom")
  
  betweenspecies<-ggplot(tmp,aes(x=1-dist)) + 
    geom_histogram(binwidth = 0.0001,aes(fill=species),position="dodge",alpha=0.75) + 
    facet_grid(within~.) +
    xlim(0.98,1) +
    theme_bw() + scale_fill_few(guide=FALSE) +
    labs(x="Jaccard similarity",y="Frequency")
  
  pdf("Histogram3.pdf",height=6,width=6)
  cowplot::plot_grid(withinspecies,betweenspecies,nrow=2,rel_heights = c(1.2,1))
  dev.off()
  
  
  
  
  a.d<-dist(t(allele[,-c(1:4)]),method = "binary")

  pheno2[pheno2$Cluster==1,6]<-"japonica"
  pheno2[pheno2$Cluster==2,6]<-"bohemica"
  pheno2 <- pheno2[order(pheno2$species),] 
  ordered<-unique(pheno2$sampleID)
  
  library(reshape2)
  toplot <- melt(data.matrix(a.d))
  head(toplot)
  
  toplot$Var1<-factor(toplot$Var1,levels = ordered)
  toplot$Var2<-factor(toplot$Var2,levels = ordered)
  
  
  ggplot(data = toplot, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  
  
  
  colnames(toplot)[1]<-"sampleID"
  tp2<-merge(toplot,pheno2,by="sampleID",all=TRUE)
  tp2<-tp2[complete.cases(tp2),]
  
  ph3<-pheno2[,c(1,7)]
  ph3<-ph3[!duplicated(ph3),]
  
  a <- ifelse(ph3$Cluster == 1, "red", "blue")
  
  tp2$value<-1-tp2$value
  
  pdf("Heatmap.jaccard.pdf",height=7,width=7)
  ggplot(tp2,aes(x=sampleID,y=Var2)) + geom_tile(aes(fill=(value))) +
    scale_fill_continuous_tableau(name="Jaccard similarity") +
    labs(x="Sample") +
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = a)) +
    theme(axis.text.y = element_text(colour = a)) +
    theme(axis.title.y=element_blank())
  dev.off()
}







#### Analyses

# This function contains the analyses performed for the genetic information: RDA, outlier tests, and principal components analysis
GeneticAnalysis<-function(allele=allele){
  
  # PCA for visualization
  
  library(flashpcaR)
  
  Y <- allele[,-c(1:4)]
  pc <- prcomp(t(Y))
  
  pheno$sampleID<-gsub("FALLOPIA_","",pheno$sampleID)
  
  
  library(vegan)
  PCdata<-data.frame(pheno,pc$x)
  library(ggplot2)
  library(ggthemes)
  library(ggrepel)
  
  PCA<-ggplot(PCdata,aes(x=PC1,y=PC2,color=treatment,label=sampleID)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) + labs(title="A") +
    scale_color_fivethirtyeight() + 
    theme_bw() +
    guides(color=guide_legend(title="Habitat"))
  
  # IBD
  Y <- data.frame(t(allele[,-c(1:4)]))
  
  fit<-rda(Y~lat+long,data=pheno)
  anova(fit,parallel=6,permutations = 9999,by="marg")
  #Lat varExp=2.621 F=1.2960 P=0.1259 #Japonica: varExp=6.635 F=1.2276 P=0.2595 #Bohemica: varExp=15.355 F=3.3433 P=0.0002
  #Long varExp=2.807 F=1.3878 P=0.0906 #Japonica: varExp=13.667 F=2.5285 P=0.0266 #Bohemica: varExp=7.699 F=1.6764 P=0.0192
  sc<-scores(fit,choices=1:5,display="sites")
  scdata<-data.frame(pheno,RDA1=sc[,1],RDA2=sc[,2],pc$x)
  
  IBD<-ggplot(scdata,aes(x=lat,y=long,color=treatment,label=sampleID)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) + labs(title="B") +
    scale_color_fivethirtyeight() + theme_bw() + guides(shape=FALSE,color=FALSE)
  IBD_RDA<-ggplot(scdata,aes(x=RDA1,y=RDA2,color=treatment,label=sampleID)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) + labs(title="B") +
    scale_color_fivethirtyeight() + theme_bw() + guides(shape=FALSE,color=FALSE)
  
  library(cowplot)
  pdf("IBD.pdf",height=5,width=10)
  plot_grid(IBD,IBD_RDA,ncol=2)
  dev.off()
  
  # RDA as both a test for differentiation and a genome scan
  
 
  fit<-rda(Y~treatment+Condition(LF.allele),data=pheno)
  sc<-scores(fit,choices=1:5,display="sites")
  scdata<-data.frame(pheno,RDA1=sc[,1],RDA2=sc[,2],pc$x)
  
  # Permutation test - takes a few min
  anova(fit,parallel=6,permutations = 9999,by="term") # varExp = 3.2%, F=1.0953,P=0.0511 # japonica: varExp = 8.2643%, F=1.0078,P=0.4781 # bohemica: varExp = 10.361%, F=1.1817,P=0.1307
  anova(fit,parallel=6,permutations = 9999,by="axis") 
  # RDA1: varExp = 1.7%, F=1.2049,P=0.0434 # japonica: varExp = 4.4019%, F=1.0736,P=0.7724  # bohemica: varExp = 5.636%, F=1.2856,P=0.2161 
  # RDA2: varExp = 1.45%, F=0.9858,P=0.57 # japonica: varExp = 3.8624%, F=1.0078,P=0.5813 # bohemica: varExp = 4.725%, F=1.0779,P=0.4016 
  
  # differences in dispersion
  tempdist<-dist((Y))
  #colnames(genomacau)<-colnames(allele[,-c(1:4)])
  beta.test<-betadisper(d=tempdist,group=as.factor(pheno$treatment),bias.adjust = TRUE)
  anova(beta.test)
  
  #pdf("RDA.genetic.pdf",height=6,width=6)
  RDA1<-ggplot(scdata,aes(y=RDA1,x=PC1,color=treatment,label=sampleID)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) + labs(title="B") +
    scale_color_fivethirtyeight() + theme_bw() + guides(shape=FALSE,color=FALSE)
  RDA2<-ggplot(scdata,aes(y=RDA2,x=PC1,color=treatment,label=sampleID)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) + labs(title="C") +
    scale_color_fivethirtyeight() + theme_bw() + guides(shape=FALSE,color=FALSE)
  RDA_12<-ggplot(scdata,aes(x=RDA1,y=RDA2,color=treatment,label=sampleID)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) + labs(title="D") +
    scale_color_fivethirtyeight() + theme_bw() + guides(shape=FALSE,color=FALSE)
  
  #dev.off()
  
  
  # LFMM outliers - this test is overly conservative from P-values histograms, do not use
  
  library(lfmm)
  modmat<-model.matrix(~treatment,data=pheno)
  mod.lfmm <- lfmm_ridge(Y = (Y), 
                         X = modmat, 
                         K = 3)
  pv <- lfmm_test(Y = Y, 
                  X = modmat, 
                  lfmm = mod.lfmm, 
                  calibrate = "NULL")
  
  pvalues <- pv$pvalue 
  qqplot(rexp(length(pvalues), rate = log(10)),
         -log10(pvalues), xlab = "Expected quantile",
         pch = 19, cex = .4)
  abline(0,1)
  
  #pv$B
  pvalues <- data.frame(allele[,1:2],pv$pvalue)
  colnames(pvalues)[4:5]<-c("beach_marsh","beach_road")
  
  tmp<-factor(pheno$treatment)
  tmp<-relevel(tmp,"marsh")
  modmat<-model.matrix(~tmp,data=pheno)#[,-1]
  mod.lfmm <- lfmm_ridge(Y = Y, 
                         X = modmat, 
                         K = 3)
  pv <- lfmm_test(Y = Y, 
                  X = modmat, 
                  lfmm = mod.lfmm, 
                  calibrate = "NULL")
  #colnames(pvalues)[4:5]<-c("beach_marsh","beach_road")
  pvalues$marsh_road<-pv$pvalue[,3]

  hist(pvalues$beach_marsh)
  hist(pvalues$beach_road)
  hist(pvalues$marsh_road)
  
  
  library(qvalue)
  pvalues$BM_q<-qvalue(pvalues$beach_marsh)$qvalues
  pvalues$BR_q<-qvalue(pvalues$beach_road)$qvalues
  pvalues$MR_q<-qvalue(pvalues$marsh_road)$qvalues
  
  pvalues$BM_sig<-pvalues$BM_q<=0.05
  pvalues$BR_sig<-pvalues$BR_q<=0.05
  pvalues$MR_sig<-pvalues$MR_q<=0.05
  
  
  table(pvalues$BM_sig)
  table(pvalues$BR_sig)
  table(pvalues$MR_sig)
  
  Sig2<-reshape2::melt(pvalues,id.vars=colnames(pvalues)[-c(4:6)])
  Sig2$ID<-paste(Sig2$Chr,Sig2$Pos,sep=":")
  
  #SigLoci<-dplyr::filter(pvalues,MR_sig==TRUE)
  fwrite(pvalues,"Genetic.sig.txt",quote=FALSE,row.names = FALSE,sep="\t")
  
  
  # Analysis of outliers - but there are none
  
  SigLoci<-pvalues#[,-c(3:6)]
  
  # get change through simple regression
  GetChange<-function(){
    library(car)
    atemp<-data.frame(t(allele[,-c(1:4)]))
    atemp<-data.frame(sampleID=rownames(atemp),atemp)
    atemp$sampleID<-gsub("FALLOPIA_","",atemp$sample)
    pheno$sampleID<-gsub("FALLOPIA_","",pheno$sample)
    atemp<-merge(pheno,atemp,by="sampleID")
    atemp<-atemp[,-c(1)]
    X=model.matrix(~atemp[,2])
    Loop<-function(loc){
      tmpframe<-data.frame(atemp[,c(2,loc)])
      tmpframe[,1]<-as.factor(tmpframe[,1])
      lmfit<-lm(atemp[,loc]~treatment,data=tmpframe)
      #coef(lmfit)
      tmp<-data.frame(ChrPos=paste(allele[loc-2,1],allele[loc-2,2],sep=":"),Beach_Marsh=coef(lmfit)[2],Beach_Road=coef(lmfit)[3])
      tmpframe[,1]<-relevel(tmpframe[,1],ref="road")
      lmfit<-lm(atemp[,loc]~treatment,data=tmpframe)
      #coef(lmfit)
      tmp$Road_Marsh<-coef(lmfit)[3]
      #lmfit<-.lm.fit(x=X,y=atemp[,loc])
      #lmfit$coefficients
      #lmfit2$coefficients
      #cor(atemp[,i],as.numeric(as.factor(atemp$treatment)))
      #lmfit<-lm(atemp[,loc]~treatment,data=atemp)
      #tmp<-Anova(lmfit)
      return(tmp)
    }
    library(parallel)
    #Avg<-mcmapply(X=3:ncol(atemp),FUN=Loop,mc.cores = 6)
    Avg<-mclapply(6:ncol(atemp),Loop,mc.cores = 6)
    Avg<-rbindlist(Avg)
    colnames(Avg)[1]<-"ID"
    #atemp<-data.frame(allele,Change=Avg)
    #atemp<-atemp[,c(1:2,ncol(atemp))]
    return(Avg)
    
  }
  atemp<-GetChange()
  
  
  SigLoci$ID<-paste(SigLoci$Chr,SigLoci$Pos,sep=":")
  #pvalues$ID<-paste(pvalues$Chr,pvalues$Pos,sep=":")
  #atemp$ID<-paste(atemp$Chr,atemp$Pos,sep=":")
  
  SigLoci<-merge(SigLoci,atemp,by="ID",all.x=TRUE)
  #pvals2<-merge(pvalues,atemp,by="ID")
  #pvals2<-pvals2[,c(1:3,8:13)]

  
  # just counting stuff
  #sum(abs(Sig2$value)>0.05 & Sig2$sig == TRUE)
  #sum(abs(Sig2$value)>0.2 & Sig2$sig == TRUE)
  #sum(Sig2$value>0 & Sig2$sig == TRUE)
  #sum(Sig2$value<0 & Sig2$sig == TRUE)
  #sum(Sig2$value>0 & Sig2$sig == TRUE & abs(Sig2$value)>0.05)
  #sum(Sig2$value<0 & Sig2$sig == TRUE & abs(Sig2$value)>0.05)
  #sum(Sig2$value>0 & Sig2$sig == TRUE & abs(Sig2$value)>0.2)
  #sum(Sig2$value<0 & Sig2$sig == TRUE & abs(Sig2$value)>0.2)
  
  library(cowplot)
  library(gridExtra)
  pdf("Genetic.plots.pdf",height=5,width=10)
  plot_grid(PCA,RDA1,rel_widths = c(1.2,1))
  #grid.arrange(arrangeGrob(PCA, RDA1, ncol = 2,widths=c(1.2,1)), # Second row with 2 plots in 2 different columns
  #             GenPlot,                             # First row with one plot spaning over 2 columns
  #             nrow = 2)                       # Number of rows
  dev.off()
  
  #pdf("Genetic.plots.pdf",height=7,width=10)
  #plot_grid(PCA, RDA1, RDA2,RDA_12,rel_widths = c(1.2,1))
  #grid.arrange(arrangeGrob(PCA, RDA1, RDA2,RDA_12, ncol = 2,widths=c(1.2,1)), # Second row with 2 plots in 2 different columns
  #             GenPlot,                             # First row with one plot spaning over 2 columns
  #             nrow = 2)                       # Number of rows
  #dev.off()
  
}
GeneticAnalysis()

# This function contains the analyses performed for the methylation information: RDA (both conditional on genetic structure and unconditional) and principal components analysis
MethylAnalysis<-function(meth=meth){
  
  library(vegan)
  library(ggplot2)
  library(ggthemes)
  library(ggrepel)
  
  Y <- meth[,-c(1:4)]
  pc <- prcomp(t(Y))
  PCdata<-data.frame(pheno,pc$x)
  
  #pdf("PCA.methylation.pdf",height=6,width=6)
  PCA<-ggplot(PCdata,aes(x=PC1,y=PC2,color=treatment,label=sampleID)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) +
    scale_color_fivethirtyeight() + theme_bw() + 
    guides(color=guide_legend(title="Habitat"))
  #dev.off()
  
  genpc<-prcomp(t(allele[,-c(1:4)]))
  genpc<-genpc$x
  
  Y <- data.frame(t(meth[,-c(1:4)]))
  fit.uncontrolled<-rda(Y~treatment+Condition(LF.meth),data=pheno)
  fit<-rda(Y~treatment+Condition(LF.meth+genpc[,1:5]),data=pheno)
  sc<-scores(fit,choices=1:5,display="sites")
  scdata<-data.frame(pheno,sc)
  
  # Memory intensive
  anova(fit.uncontrolled,permutations = 9999,by="terms",parallel=6) # varExp = 3.377, F=1.0851, P=0.057 
  #Japonica: varExp = 11.190, F=1.1849, P=0.2359
  #Bohemica: varExp = 10.163, F=1.1879, P=0.1692
  anova(fit,permutations = 9999,by="terms",parallel=6) # varExp = 3.149, F=1.0293, P=0.3922
  # Japonica: varExp = 8.6865, F=0.9091, P=0.5674
  # Bohemica: varExp = 8.5572, F=1.0356, P=0.4656
  anova(fit,permutations = 9999,by="axis",parallel=6) 
  # RDA1: varExp = 1.646%, F=1.0761,P=0.6596 # Japonica: varExp = 5.0624%, F=1.0596,P=0.7117 # Bohemica: varExp = 4.5955%, F=1.1123,P=0.7043
  # RDA2: varExp = 1.503%, F=0.9825,P=0.5710 # Japonica: varExp = 3.6242%, F=0.7586,P=0.5257 # Bohemica: varExp = 3.9617%, F=0.9589,P=0.4642
  
  #pdf("RDA.meth.pdf",height=6,width=6)
  RDA<-ggplot(scdata,aes(y=RDA1,x=PC1,color=treatment,label=sampleID)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) +
    scale_color_fivethirtyeight() + theme_bw()+ guides(shape=FALSE,color=FALSE)
  #dev.off()
  
  library(cowplot)
  pdf("Methyl.plots.pdf",height=5,width=10)
  plot_grid(PCA,RDA,rel_widths = c(1.2,1))
  dev.off()
  
  #sc<-scores(fit,choices=1:5,display="species")
  #scdata<-data.frame(Locus=1:nrow(meth),meth[,c(1:4)],sc)
  #std<-sd(scdata$RDA1)
  #scdata$sig<-abs(scdata$RDA1)>(std*3)
  
  #SigLoci<-dplyr::filter(scdata, sig==TRUE)
  fwrite(scdata,"Methylation.RDA.csv",quote=FALSE,row.names = FALSE)
  
  library(lfmm)
  modmat<-model.matrix(~treatment,data=pheno)
  mod.lfmm <- lfmm_ridge(Y = Y, 
                         X = modmat, 
                         K = 4)
  pv <- lfmm_test(Y = Y, 
                  X = modmat, 
                  lfmm = mod.lfmm, 
                  calibrate = "NULL")
  #pv$B
  pvalues <- data.frame(meth[,1:2],pv$pvalue)
  
  tmp<-factor(pheno$treatment)
  tmp<-relevel(tmp,"marsh")
  modmat<-model.matrix(~tmp,data=pheno)#[,-1]
  mod.lfmm <- lfmm_ridge(Y = Y, 
                         X = modmat, 
                         K = 4)
  pv <- lfmm_test(Y = Y, 
                  X = modmat, 
                  lfmm = mod.lfmm, 
                  calibrate = "NULL")
  colnames(pvalues)[4:5]<-c("beach_marsh","beach_road")
  pvalues$marsh_road<-pv$pvalue[,3]
  library(qvalue)
  
  pvalues$BM_q<-qvalue(pvalues$beach_marsh)$qvalues
  pvalues$BR_q<-qvalue(pvalues$beach_road)$qvalues
  pvalues$MR_q<-qvalue(pvalues$marsh_road)$qvalues
  
  pvalues$BM_sig<-pvalues$BM_q<=0.05
  pvalues$BR_sig<-pvalues$BR_q<=0.05
  pvalues$MR_sig<-pvalues$MR_q<=0.05
  
  
  #SigLoci<-dplyr::filter(pvalues,MR_sig==TRUE)
  fwrite(pvalues,"Methylation.sig.txt",quote=FALSE,row.names = FALSE,sep="\t")
  
  # look for context levels
  
  CG<-dplyr::filter(meth,context=="CG")
  CHG<-dplyr::filter(meth,context=="CHG")
  CHH<-dplyr::filter(meth,context=="CHH")
  
  MethCounter<-function(dat){
    sums<-colSums(dat[,-c(1:4)])
    sums<-sums/nrow(dat)
    sums<-data.frame(Context=unique(dat[,3]),t(sums))
    return(sums)
  }
  CG.s<-MethCounter(CG)
  CHG.s<-MethCounter(CHG)
  CHH.s<-MethCounter(CHH)
  
  MethSums<-rbind(CG.s,CHG.s,CHH.s)
  MethSums2<-data.frame(t(MethSums))
  colnames(MethSums2)<-c("CG","CHG","CHH")
  MethSums2<-MethSums2[-1,]
  MethSums<-data.frame(sampleID=rownames(MethSums2),MethSums2)
  rm(MethSums2)
  MethSums<-merge(pheno,MethSums,by="sampleID")
  MethSums[7:9] <- lapply(MethSums[7:9], as.character)
  MethSums[7:9] <- lapply(MethSums[7:9], as.numeric)
  MethSums<-data.frame(MethSums,genpc[,1:5],LF.meth)
  
  library(car)
  hist(MethSums$CHH)
  modmat<-model.matrix(~treatment+lat+long+genpc[,1:5]+LF.meth,data=MethSums)
  fit.cg<-lm(CG~treatment+PC1+PC2+PC4+PC5+X1+X2+X3+X4,data=MethSums)
  fit.chg<-lm(CHG~treatment+PC1+PC2+PC4+PC5+X1+X2+X3+X4,data=MethSums)
  fit.chh<-lm(CHH~treatment+PC1+PC2+PC4+PC5+X1+X2+X3+X4,data=MethSums)
  
  Anova(fit.cg)
  Anova(fit.chg)
  Anova(fit.chh)
  
  MethSums<-reshape2::melt(MethSums,id.vars=colnames(MethSums[,-c(7:9)]))
  
  library(ggthemes)
  library(ggsci)
  
  ggplot(MethSums,aes(x=variable,y=value)) + 
    geom_boxplot(aes(fill=treatment)) +
    theme_bw() + scale_fill_fivethirtyeight(name="Habitat") +
    labs(x="Methylation context",y="Methylation frequency") +
    theme(legend.position = "bottom")
  
  
  # look for context overrepresentation
  
  #M <- data.frame(table(meth$context))
  #M<-cbind(M,Obs=data.frame(table(SigLoci$))[,2])
  #M$Freq<-as.numeric(as.character(M$Freq))
  #M$Obs<-as.numeric(as.character(M$Obs))
  #M2<-data.frame(t(M))
  #colnames(M2)<-c("Unknown","CG","CHG","CHH")
  #M2<-M2[-1,]
  #for(i in 1:ncol(M2)){
  #  M2[,i]<-as.numeric(as.character(M2[,i]))
  #}
  #M <- as.table(as.numeric(M2[1,]),as.numeric(M2[2,]))
  #rownames(M2)<-c("Expected","Observed")
  #chisq.test((M))
}
MethylAnalysis

# This checks the epiGBS reads vs the assembly
FragmentAnalysis(){
  #setwd("/Volumes/DonohueLab1/MFA/2016Fallopia_epiGBS")
  frags<-fread("Consensus_blast.txt")
  colnames(frags)<-c("queryID","seqID","pctIdent","alignmentLength",
                          "mismatch","gapopen","queryStart","queryEnd","seqStart","seqEnd","evalue","bitscore")
  frags<-dplyr::filter(frags,evalue<1e-30)
  frags<- group_by(frags, queryID)
  frags<- dplyr::slice(frags,which.min(evalue))
  length(unique(frags$queryID))
  length(unique(frags$seqID))
  coverage<-data.frame(table(frags$seqID))
  nrow(dplyr::filter(coverage,Freq>1))
  nrow(dplyr::filter(frags,seqStart<=200))
  fwrite(frags,"frags.filtered.csv",quote=FALSE,row.names = FALSE)
  
}


# This function contains the methylation figures
Figures<-function(){
  library(ggthemes)
  library(ggsci)
  ### Methylation
  
  MethSig<-fread("methylation.all.csv",data.table = FALSE)
  MethRDA<-fread("methylation.rda.csv",data.table = FALSE)
  
  library(flashpcaR)
  PCmeth<-flashpca(t(meth[,-c(1:4)]),stand="none")
  PCmeth$vectors
  PCdata<-data.frame(pheno,PCmeth$vectors)
  
  FixPoly<-function(){
    temp<-dplyr::filter(MethSig,!context %in% c("CG","CHG","CHH"))
    temp[,4]<-"Polymorphic"
    temp2<-dplyr::filter(MethSig,context %in% c("CG","CHG","CHH"))
    temp<-rbind(temp,temp2)
    return(temp)
  }
  MethSig<-FixPoly()
  
  scdata<-data.frame(pheno,RDA1=MethRDA$RDA1,PCmeth$vectors)
  
  library(ggrepel)
  PCA<-ggplot(PCdata,aes(x=X1,y=X2,color=treatment,label=sampleID, shape=site)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) + labs(x="PC1",y="PC2") +
    scale_color_fivethirtyeight() + theme_bw() + 
    guides(color=guide_legend(title="Exposure"),shape=guide_legend(title="Site"))
  
  RDA<-ggplot(scdata,aes(x=X1,y=RDA1,color=treatment,label=sampleID, shape=site)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) + labs(x="PC1",y="RDA1") +
    scale_color_fivethirtyeight() + theme_bw()+ guides(shape=FALSE,color=FALSE)
  
  library(ggplot2)
  
  contexts_frags<-data.frame(table(MethSig$chr,MethSig$context))
  contexts_frags$Var2<-as.character(contexts_frags$Var2)
  #contexts_frags[1:6279,2]<-"Polymorphic"
  #contexts_frags$Var2<-gsub(".","Polymorphic",contexts_frags,fixed = TRUE)
  #contexts_frags<-dplyr::filter(contexts_frags, Var2 %in% c("CG","CHG","CHH"))
  
  # context distribution
  Methy_cxt<-ggplot(contexts_frags,aes(y=as.numeric(as.factor(Var1)),x=1,fill=Freq)) +
    geom_raster() + scale_fill_distiller(type="seq",palette = 7,direction=1,name="Frequency") +
    facet_grid(~as.factor(Var2)) +
    theme_bw() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()
    ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "bottom") +
    labs(y="Consensus fragment number",x="Trinucleotide context",title="A")
  
  # Methylation outliers
  MethSig$sig<-MethSig$qval<=0.05
  MethyOut<-ggplot(MethSig,aes(y=-log10(pvalue),x=Change,color=sig)) + 
    geom_point(alpha=0.5) + theme_bw() +
    facet_wrap(~context) + 
    geom_vline(xintercept = 0.2,linetype="dashed",alpha=0.5) + 
    geom_vline(xintercept = -0.2,linetype="dashed",alpha=0.5) +
    #geom_vline(xintercept = 0) +
    labs(x="Change in methylation frequency",y="Negative log10 P-value")+
    theme(legend.position = "bottom") +
    scale_color_hc(name="Significant?")
  
  library(gridExtra)
  library(cowplot)
  
  pdf("Methylation.fig.pdf",height=9,width=9)
  #plot_grid(Methy_cxt,MethyOut,align="v",nrow=2)
  plot1<-plot_grid(PCA,RDA,ncol=2,labels = "AUTO", rel_widths = c(1.2,1))
  plot2<-plot_grid(plot1,MethyOut,nrow=2,labels=c("","C"))
  print(plot2)
  dev.off()
  
}


#### Subsets

setwd("/Volumes/ANALYSIS3/2018Knotweed/japonicasubset")

phenoJ<-dplyr::filter(pheno,species=="japonica")
alleleJ<-dplyr::select(allele, one_of(phenoJ$sampleID))
alleleJ<-data.frame(allele[,c(1:4)],alleleJ)
methJ<-dplyr::select(meth, one_of(phenoJ$sampleID))
methJ<-data.frame(meth[,c(1:4)],methJ)

allele<-alleleJ
meth<-methJ
pheno<-phenoJ

LF.allele<-GetLFs(df=alleleJ,K=3) 
LF.meth<-GetLFs(df=methJ,K=4)


setwd("/Volumes/ANALYSIS3/2018Knotweed")
library(tidyverse)
options(contrasts = c("contr.treatment","contr.poly"))  # The default


# Bohemica subset


setwd("/Volumes/ANALYSIS3/2018Knotweed/bohemicasubset")

phenoB<-dplyr::filter(pheno,species=="bohemica")
alleleB<-dplyr::select(allele, one_of(phenoB$sampleID))
alleleB<-data.frame(allele[,c(1:4)],alleleB)
methB<-dplyr::select(meth, one_of(phenoB$sampleID))
methB<-data.frame(meth[,c(1:4)],methB)

allele<-alleleB
meth<-methB
pheno<-phenoB

LF.allele<-GetLFs(df=alleleB,K=3) 
LF.meth<-GetLFs(df=methB,K=4)
