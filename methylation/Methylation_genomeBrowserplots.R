
library(dplyr)
library(Gviz)
library(plotgardener)
library("foreach")  # %do%
library(GenomicRanges)
library(ggplot2)

####################################
# bigwig_to_GR 
#      Converts and merges bigwig files into a single genomic ranges object
bigwig_to_GR <- function(fl, chromosome = 'chr7', start = 25897559 , end = 25899313,
                        bigWig_dir="/g/data/zk16/projects/tammela/methylation/Tracks/")
{

    covTracksBW <- list()
    
    for(i in 1:length(fl))
#    covTracksBW <- foreach(i=1:length(fl)) %do%
    {   singleBW <- plotgardener::readBigwig(file = paste0(bigWig_dir,fl[i]),
                            chrom = chromosome,
                            chromstart = start,
                            chromend = end)
 
#        GenomicRanges::makeGRangesFromDataFrame(df = singleBW, keep.extra.columns=TRUE)
        covTracksBW[i] <- GenomicRanges::makeGRangesFromDataFrame(df = singleBW, keep.extra.columns=TRUE)
    }

    names(covTracksBW) <- fl
    covTracksBW.gr <- BRGenomics::mergeGRangesData(covTracksBW,multiplex = TRUE)
}


levelstoIDX <- function(factorObject)
{
    allLevels <- levels(factorObject)
    idx <- {}
    for(i in 1:length(allLevels))
    {
        idx <- c(idx, which(factorObject %in% allLevels[i]))
    }
    return(idx)
}


gvizGene <- function(chromosome='chr19',start=9083636, end=9087958, genome='mm10', 
                     bigWig_dir="/g/data/zk16/projects/tammela/methylation/Tracks/", 
                     groups=rep(c("adult", "young"), each = 4), saveDataPrefix='', 
                     flank=1000, Fig_title = NULL,
                     granges.HT = NULL, # granges object used as a highlight track
                     verbose=FALSE
                    )
{
    message(bigWig_dir)
    idx.levels <- {}
    if (is.factor(groups))
    {   idx.levels <- levelstoIDX(groups)
        l_Order <- levels(groups)
        groups <- factor(groups[idx.levels], levels=l_Order)  
    }
    if (verbose) 
    {    message(groups)}
    
    gtrack <- Gviz::GenomeAxisTrack(littleTicks=TRUE)
  #  itrack <- Gviz::IdeogramTrack(genome = genome, chromosome = chromosome)

    # Load UCSC gene model
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene

    UCSCgrtrack <- Gviz::GeneRegionTrack(txdb, genome = genome, chromosome = chromosome, 
                                     transcriptAnnotation = "symbol",
                                     #symbol=TRUE,
                                     name = "UCSC")

    library("AnnotationDbi")
    library("org.Mm.eg.db")


    symbols <- unlist(mapIds(org.Mm.eg.db, gene(UCSCgrtrack), "SYMBOL", "ENTREZID", multiVals = "first"))
    symbol(UCSCgrtrack) <- symbols[gene(UCSCgrtrack)]
    
    # Load data
    setwd(bigWig_dir)

    fl <- list.files(pattern = "*Cov.bw$",recursive = FALSE)
    fl <- grep(pattern = "8_25",invert = TRUE, value=TRUE,x=fl)  # Remove outlier sample
    if (is.factor(groups))
    {  fl <- fl[idx.levels]  }
    if (verbose)
    {    message(fl)}
        
        
    covTracksBW.gr <- bigwig_to_GR(fl, chromosome = chromosome, start = start - 3000 , end = end + 3000,bigWig_dir = bigWig_dir)
    
    fl <- list.files(pattern = "*methy.bw$",recursive = FALSE)
    fl <- grep(pattern = "8_25",invert = TRUE, value=TRUE,x=fl)  # Remove outlier sample
    if (is.factor(groups))
    {  fl <- fl[idx.levels]  }
    methyTracksBW.gr <- bigwig_to_GR(fl, chromosome = chromosome, start = start - 3000 , end = end + 3000,bigWig_dir=bigWig_dir)

    # "#00800f" = dark green
    cov_all <- DataTrack(range= covTracksBW.gr,   col = c("#00800f","#6fb7f8", "#FF08E8"),
                    genome=genome, name="Coverage", # col.histogram=c('blue'),
                    chromosome=chromosome)

    methyl_all <- DataTrack(range= methyTracksBW.gr, col = c("#00800f","#6fb7f8", "#FF08E8"), 
                    genome=genome, name="Methylation",# col.histogram=c('blue'),
                    chromosome=chromosome)
    
    if (saveDataPrefix != '')
    {   saveDir <- '/g/data/zk16/genomicsCore/jupyter/results/temp/'
        filename <- paste0(saveDir, saveDataPrefix,"_methylationGVIZ_Plotdata.RData")
        save(#itrack, 
            gtrack, UCSCgrtrack, cov_all, methyl_all, start, end, groups, file=filename)
    }
    
    allTracks <- list(#itrack,
        gtrack, UCSCgrtrack, cov_all, methyl_all)
    if (! is.null(granges.HT))
    {      allTracks <- HighlightTrack(trackList = list(#itrack,
                                        gtrack, UCSCgrtrack, cov_all, methyl_all),
                     granges.HT,
                     inBackground=FALSE, alpha=0.4)
     }
    
    plotTracks(allTracks, main=Fig_title,
           groups = as.character(groups),             
           type = c("a", "p", "confint"),
           aggregateGroups = TRUE,
           legend = TRUE,
          from = start-flank, 
           to = end+flank,
           fontsize=14)


 #   plotTracks(list(itrack, gtrack, UCSCgrtrack, cov_all, methyl_all), 
 #          groups = groups,
 #          type = c("a", "p", "confint"),
 #          aggregateGroups = TRUE,
 #          legend = TRUE,
 #         from = start-2000, 
 #          to = end+2000,
 #          fontsize=14)
}







# Nupr1

run1run2_dir <- '/g/data/zk16/projects/tammela/methylation/Project_13012_B/tracks_consolidated/'
gvizGene(chromosome='chr7',start=126623249, end=126625676, bigWig_dir = run1run2_dir, saveDataPrefix='AT2_Nupr1')

Lcn2

run1run2_dir <- '/g/data/zk16/projects/tammela/methylation/Project_13012_B/tracks_consolidated/'
gvizGene(chromosome='chr2',start=32384503, end=32393857, bigWig_dir = run1run2_dir, saveDataPrefix='AT2_Lcn2')  # Run 1 + 2 coverage