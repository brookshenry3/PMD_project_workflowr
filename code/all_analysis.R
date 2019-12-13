
#This script contains another analysis of the regions, however this time I am using a GRanges/DF that 
#contains all regions of the genomes broken up into either PMD or non-PMD regions

#Important to know that ALL of this analysis is done in PRAD

####Importing PMDs and subsetting####

PMDs <- import.bedGraph("data/raw/PMDs/PMD_coordinates_hg19.bed.gz", genome = "hg19")

PMDs <- PMDs[, "NA.1"]
PMDs <- PMDs[ PMDs$NA.1 == "commonPMD"]
PMDs <- keepStandardChromosomes(PMDs)
PMDs.r <- reduce(PMDs, min.gapwidth=3, drop.empty.ranges = TRUE)
df.PMDs.r <- data.frame(PMDs.r)

hg19 <- GRangesForBSGenome(genome = 'hg19', chrom = NULL, ranges = NULL)
df.hg19 <- data.frame(hg19)
df.hg19 <- df.hg19[c(1:22), ]
hg19 <- makeGRangesFromDataFrame(df.hg19, keep.extra.columns = FALSE, seqnames.field = "seqnames",
                                 start.field = "start", end.field = "end")

hits <- findOverlaps(hg19, PMDs.r)
grl <- extractList(PMDs.r, as(hits, "List"))
nonPMDs <- psetdiff(hg19, grl)
df.nonPMDs <- data.frame(nonPMDs)


df.nonPMDs$region <- 'non-PMD'
df.nonPMDs <- df.nonPMDs[, c(3, 4, 5, 6, 8)]
nonPMDs <- makeGRangesFromDataFrame(df.nonPMDs, 
                                    keep.extra.columns = TRUE,
                                    seqnames.field = 'seqnames',
                                    start.field = 'start',
                                    end.field = 'end')

df.PMDs.r$region <- 'PMD'
df.PMDs <- df.PMDs.r[, c(1, 2, 3, 4, 6)]

PMDs <- makeGRangesFromDataFrame(df.PMDs, 
                                 keep.extra.columns = TRUE,
                                 seqnames.field = 'seqnames',
                                 start.field = 'start',
                                 end.field = 'end')

#The GRanges and DF below have all regions of the hg19 genome split up into either PMD or non-PMD regions, let's try running the analysis with this
#Analysis will be done over in a new script called all_analysis
all.regions <- c(PMDs, nonPMDs)
df.all <- data.frame(all.regions)

#Cleaning things up 
rm(df.hg19, df.nonPMDs, df.PMDs, df.PMDs.r, grl, hg19, hits, nonPMDs, PMDs, PMDs.r)

####Loading PRAD & making matrices######

load(file = "~/Desktop/KBH_PMD_Project_master/data/output/MAEs/PRAD-MAE")
PRAD.mae <- mae
rm(mae)

PRADmeth <- assay(PRAD.mae[, , 1])
PRADexp <- assay(PRAD.mae[, , 2])

summary(as.factor(PRAD.mae@colData@listData[["definition"]]))
summary(as.factor(PRAD.mae@colData@listData[['ajcc_clinical_t']]))
summary(as.factor(PRAD.mae@colData@listData[['ajcc_clinical_m']]))
summary(as.factor(PRAD.mae@colData@listData[['primary_gleason_grade']]))
summary(as.factor(PRAD.mae@colData@listData[['TN']]))


####Promoter probes and genes########

#This time I am going to try working exclusively with RnBeads for annotation/data processing, etc.
#Getting 450k probes, hg19 genes & probes using RnBeads

probes.450k <- rnb.annotation2data.frame(rnb.get.annotation(type = 'probes450'))
probes.450k <- probes.450k[, c(1, 2, 3, 15, 16, 22)]
probes.gr <- makeGRangesFromDataFrame(probes.450k, 
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE,
                                      seqnames.field = 'Chromosome',
                                      start.field = 'Start',
                                      end.field = 'End')

genes.hg19 <- rnb.annotation2data.frame(rnb.get.annotation(type = 'genes', assembly = 'hg19'))
genes.hg19 <- genes.hg19[, c(1, 2, 3, 5, 7, 11)]
genes.gr <- makeGRangesFromDataFrame(genes.hg19, 
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE,
                                      seqnames.field = 'Chromosome',
                                      start.field = 'Start',
                                      end.field = 'End')

promoters <- rnb.annotation2data.frame(rnb.get.annotation(type = 'promoters', assembly = 'hg19'))
promoters <- promoters[, c(1, 2, 3, 5, 7, 11)]
promoters.gr <- makeGRangesFromDataFrame(promoters,
                                         keep.extra.columns = TRUE,
                                         seqnames.field = 'Chromosome',
                                         start.field = 'Start',
                                         end.field = 'End')

#Now to filter the genes and probes to get only the ones that are in the regions objects

promoter.probes.gr <- subsetByOverlaps(probes.gr, promoters.gr)
region.prom.probes.gr <- subsetByOverlaps(promoter.probes.gr, all.regions)
region.prom.probes.anno <- splicejam::annotateGRfromGR(GR1 = region.prom.probes.gr,
                                                       GR2 = all.regions)

df.region.prom.probes <- as.data.frame(region.prom.probes.anno)
#From the above code I now have a list of the 450k probes that are within promoters and with their region (PMD vs. non-PMD) defined

df.region.prom.probes$region <- as.factor(df.region.prom.probes$region)

#Now doing the above filtering with the genes

region.genes.gr <- subsetByOverlaps(genes.gr, all.regions)
region.genes.anno <- splicejam::annotateGRfromGR(GR1 = region.genes.gr,
                                                 GR2 = all.regions)

df.region.genes <- as.data.frame(region.genes.anno)
df.region.genes$region <- as.factor(df.region.genes$region)
#Removing genes that are in overlapping regions
df.region.genes <- subset(df.region.genes, df.region.genes$region != 'non-PMD,PMD')

#Ok so now all the probes and genes are annotated to show which region they belong in

#Lastly restoring the GRanges objects after all the filtering

df.region.genes <- df.region.genes[, c(1, 2, 3, 6, 7, 8, 9)]
df.region.prom.probes <- df.region.prom.probes[, c(1, 2, 3, 6, 7, 8, 9)]
row.names(df.region.prom.probes) <- c(1:168835)

region.prom.probes.gr <- makeGRangesFromDataFrame(df.region.prom.probes, 
                                                  keep.extra.columns = TRUE,
                                                  ignore.strand = TRUE,
                                                  seqnames.field = 'seqnames',
                                                  start.field = 'start',
                                                  end.field = 'end')

region.genes.gr <- makeGRangesFromDataFrame(df.region.genes, 
                                            keep.extra.columns = TRUE,
                                            ignore.strand = TRUE,
                                            seqnames.field = 'seqnames',
                                            start.field = 'start',
                                            end.field = 'end')
  

####Basic analysis#####

#Now I can see how many genes are in each region

summary(df.region.genes$region)

#So there are 38234 gene annotations outside of PMDs and 14827 annotations inside of PMDs

#As a better method let's see gene density inside vs outside PMDs

human.genes = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
gene.density <- all.regions
gene.density$totgenes <- countOverlaps(gene.density, human.genes)

df.gene.density <- as.data.frame(gene.density)

ggplot(data=df.gene.density, aes(x=totgenes, group=region, fill=region)) +
  geom_density(adjust=1.5, alpha=.4)

#So gene density seems to be lower inside PMDs than outside
#Checking the distribution with a KS test - tests differences between distributions 

PMD.gene.density <- subset(df.gene.density, df.gene.density$region == 'PMD')
nonPMD.gene.density <- subset(df.gene.density, df.gene.density$region == 'non-PMD')

ks.test(PMD.gene.density$totgenes, nonPMD.gene.density$totgenes)

#And how many of each CGI relation are in each region

summary.regions <- group_by(df.region.prom.probes, CGI.Relation, region) %>% dplyr::summarise(count = n())
summary.regions$percent <- summary.regions$count / rep(c(141693, 27142), 6) * 100

ggplot(summary.regions, aes(fill=region, y=percent, x=CGI.Relation)) + 
  geom_bar(position="dodge", stat="identity")

PRADmeth <- subset(PRADmeth, rownames(PRADmeth) %in% df.region.prom.probes$ID)
df.region.prom.probes <- subset(df.region.prom.probes, df.region.prom.probes$ID %in% rownames(PRADmeth))

summary(df.region.prom.probes$region)

#Ordering the PRADmeth object so that the probes are in the same order as df.region.prom.probes

PRADmeth <- PRADmeth[order(match(rownames(PRADmeth), df.region.prom.probes$ID)), , drop = FALSE]

#Plotting the beta value (methylation) deviation for the PRAD data set based on region
#For some reason the below doesn't work, but if I leave out the c.values and c.legend list it does
deviation.plot.beta(betas = PRADmeth, c.values = region.list, c.legend = 2)

####Finding probe-gene pairs#######

annots <- c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
            'hg19_genes_intronexonboundaries')

annotations <- build_annotations(genome = 'hg19', annotations = annots)

region.annotated.gr <- annotate_regions(regions = region.prom.probes.gr,
                                       annotations = annotations,
                                       ignore.strand = TRUE,
                                       quiet = FALSE)

df.region.anno <- as.data.frame(region.annotated.gr)
df.region.anno <- df.region.anno[, c(1, 2, 3, 6, 8, 9, 17, 18, 19)]

#Now to map the gene symbols given above to ENSEMBL IDs

df.region.anno$ENSEMBL <- NA

gene.map <- df.region.anno[, c(8, 10)]

gene.map$ENSEMBL <- mapIds(org.Hs.eg.db,
                            keys=gene.map$annot.symbol,
                           column = "ENSEMBL",
                            keytype="SYMBOL",
                            multiVals="first")

ensembl.ids <- as.factor(as.character(gene.map$ENSEMBL))

df.region.anno$ENSEMBL <- ensembl.ids

df.region.anno$annot.type <- as.factor(df.region.anno$annot.type)

region.annotated.gr <- makeGRangesFromDataFrame(df.region.anno,
                                                keep.extra.columns = TRUE,
                                                ignore.strand = TRUE,
                                                seqnames.field = 'seqnames',
                                                start.field = 'start',
                                                end.field = 'end')

#Filtering the probe gene pairs so that they are only the ones in the PRAD.mae

df.region.pairs.FINAL <- subset(df.region.anno, df.region.anno$ID %in% rownames(PRADmeth))

#Although the above annotations are "cleaner' they also lose the ENSEMBL ID, will need to figure out a way to add those in
df.region.pairs.FINAL <- subset(df.region.pairs.FINAL, df.region.pairs.FINAL$ENSEMBL %in% rownames(PRADexp))

#Changing the gene annotations so that they are in more standardized categories 

df.region.pairs.FINAL$annot.type <- as.factor(df.region.pairs.FINAL$annot.type)

rm(ensembl.ids, gene.map)

####Correlation Test#####

df.region.pairs.FINAL$r <- NA
df.region.pairs.FINAL$p.value <- NA

pb <- txtProgressBar(min = 0, max = nrow(df.region.pairs.FINAL), style = 3)

for (i in seq_len(nrow(df.region.pairs.FINAL))) {
  test <- cor.test(PRADmeth[df.region.pairs.FINAL$ID[i], ], PRADexp[df.region.pairs.FINAL$ENSEMBL[i], ])
  df.region.pairs.FINAL$p.value[i] <- test$p.value
  df.region.pairs.FINAL$r[i] <- test$estimate
  setTxtProgressBar(pb, i)
}

#save(df.region.pairs.FINAL, file = 'PRAD_all_regions_corr_test')
load(file = '~/Desktop/KBH_PMD_Project_master/data/output/Correlation_tests/PRAD_all_regions_corr_test')

#Just playing around looking at how different regions influence the correlation coefficient  

ggplot(df.region.pairs.FINAL, aes(x=region, y=r, fill=annot.type)) + 
  geom_violin()


#Ok the above figure is interesting and possibly meaningful maybe?
#Now getting the most significantly hypomethylated pairs

df.region.sig <- subset(df.region.pairs.FINAL, p.value < 0.05 & r < (-0.5))

summary(df.region.sig$region)

#So for the significant regions 1451 of them are in non-PMD regions and 176 of them are in PMDs
#For now I will focus on the ones in PMDs

df.PMD.sig <- subset(df.region.sig, df.region.sig$region == 'PMD')

####Differentially expressed gene analysis######

sig.diff <- get.diff.meth(data = PRAD.mae, 
                          group.col = "definition",
                          group1 =  "Primary solid Tumor",
                          group2 = "Solid Tissue Normal",
                          minSubgroupFrac = 0.2, # if supervised mode set to 1
                          sig.dif = 0.3,
                          diff.dir = "hypo", # Search for hypomethylated probes in group 1
                          cores = 6, 
                          dir.out ="141119_sig_diff_PRAD", 
                          pvalue = 0.01)

sig.diff <- read.csv('~/Desktop/KBH_PMD_Project_master/data/output/Methylation_Diff/141119_sig_diff_PRAD/getMethdiff.hypo.probes.significant.csv')

sig.diff.in.regions <- subset(sig.diff, sig.diff$probe %in% df.region.prom.probes$ID)

sig.diff.in.regions$region <- subset(df.region.prom.probes$region, df.region.prom.probes$ID %in% sig.diff.in.regions$probe)

#Below doesn't work :(
#sig.diff.in.regions$annot <- subset(df.region.pairs.FINAL$annot.type, df.region.pairs.FINAL$ID %in% sig.diff.in.regions$probe)

regions.sig.meth.diff.anno <- subset(df.region.pairs.FINAL, df.region.pairs.FINAL$ID %in% sig.diff.in.regions$probe)

ggplot(regions.sig.meth.diff.anno, aes(x=region, y=r, fill=annot.type)) +
  geom_violin()

#The above GG plot is the distribution of the correlation coefficients for significantly hypomethylated PRAD promoter gene-probe pairs

####ChromHMM annotations####

#Starting with PC3 (cancer cell line)

PC3_chrom_annots <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/Cell_Lines/GSE57498_PC3_ChromHMM.bed")

df.PC3.annots <- data.frame(PC3_chrom_annots$thick, PC3_chrom_annots$name)

df.PC3.annots$chr <- NA
df.PC3.annots <- df.PC3.annots[, c(5, 1, 2, 3, 4)]
#To find where the ends/beginnings of the chromosomes are:
subset(df.PC3.annots, df.PC3.annots$start == 1)
#Adding names to each of the regions:
df.PC3.annots$chr <- c(rep('chr10', 15153), rep('chr11', 18014), rep('chr12', 15980), rep('chr13', 8843), rep('chr14', 9617),
                       rep('chr15', 11256), rep('chr16', 10779), rep('chr17', 15788), rep('chr18', 6288), rep('chr19', 13102),
                       rep('chr1', 30894), rep('chr20', 10765), rep('chr21', 3759), rep('chr22', 6335), rep('chr2', 25002),
                       rep('chr3', 19347), rep('chr4', 13371), rep('chr5', 16650), rep('chr6', 16512), rep('chr7', 17503),
                       rep('chr8', 14711), rep('chr9', 12774), rep('chrM', 23), rep('chrX', 8318), rep('chrY', 152))

#Making a GRanges object from the above data frame
PC3.annots.gr <- makeGRangesFromDataFrame(df.PC3.annots, 
                                          keep.extra.columns = TRUE, 
                                          seqnames.field = "chr",
                                          start.field = "start", 
                                          end.field = "end")

#Now to do the same with the PrEC normal prostate cell line

PrEC_chrom_annots <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/Cell_Lines/GSE57498_PrEC_ChromHMM.bed")

df.PrEC.annots <- data.frame(PrEC_chrom_annots@ranges@start, PrEC_chrom_annots@ranges@width, PrEC_chrom_annots$name)
#Ok annoyingly the PrEC data is not in the same format as the PC3 data so I will have to re-organize and clean it up a bit 

df.PrEC.annots$chr <- NA
df.PrEC.annots$end <- NA
df.PrEC.annots <- df.PrEC.annots[, c(4, 1, 5, 2, 3)]
df.PrEC.annots$PrEC_chrom_annots.type <- gsub('hg19_custom_', '', df.PrEC.annots$PrEC_chrom_annots.type)
df.PrEC.annots$end <- df.PrEC.annots$PrEC_chrom_annots.ranges.start + df.PrEC.annots$PrEC_chrom_annots.ranges.width
names(df.PrEC.annots)[2] <- 'start'
names(df.PrEC.annots)[4] <- 'width'
names(df.PrEC.annots)[5] <- 'annots.type'

#To find where the ends/beginnings of the chromosomes are:
subset(df.PrEC.annots, df.PrEC.annots$PrEC_chrom_annots.ranges.start == 1)

#Adding names to each of the regions:
df.PrEC.annots$chr <- c(rep('chr10', 17898), rep('chr11', 19786), rep('chr12', 19065), rep('chr13', 9817), rep('chr14', 12644),
                        rep('chr15', 12851), rep('chr16', 13171), rep('chr17', 18042), rep('chr18', 8637), rep('chr19', 13516),
                        rep('chr1', 37335), rep('chr20', 10964), rep('chr21', 5004), rep('chr22', 8100), rep('chr2', 31022),
                        rep('chr3', 24901), rep('chr4', 18881), rep('chr5', 20390), rep('chr6', 22454), rep('chr7', 19577),
                        rep('chr8', 17142), rep('chr9', 16559), rep('chrM', 23), rep('chrX', 7809), rep('chrY', 707))

#Making a GRanges object from the above data frame
PrEC.annots.gr <- makeGRangesFromDataFrame(df.PrEC.annots, 
                                           keep.extra.columns = TRUE, 
                                           seqnames.field = "chr",
                                           start.field = "start", 
                                           end.field = "end")

rm(PC3_chrom_annots, PrEC_chrom_annots)

#Now to actually overlay the ChromHMM chromatin state annotations with the regions

#Have to get rid of sex chromosomes

chromosomes <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
                 'chr9', 'chr10', 'chr11', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22')

#Starting with PC3 annotations

region.PC3.annots <- splicejam::annotateGRfromGR(GR1 = PC3.annots.gr,
                                                 GR2 = all.regions)
region.PC3.annots@elementMetadata@listData$region <- as.factor(region.PC3.annots@elementMetadata@listData$region)

df.region.PC3.annots <- data.frame(region.PC3.annots)

df.region.PC3.annots <- subset(df.region.PC3.annots, df.region.PC3.annots$seqnames %in% chromosomes)

#Now for PrEC

region.PrEC.annots <- splicejam::annotateGRfromGR(GR1 = PrEC.annots.gr,
                                                 GR2 = all.regions)
region.PrEC.annots@elementMetadata@listData$region <- as.factor(region.PrEC.annots@elementMetadata@listData$region)

df.region.PrEC.annots <- data.frame(region.PrEC.annots)

df.region.PrEC.annots <- subset(df.region.PrEC.annots, df.region.PrEC.annots$seqnames %in% chromosomes)


#Summarizing regions in each categories 

PC3.annots.sum <- table(df.region.PC3.annots$PC3_chrom_annots.name, df.region.PC3.annots$region)
df.PC3.annots.sum <- data.frame(PC3.annots.sum)
df.PC3.annots.sum$cell_line <- 'PC3'

PrEC.annots.sum <- table(df.region.PrEC.annots$annots.type, df.region.PrEC.annots$region)
df.PrEC.annots.sum <- data.frame(PrEC.annots.sum)
df.PrEC.annots.sum$cell_line <- 'PrEC'

chrom.annots.sum <- rbind(df.PC3.annots.sum, df.PrEC.annots.sum)
chrom.annots.sum$total.annots <- c(rep(248860, 10), rep(3489, 10), rep(60094, 10), rep(293943, 10), rep(3556, 10), rep(80257, 10))
chrom.annots.sum$percent <- chrom.annots.sum$Freq / chrom.annots.sum$total.annots * 100
chrom.annots.sum <- subset(chrom.annots.sum, chrom.annots.sum$Var2 != 'non-PMD,PMD')

#Plotting the above

ggplot(chrom.annots.sum, aes(fill=Var2, y=percent, x=Var1)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Chromatin annotations in PMD and non-PMD regions for PC3 and PrEC cells") +
  facet_wrap(~cell_line) +
  theme(legend.position="top") +
  xlab("Annotation") +
  ylab('Percent total annotations') +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue")) +
  coord_flip()

####Generating random regions and making enrichment plot####

PMDs <- subset(region.annotated.gr, region.annotated.gr$region == 'PMD')
nonPMDs <- subset(region.annotated.gr, region.annotated.gr$region == 'non-PMD')
rnd.regions <-  randomizeRegions(PMDs, genome = 'hg19')

rnd.annots <- annotate_regions(regions = rnd.regions, 
                               annotations = annotations,
                               ignore.strand = TRUE)

PMD.annsum <- summarize_annotations(annotated_regions = PMDs,
                                    quiet = TRUE)

nonPMD.annsum <- summarize_annotations(annotated_regions = nonPMDs,
                                       quiet = TRUE)

rnd.annsum <- summarize_annotations(annotated_regions = rnd.annots,
                                    quiet = TRUE)

PMD.annsum$region <- 'PMD'
nonPMD.annsum$region <- 'non-PMD'
rnd.annsum$region <- 'random'

annots.sum <- rbind(PMD.annsum, nonPMD.annsum, rnd.annsum)
annots.sum$annot.type <- gsub("hg19_", "", annots.sum$annot.type)
annots.sum$total.annots <- c(rep(73093, 12), rep(425371, 12), rep(252304, 12))
annots.sum$percent <- annots.sum$n / annots.sum$total.annots * 100

ggplot(annots.sum, aes(fill=region, y=percent, x=annot.type)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Annotation enrichment in PRAD for PMD and non-PMD regions") +
  theme(legend.position="top") +
  xlab("Annotation") +
  ylab('Percent total annotations') +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue")) +
  coord_flip()

#Ok so inclusion of the random regions doesn't seem to do much for the annotation enrichment

####Heatmap stuff####

PRADmeth.r <- subset(PRADmeth, rownames(PRADmeth) %in% sig.diff.in.regions$probe)
PRADmeth.r <- PRADmeth.r[order(match(rownames(PRADmeth.r), sig.diff.in.regions$probe)), , drop = FALSE]

col.metadata <- data.frame(Definition = as.factor(PRAD.mae@colData@listData[['definition']]),
                           Stage = as.factor(PRAD.mae@colData@listData[['ajcc_clinical_t']]),
                           Primary_Gleason_Grade = as.factor(PRAD.mae@colData@listData[['primary_gleason_grade']]),
                      row.names = colnames(PRADmeth.r))

row.metadata <- data.frame(region = sig.diff.in.regions$region,
                           row.names = sig.diff.in.regions$probe)

pheatmap(PRADmeth.r,
         cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = col.metadata,
         annotation_row = row.metadata)

#The above heatmap shows methylation beta values for the top 1912 most significantly differentially methylated probes, and groups them according to region
#Patients are grouped according to primary gleason grade, stage, and whether they are tumor or normal 


#Trying a different subsetting 



PRADmeth.r2 <- subset(PRADmeth, rownames(PRADmeth) %in% df.region.sig$ID)
PRADmeth.r2 <- PRADmeth.r2[order(match(rownames(PRADmeth.r2), df.region.sig$ID)), , drop = FALSE]

row.metadata2 <- data.frame(region = df.region.sig$region,
                            CGI.relation = df.region.sig$CGI.Relation,
                            Annotation = df.region.sig$annot.type,
                            row.names = df.region.sig$ID)
#The above won't work because there are multiple copies of the same probe, so they can't be coerced into row names of the data frame
#The heatmap can be made without using the row metadata2

pheatmap(PRADmeth.r2,
         cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = col.metadata)





####Parallel####

num.cores <- 6
parallel.setup(num.cores)

parallel.isEnabled()

parallel.teardown()

#!!!!!!Remember to stop parallel processing using parallel.teardown() once I am done in today's session

####Packages####

#Setting the wd to the file on the desktop that contains the scripts and data

setwd("~/Desktop/KBH_PMD_Project_master")

#Going to try and keep all the libraries I need up here in order for clarity's sake, will update as I add more

library(rtracklayer)
library(GenomicRanges)
library(AnnotationHub)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(annotatr)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(limma)
library(GO.db)
library(regioneR)
library(dplyr)
library(RnBeads)
library(TCGA2STAT)
library(data.table)
library(stringr)
library(magrittr)
library(pheatmap)
library(CCA)
library(ChIPseeker)
library(ggbio)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(ELMER)
library(arsenal)
library(TCGAbiolinks)
library(maftools)
library(EnhancedVolcano)
library(MultiAssayExperiment)
library(EnsDb.Hsapiens.v75)
library(biomaRt)
library(doParallel)
library(foreach)
library(RColorBrewer)



