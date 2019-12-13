
#This is a workflow for re-running most of the PMD/PRAD workflow, however this time I am defining 3 regions (PMD, HMD, and neither) using the original BED file
#Important to note that this workflow uses the PRAD data set!!


####Importing PMD and subsetting####
#Importing the bed file again but NOT filtering it out for the PMDs

regions <- import.bedGraph("data/raw/PMDs/PMD_coordinates_hg19.bed.gz", genome = "hg19")

df.regions <- as.data.frame(regions)
df.regions <- df.regions[, c(1, 2, 3, 8)]
names(df.regions)[4] <- 'region'

region.length <- sum(width(regions))


hg19 <- GRangesForBSGenome(genome = 'hg19', chrom = NULL, ranges = NULL)
df.hg19 <- data.frame(hg19[c(1:22), ])
hg19 <- makeGRangesFromDataFrame(df.hg19, 
                                 keep.extra.columns = TRUE,
                                 seqnames.field = 'seqnames',
                                 start.field = 'start',
                                 end.field = 'end')

hg19length <- sum(width(hg19))
region.length / hg19length * 100

rm(df.hg19, hg19, hg19length, region.length)

regions <- makeGRangesFromDataFrame(df.regions, 
                                 keep.extra.columns = TRUE,
                                 seqnames.field = 'seqnames',
                                 start.field = 'start',
                                 end.field = 'end')

#From the above line it looks like the regions (PMDs, HMDs, and 'neither') cover 93.7% of the hg19 genome

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

#Removing SNP probes & probes in "blacklist" that Reka sent me

probe.names <- data.frame(probe = rownames(PRADmeth))

PRADmeth <- PRADmeth[c(1:394289), ]

rm(probe.names)

####Promoters, probes, and genes########

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

promoters <- rnb.annotation2data.frame(rnb.get.annotation(type = 'promoters', assembly = 'hg19'))
promoters <- promoters[, c(1, 2, 3, 5, 7, 11)]
promoters.gr <- makeGRangesFromDataFrame(promoters,
                                         keep.extra.columns = TRUE,
                                         seqnames.field = 'Chromosome',
                                         start.field = 'Start',
                                         end.field = 'End')

probes.regions <- splicejam::annotateGRfromGR(probes.gr,
                                              regions)

df.probes.regions <- as.data.frame(probes.regions)

probes.regions.genes <- splicejam::annotateGRfromGR(probes.regions,
                                                    promoters.gr)

df.probes.regions.genes <- as.data.frame(probes.regions.genes)
df.probes.regions.genes <- df.probes.regions.genes[ , c(1, 2, 3, 6, 8, 9, 10, 11)]
names(df.probes.regions.genes)[8] <- 'ENSEMBL'

df.probes.regions.genes <- subset(df.probes.regions.genes, df.probes.regions.genes$ID %in% rownames(PRADmeth) & df.probes.regions.genes$ENSEMBL %in% rownames(PRADexp))
df.probes.regions.genes <- na.omit(df.probes.regions.genes)
df.probes.regions.genes$region <- as.factor(df.probes.regions.genes$region)

#I am switching all of the subsequent analysis to the above data frame, both because it's smaller and because there are no duplicates, so the distributions won't be as skewed
#!!!!Importantly though the above are just promoter probe-gene pairs!

rm(probes.450k, probes.gr, promoters, promoters.gr, probes.regions)

####Basic analysis#####

#Now I can see how many genes are in each region

#summary(df.region.genes$region)

#So there are 38234 gene annotations outside of PMDs and 14827 annotations inside of PMDs

#As a better method let's see gene density inside vs outside PMDs

human.genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
gene.density <- regions
gene.density$totgenes <- countOverlaps(gene.density, human.genes)

df.gene.density <- as.data.frame(gene.density)

ggplot(data=df.gene.density, aes(x=totgenes, group=region, fill=region)) +
  geom_density(adjust=1.5, alpha=.4)

#So gene density seems to be lower inside PMDs than outside
#Checking the distribution with a KS test - tests differences between distributions 

PMD.gene.density <- subset(df.gene.density, df.gene.density$region == 'commonPMD')
nonPMD.gene.density <- subset(df.gene.density, df.gene.density$region == c('neither','commonHMD'))

ks.test(PMD.gene.density$totgenes, nonPMD.gene.density$totgenes)

rm(human.genes, gene.density, df.gene.density, PMD.gene.density, nonPMD.gene.density)

#And how many of each CGI relation are in each region

summary.regions <- group_by(df.probes.regions.genes, CGI.Relation, region) %>% dplyr::summarise(count = n())
summary.regions$percent <- summary.regions$count / rep(c(47109, 14339, 35007), 6) * 100

ggplot(summary.regions, aes(fill=region, y=percent, x=CGI.Relation)) + 
  geom_bar(position="dodge", stat="identity")



#Plotting the beta value (methylation) deviation for the PRAD data set based on region
#For some reason the below doesn't work, but if I leave out the c.values and c.legend list it does
#deviation.plot.beta(betas = PRADmeth, c.values = region.list, c.legend = 2)

rm(summary.regions)


#######Correlation Test########

df.probes.regions.genes$r <- NA
df.probes.regions.genes$p.value <- NA

pb <- txtProgressBar(min = 0, max = nrow(df.probes.regions.genes), style = 3)

for (i in seq_len(nrow(df.probes.regions.genes))) {
  test <- cor.test(PRADmeth[df.probes.regions.genes$ID[i], ], PRADexp[df.probes.regions.genes$ENSEMBL[i], ])
  df.probes.regions.genes$p.value[i] <- test$p.value
  df.probes.regions.genes$r[i] <- test$estimate
  setTxtProgressBar(pb, i)
}

#save(df.probes.regions.genes, file = 'promoterprobe_pairs_v2')
load(file = "~/Desktop/KBH_PMD_Project_master/data/output/Correlation_tests/promoterprobe_pairs_v2")

df.probes.regions.genes <- na.omit(df.probes.regions.genes)
df.probes.regions.genes$region <- as.factor(df.probes.regions.genes$region)

ggplot(df.probes.regions.genes, aes(x=region, y=r, fill=region)) + 
  ylim(-1, 1) +
  geom_violin()

df.pairs.FINAL <- df.probes.regions.genes


####Classifying Correlation####

#21.11.19: After meeting with Reka, moving forward I will break up genes in PMD regions into those that are not expressed, 
#negatively correlated, positively correlated, and not correlated and use annotatr to look at TF binding sites, ChromHMM annots, 
#and other thigns
#The below ranges were found using: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3576830/

Correlation <- data.frame(correlation = c('very high negative', 'high negative', 'moderate negative', 'low negative', 'none', 'low positive', 'moderate positive', 'high positive', 'very high positive'),
                          lower = c(-1, -.9, -.7, -.5, -.3, .3, .5, .7, .9),
                          upper = c(-.9, -.7, -.5, -.3, .3, .5, .7, .9, 1))

df.pairs.FINAL <- data.frame(df.pairs.FINAL,
                      'correlation' = cut(df.pairs.FINAL$r,
                                          breaks = Correlation$lower,
                                          right = T,
                                          include.lowest = T))
#unique(df.pairs.FINAL$correlation)
df.pairs.FINAL$correlation <- revalue(df.pairs.FINAL$correlation,
                               c('(-0.3,0.3]'='none', '(-0.5,-0.3]'='low negative',
                                 '(0.3,0.5]'='low positive', '(-0.7,-0.5]'='moderate negative',
                                 '(0.5,0.7]'='moderate positive', '(-0.9,-0.7]'='high negative',
                                 '[-1,-0.9]'='very high negative', '(0.7,0.9]'='high positive'))
summary(df.pairs.FINAL$correlation)



#Now subsetting to just look at PMDs

df.PMDs <- subset(df.pairs.FINAL, df.pairs.FINAL$region == 'commonPMD' & df.pairs.FINAL$p.value < (0.05))

df.hi.PMD.corr <- subset(df.PMDs, df.PMDs$correlation == "high negative")

PMD.hi.corr.meth <- subset(PRADmeth, rownames(PRADmeth) %in% df.hi.PMD.corr$ID)


Heatmap(PMD.hi.corr.meth,
        top_annotation = ha.col,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_title = 'Probe',
        column_title = 'Sample')



rm(df.probes.regions.genes, Correlation, df.PMDs, df.hi.PMD.corr, PMD.hi.corr.meth)

####Looking at highly-correlated PMD regions####

PMDs <- subset(df.pairs.FINAL, df.pairs.FINAL$region == 'commonPMD' & df.pairs.FINAL$correlation == 'high negative')

scatter.plot(data = PRAD.mae,
             byPair = list(probe = c("cg19241352"), gene = c("ENSG00000143194")), 
             category = "definition", save = FALSE, lm_line = TRUE) 
#Notes, at this point, looking only at highly correlated PMDs, it seems like OR51E2 or TMEM26 could be potential targets 
#hypomethylation of OR51E2 promoter leads to proliferation(?) https://www.ncbi.nlm.nih.gov/pubmed/28032594

PMD.sig.exp <- subset(PRADexp, rownames(PRADexp) %in% PMDs$ENSEMBL)
PMD.sig.meth <- subset(PRADmeth, rownames(PRADmeth) %in% PMDs$ID)

PMD.sig.exp.hm <- Heatmap(PMD.sig.exp,
                          top_annotation = ha.col,
                          show_row_names = TRUE,
                          show_column_names = FALSE,
                          row_title = 'Gene',
                          column_title = 'Sample')

PMD.sig.exp.hm

PMD.sig.meth.hm <- Heatmap(PMD.sig.meth,
                           top_annotation = ha.col,
                           show_row_names = TRUE,
                           show_column_names = FALSE,
                           row_title = 'Probe',
                           column_title = 'Sample')

PMD.sig.meth.hm

#Trying to make the heatmaps again using pheatmap so that I can scale the rows 

pheatmap(PMD.sig.meth,
         cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = FALSE,
         annotation_col = all.col.metadata)

pheatmap(PMD.sig.exp,
         cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = FALSE,
         annotation_col = all.col.metadata)


#####Methylation SD/heatmaps####

meth.sd <- rowSds(PRADmeth, na.rm = TRUE)

#Finding the top 10,000 probes with the highest SD

top10k.meth.sd <- head(PRADmeth[order(rowSds(PRADmeth), decreasing = TRUE), ], 10000)
top10k.meth.sd <- subset(top10k.meth.sd, rownames(top10k.meth.sd) %in% df.probes.regions$ID)

#When the above is done there are 9326 probes, so some of the PRADmeth probes are not in the "regions" object

PRAD.tp <- subset(Tumor.purity, Tumor.purity$Sample.ID %in% PRAD.mae@sampleMap@listData[['primary']])
PRAD.tp$LUMP <- as.numeric(sub(",", ".", sub(".", "", PRAD.tp$LUMP, fixed=TRUE), fixed=TRUE))
row.names(PRAD.tp) <- PRAD.tp$Sample.ID
PRAD.tp <- PRAD.tp[ , c(1, 5)]
names(PRAD.tp)[1] <- 'primary'

sampleMap <- data.frame(PRAD.mae@sampleMap)
sampleMap <- subset(sampleMap, sampleMap$assay == "DNA methylation")
row.names(sampleMap) <- sampleMap$primary
sampleMap$LUMP <- NA

sampleMap$LUMP <- PRAD.tp$LUMP[match(sampleMap$primary, PRAD.tp$primary)]
sampleMap[sampleMap == 'NaN'] <- NA

all.col.metadata <- data.frame(row.names = colnames(PRADmeth),
                               Gleason = PRAD.mae@colData@listData[['primary_gleason_grade']],
                               Stage = PRAD.mae@colData@listData[['ajcc_clinical_t']],
                               Definition = PRAD.mae@colData@listData[['TN']],
                               LUMP = sampleMap$LUMP)

row.annots <- data.frame(row.names = rownames(top10k.meth.sd),
                     region = subset(df.probes.regions$region, df.probes.regions$ID %in% rownames(top10k.meth.sd)))

ha.row <- HeatmapAnnotation(
  region = row.annots$region,
  annotation_name_side = 'bottom',
  show_legend = TRUE,
  which = 'row')

ha.col <- HeatmapAnnotation(
  Definition = all.col.metadata$Definition,
  Stage = all.col.metadata$Stage,
  Gleason = all.col.metadata$Gleason,
  LUMP = all.col.metadata$LUMP,
  na_col = 'black',
  show_legend = TRUE)

top10k.hm <- Heatmap(top10k.meth.sd,
        top_annotation = ha.col,
        show_row_names = FALSE,
        show_column_names = FALSE,
        use_raster = TRUE,
        raster_device = 'png',
        row_title = 'Probe',
        column_title = 'Sample')

draw(top10k.hm + ha.row)

PMD_annots <- subset(row.annots, row.annots$region == 'commonPMD')
top10kpmd <- subset(top10k.meth.sd, rownames(top10k.meth.sd) %in% rownames(PMD_annots))

top10k.PMD.hm <- Heatmap(top10kpmd,
                     top_annotation = ha.col,
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     use_raster = TRUE,
                     raster_device = 'png',
                     row_title = 'Probe',
                     column_title = 'Sample')

top10k.PMD.hm



#Let's see what the sds look like for promoter probes (as opposed to just all probes ^)

df.pairs.FINAL$meth.sd <- rowSds(PRADmeth[df.pairs.FINAL$ID, ], na.rm = TRUE)

ggplot(df.pairs.FINAL, aes(x=region, y=meth.sd, fill=region)) + 
  geom_violin()

#Now trying to plot the above again but this time with only the probes correlated with gene expression

df.pairs.hi.corr <- subset(df.pairs.FINAL, df.pairs.FINAL$r <(-0.5) | df.pairs.FINAL$r >(0.5))

ggplot(df.pairs.hi.corr, aes(x=region, y=meth.sd, fill=region)) + 
  geom_violin()

#Ok so from the above, even highly probes highly correlated with gene expression don't appear to show any difference in the variance of the mehtylation sd

#Now to look at methylation in only PMD regions in tumor only samples

sample.map <- data.frame(sample.type = rep(as.factor(PRAD.mae@colData@listData[['TN']]), 2),
                         patient.id = PRAD.mae@sampleMap@listData[['colname']])
sample.map <- sample.map[c(1:537), ]
sample.map <- subset(sample.map, sample.map$sample.type == 'Tumor')
sample.map$patient.id <- factor(sample.map$patient.id)
sample.map$sample.type <- factor(sample.map$sample.type)
row.names(sample.map) <- sample.map$patient.id

#Finding PMDs with high correlation values

PMD.meth.corr <- subset(df.pairs.FINAL, df.pairs.FINAL$region == 'commonPMD' & df.pairs.FINAL$r <(-0.5))


PRADmeth.r <- subset(PRADmeth, colnames(PRADmeth) %in% rownames(sample.map))
PRADmeth.r <- subset(PRADmeth.r, row.names(PRADmeth.r) %in% PMD.meth.corr$ID)

#Now making the heatmap 

pheatmap(PRADmeth.r,
         cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE)

#The above heatmap shows the most highly negatively-correlated probes in PMD tumor samples 

rm(all.col.metadata, ha.row, PRAD.tp, row.annots, sample.id, top10k.meth.sd, Tumor.purity,
   sampleMap, PMD_annots, top10k.hm, top10kpmd, top10k.PMD.hm)


######ChromHMM###########

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

#Getting the chromatin regions onto the df.region.pairs.FINAL

pairs.FINAL.gr <- makeGRangesFromDataFrame(df.pairs.FINAL,
                                                  keep.extra.columns = TRUE,
                                                  ignore.strand = TRUE,
                                                  seqnames.field = 'seqnames',
                                                  start.field = 'start',
                                                  end.field = 'end')

pairs.FINAL.gr <- splicejam::annotateGRfromGR(GR1 = pairs.FINAL.gr,
                                                     GR2 = PC3.annots.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(GR1 = pairs.FINAL.gr,
                                                     GR2 = PrEC.annots.gr)

df.pairs.FINAL <- as.data.frame(pairs.FINAL.gr)

names(df.pairs.FINAL)[15] <- 'PC3.annots'
names(df.pairs.FINAL)[16] <- 'PrEC.annots'


#Last step is to get rid of overlapping regions

chrom.regions <- c('CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repressed', 'Transcribed')

df.pairs.FINAL <- df.pairs.FINAL[df.pairs.FINAL$PC3.annots %in% chrom.regions & df.pairs.FINAL$PrEC.annots %in% chrom.regions, ]

rm(chrom.regions)

#Now with the above the chromatin annotations are on the regions
#!!!!just like with the correlation test, the fact that there are multiple of the same probe means that the 
#annotation numbers/figurs might be skewed

#Plotting the above 
#Ok so the df.region.pairs.FINAL data frame contains PROMOTER probe & gene pairs for all of the regions described in the PMD bed file

PMDs <- dplyr::filter(df.pairs.FINAL, region=='commonPMD')
HMDs <- dplyr::filter(df.pairs.FINAL, region=='commonHMD')
Neither <- dplyr::filter(df.pairs.FINAL, region=='neither')

PMD.pc3.annots <- table(PMDs$PC3.annots) / nrow(PMDs) * 100
HMD.pc3.annots <- table(HMDs$PC3.annots) / nrow(HMDs) * 100
Neither.pc3.annots <- table(Neither$PC3.annots) / nrow(Neither) * 100

#Unrelated thought:
#CTCF are located at the border of lamina-associated domains, makes sense then that they are enriched in PMDs

res.aov <- aov(r ~ PC3.annots, data = df.pairs.FINAL)
summary(res.aov)
TukeyHSD(res.aov)

res.aov2 <- aov(r ~ PC3.annots + region, data = df.pairs.FINAL)
summary(res.aov2)

#Ok from the ANOVA done above it looks like region and annotation do have some influence on gene expression/promoter methylation correlation 


summary.regions <- group_by(df.pairs.FINAL, PC3.annots, region) %>%
                            dplyr::summarise(
                            count = n(),
                            mean = mean(r, na.rm = TRUE),
                            sd = sd(r, na.rm = TRUE))

sum(summary.regions[summary.regions$region=='commonHMD',]$count)
sum(summary.regions[summary.regions$region=='commonPMD',]$count)
sum(summary.regions[summary.regions$region=='neither',]$count)

summary.regions$total <- rep(c(46854, 14105, 34817), 9)
summary.regions$percent <- summary.regions$count / summary.regions$total * 100

ggplot(summary.regions, aes(fill=region, y=percent, x=PC3.annots)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Chromatin annotations in PMDs, HMDs, and neither regions for PC3 cells") +
  theme(legend.position="top") +
  xlab("Annotation") +
  ylab('Percent total annotations') +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue")) +
  coord_flip()

#Ok from the above figure it confirms what I've been seeing, PMDs are enriched in heterochromatin and CTCFs

df.sig.corr <- subset(df.pairs.FINAL, p.value < 0.05 & r < (-0.75))

#Doing the same in PrEC cells

PMD.prec.annots <- table(PMDs$PrEC.annots) / nrow(PMDs) * 100
HMD.prec.annots <- table(HMDs$PrEC.annots) / nrow(HMDs) * 100
Neither.prec.annots <- table(Neither$PrEC.annots) / nrow(Neither) * 100

summary.regions2 <- group_by(df.pairs.FINAL, PrEC.annots, region) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(r, na.rm = TRUE),
    sd = sd(r, na.rm = TRUE))

sum(summary.regions2[summary.regions2$region=='commonHMD',]$count)
sum(summary.regions2[summary.regions2$region=='commonPMD',]$count)
sum(summary.regions2[summary.regions2$region=='neither',]$count)

summary.regions2$total <- rep(c(46854, 14105, 34817), 9)
summary.regions2$percent <- summary.regions2$count / summary.regions2$total * 100

ggplot(summary.regions2, aes(fill=region, y=percent, x=PrEC.annots)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Chromatin annotations in PMDs, HMDs, and neither regions for PrEC cells") +
  theme(legend.position="top") +
  xlab("Annotation") +
  ylab('Percent total annotations') +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue")) +
  coord_flip()

#From comparing the PC3 and PrEC enrichment graphs, it maybe seems like there is a decrease in repressed annotations in PMDs from the normal cell line (PrEC) to PC3




rm(df.PC3.annots, df.PrEC.annots, PC3.annots.gr, PrEC.annots.gr, chrom.regions, HMDs, Neither, PMDs,
   res.aov, res.aov2, summary.regions, HMD.pc3.annots, Neither.pc3.annots, PMD.pc3.annots,
   HMD.prec.annots, Neither.prec.annots, PMD.prec.annots, summary.regions2)


####LUMP and GSTP1####

#LUMP (leukocytes unmethylation for purity) and GSTP1 are both ways to measure the purity of the tumor samples in the PRAD data set
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4671203/
#LUMP scores can be gotten directly from TCGAbiolinks

Tumor.purity <- Tumor.purity

PRAD.tp <- subset(Tumor.purity, Tumor.purity$Sample.ID %in% PRAD.mae@sampleMap@listData[['primary']])
PRAD.tp$LUMP <- as.numeric(sub(",", ".", sub(".", "", PRAD.tp$LUMP, fixed=TRUE), fixed=TRUE))
PRAD.tp$CPE <- as.numeric(sub(",", ".", sub(".", "", PRAD.tp$CPE, fixed=TRUE), fixed=TRUE))

ggplot(PRAD.tp, aes(y=LUMP, x=CPE)) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggtitle('PRAD CPE vs LUMP scores') +
  xlab("CPE")

#Now to look at GSTP1
#Actually don't need to do anything, can just search it up in the df.region.pairs.FINAL

GSTP1 <- subset(df.pairs.FINAL, df.pairs.FINAL$symbol == 'GSTP1')

#Can plot probe-gene correlation using ELMER

scatter.plot(data = PRAD.mae,
             byPair = list(probe = c("cg06928838"), gene = c("ENSG00000084207")), 
             category = "definition", save = FALSE, lm_line = TRUE) 

#The probe that I chose to plot against showed the strongest negative correlation, however I don't know if it is truly the best, probe-gene pair
#Ok so based on some background research the probe-gene pairs I have identified appear to be right


rm(PRAD.tp, Tumor.purity, GSTP1)

####Loading up pairs.FINAL data frame and GR####

#save(df.pairs.FINAL, file = 'df_pairs_FINAL_with_ChromHMM')
#load(file = '~/Desktop/KBH_PMD_Project_master/data/output/df_pairs_FINAL_with_ChromHMM')

#save(pairs.FINAL.gr, file = 'GR_pairs_FINAL_with_ChromHMM')
load(file = '~/Desktop/KBH_PMD_Project_master/data/output/GR_pairs_FINAL_with_ChromHMM')


####Annotating with TF binding sites####

#I am now adding transcription factor binding sites to my "master" file
#TF binding sites are coming from: http://tagc.univ-mrs.fr/remap/index.php?page=download

TF.CTCF.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_CTCF_nr_macs2_hg19_v1_2.bed.gz",
                         genome = 'hg19')

TF.AR.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_AR_nr_macs2_hg19_v1_2.bed.gz",
                         genome = 'hg19')

TF.ETS1.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_ETS1_nr_macs2_hg19_v1_2.bed.gz",
                       genome = 'hg19')

TF.TP53.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_TP53_nr_macs2_hg19_v1_2.bed.gz",
                         genome = 'hg19')

TF.SP1.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_SP1_nr_macs2_hg19_v1_2.bed.gz",
                         genome = 'hg19')

TF.FOXA1.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_FOXA1_nr_macs2_hg19_v1_2.bed.gz",
                        genome = 'hg19')

TF.RB1.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_RB1_nr_macs2_hg19_v1_2.bed.gz",
                        genome = 'hg19')

TF.STAT1.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_STAT1_nr_macs2_hg19_v1_2.bed.gz",
                        genome = 'hg19')

TF.FOXJ2.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_FOXJ2_nr_macs2_hg19_v1_2.bed.gz",
                        genome = 'hg19')

TF.TAF1.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_TAF1_nr_macs2_hg19_v1_2.bed.gz",
                         genome = 'hg19')

TF.NRF1.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_NRF1_nr_macs2_hg19_v1_2.bed.gz",
                         genome = 'hg19')

TF.REST.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_REST_nr_macs2_hg19_v1_2.bed.gz",
                         genome = 'hg19')

TF.E2F6.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_E2F6_nr_macs2_hg19_v1_2.bed.gz",
                         genome = 'hg19')

TF.CEBPA.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_CEBPA_nr_macs2_hg19_v1_2.bed.gz",
                         genome = 'hg19')

TF.SPI1.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/TFs/remap2018_SPI1_nr_macs2_hg19_v1_2.bed.gz",
                          genome = 'hg19')


#Don't actually need to convert them to data frames
#df.CTCF <- as.data.frame(TF.CTCF.gr)
#df.AR <- as.data.frame(TF.AR.gr)
#df.ETS1 <- as.data.frame(TF.ETS1.gr)

pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.AR.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.CTCF.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.ETS1.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.SP1.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.TP53.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.FOXA1.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.RB1.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.STAT1.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.FOXJ2.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.NRF1.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.TAF1.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.CEBPA.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.E2F6.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.REST.gr)
pairs.FINAL.gr <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                              TF.SPI1.gr)


df.pairs.FINAL <- as.data.frame(pairs.FINAL.gr)


df.pairs.FINAL <- df.pairs.FINAL[ , c(1:17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59)]
names(df.pairs.FINAL)[17] <- 'AR.binding.sites'
names(df.pairs.FINAL)[18] <- 'CTCF.binding.sites'
names(df.pairs.FINAL)[19] <- 'ETS1.binding.sites'
names(df.pairs.FINAL)[20] <- 'SP1.binding.sites'
names(df.pairs.FINAL)[21] <- 'TP53.binding.sites'
names(df.pairs.FINAL)[22] <- 'FOXA1.binding.sites'
names(df.pairs.FINAL)[23] <- 'RB1.binding.sites'
names(df.pairs.FINAL)[24] <- 'STAT1.binding.sites'
names(df.pairs.FINAL)[25] <- 'FOXJ2.binding.sites'
names(df.pairs.FINAL)[26] <- 'NRF1.binding.sites'
names(df.pairs.FINAL)[27] <- 'TAF1.binding.sites'
names(df.pairs.FINAL)[28] <- 'CEBPA.binding.sites'
names(df.pairs.FINAL)[29] <- 'E2F6.binding.sites'
names(df.pairs.FINAL)[30] <- 'REST.binding.sites'
names(df.pairs.FINAL)[31] <- 'SPI1.binding.sites'

df.pairs.FINAL[, 17:31][is.na(df.pairs.FINAL[, 17:31])] <- 'none'

#Have to convert them to factors after this step ^

cols <- colnames(df.pairs.FINAL[, 15:31])

df.pairs.FINAL[cols] <- lapply(df.pairs.FINAL[cols], factor)

#save(df.pairs.FINAL, file = 'df_pairs_FINAL_with_TFs')
load(file = '~/Desktop/KBH_PMD_Project_master/data/output/df_pairs_FINAL_with_TFs')

#Now looking at TF factor enrichment:

####TF binding site enrichment analysis####

#Need to do enrichment analysis on TF binding sites to see if PMDs are enriched in any in particular (maybe AR or CTCF?)

grouped.data <- data.frame(binding.site = rep(colnames(df.pairs.FINAL[17:31]), 3),
                           region = c(rep('commonHMD', 15), rep('commonPMD', 15), rep('neither', 15)),
                           count = NA,
                           total = c(rep(46854, 15), rep(14105, 15), rep(34817, 15)), 
                           percent = NA)

for (i in 17:31) {
  sum <- group_by(df.pairs.FINAL, df.pairs.FINAL[, i], region) %>%dplyr::summarise(count = n())
  grouped.data$count[c(i-16, i-1, i+14)] <- sum$count[c(1:3)]
  grouped.data$percent[c(i-16, i-1, i+14)] <- grouped.data$count[c(i-16, i-1, i+14)] / grouped.data$total[c(i-16, i-1, i+14)] * 100
}

ggplot(grouped.data, aes(fill=region, y=percent, x=binding.site)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Interesting, so the expected TFs don't show really any enrichment in PMDs
#TAF1 & NRF1 are known to only bind to non-methylated DNA, so it's interesting that they're enriched in PMDs
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5559737/
#The other interesting thing about the figure above is how the pattern of binding sites seems to fall into one of two categories, very weird


rm(cols, grouped.data, sum, i, TF.AR.gr, TF.CTCF.gr, TF.ETS1.gr, TF.TP53.gr, TF.SP1.gr, TF.FOXA1.gr, TF.RB1.gr, TF.STAT1.gr, 
   TF.FOXJ2.gr, TF.NRF1.gr, TF.TAF1.gr, TF.CEBPA.gr, TF.E2F6.gr, TF.REST.gr, TF.SPI1.gr)

####GeneWalk Prep####

#Below I am preparing lists of genes to be fed into GeneWalk (run in a Python VE from the terminal) which will perform functional analysis 

#Prep for PMDs

PMDs.all <- subset(df.pairs.FINAL, df.pairs.FINAL$region == 'commonPMD')
PMDs.all <- PMDs.all %>% arrange(correlation)
PMDs.all <- unique(PMDs.all$symbol)
PMDs.all <- data.frame(gene = PMDs.all)

write.csv(PMDs.all, "pmd_genes.csv",
          row.names = FALSE)

#Prep for HMDs

HMDs.all <- subset(df.pairs.FINAL, df.pairs.FINAL$region == 'commonHMD')
HMDs.all <- HMDs.all %>% arrange(correlation)
HMDs.all <- unique(HMDs.all$symbol)
HMDs.all <- data.frame(gene = HMDs.all)

write.csv(HMDs.all, "hmd_genes.csv",
          row.names = FALSE)

#Prep for Neither

neither.all <- subset(df.pairs.FINAL, df.pairs.FINAL$region == 'neither')
neither.all <- neither.all %>% arrange(correlation)
neither.all <- unique(neither.all$symbol)
neither.all <- data.frame(gene = neither.all)

write.csv(neither.all, "neither_genes.csv",
          row.names = FALSE)

#GeneWalk BASH code:

#1: Setting up VE:
#source ~/Desktop/VE/VE-env/bin/activate

#2: running GW
#genewalk --project hmds --genes hmd_genes.csv --id_type hgnc_symbol --nproc 8
#Change project ^,          gene list ^,          symbol type ^,      cores ^   to desired names/specs


####Reading in and visualizing GeneWalk results####

#PMD results

PMD.gw <- read.csv(file = '~/Desktop/KBH_PMD_Project_master/data/output/GeneWalk/pmds/genewalk_results.csv')

PMD.gw$mlog10padj <- -log10(PMD.gw$mean_padj)
PMD.gw$mlog10padj_error <- -log10(PMD.gw$cilow_pval) - PMD.gw$mlog10padj

#HMD results

HMD.gw <- read.csv(file = '~/Desktop/KBH_PMD_Project_master/data/output/GeneWalk/hmds/genewalk_results.csv')

HMD.gw$mlog10padj <- -log10(HMD.gw$mean_padj)
HMD.gw$mlog10padj_error <- -log10(HMD.gw$cilow_pval) - HMD.gw$mlog10padj

#Neither results 

neither.gw <- read.csv(file = '~/Desktop/KBH_PMD_Project_master/data/output/GeneWalk/neither/genewalk_results.csv')

neither.gw$mlog10padj <- -log10(neither.gw$mean_padj)
neither.gw$mlog10padj_error <- -log10(neither.gw$cilow_pval) - neither.gw$mlog10padj





####!Testing!####

PMDs <- subset(df.pairs.FINAL, df.pairs.FINAL$region == 'commonPMD' & df.pairs.FINAL$correlation == 'high negative')


scatter.plot(data = PRAD.mae,
             byPair = list(probe = c("cg19241352"), gene = c("ENSG00000143194")), 
             category = "definition", save = FALSE, lm_line = TRUE) 


PMDmeth <- subset(PRADmeth, rownames(PRADmeth) %in% PMDs$ID)
PMDexp <- subset(PRADexp, rownames(PRADexp) %in% PMDs$ENSEMBL)
PMDexp <- data.frame(PMDexp)
PMDmeth <- data.frame(PMDmeth)

PMDexp <- PMDexp[,!(names(PMDexp) %in% drop)]
PMDmeth <- PMDmeth[, !(names(PMDmeth) %in% drop)]

drop <- "TCGA.EJ.7782.01A.11R.2118.07"







####Workflowr####





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
library(ComplexHeatmap)
library(workflowr)





