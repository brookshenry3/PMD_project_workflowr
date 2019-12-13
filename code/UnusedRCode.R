df.PRADprimaryexp <- tibble::rownames_to_column(df.PRADprimaryexp, var = "Gene")

#This next line will give me just the PRAD primary genes that are within the PMDs using the %in% operator 

df.PRADprimaryexpinPMDS <- df.PRADprimaryexp[df.PRADprimaryexp$Gene %in% PMD_gene_symbols, ]

#Nice, so there are 3493 genes that are located within PMDs that are present in the PRAD expression data set 

#

prad.rnaseq2.tum.norm <- TumorNormalMatch(PRADexp$dat)

df.prad.tum.norm <- as.data.frame(prad.rnaseq2.tum.norm)
df.prad.tum.norm <- tibble::rownames_to_column(df.prad.tum.norm, var = "Gene")
df.prad.matched.pmd <- df.prad.tum.norm[df.prad.tum.norm$Gene %in% PMD_gene_symbols, ]
df.prad.matched.pmd <- tibble::column_to_rownames(df.prad.matched.pmd, var = "Gene")

hist(df.prad.matched.pmd)

#TRying to plot PMD coverage using covplot function (didn't work)
arbitrary <- rep.int(500, 1858)
cov.values <- matrix(data = arbitrary, nrow = length(PMDs.r), ncol = 1)
df.PMDs.r <- as.data.frame(PMDs.r)
df.PMDs.arb <- cbind(df.PMDs.r, cov.values)
values(PMDs.r) <- df.PMDs.arb
head(PMDs.r)


covplot(PMDs.r, weightCol = "cov.values", xlab = "Chromosome Size (bp", title = "PMD coverage in hg19", lower = 450)




keepSeqlevels(PMDs.r, c("chr1", "chr2", "chr3", "chr4", "chr5", 
                        "chr6", "chr7", "chr8", "chr9", "chr10", 
                        "chr11", "chr12", "chr13", "chr14", "chr15",
                        "chr16", "chr17", "chr18", "chr19", "chr20",
                        "chr21", "chr22", "chrX", "chrY"))

#######################

#Code to get the probe annotations using RnBeads, decided to go with a different package instead 


probes <- rnb.get.annotation(type = 'probes450', assembly = 'hg19')
df.probes <- as.data.frame(probes)
dd
#THe following code only gets the probes that are within the PMD regions, after this is done 97.5% of probes are retained, so most of the probes fall within PMD regions 

PMD_with_probes <- subsetByOverlaps(probes, PMDs_annotated)

df_PMDs_probes <- as.data.frame(PMD_with_probes)

#This next line was just to see what  percent of the PMD regions are covered by probes. Obviously it is very small because the probes are only 1 bp long/ 

sum(width(PMD_with_probes)) / sum(width(PMDs.r)) * 100

#Subsetting the probes data frame so that I just get chromosome number, position, and CpG island relation

df_PMDs_probes <- df_PMDs_probes[, c(2, 4, 5, 18)]


##########################################3

##The following code is my first attempt at importing data using TCGA2STAT and then trying to create a MAE for ELMER, it didn't work so now I am trying to use TCGAbiolinks

#Was originally going to do this with TCGA2STAT, but ELMER has built in functions to retreive the data that seems easier to use/nicer
#EDIT: got error using ELMER download function, going back to TCGA2STAT -___-

PRADexp <- TCGA2STAT::getTCGA(disease = 'PRAD', data.type = 'RNASeq2', type = 'RPKM')
#PRADmeth <- TCGA2STAT::getTCGA(disease = 'PRAD', data.type = 'Methylation', type = '450K')

#Getting the methylation data takes a really long time, so to expedite it in the future I am saving the "PRADmeth" data object to the folder with all my scripts and then I can load it in the future
#save(PRADmeth, file = "PRADmeth.RData")
load(file = "PRADmeth.RData")

#Subsequently getting the tumor/normal matched data for both data sets

PRADexp.matched <- TumorNormalMatch(PRADexp$dat)
PRADmeth.matched <- TumorNormalMatch(PRADmeth$dat)

#Combining the above data sets using OMICSbind, after this I will somehow need to get samples that are the same

PRAD.combined <- OMICSBind(dat1 = PRADexp$dat, dat2 = PRADmeth$dat)


#The next lbig chunk of code is the painfully long way to get two matrices, one for expression and one for methylation data, that can be fed into createMAE in an acceptable format

Patient.ids <- colnames(PRAD.combined$X)

genes <- rownames(PRADexp$dat)

#Importantly, converting the gene symbols to Ensembl ids (as is required in create MAE) will give some NAs, that is there are some gene symbols that don't have a direct ensembl id conversion. In order to assign the 
#ids to the expression data I therefore have to remove genes that don't have an ensembl id (this seems kinda weird right?)

ensemblid <- mapIds(org.Hs.eg.db, genes, 'ENSEMBL', 'SYMBOL')

gene.map <- na.omit(data.frame(genes, ensemblid))

#I got the following ensembl names that were assigned to multiple gene symbols when I tried to assign gene map to the row names of the df.PRAD.exp data frame

ensembl.duplicates <- c('ENSG00000011454', 'ENSG00000091592', 'ENSG00000127603', 
                        'ENSG00000130035', 'ENSG00000133816', 'ENSG00000143226', 
                        'ENSG00000159674', 'ENSG00000172613', 'ENSG00000174640', 
                        'ENSG00000187510', 'ENSG00000189064', 'ENSG00000189195', 
                        'ENSG00000204131', 'ENSG00000213694', 'ENSG00000225830', 
                        'ENSG00000236362', 'ENSG00000254911', 'ENSG00000273032' )

#Some function will go here that will remove rows from the gene.map data frame that match the above ensembl ids

gene.map.final <- gene.map[ ! gene.map$ensemblid %in% ensembl.duplicates, ]


ensembl.final <- gene.map.final[,2]
gene.final <- gene.map.final[,1]

probe.ids <- rownames(PRADmeth$dat)

df.PRAD.exp <- as.data.frame(PRAD.combined$X)
df.PRAD.meth <- as.data.frame(PRAD.combined$Y)

row.names(df.PRAD.exp) <- genes
df.PRAD.exp.r <- subset(df.PRAD.exp, row.names(df.PRAD.exp) %in% gene.final)

##OK now something is going wrong when I try to do the line of code above
#Found the problem, the gene names in df.PRAD.exp have "d1." before each of their names, so of course there are no matches between gene.final and the gene names 
#Solving the problem was really easy actually, I first assign the row names of df.PRAD.exp from the "genes" character vector object, and then I can do the following line of code 

row.names(df.PRAD.exp.r) <- ensembl.final
row.names(df.PRAD.meth) <- probe.ids

#Finally converting the data frames for the methylation and expression data into matrices to feed into createMAE

mat.PRAD.exp <- as.matrix.data.frame(df.PRAD.exp.r) 
mat.PRAD.meth <- as.matrix.data.frame(df.PRAD.meth)

#Creating a data frame to read into the colData arguement of createMAE

phenotype.data <- data.frame(row.names = Patient.ids,
                             samples = Patient.ids
)


#Now to actually create the MAE

mae.PRAD <- createMAE(exp = mat.PRAD.exp,
                      met = mat.PRAD.meth,
                      met.platform = '450K',
                      genome = 'hg19',
                      save = TRUE,
                      save.filename = "PRADmae.rda",
                      TCGA = FALSE
)



samples <- c('TCGA-G9-6371', 'TCGA-HC-7233', 'TCGA-KK-A6E7',
             'TCGA-XK-AAIV', 'TCGA-FC-A8O0', 'TCGA-VN-A88R')




expresults <- getResults(query.exp.PRAD)

expresults <- expresults[-c(30, 141, 269, 493, 139, 188, 249, 331), ]


with(sig.diff.meth.PRAD, plot(Primary.solid.Tumor_Minus_Solid.Tissue.Normal, -log10(adjust.p), pch=20, 
                              main='Probes hypomethylated in Primary solid tumor vs Solid tissue normal',
                              xlim = c(-0.6, 0)))



scatter.plot(data = PRAD.mae,
             byProbe = list(probe = c('cg00012148'), numFlankingGenes = 20),
             category = 'definition',
             lm = TRUE,
             save = FALSE)


#Temporarily moving the PRAD analysis here to clear up the main script 



sig.diff.meth.PRAD <- get.diff.meth(data = PRAD.mae,
                                    group.col = 'definition',
                                    group1 = "Primary solid Tumor",
                                    group2 = 'Solid Tissue Normal',
                                    minSubgroupFrac = 0.2,
                                    sig.dif = 0.3,
                                    diff.dir = 'hypo'
)

head(sig.diff.meth.PRAD)

#Wow everything above actually worked which blows my mind, now I can begin looking at the analysis
#Starting with a volcano plot of significantly hypomethylated probes

PRADdiffmethprobe <- read.csv('getMethdiff.hypo.probes.csv', header = TRUE)
PRADdiffmethprobesig <- read.csv('getMethdiff.hypo.probes.significant.csv', header = TRUE)


with(PRADdiffmethprobe, plot(Primary.solid.Tumor_Minus_Solid.Tissue.Normal, -log10(adjust.p), pch=20, 
                             main='Probes hypomethylated in PRAD Primary solid tumor vs Solid tissue normal',
                             xlab = 'DNA Methylation difference (Beta values)',
                             ylab = '-Log10 (FDR corrected P-values) [one tailed test]',
                             xlim = c(-0.6, 0)))
with(subset(PRADdiffmethprobe, Primary.solid.Tumor_Minus_Solid.Tissue.Normal< -0.3), points(Primary.solid.Tumor_Minus_Solid.Tissue.Normal, -log10(adjust.p), pch=20, col="red"))


#Ok the above volcano plot is pretty nice, saving it as 241019_Hypomethprobes

#Now looking for putative probe-gene pairs 

nearGenes <- GetNearGenes(data = PRAD.mae,
                          probes = sig.diff.meth.PRAD$probe
)

hypo.pair <- get.pair(data = PRAD.mae,
                      group.col = 'definition',
                      group1 = 'Primary solid Tumor',
                      group2 = 'Solid Tissue Normal',
                      nearGenes = nearGenes,
                      mode = 'unsupervised',
                      diff.dir = 'hypo')

#saving the file generated above, even though it's small it took ~3 hours to generate, so best not run the code above again

save(hypo.pair, file = "PRAD.hypo.pair")

#Just generating an example scatter plot from the above data, interesting but not sure what to make of it. 

scatter.plot(data = PRAD.mae, byPair = list(probe = c("cg25918833"), 
                                            gene = c("ENSG00000167332")), 
             category = 'definition', lm_line = TRUE)


#######pretty sure this is just a repeat of the workflow done in TCGA_download_and_MAE_workflow but thought that I would save it any way just to be safe 


#####################DONT RUN AGAIN####################


query.exp.PRAD <- GDCquery(project = 'TCGA-PRAD',
                           data.category = 'Gene expression',
                           data.type = 'Gene expression quantification',
                           platform = 'Illumina HiSeq',
                           file.type = 'normalized_results',
                           legacy = TRUE
)

getResults(query.exp.PRAD)

query.exp.PRAD[[1]][[1]] <- query.exp.PRAD[[1]][[1]][!duplicated(query.exp.PRAD[[1]][[1]]$cases),]

GDCdownload(query.exp.PRAD)
PRADexp <- GDCprepare(query.exp.PRAD)

rownames(PRADexp) <- values(PRADexp)$ensembl_gene_id

#now downloading the methylation data

query.meth.PRAD <- GDCquery(project = 'TCGA-PRAD',
                            data.category = 'DNA methylation',
                            platform = 'Illumina Human Methylation 450',
                            legacy = TRUE
)

GDCdownload(query.meth.PRAD, method = 'api')
PRADmeth <- GDCprepare(query.meth.PRAD)

#save(PRADmeth, file = "TCGA-PRAD-meth.RData")
load(file = "~/Desktop/KBH_PMD_Project_master/MAEs/PRADmeth.RData")


#Now on to creating the MAE for ELMER

PRAD.mae <- createMAE(exp = PRADexp,
                      met = PRADmeth,
                      TCGA = TRUE,
                      genome = 'hg19',
                      filter.probes = PMD.probe.promoters)

save(PRAD.mae, file = "TCGA-PRAD-MAE")

#############################START RUNNING AGAIN HERE#######################################3


load(file = "~/Desktop/KBH_PMD_Project_master/MAEs/TCGA-PRAD-MAE")

#once the MAE is created, I can start the actual analysis - in this case I am assessing differential (hypo-)methylation between normal and tumor tissue 



######################################################TCGA-COAD

#Now to try downloading another cancer data set, at least for the meantime - I'll go with colon adenocarcinoma (COAD)

query.exp.COAD <- GDCquery(project = 'TCGA-COAD',
                           data.category = 'Gene expression',
                           data.type = 'Gene expression quantification',
                           platform = 'Illumina HiSeq',
                           file.type = 'normalized_results',
                           legacy = TRUE
)



query.exp.COAD[[1]][[1]] <- query.exp.COAD[[1]][[1]][!duplicated(query.exp.COAD[[1]][[1]]$cases),]

GDCdownload(query.exp.COAD)
COADexp <- GDCprepare(query.exp.COAD)

rownames(COADexp) <- values(COADexp)$ensembl_gene_id

#now downloading the methylation data

query.meth.COAD <- GDCquery(project = 'TCGA-COAD',
                            data.category = 'DNA methylation',
                            platform = 'Illumina Human Methylation 450',
                            legacy = TRUE
)

GDCdownload(query.meth.COAD, method = 'api')
COADmeth <- GDCprepare(query.meth.COAD)

#save(COADmeth, file = "TCGA-COAD-meth.RData")

load(file = "~/Desktop/KBH_PMD_Project_master/MAEs/TCGA-COAD-meth.RData")

#Now on to creating the MAE for ELMER

COAD.mae <- createMAE(exp = COADexp,
                      met = COADmeth,
                      TCGA = TRUE,
                      genome = 'hg19',
                      filter.probes = PMD.probe.promoters)

save(COAD.mae, file = "TCGA-COAD-MAE")

#I am going to repeat creating MAEs for different cancer type in the following lines 

#Starting with Pancreatic cancer (PAAD)

query.exp.PAAD <- GDCquery(project = 'TCGA-PAAD',
                           data.category = 'Gene expression',
                           data.type = 'Gene expression quantification',
                           platform = 'Illumina HiSeq',
                           file.type = 'normalized_results',
                           legacy = TRUE
)



query.exp.PAAD[[1]][[1]] <- query.exp.PAAD[[1]][[1]][!duplicated(query.exp.PAAD[[1]][[1]]$cases),]

GDCdownload(query.exp.PAAD)
PAADexp <- GDCprepare(query.exp.PAAD)

rownames(PAADexp) <- values(PAADexp)$ensembl_gene_id

#now downloading the methylation data

query.meth.PAAD <- GDCquery(project = 'TCGA-PAAD',
                            data.category = 'DNA methylation',
                            platform = 'Illumina Human Methylation 450',
                            legacy = TRUE
)

GDCdownload(query.meth.PAAD, method = 'api')
PAADmeth <- GDCprepare(query.meth.PAAD)

save(PAADmeth, file = "TCGA-PAAD-meth.RData")
#Load() will go here upon subsequent reloadings 

#Now on to creating the MAE for ELMER

PAAD.mae <- createMAE(exp = PAADexp,
                      met = PAADmeth,
                      TCGA = TRUE,
                      genome = 'hg19',
                      filter.probes = PMD.probe.promoters)

save(PAAD.mae, file = "TCGA-PAAD-MAE")


#Now on to stomach adenocarcinoma (STAD)


query.exp.STAD <- GDCquery(project = 'TCGA-STAD',
                           data.category = 'Gene expression',
                           data.type = 'Gene expression quantification',
                           platform = 'Illumina HiSeq',
                           file.type = 'normalized_results',
                           legacy = TRUE
)



query.exp.STAD[[1]][[1]] <- query.exp.STAD[[1]][[1]][!duplicated(query.exp.STAD[[1]][[1]]$cases),]

GDCdownload(query.exp.STAD)
STADexp <- GDCprepare(query.exp.STAD)

rownames(STADexp) <- values(STADexp)$ensembl_gene_id

#now downloading the methylation data

query.meth.STAD <- GDCquery(project = 'TCGA-STAD',
                            data.category = 'DNA methylation',
                            platform = 'Illumina Human Methylation 450',
                            legacy = TRUE
)

GDCdownload(query.meth.STAD, method = 'api')
STADmeth <- GDCprepare(query.meth.STAD)

save(STADmeth, file = "TCGA-STAD-meth.RData")
#Load() will go here upon subsequent reloadings 

#Now on to creating the MAE for ELMER

STAD.mae <- createMAE(exp = STADexp,
                      met = STADmeth,
                      TCGA = TRUE,
                      genome = 'hg19',
                      filter.probes = PMD.probe.promoters)

save(STAD.mae, file = "TCGA-STAD-MAE")



probe.gene.pairs <- read.csv(file = "~/Desktop/KBH_PMD_Project_master/ELMER_Results/PRAD/getPair..all.pairs.statistic.csv")

probe.gene.pairs.sub <- probe.gene.pairs[1:100,]

experiments(PRAD.mae)

PRADmeth.sub <- PRAD.mae[probe.gene.pairs$Probe, ,1]
PRADexp.sub <- PRAD.mae[probe.gene.pairs$GeneID, ,2]

merged.PRAD.mae <- mergeReplicates(intersectColumns(PRAD.mae))

merged.PRADmeth.df <- longFormat(merged.PRAD.mae[, ,1])

scatter.plot(data = PRAD.mae,
             byProbe = list(probe = c("cg25918833"), numFlankingGenes = 20),
             category = 'definition',
             lm = TRUE,
             save = FALSE)


scatter.plot(data = PRAD.nonPMD.mae,
             byProbe = list(probe = c("cg05480730"), numFlankingGenes = 20),
             category = 'definition',
             lm = TRUE,
             save = FALSE)


nearPMDGenes <- GetNearGenes(data = PRAD.mae,
                             probes = df.PMD.probe.promoters$Name
)

#Just seeing if what I want to do will actually work below

probe.gene.pairs <- read.csv(file = "~/Desktop/KBH_PMD_Project_master/ELMER_Results/PRAD/getPair..all.pairs.statistic.csv")

probe.gene.pairs.sub <- probe.gene.pairs[1:100,]

experiments(PRAD.mae)

PRADmeth.sub <- PRAD.mae[probe.gene.pairs$Probe, ,1]
PRADexp.sub <- PRAD.mae[probe.gene.pairs$GeneID, ,2]

merged.PRAD.mae <- mergeReplicates(intersectColumns(PRAD.mae))

merged.PRADmeth.df <- longFormat(PRADmeth.sub[, ,1])

sig.diff.meth.non.PMD <- get.diff.meth(data = PRAD.nonPMD.mae,
                                       group.col = 'definition',
                                       group1 = 'Primary solid Tumor',
                                       group2 = 'Solid Tissue Normal',
                                       diff.dir = 'hypo')


nearGenes <- GetNearGenes(data = PRAD.nonPMD.mae,
                          probes = sig.diff.meth.non.PMD$probe
)


scatter.plot(data = PRAD.nonPMD.mae,
             byProbe = list(probe = c("cg05115308"), numFlankingGenes = 100),
             category = 'definition',
             lm = TRUE,
             save = FALSE)


scatter.plot(data = PRAD.mae,
             byProbe = list(probe = c("cg02328010"), numFlankingGenes = 100),
             category = 'definition',
             lm = TRUE,
             save = FALSE)

#Code below doesn't work, needs a MEE object (ELMER v1) not a MAE object (ELMER v2)

PMD.promoter.meth <- promoterMeth(PRAD.mae,
                                  sig.pvalue = 0.01,
                                  minSubgroupFrac = 0.4,
                                  upstream = 200,
                                  downstream = 2000,
                                  save = FALSE)










probe.gene.PMD.pairs <- read.csv(file = "~/Desktop/KBH_PMD_Project_master/ELMER_Results/PRAD/getPair..all.pairs.statistic.csv")


sig.diff.meth.nonPMD.PRAD <- get.diff.meth(data = PRAD.nonPMD.mae,
                                           group.col = 'definition',
                                           group1 = "Primary solid Tumor",
                                           group2 = 'Solid Tissue Normal',
                                           minSubgroupFrac = 0.2,
                                           sig.dif = 0.3,
                                           diff.dir = 'hypo')


near.non.PMD.Genes <- GetNearGenes(data = PRAD.nonPMD.mae,
                                   probes = sig.diff.meth.nonPMD.PRAD$probe)

hypo.pair.nonPMD <- get.pair(data = PRAD.nonPMD.mae,
                             group.col = 'definition',
                             group1 = 'Primary solid Tumor',
                             group2 = 'Solid Tissue Normal',
                             nearGenes = near.non.PMD.Genes,
                             mode = 'unsupervised',
                             diff.dir = 'hypo')

save(hypo.pair.nonPMD, file = "PRAD.nonPMD.hypo.pair")



PMDpromoter.gene.pairs <- GetNearGenes(data = PRAD.mae,
                                       probes = df.PMD.probe.promoters$Name,
                                       numFlankingGenes = 2)

nonPMDpromoter.gene.pairs <- GetNearGenes(data = PRAD.nonPMD.mae,
                                          probes = df.not.PMD.promoter.probes$Name,
                                          numFlankingGenes = 2)



load(file = "~/Desktop/KBH_PMD_Project_master/MAEs/Multi_Assay_Experiments/TCGA-PRAD-MAE")
load(file = "~/Desktop/KBH_PMD_Project_master/MAEs/TCGA-COAD-MAE")
load(file = "~/Desktop/KBH_PMD_Project_master/MAEs/TCGA-PAAD-MAE")
load(file = "~/Desktop/KBH_PMD_Project_master/MAEs/TCGA-STAD-MAE")
load(file = "~/Desktop/KBH_PMD_Project_master/MAEs/TCGA-SKCM-MAE")
load(file = "~/Desktop/KBH_PMD_Project_master/MAEs/TCGA-BRCA-MAE") 

load(file = "~/Desktop/KBH_PMD_Project_master/MAEs/Multi_Assay_Experiments/TCGA-PRAD-MAE")
load(file = "~/Desktop/KBH_PMD_Project_master/MAEs/Multi_Assay_Experiments/TCGA-PRAD-nonPMD-MAE")


PRADexp.data <- assays(PRADexp)$normalized_count
PRADmeth.data <- assay(PRADmeth)

scatter.plot(data = mae.PRAD,
             byPair = list(probe = c('cg20846447'), gene = c('ENSG00000198797')),
             category = 'definition', 
             lm_line = TRUE)



PMD.gene.probe.pairs <- getNearestGene(PMD.probe.promoters)


experiments(PRAD.mae)

XX.exp <- PRAD.mae[["Gene expression"]]
XX.mat.exp <- assay(XX.exp)




XX.PMD.ensembl <- df.PMD.pairs.r$ENSEMBL
XX.PMD.probe <- df.PMD.pairs.r$Name

XX.mat.exp.r <- subset(XX.mat.exp, rownames(XX.mat.exp) %in% XX.PMD.ensembl)

XX.meth <- PRAD.mae[['DNA methylation']]
XX.mat.meth <- assay(XX.meth)

XX.mat.meth.r <- subset(XX.mat.meth, rownames(XX.mat.meth) %in% XX.PMD.probe)


pb <- txtProgressBar(min = 0, max = nrow(df.nonPMD.pairs.r), style = 3)

#Now to iterate the correlation test over every probe gene pair within PRAD.mae, sstarting with the probe gene pairs from the PMD:

for (i in seq_len(nrow(df.nonPMD.pairs.r))){
  test <- cor.test(assay(PRAD.mae[df.PMD.pairs.r$Name[i], , drop = T]), assay(PRAD.mae[df.PMD.pairs.r$ENSEMBL[i], , drop = T]))
  df.PMD.pairs.r$p.value[i] <- test$p.value
  df.PMD.pairs.r$r[i] <- test$estimate
  setTxtProgressBar(pb, i)
}
close(pb)


##The good, clean, original one

for (i in seq_len(nrow(df.PMD.pairs.r))){
  test <- cor.test(assay(PRAD.mae[df.PMD.pairs.r$Name[i], , drop = T]), assay(PRAD.mae[df.PMD.pairs.r$ENSEMBL[i], , drop = T]))
  df.PMD.pairs.r$p.value[i] <- test$p.value
  df.PMD.pairs.r$r[i] <- test$estimate
  setTxtProgressBar(pb, i)
}
close(pb)


PRADexp <- assay(PRAD.mae[df.PMD.pairs.r$ENSEMBL, , 2])
PRADmeth <- assay(PRAD.mae[df.PMD.pairs.r$Name, , 1])

PRADmeth.r <- na.omit(PRADmeth)

PRADmeth.r.rownames <- rownames(PRADmeth.r)

df.PMD.pairs.r2 <- subset(df.PMD.pairs.r, df.PMD.pairs.r$Name %in% PRADmeth.r.rownames)


for (i in seq_len(nrow(df.PMD.pairs.r))) {
  if (assay(PRAD.mae[df.PMD.pairs.r$Name[i], , drop = T]) > 0 & assay(PRAD.mae[df.PMD.pairs.r$ENSEMBL[i], , drop = T]) > 0 ) {
  test <- cor.test(assay(PRAD.mae[df.PMD.pairs.r$Name[i], , drop = T]), assay(PRAD.mae[df.PMD.pairs.r$ENSEMBL[i], , drop = T]))
  df.PMD.pairs.r$p.value[i] <- test$p.value
  df.PMD.pairs.r$r[i] <- test$estimate
  setTxtProgressBar(pb, i)
  }
  else print('NA')
}
close(pb)


PRADmeth.r <- subset(PRADmeth, rownames(PRADmeth) %in% df.PMD.pairs.r$Name)



experiments(PRAD.mae)


PRADmeth.r <- na.omit(PRADmeth)
PRADexp.r <- na.omit(PRADexp)

dim(PRADmeth)
dim(PRADmeth.r)

dim(PRADexp)
dim(PRADexp.r)


for (i in seq_len(nrow(df.PMD.pairs.r))){
  test <- cor.test(PRADmeth.r[,df.PMD.pairs.r$Name[i] ], PRADexp.r[, df.PMD.pairs.r$ENSEMBL[i]])
  df.PMD.pairs.r$p.value[i] <- test$p.value
  df.PMD.pairs.r$r[i] <- test$estimate
  setTxtProgressBar(pb, i)
}
close(pb)




X.PAIRS <- subset(df.PMD.pairs.r, df.PMD.pairs.r$Name %in% rownames(PRADmeth.r))

rownames(X.PAIRS) <- seq(length = nrow(X.PAIRS))


for (i in seq_len(nrow(X.PAIRS))) {
  test <- cor.test(assay(PRAD.mae[X.PAIRS$Name[i], , drop = T]), assay(PRAD.mae[X.PAIRS$ENSEMBL[i], , drop = T]))
  X.PAIRS$p.value[i] <- test$p.value
  X.PAIRS$r[i] <- test$estimate
  setTxtProgressBar(pb, i)
}

close(pb)




X.meth.test <- PRADmeth['cg00091285', ]

X.meth.test


ENSG00000284247


X.exp.test <- PRADexp['ENSG00000284247', ]

X.PAIRS <- subset(X.PAIRS, X.PAIRS$ENSEMBL %in% rownames(PRADexp))


for (i in seq_len(nrow(X.PAIRS))) {
  test <- cor.test(assay(PRAD.mae[X.PAIRS$Name[i], , drop = T]), assay(PRAD.mae[X.PAIRS$ENSEMBL[i], , drop = T]))
  X.PAIRS$p.value[i] <- test$p.value
  X.PAIRS$r[i] <- test$estimate
  setTxtProgressBar(pb, i)
}

close(pb)

mcapply



rownames(assay(PRAD.mae[ , ,2]))


rownames(df.nonPMD.pairs.r) <- seq(length = nrow(df.nonPMD.pairs.r))



?mcapply

#The code below is to even add a progress bar to the function so that I know how long it will take

pb <- txtProgressBar(min = 0, max = nrow(df.nonPMD.pairs.r), style = 3)

setTxtProgressBar(pb, i)



cl <- makeCluster(8)
registerDoParallel(cl)

#Now to iterate the correlation test over every probe gene pair within PRAD.mae, starting with the probe gene pairs from the PMD:

foreach(i=1:5, .errorhandling = 'pass') %dopar% 
{
  test <- cor.test(assay(PRAD.mae[df.PMD.pairs.FINAL$Name[i], , drop = T]), assay(PRAD.mae[df.PMD.pairs.r$ENSEMBL[i], , drop = T]))
  df.PMD.pairs.r$p.value[i] <- test$p.value
  df.PMD.pairs.r$r[i] <- test$estimate
}


df.PMD.pairs.FINAL <- subset(df.PMD.pairs.r, df.PMD.pairs.r$ENSEMBL %in% rownames(assay(PRAD.mae[ , ,2])))
df.PMD.pairs.FINAL <- subset(df.PMD.pairs.FINAL, df.PMD.pairs.FINAL$Name %in% rownames(assay(PRAD.mae[ , ,1])))
rownames(df.PMD.pairs.FINAL) <- seq(length = nrow(df.PMD.pairs.FINAL))





pb <- txtProgressBar(min = 0, max = nrow(df.PMD.pairs.FINAL), style = 3)

for (i in seq_len(nrow(df.PMD.pairs.FINAL))) {
  test <- cor.test(assay(PRAD.mae[df.PMD.pairs.FINAL$Name[i], , drop = T]), assay(PRAD.mae[df.PMD.pairs.FINAL$ENSEMBL[i], , drop = T]))
  df.PMD.pairs.FINAL$p.value[i] <- test$p.value
  df.PMD.pairs.FINAL$r[i] <- test$estimate
  setTxtProgressBar(pb, i)
}

close(pb)

#Ok for some reason it seems to have worked on the "X.PAIRS" data set, I will have to go back through and see what the difference with that one was
#For now to at least see if I can see what the figure looks like from the X.PAIRS data
df.PMD.pairs.FINAL <- X.PAIRS

save(df.PMD.pairs.FINAL, file = "PRAD_PMD_Exp-meth_corr_test")





#First to set up and register the cores
cl <- makeCluster(6)
registerDoParallel(cl)

#Now to try running the foreach loop in parallel

foreach(i=seq_len(nrow(df.nonPMD.pairs.FINAL)), .errorhandling = 'pass') %dopar% 
{
  test2 <- cor.test(assay(PRAD.mae[df.nonPMD.pairs.FINAL$Name[i], , drop = T]), assay(PRAD.mae[df.nonPMD.pairs.FINAL$ENSEMBL[i], , drop = T]))
  df.nonPMD.pairs.FINAL$p.value[i] <- test2$p.value
  df.nonPMD.pairs.FINAL$r[i] <- test2$estimate
}

######################################################################################33

test.data.frame <- df.nonPMD.pairs.FINAL[1:10, ]
test.data.frame$p.value <- NA
test.data.frame$r <- NA





c.test <- function(MAE, df, i) {
  test2 <- cor.test(assay(MAE[df$Name[i], , drop = T]), assay(MAE[df$ENSEMBL[i], , drop = T]))
  return(test2)
}




corr <- sapply(1:500, function(i){
              result = cor.test(assay(PRAD.mae[df.nonPMD.pairs.FINAL$Name[i], , drop = T]), assay(PRAD.mae[df.nonPMD.pairs.FINAL$ENSEMBL[i], , drop = T]))
              return(c(result$estimate, result$p.value))
              })


corr2 <- mclapply(1:500, function(i){
  result = cor.test(assay(PRAD.mae[df.nonPMD.pairs.FINAL$Name[i], , drop = T]), assay(PRAD.mae[df.nonPMD.pairs.FINAL$ENSEMBL[i], , drop = T]))
  return(c(result$estimate, result$p.value))
}, mc.preschedule = TRUE)

###########################################################


PRADexp2 <- assay(PRAD.mae[df.nonPMD.pairs.r$ENSEMBL, , 2])
PRADmeth2 <- assay(PRAD.mae[df.nonPMD.pairs.r$Name, , 1])

PRADmeth.r2 <- na.omit(PRADmeth2)

df.nonPMD.pairs.FINAL <- subset(df.nonPMD.pairs.r, df.nonPMD.pairs.r$Name %in% rownames(PRADmeth.r2))
rownames(df.nonPMD.pairs.FINAL) <- seq(length = nrow(df.nonPMD.pairs.FINAL))
df.nonPMD.pairs.FINAL <- subset(df.nonPMD.pairs.FINAL, df.nonPMD.pairs.FINAL$ENSEMBL %in% rownames(PRADexp2))
rownames(df.nonPMD.pairs.FINAL) <- seq(length = nrow(df.nonPMD.pairs.FINAL))


#Now to actually run the correlation tests, this time I will try to run things in parallel using an apply function because for the nonPMD data I have roughly 6x the amount of samples to go through :/

#For writing a function, I would like to enter the MAE, specify which list Exp and Meth are in (i.e. 1 and 2)
#And whether i want to run the correlations on PMD or non-PMD regions. The results of the correlation test should
#be output into a data frame containing the gene-probe pairs, correlation coefficients, and p values. 

c.test.fun <- function(MAE, Exp, Meth, PMD=TRUE, Output){
  
  if PMD==TRUE return (Meth = assay(MAE[df.PMD.pairs.FINAL$Name[i], , drop = T])) & Exp = assay(MAE[df.PMD.pairs.FINAL$ENSEMBL[i], , drop = T])))
if PMD==FALSE return (Meth = assay(MAE[df.nonPMD.pairs.FINAL$Name[i], , drop = T]) & Exp = assay(MAE[df.nonPMD.pairs.FINAL$ENSEMBL[i], , drop = T])))
test <- cor.test(Exp, Meth)
Output$p.value[i] <- test$p.value
Output$r[i] <- test$estimate
}

#Putting the above function on hold for right now, trying to get mcapply to work 


#creating a function for the correlation test


#########################################################

c.test <- function(MAE, probeList, exp, meth, cores) {
  
  #cl <- parallel::makeCluster(cores)
  #doParallel::registerDoParallel(cl)
  
  
  
  df <- probeList
  
  foreach (i = probeList) %do% {
    browser()
    i
    test <- cor.test(assay(MAE[i$Name, , meth]), assay(MAE[i$ENSEMBL, , exp]))
    df$p.value[df$Name == i$Name] <- test$p.value
    df$r[df$Name == i$Name] <- test$estimate
  }
  
  return(df)
}




df[1,]$r <- 5
df$Name
probeList$Name


probeList <- probeList[1,]

probeList <- df.nonPMD.pairs.FINAL
df <- XX.test.data[1,]
df$r[df$Name == probeList$Name] <- probeList$r
df$r <- 5

c.test(PRAD.mae, XX.test.data, 2,1,1)

XX.test.data <- df.nonPMD.pairs.FINAL[c(1:5), ]
XX.test.data$p.value <- NA
XX.test.data$r <- NA

close(pb)

pb <- txtProgressBar(min = 0, max = nrow(df.nonPMD.pairs.FINAL), style = 3)
index = 0


pmap(XX.test.data,c.test)
# lapply(df.nonPMD.pairs.FINAL,c.test)

c.test <- function(thing) {
  thing
  # thing[,1]
  # thing[,2]
  
  # test <- cor.test(assay(PRAD.mae[i[,'Name'], , drop = T]), assay(PRAD.mae[i[,'ENSEMBL'], , drop = T]))
  # i[,"p.value"] <- test$p.value
  # i[,"r"] <- test$estimate
  # index <- index+1
  # setTxtProgressBar(pb, index)
}

close(pb)




#My last ditch backup idea for getting the data 


# pb <- txtProgressBar(min = 0, max = nrow(df.nonPMD.pairs.FINAL), style = 3)
# 
# for (i in seq_len(nrow(df.nonPMD.pairs.FINAL))) {
#   test <- cor.test(assay(PRAD.mae[df.nonPMD.pairs.FINAL$Name[i], , drop = T]), assay(PRAD.mae[df.nonPMD.pairs.FINAL$ENSEMBL[i], , drop = T]))
#   df.nonPMD.pairs.FINAL$p.value[i] <- test$p.value
#   df.nonPMD.pairs.FINAL$r[i] <- test$estimate
#   setTxtProgressBar(pb, i)
# }


subset(df.PC3.annots, df.PC3.annots$start == 1)

df.PC3.chrom.annots <- subset(df.chrom.annots, df.chrom.annots$Cell_Line == 'PC3') 
df.PrEC.chrom.annots <- subset(df.chrom.annots, df.chrom.annots$Cell_Line == 'PrEC') 

df.PC3.chrom.annots <- subset(df.chrom.annots, df.chrom.annots$Cell_Line == 'PC3') 
df.PrEC.chrom.annots <- subset(df.chrom.annots, df.chrom.annots$Cell_Line == 'PrEC') 


ggplot(df.PC3.chrom.annots, aes(Annot.Type, percent, fill=Region)) + 
  geom_bar(position = 'dodge', stat = "identity", width=0.5) +
  coord_flip()

ggplot(df.PC3.chrom.annots, aes(Annot.Type, percent, fill=Region)) + 
  geom_bar(position = 'dodge', stat = "identity", width=0.5) +
  coord_flip()


boxplot(nonPMD.corr.test$r, df.PMD.pairs.FINAL$r,
        main = 'Distribution of Correlation Coefficients',
        ylab = 'Correlation Coefficient (r)',
        xlab = 'Group',
        names = c('Non-PMD regions', 'PMD regions'))



######Time for Motif enrichment##########

df.BRCA.PMD.corr.test.sig <- subset(df.BRCA.PMD.corr.test, df.PRAD.PMD.corr.test$p.value < (0.05))
df.BRCA.nonPMD.corr.test.sig <- subset(df.PRAD.nonPMD.corr.test, df.PRAD.nonPMD.corr.test$p.value < (0.05))


BRCA.PMD.EM.sig.pairs <- get.enriched.motif(data = BRCA.mae,
                                            probes = df.PRAD.PMD.co$Name, 
                                            label = "hypo",
                                            min.incidence = 10,
                                            lower.OR = 1.1)




covplot(region.annots.gr, weightCol = 'r', 
        xlab = 'Chromosome Size (bp)',
        ylab = 'Correlation Coefficient (r)',
        title = 'R Values of Probe-Gene Pairs',
        chrs = NULL,
        lower = 0)



#Below I am trying a different way of organizing the big data frame, based around the probe rather than the regions

df.regions <- df.regions[ , c(1, 2, 3, 8)]
names(df.regions)[4] <- 'region'

promoter.probes <- get.feature.probe(feature = NULL, genome = 'hg19', met.platform = "450K", promoter = TRUE)
region.promoter.probes <- subsetByOverlaps(promoter.probes, regions)
df.region.prom.probes <- as.data.frame(region.promoter.probes)

df.region.prom.probes <- df.region.prom.probes[, c(1, 2, 3)]
region.promoter.probes <- splicejam::annotateGRfromGR(GR1 = region.promoter.probes, GR2 = regions)

df.region.prom.probes <- as.data.frame(region.promoter.probes)
df.region.prom.probes <- df.region.prom.probes[, c(1, 2, 3, 53)]
names(df.region.prom.probes)[4] <- 'region'
df.region.prom.probes$probe <- rownames(df.region.prom.probes)

regions.gr <- makeGRangesFromDataFrame(df.region.prom.probes, 
                                       keep.extra.columns = TRUE, 
                                       seqnames.field = "seqnames",
                                       start.field = "start", 
                                       end.field = "end")

Region.new.annotated <- annotate_regions(regions = regions.gr,
                                         annotations = annotations,
                                         ignore.strand = TRUE,
                                         quiet = FALSE)

probes.regions.annots <- data.frame(Region.new.annotated)
probes.regions.annots <- probes.regions.annots[, c(1, 2, 3, 6, 7, 15, 16, 17)]

str(probes.regions.annots)

probes.regions.annots$region <- as.factor(probes.regions.annots$region)
probes.regions.annots$annot.type <- as.factor(probes.regions.annots$annot.type)


summary(probes.regions.annots$annot.type)


#See if highly correlated genes are enriched in certain regions of PMDs - or associated with certain chromatin marks 




col <- colorRampPalette(brewer.pal(8, "RdYlGn"))(25)

heatmap(PRADmeth.r, 
        scale = 'column',
        Colv = NA,
        Rowv = NA,
        col = col,
        xlab = 'Patient',
        ylab = 'Probe')


meth.samples <- PRAD.mae@sampleMap@listData[['colname']]
meth.samples <- meth.samples[1:537]

#Trying to make the heatmap on a different subset of the data

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



#Now to try looking at the ChromHMM annotations


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

#getting and formatting the PrEC data:

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


#####Getting annotations for both PC3 and PrEC onto the regions#####

df.region.pairs.FINAL$PC3.annots <- NA
df.region.pairs.FINAL$PrEC.annots <- NA

X.PC3.annots <- splicejam::annotateGRfromGR(GR1 = region.annotated.gr, GR2 = PC3.annots.gr)
X.df.PC3 <- as.data.frame(X.PC3.annots)

X.df.PC3 <- subset(X.df.PC3, X.df.PC3$ID %in% rownames(PRADmeth))
X.df.PC3 <- subset(X.df.PC3, X.df.PC3$ENSEMBL %in% rownames(PRADexp))

X.PrEC.annots <- splicejam::annotateGRfromGR(GR1 = region.annotated.gr, GR2 = PrEC.annots.gr)
X.df.PrEC <- as.data.frame(X.PrEC.annots)

X.df.PrEC <- subset(X.df.PrEC, X.df.PrEC$ID %in% rownames(PRADmeth))
X.df.PrEC <- subset(X.df.PrEC, X.df.PrEC$ENSEMBL %in% rownames(PRADexp))

df.region.pairs.FINAL$PC3.annots <- X.df.PC3$PC3_chrom_annots.name
df.region.pairs.FINAL$PrEC.annots <- X.df.PrEC$annots.type

#Cleaning up things
rm(X.PC3.annots, X.df.PC3, X.df.PrEC, X.PrEC.annots)

#Making df.region.pairs.FINAL into a granges object and removing NA values
df.region.pairs.FINAL <- na.omit(df.region.pairs.FINAL)

region.annots.FINAL.gr <- makeGRangesFromDataFrame(df.region.pairs.FINAL, 
                                                   keep.extra.columns = TRUE, 
                                                   seqnames.field = "seqnames",
                                                   start.field = "start", 
                                                   end.field = "end")


#Have to get rid of sex chromosomes

chromosomes <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
                 'chr9', 'chr10', 'chr11', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22')

#Starting with PC3 annotations

region.PC3.annots <- splicejam::annotateGRfromGR(GR1 = PC3.annots.gr,
                                                 GR2 = regions)
region.PC3.annots@elementMetadata@listData$region <- as.factor(region.PC3.annots@elementMetadata@listData$region)

df.region.PC3.annots <- data.frame(region.PC3.annots)

df.region.PC3.annots <- subset(df.region.PC3.annots, df.region.PC3.annots$seqnames %in% chromosomes)

#Now for PrEC

region.PrEC.annots <- splicejam::annotateGRfromGR(GR1 = PrEC.annots.gr,
                                                  GR2 = regions)
region.PrEC.annots@elementMetadata@listData$region <- as.factor(region.PrEC.annots@elementMetadata@listData$region)

df.region.PrEC.annots <- data.frame(region.PrEC.annots)

df.region.PrEC.annots <- subset(df.region.PrEC.annots, df.region.PrEC.annots$seqnames %in% chromosomes)


#Trying with heatmap2 to see if I can plot more#

#top30k.methsd <- probe.meth.sd %>% top_n(30000, sd)
#top30k.methsd <- subset(top30k.methsd, top30k.methsd$probe %in% df.probes.regions$ID)
#top30k.methsd$region <- subset(df.probes.regions$region, df.probes.regions$ID %in% top30k.methsd$probe)

#PRADmeth.r2 <- subset(PRADmeth, row.names(PRADmeth) %in% top30k.methsd$probe)

#row.metadata <- data.frame(region = top30k.methsd$region,
#                           row.names = top30k.methsd$probe)

#heatmap.2(PRADmeth.r2,
#          labRow = row.metadata$region)


meth.sd <- rowMads(PRADmeth, na.rm = TRUE)
probe.meth.sd <- data.frame(probe = row.names(PRADmeth), sd = meth.sd)
top2k.methsd <- probe.meth.sd %>% top_n(2000, sd)
top2k.methsd <- subset(top2k.methsd, top2k.methsd$probe %in% df.probes.regions$ID)
top2k.methsd$region <- subset(df.probes.regions$region, df.probes.regions$ID %in% top2k.methsd$probe)

PRADmeth.r <- subset(PRADmeth, row.names(PRADmeth) %in% top2k.methsd$probe)

#Now making the heatmap 
row.metadata <- data.frame(region = top2k.methsd$region,
                           row.names = top2k.methsd$probe)

pheatmap(PRADmeth.r,
         cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE, show_colnames = FALSE,
         annotation_row = row.metadata)

#sample.id <- data.frame(ID = unique(PRAD.mae@sampleMap@listData[['primary']]))

PRAD.tp <- rbind.fill(PRAD.tp, sample.id)





#genes.hg19 <- rnb.annotation2data.frame(rnb.get.annotation(type = 'genes', assembly = 'hg19'))
#genes.hg19 <- genes.hg19[, c(1, 2, 3, 5, 7, 11)]
#genes.gr <- makeGRangesFromDataFrame(genes.hg19, 
#                                     keep.extra.columns = TRUE,
#                                     ignore.strand = TRUE,
#                                     seqnames.field = 'Chromosome',
#                                     start.field = 'Start',
#                                     end.field = 'End')


#Now to filter the genes and probes to get only the ones that are in the regions objects

probes.regions <- subsetByOverlaps(probes.gr, regions)
probes.regions <- splicejam::annotateGRfromGR(GR1 = probes.regions, 
                                              GR2 = regions)
df.probes.regions <- as.data.frame(probes.regions)

#Now getting just promoter probes

promoter.probes.gr <- subsetByOverlaps(probes.gr, promoters.gr)
region.prom.probes.gr <- subsetByOverlaps(promoter.probes.gr, regions)
region.prom.probes.anno <- splicejam::annotateGRfromGR(GR1 = region.prom.probes.gr,
                                                       GR2 = regions)

df.region.prom.probes <- as.data.frame(region.prom.probes.anno)
#From the above code I now have a list of the 450k probes that are within promoters and with their region (PMD vs. non-PMD) defined

df.region.prom.probes$region <- as.factor(df.region.prom.probes$region)

#Now doing the above filtering with the genes

#region.genes.gr <- subsetByOverlaps(genes.gr, regions)
#region.genes.anno <- splicejam::annotateGRfromGR(GR1 = region.genes.gr,
#                                                 GR2 = regions)

#df.region.genes <- as.data.frame(region.genes.anno)
#df.region.genes$region <- as.factor(df.region.genes$region)
#Removing genes that are in overlapping regions
#df.region.genes <- subset(df.region.genes, df.region.genes$region != 'non-PMD,PMD')

#Ok so now all the probes and genes are annotated to show which region they belong in

#Lastly restoring the GRanges objects after all the filtering

#df.region.genes <- df.region.genes[, c(1, 2, 3, 6, 7, 8, 9)]
df.region.prom.probes <- df.region.prom.probes[, c(1, 2, 3, 6, 7, 8, 9)]
row.names(df.region.prom.probes) <- c(1:168835)

region.prom.probes.gr <- makeGRangesFromDataFrame(df.region.prom.probes, 
                                                  keep.extra.columns = TRUE,
                                                  ignore.strand = TRUE,
                                                  seqnames.field = 'seqnames',
                                                  start.field = 'start',
                                                  end.field = 'end')

#region.genes.gr <- makeGRangesFromDataFrame(df.region.genes, 
#                                            keep.extra.columns = TRUE,
#                                            ignore.strand = TRUE,
#                                            seqnames.field = 'seqnames',
#                                            start.field = 'start',
#                                            end.field = 'end')



PRADmeth.r1 <- subset(PRADmeth, rownames(PRADmeth) %in% df.region.prom.probes$ID)

df.region.prom.probes <- subset(df.region.prom.probes, df.region.prom.probes$ID %in% rownames(PRADmeth.r1))

summary(df.probes.regions.genes$region)

#Ordering the PRADmeth object so that the probes are in the same order as df.region.prom.probes

PRADmeth.r1 <- PRADmeth.r1[order(match(rownames(PRADmeth.r1), df.region.prom.probes$ID)), , drop = FALSE]

####Finding probe-gene pairs#######

#using annotatr to annotate regions 

annots <- c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
            'hg19_genes_intronexonboundaries')

annotations <- build_annotations(genome = 'hg19', annotations = annots)

####

#XXX.region.anno.test <- splicejam::annotateGRfromGR(GR1 = region.prom.probes.gr,
#                                                    GR2 = annotations)
#XXX.df <- as.data.frame(XXX.region.anno.test)

###

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

df.region.pairs.FINAL <- subset(df.region.anno, df.region.anno$ID %in% rownames(PRADmeth.r1))

#Although the above annotations are "cleaner' they also lose the ENSEMBL ID, will need to figure out a way to add those in
df.region.pairs.FINAL <- subset(df.region.pairs.FINAL, df.region.pairs.FINAL$ENSEMBL %in% rownames(PRADexp))

#Changing the gene annotations so that they are in more standardized categories 

df.region.pairs.FINAL$annot.type <- as.factor(df.region.pairs.FINAL$annot.type)

rm(ensembl.ids, gene.map)

#Just trying a different method of getting probe-gene pairs# 

#nearGenes <- GetNearGenes(data = PRAD.mae, 
#                          probes = df.region.prom.probes$ID, 
#                          numFlankingGenes = 10)

#Hypo.pair <- get.pair(data = PRAD.mae,
#                      group.col = "definition",
#                      group1 =  "Primary solid Tumor",
#                      group2 = "Solid Tissue Normal",
#                      diff.dir = 'hypo',
#                      nearGenes = nearGenes,
#                      mode = "unsupervised",
#                      permu.size = 100000, # Please set to 100000 to get significant results
#                      raw.pvalue = 0.05,   
#                      Pe = 0.001, # Please set to 0.001 to get significant results
#                      dir.out = "PRAD_hypo_pairs",
#                      cores = 6)


#rm(nearGenes)

#Now I will re-run the correlation test

df.region.pairs.FINAL <- subset(df.region.pairs.FINAL, df.region.pairs.FINAL$ID %in% rownames(PRADmeth))
df.region.pairs.FINAL <- subset(df.region.pairs.FINAL, df.region.pairs.FINAL$ENSEMBL %in% rownames(PRADexp))

df.region.pairs.FINAL$r <- NA
df.region.pairs.FINAL$p.value <- NA

pb <- txtProgressBar(min = 0, max = nrow(df.region.pairs.FINAL), style = 3)

for (i in seq_len(nrow(df.region.pairs.FINAL))) {
  test <- cor.test(PRADmeth[df.region.pairs.FINAL$ID[i], ], PRADexp[df.region.pairs.FINAL$ENSEMBL[i], ])
  df.region.pairs.FINAL$p.value[i] <- test$p.value
  df.region.pairs.FINAL$r[i] <- test$estimate
  setTxtProgressBar(pb, i)
}

#save(df.region.pairs.FINAL, file = 'all_region_corr_test_annotatr')

####Loading up the Correlation Test data####

load(file = "~/Desktop/KBH_PMD_Project_master/data/output/Correlation_tests/all_region_corr_test_annotatr")

#Revelation, the above for loop runs WAY faster when I use the Matrices for PRADmeth and PRADexp, as opposed to using assay() on the MAE, for the future this is the way to do it

#Now I can look at the correlation in the different groups, and then do ChromHMM annotations, etc. 

ggplot(df.region.pairs.FINAL, aes(x=region, y=r, fill=region)) + 
  ylim(-1, 1) +
  geom_violin()

#!!!!!!!Something that I hadn't really thought about, the distributions will be skewed because there are multiple values for each probe-gene pair

####DOCK8####

#Curious about DOCK8, looks like there is some work suggesting it might be a T.S.G.
#https://www.ncbi.nlm.nih.gov/pubmed/16391785

DOCK8 <- subset(df.region.pairs.FINAL, df.region.pairs.FINAL$annot.symbol == 'DOCK8')

#Can plot probe-gene correlation using ELMER

scatter.plot(data = PRAD.mae,
             byPair = list(probe = c("cg19252956"), gene = c("ENSG00000107099")), 
             category = "definition", save = FALSE, lm_line = TRUE) 

scatter.plot(data = PRAD.mae,
             byPair = list(probe = c("cg13470032"), gene = c("ENSG00000107099")), 
             category = "definition", save = FALSE, lm_line = TRUE) 

scatter.plot(data = PRAD.mae,
             byPair = list(probe = c("cg19201719"), gene = c("ENSG00000107099")), 
             category = "definition", save = FALSE, lm_line = TRUE) 



rm(meth.sd, probe.meth.sd, top2k.methsd, PRADmeth.r, row.metadata, df.region.pairs.highly.corr, 
   sample.map, DOCK8)




#################################################################
####ChIP-seq####

extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

#The data imported below is from: https://www.encodeproject.org/experiments/ENCSR946MNG/

TF.CTCF.gr <- import.bed("~/Desktop/KBH_PMD_Project_master/data/raw/Chip-seq/ENCFF165JNG.bed.gz",
                         extraCols = extraCols_narrowPeak,
                         genome = 'hg19')




df.TF.CTCF <- as.data.frame(TF.CTCF.gr)
df.TF.CTCF$type <- 'CTCF'
TF.CTCF.gr <- makeGRangesFromDataFrame(df.TF.CTCF, 
                                       keep.extra.columns = TRUE,
                                       seqnames.field = 'seqnames',
                                       start.field = 'start',
                                       end.field = 'end')

XXX.test <- splicejam::annotateGRfromGR(pairs.FINAL.gr,
                                        TF.CTCF.gr)
XXX.df <- as.data.frame(XXX.test)


#Here I think I can use ELMER's tool for identifying regulatory TF's, I'll first subset the data so that I am only looking at probes within PMDs that have significant negative correlation
####ELMER motif enrichment####

PMD.enriched.motif <- get.enriched.motif(data = PRAD.mae,
                                         probes = df.PMDs$ID, 
                                         label = "hypo",
                                         min.incidence = 10,
                                         lower.OR = 1.1)

PMD.TFs <- get.TFs(data = PRAD.mae, 
                   group.col = "definition",
                   group1 =  "Primary solid Tumor",
                   group2 = "Solid Tissue Normal",
                   mode = "unsupervised",
                   enriched.motif = PMD.enriched.motif,
                   cores = 6, 
                   label = "hypo")

load("~/Desktop/KBH_PMD_Project_master/data/output/PMD.getTF.hypo.TFs.with.motif.pvalue.rda")
motif <- colnames(TF.meth.cor)[1]

#TF ranking plot: For a given enriched motif, all human TF are ranked by the statistical 
#− log10(P − value) assessing the anti-correlation level of candidate Master Regulator 
#TF expression with average DNA methylation level for sites with the given motif.
#As a result, the most anti-correlated TFs will be ranked in the first positions.
#By default, the top 3 most anti-correlated TFs, and all TF classified by TFClass database 
#in the same family and subfamily are highlighted with colors blue, red and orange, respectively. 

TF.rank.plot(motif.pvalue = TF.meth.cor, 
             motif = motif,
             save = FALSE) 

#################################################################

ARsum <- group_by(df.pairs.FINAL, AR.binding.sites, region) %>%
  dplyr::summarise(
    count = n())

CTCFsum <- group_by(df.pairs.FINAL, CTCF.binding.sites, region) %>%
  dplyr::summarise(
    count = n())

ETS1sum <- group_by(df.pairs.FINAL, ETS1.binding.sites, region) %>%
  dplyr::summarise(
    count = n())

SP1sum <- group_by(df.pairs.FINAL, SP1.binding.sites, region) %>%
  dplyr::summarise(
    count = n())

TP53sum <- group_by(df.pairs.FINAL, TP53.binding.sites, region) %>%
  dplyr::summarise(
    count = n())

FOXA1sum <- group_by(df.pairs.FINAL, FOXA1.binding.sites, region) %>%
  dplyr::summarise(
    count = n())

RB1sum <- group_by(df.pairs.FINAL, RB1.binding.sites, region) %>%
  dplyr::summarise(
    count = n())


FOXA1sum$total <- rep(c(46854, 14105, 34817), 2)
FOXA1sum$percent <- FOXA1sum$count / FOXA1sum$total * 100

ggplot(FOXA1sum, aes(fill=region, y=percent, x=FOXA1.binding.sites)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("FOXA1 binding sites in different regions") +
  theme(legend.position="top") +
  xlab("Annotation") +
  ylab('Percent total annotations') +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"))
#+
#  coord_flip()


genes.hg19 <- rnb.annotation2data.frame(rnb.get.annotation(type = 'genes', assembly = 'hg19'))
genes.hg19 <- genes.hg19[, c(1, 2, 3, 5, 7, 11)]
genes.gr <- makeGRangesFromDataFrame(genes.hg19, 
                                     keep.extra.columns = TRUE,
                                     ignore.strand = TRUE,
                                     seqnames.field = 'Chromosome',
                                     start.field = 'Start',
                                     end.field = 'End')

probe.genes <- splicejam::annotateGRfromGR(probes.gr,
                                           genes.gr)


df.probes.genes <- as.data.frame(probe.genes)

OR51E2 <- subset(df.probes.genes, df.probes.genes$symbol == 'OR51E2')


#f   

gonames = as.character(dplot$go_name)
lchar = sapply(gonames, nchar)
w0 = max(lchar)

fname = paste('pmd_mael_OR51E2_x_mlog10padj_y_GO_',gene,'.',sep="") 
ggsave(paste(path,fname,'pdf',sep=""),
       f, 
       width = (ymax*0.3+w0*0.1), 
       height=(0.2+0.25*length(dplot$go_name)))#units:inch
ggsave(paste(path,fname,'png',sep=""),
       f, 
       width = (ymax*0.3+w0*0.1), 
       height=(0.2+0.25*length(dplot$go_name)))


alpha_FDR_plot = 0.05
GENES = c('MAEL','OR51E2')

for (gene in GENES){
  
  dplot = PMD.gw[which(PMD.gw$hgnc_symbol==gene),]    
  #for plotting: order go terms (from all GO domains) according to mean_padj (and if equal then sort by mean_sim)
  dplot$go_name = factor(dplot$go_name, 
                         levels = dplot$go_name[with(dplot,order(-mean_padj,mean_sim))])
  
  font_sz=12
  ymax=max(dplot$mlog10padj + dplot$mlog10padj_err) + 0.3
  
  f = ggplot(dplot, aes(go_name,mlog10padj)) + 
    theme_light(base_size=font_sz) +
    geom_bar(stat='identity',fill='royalblue',width=0.8) +
    geom_errorbar(aes(ymin=mlog10padj-mlog10padj_error, ymax=mlog10padj+mlog10padj_error),
                  width=0.2, color='grey') +
    geom_hline(yintercept=-log10(alpha_FDR_plot), color='red', linetype='dashed') +
    xlab('') +
    ylab('-log10(p-adjust)') +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_text(size=font_sz*0.8)) +
    scale_y_continuous(breaks=c(0,1,2,3,4),
                       limits = c(0, ymax)) +
    coord_flip()
  
  
}

f





```{r}
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

df.region.pairs.FINAL <- subset(df.region.anno, df.region.anno$ID %in% rownames(PRADmeth.r1))

#Although the above annotations are "cleaner' they also lose the ENSEMBL ID, will need to figure out a way to add those in
df.region.pairs.FINAL <- subset(df.region.pairs.FINAL, df.region.pairs.FINAL$ENSEMBL %in% rownames(PRADexp))

#Changing the gene annotations so that they are in more standardized categories 

df.region.pairs.FINAL$annot.type <- as.factor(df.region.pairs.FINAL$annot.type)

```



df.PMDs <- subset(df.pairs.FINAL, df.pairs.FINAL$region == 'commonPMD' & df.pairs.FINAL$p.value < (0.05))


Correlation <- data.frame(correlation = c('very high negative', 'high negative', 'moderate negative', 'low negative', 'none', 'low positive', 'moderate positive', 'high positive', 'very high positive'),
                          lower = c(-1, -.9, -.7, -.5, -.3, .3, .5, .7, .9),
                          upper = c(-.9, -.7, -.5, -.3, .3, .5, .7, .9, 1))

df.PMDs <- data.frame(df.PMDs,
                      'correlation' = cut(df.PMDs$r,
                                          breaks = Correlation$lower,
                                          right = T,
                                          include.lowest = T))
unique(df.PMDs$correlation)
df.PMDs$correlation <- revalue(df.PMDs$correlation,
                               c('(-0.3,0.3]'='none', '(-0.5,-0.3]'='low negative',
                                 '(0.3,0.5]'='low positive', '(-0.7,-0.5]'='moderate negative',
                                 '(0.5,0.7]'='moderate positive', '(-0.9,-0.7]'='high negative',
                                 '[-1,-0.9]'='very high negative', '(0.7,0.9]'='high positive'))
summary(df.PMDs$correlation)



#######################################


#Basic analysis
**4: Basic Analysis**
  
  Now with regions loaded up and annotations added, I can start looking at what sorts of annotations are present in different regions. Below I am looking at the different CpG island relation annotations present in the different regions

```{r}
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

#And how many of each CGI relation are in each region

summary.regions <- group_by(df.probes.regions.genes, CGI.Relation, region) %>% dplyr::summarise(count = n())
summary.regions$percent <- summary.regions$count / rep(c(47109, 14339, 35007), 6) * 100

ggplot(summary.regions, aes(fill=region, y=percent, x=CGI.Relation)) + 
  geom_bar(position="dodge", stat="identity")


```

```{r}
rm(human.genes, gene.density, df.gene.density, PMD.gene.density, nonPMD.gene.density, summary.regions)
```

Based on the above figure, it seems like PMDs are enriched in 'open sea' areas, but this should be expected.

Some previous literature characterizing PMDs:
  


##################################


ARsum <- group_by(df.pairs.FINAL, AR.binding.sites, region) %>%
  dplyr::summarise(
    count = n())

ARsum$total <- rep(c(46854, 14105, 34817), 2)
ARsum$percent <- ARsum$count / ARsum$total * 100

ETS1sum <- group_by(df.pairs.FINAL, ETS1.binding.sites, region) %>%
  dplyr::summarise(
    count = n())

SP1sum <- group_by(df.pairs.FINAL, SP1.binding.sites, region) %>%
  dplyr::summarise(
    count = n())

TP53sum <- group_by(df.pairs.FINAL, TP53.binding.sites, region) %>%
  dplyr::summarise(
    count = n())

FOXA1sum <- group_by(df.pairs.FINAL, FOXA1.binding.sites, region) %>%
  dplyr::summarise(
    count = n())

RB1sum <- group_by(df.pairs.FINAL, RB1.binding.sites, region) %>%
  dplyr::summarise(
    count = n())

ETS1sum$total <- rep(c(46854, 14105, 34817), 2)
ETS1sum$percent <- ETS1sum$count / ETS1sum$total * 100

SP1sum$total <- rep(c(46854, 14105, 34817), 2)
SP1sum$percent <- SP1sum$count / SP1sum$total * 100

TP53sum$total <- rep(c(46854, 14105, 34817), 2)
TP53sum$percent <- SP1sum$count / SP1sum$total * 100

SP1sum$total <- rep(c(46854, 14105, 34817), 2)
SP1sum$percent <- SP1sum$count / SP1sum$total * 100



grouped.data <- data.frame(region = rep(c('commonHMD', 'commonPMD', 'neither'), 4),
                           binding.site = c(rep('AR', 3), rep('CTCF', 3), rep('ETS1', 3), rep('SP1', 3)),
                           count = c(ARsum$count[1:3], CTCFsum$count[1:3], ETS1sum$count[1:3], SP1sum$count[1:3]),
                           total = rep(c(46854, 14105, 34817), 4),
                           percent = c(ARsum$percent[1:3], CTCFsum$percent[1:3], ETS1sum$percent[1:3], SP1sum$percent[1:3]))


XXX.test <- df.pairs.FINAL[, c(17:25)]


i = 18

sum <- group_by(df.pairs.FINAL, df.pairs.FINAL[, i], region) %>%dplyr::summarise(count = n())
grouped.data$count[c(i-16, i-7, i+2)] <- sum$count[c(1:3)]
grouped.data$percent[c(i-16, i-7, i+2)] <- grouped.data$count[c(i-16, i-7, i+2)] / grouped.data$total[c(i-16, i-7, i+2)] * 100


CTCFsum <- group_by(df.pairs.FINAL, CTCF.binding.sites, region) %>%
  dplyr::summarise(
    count = n())

CTCFsum$total <- rep(c(46854, 14105, 34817), 2)
CTCFsum$percent <- CTCFsum$count / CTCFsum$total * 100

rm(grouped.data, ARsum, CTCFsum, ETS1sum, FOXA1sum, RB1sum, SP1sum, TP53sum)


#Comparing PC3 CTCF annots to CTCF TF binding site annots

PMD_ctcf_sum <- group_by(df.pairs.FINAL, PC3.annots, region) %>% dplyr::summarise(count = n())

#1287 total annots in the PC3 group for PMD regions that include CTCF


summary(df.pairs.FINAL$region)

1287/14105 * 100

#Discrepancy between PC3 annots and CTCF binding site annots might be because observed areas are much smaller than all the total possible CTCF binding sites







Lastly, for the regions I can summarize what annotations are present to see if there are any standout areas. 

```{r region_summary}

summary.regions <- group_by(df.probes.regions.genes, CGI.Relation, region) %>% dplyr::summarise(count = n())
summary.regions$percent <- summary.regions$count / rep(c(47109, 14339, 35007), 6) * 100

ggplot(summary.regions, aes(fill=region, y=percent, x=CGI.Relation)) + 
  geom_bar(position="dodge", stat="identity")
```








#Now to look at GSTP1
#Actually don't need to do anything, can just search it up in the df.region.pairs.FINAL

GSTP1 <- subset(df.pairs.FINAL, df.pairs.FINAL$annot.symbol == 'GSTP1')

#Can plot probe-gene correlation using ELMER

scatter.plot(data = PRAD.mae,
             byPair = list(probe = c("cg06928838"), gene = c("ENSG00000084207")), 
             category = "definition", save = FALSE, lm_line = TRUE) 

#The probe that I chose to plot against showed the strongest negative correlation, however I don't know if it is truly the best, probe-gene pair
#Ok so based on some background research the probe-gene pairs I have identified appear to be right



#probes.in.proms <- splicejam::annotateGRfromGR(probes.gr, promoters.gr)
#df.probes.proms <- as.data.frame(probes.in.proms)
#df.probes.proms <- df.probes.proms[, c(1, 2, 3, 8, 10)]
#names(df.probes.proms)[4] <- 'probe'
#names(df.probes.proms)[5] <- 'promoter.id'
#df.probes.proms <- na.omit(df.probes.proms)
#df.probes.proms <- subset(df.probes.proms, df.probes.proms$promoter.id %in% promoters$ID)

#df.probes.proms.FINAL <- aggregate(df.probes.proms[,4], list(df.probes.proms[,5]), function(x) paste0(unique(x)))
#names(df.probes.proms.FINAL)[1] <- 'promoter'
#names(df.probes.proms.FINAL)[2] <- 'probes'


#test2 <- subset(df.probes.proms.FINAL, df.probes.proms.FINAL$probes %in% rownames(PRADmeth))






