
#The following code is the workflow needed to access, filter, and download TCGA legacy gene expression and methylation data for different types of cancer (COAD is shown here)

#IMPORTANT: THE INITAL LOADING OF THE DATA STILL HAS TO BE DONE IN KBH_PMD_PROJECT1 BECAUSE IT HAS THE PROPER WORKING DIRECTORY

##################PRAD######################

#Finding and downloading expression data

query.exp.PRAD <- GDCquery(project = 'TCGA-PRAD',
                           data.category = 'Gene expression',
                           data.type = 'Gene expression quantification',
                           platform = 'Illumina HiSeq',
                           file.type = 'normalized_results',
                           legacy = TRUE
)


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
load(file = "~/Desktop/KBH_PMD_Project_master/MAEs/Methylation Data/TCGA-PRAD-meth.RData")

#Now on to creating the MAE for ELMER

PRAD.mae <- createMAE(exp = PRADexp,
                      met = PRADmeth,
                      TCGA = TRUE,
                      genome = 'hg19',
                      filter.probes = PMD.probe.promoters)

save(PRAD.mae, file = "TCGA-PRAD-MAE")


#Creating the MAE for PRAD using promoter probes NOT in PMDs

PRAD.nonPMD.mae <- createMAE(exp = PRADexp,
                      met = PRADmeth,
                      TCGA = TRUE,
                      genome = 'hg19',
                      filter.probes = not.PMD.promoter.probes)

save(PRAD.nonPMD.mae, file = "TCGA-PRAD-nonPMD-MAE")


##################COAD######################


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

save(COADmeth, file = "TCGA-COAD-meth.RData")
#Load() will go here upon subsequent reloadings 

#Now on to creating the MAE for ELMER

COAD.mae <- createMAE(exp = COADexp,
                      met = COADmeth,
                      TCGA = TRUE,
                      genome = 'hg19',
                      filter.probes = PMD.probe.promoters)

save(COAD.mae, file = "TCGA-COAD-MAE")

##################PAAD######################

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


##################STAD######################

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

##################SKCM######################

query.exp.SKCM <- GDCquery(project = 'TCGA-SKCM',
                           data.category = 'Gene expression',
                           data.type = 'Gene expression quantification',
                           platform = 'Illumina HiSeq',
                           file.type = 'normalized_results',
                           legacy = TRUE
)


query.exp.SKCM[[1]][[1]] <- query.exp.SKCM[[1]][[1]][!duplicated(query.exp.SKCM[[1]][[1]]$cases),]



GDCdownload(query.exp.SKCM)
SKCMexp <- GDCprepare(query.exp.SKCM)

rownames(SKCMexp) <- values(SKCMexp)$ensembl_gene_id

#now downloading the methylation data

query.meth.SKCM <- GDCquery(project = 'TCGA-SKCM',
                            data.category = 'DNA methylation',
                            platform = 'Illumina Human Methylation 450',
                            legacy = TRUE
)

GDCdownload(query.meth.SKCM, method = 'api')
SKCMmeth <- GDCprepare(query.meth.SKCM)

save(SKCMmeth, file = "TCGA-SKCM-meth.RData")


#Now on to creating the MAE for ELMER

SKCM.mae <- createMAE(exp = SKCMexp,
                      met = SKCMmeth,
                      TCGA = TRUE,
                      genome = 'hg19',
                      filter.probes = PMD.probe.promoters)

save(SKCM.mae, file = "TCGA-SKCM-MAE")


##################BRCA######################

query.exp.BRCA <- GDCquery(project = 'TCGA-BRCA',
                           data.category = 'Gene expression',
                           data.type = 'Gene expression quantification',
                           platform = 'Illumina HiSeq',
                           file.type = 'normalized_results',
                           legacy = TRUE 
)


query.exp.BRCA[[1]][[1]] <- query.exp.BRCA[[1]][[1]][!duplicated(query.exp.BRCA[[1]][[1]]$cases),]



GDCdownload(query.exp.BRCA)
BRCAexp <- GDCprepare(query.exp.BRCA)

rownames(BRCAexp) <- values(BRCAexp)$ensembl_gene_id


#now downloading the methylation data

query.meth.BRCA <- GDCquery(project = 'TCGA-BRCA',
                            data.category = 'DNA methylation',
                            platform = 'Illumina Human Methylation 450',
                            legacy = TRUE
)

GDCdownload(query.meth.BRCA, method = 'api')
BRCAmeth <- GDCprepare(query.meth.BRCA)

save(BRCAmeth, file = "TCGA-BRCA-meth.RData")


#Now on to creating the MAE for ELMER

BRCA.mae <- createMAE(exp = BRCAexp,
                      met = BRCAmeth,
                      TCGA = TRUE,
                      genome = 'hg19',
                      filter.probes = PMD.probe.promoters)

save(BRCA.mae, file = "TCGA-BRCA-MAE")

##################PRAD-real######################

query.exp.PRAD <- GDCquery(project = 'TCGA-PRAD',
                           data.category = 'Gene expression',
                           data.type = 'Gene expression quantification',
                           platform = 'Illumina HiSeq',
                           file.type = 'normalized_results',
                           legacy = TRUE)


query.exp.PRAD[[1]][[1]] <- query.exp.PRAD[[1]][[1]][!duplicated(query.exp.PRAD[[1]][[1]]$cases),]

GDCdownload(query.exp.PRAD)
PRADexp <- GDCprepare(query.exp.PRAD)

rownames(PRADexp) <- values(PRADexp)$ensembl_gene_id


save(PRADexp, file = "TCGA-PRAD-expression.RData")

#now downloading the methylation data

query.meth.PRAD <- GDCquery(project = 'TCGA-PRAD',
                            data.category = 'DNA methylation',
                            platform = 'Illumina Human Methylation 450',
                            legacy = TRUE)

GDCdownload(query.meth.PRAD, method = 'api')
PRADmeth <- GDCprepare(query.meth.PRAD)

save(PRADmeth, file = "TCGA-PRAD-methylation.RData")

#Now to create a MAE that contains ALL info (sorting out PMDs/promoters/etc can be done later)


mae.PRAD <- createMAE(exp = PRADexp,
                      met = PRADmeth,
                      met.platform = '450K',
                      genome = 'hg19',
                      save = TRUE,
                      save.filename = "PRAD-MAE",
                      TCGA = TRUE)

summary(complete.cases(mae.PRAD))

intersectColumns(mae.PRAD)

##################BRCA-real######################

query.exp.BRCA <- GDCquery(project = 'TCGA-BRCA',
                           data.category = 'Gene expression',
                           data.type = 'Gene expression quantification',
                           platform = 'Illumina HiSeq',
                           file.type = 'normalized_results',
                           legacy = TRUE 
)


query.exp.BRCA[[1]][[1]] <- query.exp.BRCA[[1]][[1]][!duplicated(query.exp.BRCA[[1]][[1]]$cases),]



GDCdownload(query.exp.BRCA)
BRCAexp <- GDCprepare(query.exp.BRCA)

rownames(BRCAexp) <- values(BRCAexp)$ensembl_gene_id


#now downloading the methylation data

query.meth.BRCA <- GDCquery(project = 'TCGA-BRCA',
                            data.category = 'DNA methylation',
                            platform = 'Illumina Human Methylation 450',
                            legacy = TRUE
)

GDCdownload(query.meth.BRCA, method = 'api')
BRCAmeth <- GDCprepare(query.meth.BRCA)

save(BRCAmeth, file = "TCGA-BRCA-meth.RData")


#Now on to creating the MAE for ELMER

BRCA.mae <- createMAE(exp = BRCAexp,
                      met = BRCAmeth,
                      met.platform = '450K',
                      TCGA = TRUE,
                      genome = 'hg19',
                      save.filename = 'BRCA-MAE')

intersectColumns(BRCA.mae)











