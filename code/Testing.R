
`%notin%` = Negate(`%in%`)


normalmeth <- subset(PRADmeth, colnames(PRADmeth) %notin% sample.map$patient.id)

normal.meth.sd <- rowSds(PRADmeth[ , sample.map$patient.id], na.rm = TRUE)
XXX <- data.frame(id = rownames(PRADmeth),
                  sd = tumor.meth.sd)
tumor.region.meth.sd <- df.pairs.FINAL
tumor.region.meth.sd$tumormethsd <- subset(XXX$sd, XXX$id %in% df.pairs.FINAL$ID)

ggplot(tumor.region.meth.sd, aes(x=region, y=tumormethsd, fill=region)) + 
  geom_violin()

df.pairs.hi.corr2 <- subset(tumor.region.meth.sd, tumor.region.meth.sd$r <(-0.5) | tumor.region.meth.sd$r >(0.5))

ggplot(df.pairs.hi.corr2, aes(x=region, y=tumormethsd, fill=region)) + 
  geom_violin()

names(PRADmeth)

colnames(PRADmeth)

colnames(PRADmeth.r)