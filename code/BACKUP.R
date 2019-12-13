df.probes.proms.FINAL$length <- NA

for (i in seq_len(nrow(df.probes.proms.FINAL))) {
  df.probes.proms.FINAL$length[i] <- length(df.probes.proms.FINAL$probe[[i]])
}


df.probes.proms.FINAL$mean.sd <- NA
probes <- apply(df.probes.proms.FINAL, 1, function(x) x$probe)


for (i in seq_len(nrow(df.probes.proms.FINAL))) {
  if (df.probes.proms.FINAL$length[i] == 1){
    df.probes.proms.FINAL$mean.sd[i] <- NA
  }
  else {
    df.probes.proms.FINAL$mean.sd[i] <- mean(colSds(PRADmeth[probes[[i]],], na.rm=T), na.rm=T)
  }
}

df.probes.proms.FINAL$region <- as.factor(df.probes.proms.FINAL$region)

ggplot(df.probes.proms.FINAL, aes(x = region, y=mean.sd)) +
  geom_boxplot(fill="slateblue", alpha=0.2)





## tumor purity
LUMP <- PRAD.tp$LUMP

## Adjust for tumor purity

beta_res <- apply(PRADmeth.tumors, 1, function(x) {
  mod<-lm(x~LUMP, na.action=na.exclude)
  pval<-pf(summary(mod)$fstatistic[1], summary(mod)$fstatistic[2], summary(mod)$fstatistic[3], lower.tail = FALSE)
  if(pval<0.01){
    residuals(mod)
  }else{
    scale(x, scale = FALSE)[,1, drop=T]
  }
})

