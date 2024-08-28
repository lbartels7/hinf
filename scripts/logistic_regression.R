  library(dplyr)
  library(data.table)
  library(arrow)
  library(IRanges)
  
  

  df <- read_feather(snakemake@input[[1]])
  df$AMP <- as.factor(df$AMP)
  df <- subset(df, select = -c(samples, Institute, HaemoSeq_species, HaemoSeq_serotype,
       SequenceType, beta_lactamase, AMP_MIC, group, Haplotype, perc_reads_mapped,
       coverage_mapped_reads, libID, Bioproject, AccessionNumber,
       BioSample))


n_variants = length(df) -1 

log.reg.res <- data.table(
  id = character(n_variants),
  pos = integer(n_variants),
  alt.allele = character(n_variants),
  pvalue = numeric(n_variants),
  effect = numeric(n_variants),
  odds.ratio = numeric(n_variants)
)



for (i in head(seq_along(df),-1))  {
  id = colnames(df)[i]
  pos <- strsplit(id,'_')[[1]][2]
  allele <- strsplit(id, '_')[[1]][3]
  model <- glm(AMP ~ get(id), data=df, family = "binomial")
  set(
    log.reg.res,
    i = i,
    j = c("id", "pos", "alt.allele", "pvalue", "effect", 'odds.ratio'),
    value = list(
      id = id,
      pos = pos,
      alt.allele = allele,
      pvalue = coef(summary(model))[2, "Pr(>|z|)"],
      effect = coef(summary(model))[2, "Estimate"],
      odds.ratio = exp(coef(summary(model))[2, "Estimate"])
    )
  )

}





# Annotate variants with gene information

genes <- read.csv(snakemake@input[[2]], sep = '\t')

gene.ranges <- with(genes, IRanges(start_adjst, end_adjst))
zone.ind <- findOverlaps(log.reg.res$pos, gene.ranges, select="arbitrary")

log.reg.res$gene.name <- genes$name[zone.ind]
log.reg.res$gene.product <- genes$product[zone.ind]
log.reg.res$gene.type <- genes$type[zone.ind]
log.reg.res$gene.id <- genes$X..ID[zone.ind]


log.reg.res$p_values.adj <- p.adjust(log.reg.res$pvalue, method = "fdr", n = length(log.reg.res$pvalue))

write.csv(log.reg.res, snakemake@output[[1]])
