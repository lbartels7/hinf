library(arrow)
library(dplyr)
library(data.table)
library(ape)
library(IRanges)
library(RNOmni)



df <- read_feather(snakemake@input[[1]])
# df["AMP_MIC"][df["AMP_MIC"] == '>8'] <- '8'
df$AMP_MIC <- as.numeric(df$AMP_MIC)
# df <- df[df$AMP_MIC < 200,]
df$AMP_MIC.logscaled <- log(df$AMP_MIC,2)
df <- subset(df, select = -c(FullID, Perc.ReadsMapped, CoverageMappedReads,
                             origin, AMP, serotype, beta_lactamase, samples, AMP_MIC))
n_variants = length(df) -1 

lin.reg.res <- data.table(
  id = character(n_variants),
  pos = integer(n_variants),
  alt.allele = character(n_variants),
  pvalue = numeric(n_variants),
  effect = numeric(n_variants),
  r.squared = numeric(n_variants),
  adj.r.squared = numeric(n_variants)
)



for (i in head(seq_along(df),-1))  {
  id = colnames(df)[i]
  pos <- strsplit(id,'_')[[1]][2]
  allele <- strsplit(id, '_')[[1]][3]
  model <- lm(AMP_MIC.logscaled ~ get(id), df)
  set(
    lin.reg.res,
    i = i,
    j = c("id", "pos", "alt.allele", "pvalue", "effect", "r.squared", "adj.r.squared"),
    value = list(
      id = id,
      pos = pos,
      alt.allele = allele,
      pvalue = coef(summary(model))[2, "Pr(>|t|)"],
      effect = coef(summary(model))[2, "Estimate"],
      r.squared = summary(model)$r.squared,
      adj.r.squared = summary(model)$adj.r.squared
    )
  )
}


# create plots ------------------------------------------------------------



genes <- read.csv(snakemake@input[[2]], sep = '\t')

gene.ranges <- with(genes, IRanges(start_adjst, end_adjst))
zone.ind <- findOverlaps(lin.reg.res$pos, gene.ranges, select="arbitrary")
# lin.reg.res$gene <- genes$attributes[zone.ind]
lin.reg.res$gene.name <- genes$name[zone.ind]
lin.reg.res$gene.product <- genes$product[zone.ind]
lin.reg.res$gene.type <- genes$type[zone.ind]

# lin.reg.res$p_values.adj <- -log(p.adjust(lin.reg.res$pvalue, method = "fdr", n = length(lin.reg.res$pvalue)),10)
lin.reg.res$p_values.adj <- p.adjust(lin.reg.res$pvalue, method = "fdr", n = length(lin.reg.res$pvalue))

write.csv(lin.reg.res, snakemake@output[[1]])
