library(arrow)
library(dplyr)
library(data.table)
library(ape)
library(IRanges)
library(RNOmni)
library(stringr)



lmm.results <- read.csv(snakemake@input[[1]], sep = '\t') %>%
                mutate(position = as.integer(str_split_i(variant, '_', i=3)))
genes <- read.csv(snakemake@input[[2]], sep = '\t')

gene.ranges <- with(genes, IRanges(start_adjst, end_adjst))
zone.ind <- IRanges::findOverlaps(lmm.results$position, gene.ranges, select="arbitrary")
# lmm.results$gene <- genes$attributes[zone.ind]
lmm.results$gene.name <- genes$name[zone.ind]
lmm.results$gene.product <- genes$product[zone.ind]
lmm.results$gene.type <- genes$type[zone.ind]


write.csv(lmm.results, snakemake@output[[1]])
