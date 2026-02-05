library(GSVA)
library(ssGSEA2)

# run_ssGSEA2(input.ds = ??,
#             output.prefix = "",
#             gene.set.databases = ??,
#             sample.norm.type = "none",    #c("rank", "log", "log.rank", "none"),
#             correl.type = "z.score",      #c("rank", "z.score", "symm.rank"),
#             statistic = "area.under.RES", #c("area.under.RES", "Kolmogorov-Smirnov"),
#             output.score.type = "NES",    #c("NES", "ES"),
#             weight = 0,
#             nperm = 1000,
#             min.overlap = 10,
#             combine.mode = c("combine.off", "combine.replace", "combine.add"),
#             output.directory = ".",
#             global.fdr = FALSE,
#             extended.output = TRUE,
#             log.file = "run.log")

# If weight==0, sample.norm.type and correl.type do not matter;
# If weight > 0, the combination of sample.norm.type and correl.type dictate how
# the gene expression values in input.ds are transformed to obtain the score, so
# use this setting with care as the transformations can skew scores towards +ve
# or -ve values)

# sample.norm.type=='rank' weights genes proportional to rank
# sample.norm.type=='log' can be used for log-transforming input data
# sample.norm.type=='none' uses actual expression values;
# combined with correl.type=='rank', genes are weighted by actual values
# correl.type=='z.score' standardizes the (normalized) input values before using
# them to calculate scores.

# If "combine.off", do not combine '*_UP' and '*_DN'versions in a single score.
# If "combine.replace", combine '*_UP' and '*_DN' versions in a single score.
# If "combine.add" combine '*_UP' and '*_DN' versions in a single score and add
# it but keeping the individual '*_UP' and '*_DN' versions.


ssgsea_par <- GSVA::ssgseaParam(exprData = x,
                                geneSets = x,
                                assay = NA_character_,
                                annotation = NA_character_,
                                minSize = 1,
                                maxSize = Inf,
                                alpha = 0.25,
                                normalize = TRUE)

ssgsea_es <- GSVA::gsva(param = ssgsea_par)
