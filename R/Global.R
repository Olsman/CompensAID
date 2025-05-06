# CalculateSSI
utils::globalVariables(c("mfi.pos", "mfi.neg", "sd.neg"))

# DensityGating
utils::globalVariables(c("cutoffs"))

# EmptyMatrixInfo
utils::globalVariables(globalVariables(c("file", "primary.marker", "primary.channel", "pretty.primary", "secondary.marker",
                                       "secondary.channel", "pretty.secondary", "segment", "secondary.cutoff", "primary.cutoff.neg",
                                       "primary.cutoff.pos", "segment.min", "segment.max", "event.count", "event.count.merge",
                                       "mfi.neg", "sd.neg", "mfi.pos", "ssi", "message", "primary.cutoff.adjusted",
                                       "primary.perc.neg.before", "primary.perc.neg.after", "primary.perc.pos.before",
                                       "primary.perc.pos.after", "secondary.cutoff.adjusted", "secondary.perc.before", "secondary.perc.after",
                                       "ff", "range.value", "mc", "co", "separation.distance", ".")))

# plotMatrix/plotDotSSI
utils::globalVariables(globalVariables(c("p", "value", "Var1", "Var2", "value2", "color")))
