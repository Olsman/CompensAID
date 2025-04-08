ff <- read.FCS(paste0("/mnt/external/it_nfs_share_immuno/CBBI_Projects/CUPS_DR196/processed/rolsman/temp/eurflowBALL/preprocessed.peacoqc/KSA-0212-Tube2.fcs"))
ff.red <- FlowSOM::AggregateFlowFrames("/mnt/external/it_nfs_share_immuno/CBBI_Projects/CUPS_DR196/processed/rolsman/temp/eurflowBALL/preprocessed.peacoqc/KSA-0212-Tube2.fcs",
                                   cTotal = 10000)
range.value <- 4
events.value <- 20
center.plot <- 2
separation.distance <- 0.25

functions <- list.files("/data/rolsman/AnalysisRosan_R/euroflow_bcp.all/Rosan/Compensation/Final/Functions/")
for (f in functions) {
  source(paste0("/data/rolsman/AnalysisRosan_R/euroflow_bcp.all/Rosan/Compensation/Final/Functions/", f))
}

rm(f, functions)

nrow(ff.red@exprs)

res <- CompensAID(ff = ff.red,
                  range.value = range.value,
                  events.value = events.value,
                  center.plot = center.plot,
                  separation.distance = separation.distance)

ggcyto::autoplot(ff.red,
                 x = "CD19",
                 y = "IgL",
                 bins = 100)

