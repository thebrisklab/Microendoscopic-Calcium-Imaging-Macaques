# Microendoscopic-Calcium-Imaging-Macaques
R code supporting the synchrony analysis in the paper "Microendoscopic calcium imaging in supplementary motor area and primary motor cortex of rhesus macaques at rest and during arm movement"

- The file "1_JaccardIndexAnalysisAllSessions.R" conducts the synchrony/co-activation analysis that appears in Figure 3C. This file needs to be run first. See comments about toggling the file to run the full analysis (lines 39 and 61, set if(1)) versus running the plots only (set if(0)). 

- The file "2_JaccardDistancePlots_v2.R" creates the plots of the normalized Jaccard index for Monkey U and the relationship between normalized Jaccard and distance, which appear in Figures 3A and B.

- The file "permfun.R" contains R functions that are used in "1_JaccardIndexAnalysisAllSessions.R"
