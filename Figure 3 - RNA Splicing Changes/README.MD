## Figure 3 - Corin induces RNA splicing changes in melanoma
Here, RNAseq (100M, 150bp) conducted on 6 melanoma cell lines (451Lu, 1205Lu, SKMEL5, SKMEL24, SKMEL28, WM793) treated +/- corin  (24h, 2.5uM) was used to determine if corin induces RNA splicing changes in melanoma. rMATs was used to detect SE, A5SS, A3SS, and MXE events while LeafCutter was used to detect the number of RI events. Significant events were called as: q < 0.05, deltaPSI >= |0.1|.  
  - `New RMATS analysis.Rmd` contains the code used to create filter raw outputes for splicing events, create the stacked bar plot, UpSet plot, heatmap, and GO dotplot for figrue 3. 
  - `Splicing_Datafile.xlsx` contains the raw data for Figure 3 (significant SE splicing counts, splicing events per cell line, commone SE events, significant GO pathways) and outputs for the R script.
  - `scripts` directory contains the rMATs and LeafCutter bash scripts run in linus
  - `raw_splicing_files` directory contains the raw outputs for rMATs and LeafCutter pre-filtering (directly outputs from scripts)
