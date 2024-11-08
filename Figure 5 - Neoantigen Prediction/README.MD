## Figure 6 - Corin produces neopeptides presented on human MHC. 
Here, we conducted neopeptide generation analysis based on the splicing changes identified from RNAseq using Splicetools and SNAF. Binding prediction values were determined by HLAthena and NetMHCpan. Candidate neopeptides were selected based on junction count, PSI value, bidning score, and IP-MS validation of binding.
  - `Neopeptide_Prediction.Rmd` contains the code used to process the neopeptide binding outputs from HLAthena and NetMHCpan and the neopeptide identification tools Splicetools and SNAF. 
  - `Neopeptide_Prediction_Datafile.xlsx` contains IP-MS data, neopeptides by cell line, HLAthena and NetMHCpan binding predictions, neopeptide canddiate summary, and a unique candidate list.
  - `scripts` directory contains the bash scripts for SNAF analysis.
