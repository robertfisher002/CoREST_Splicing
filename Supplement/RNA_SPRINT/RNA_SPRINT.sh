downloaded decoupler_environment_linux.yml from github page (https://github.com/frankligy/SNAF/blob/main/RNA-SPRINT/README.md)

module load miniconda
conda env create -f decoupler_environment_linux.yml
conda activate decoupler_env

downloaded prior.tsv from https://www.synapse.org/#!Synapse:syn53038679
made meta file to include sample names from annotation output, treatment, and cell line info
downloaded RNA_SPRINT.py from github https://github.com/frankligy/SNAF/blob/main/RNA-SPRINT/RNA_SPRINT.py

python RNA_SPRINT.py --splicing /projectnb/epigen/Robert/melanoma_splice/SNAF/sing/altanalyze_output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt --prior ./prior.tsv --outdir results --meta ./meta.txt 