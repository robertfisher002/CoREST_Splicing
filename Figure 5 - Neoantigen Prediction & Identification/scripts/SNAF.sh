Method
Follow this: https://snaf.readthedocs.io/en/latest/installation.html
Installation
Download NetMHCpan4 with wget to scc within SNAF folder (follow video tutorial on installation page of SNAF)
Download data file with wget from SNAF page
Have to use singularity since docker is not allowed on BU scc. This was a bit trickier than the tutorial explained so I found a resource through BU and followed that: http://rcs.bu.edu/examples/machine_learning/singularity/
[rjfisher@scc1 SNAF]$ export SINGULARITY_CACHEDIR=/projectnb/epigen/Robert/melanoma_splice/SNAF/.singularity/
[rjfisher@scc1 SNAF]$ export SINGULARITY_TMPDIR=/projectnb/epigen/Robert/melanoma_splice/SNAF/.singularity/
[rjfisher@scc1 SNAF]$ export SINGULARITY_LOCALCACHEDIR=/projectnb/epigen/Robert/melanoma_splice/SNAF/.singularity/
[rjfisher@scc1 SNAF]$ singularity pull docker://frankligy123/altanalyze:0.7.0.1


# assuming you are at a empty folder /path/to/my/project
mkdir sing
cd sing
mkdir .singularity
export SINGULARITY_CACHEDIR=projectnb/epigen/Robert/melanoma_splice/SNAF/sing/.singularity/
export SINGULARITY_TMPDIR=projectnb/epigen/Robert/melanoma_splice/SNAF/sing/.singularity/
export SINGULARITY_LOCALCACHEDIR=projectnb/epigen/Robert/melanoma_splice/SNAF/sing/.singularity/

singularity build --sandbox altanalyze/ docker://frankligy123/altanalyze:0.7.0.1
cd altanalyze
mkdir usr3
Run
#!/bin/bash -l
#$ -pe omp 24
#$ -P epigen

module load python2
singularity run -B $PWD:/mnt --writable altanalyze/ identify bam 24
*will need to account for replicates and treatment. So make sure eventual neopeptide shows up in both corin replicates and with greater PSI than dmso

HLA typing - use SpecHLA: https://github.com/deepomicslab/SpecHLA. Then check with cellosaurus for validation. Ultimately, pretty much followed what cellosaurus had but it was a good check for the cells/RNAseq data which aligned nicely with what HLAs have been reported for this cell line.

# Identify columns to modify: need to be a speciifc format - see SNAF tutorial for example
hla_columns <- grep("^hla", names(hla), value = TRUE)

# Loop through the identified columns
for (col in hla_columns) {
  # Use sub to replace everything after and including the second ":"
  hla[[col]] <- sub("^(.*?:.*?):.*$", "\\1", hla[[col]])
  hla[[col]] <- paste0("HLA-", hla[[col]])
}


hla_combined <- hla %>%
  unite(hla, starts_with("hla"), sep = ",", na.rm = TRUE)

# Print the result
print(hla_combined)

First, create the env with conda, and activate the env.

git clone https://github.com/deepomicslab/SpecHLA.git --depth 1
cd SpecHLA/
conda env create --prefix=./spechla_env -f environment.yml
conda activate ./spechla_env

Second, make the softwares in bin/ executable.
chmod +x -R bin/*

Third, index the database and install the packages.
unset LD_LIBRARY_PATH && unset LIBRARY_PATH 
bash index.sh

Extract HLA related reads
bash script/ExtractHLAread.sh -s SKMEL5_1 -b /projectnb/epigen/Robert/melanoma_splice/SNAF/bam/D1.bam -r hg38 -o ../SKMEL5
bash script/ExtractHLAread.sh -s SKMEL5_2 -b /projectnb/epigen/Robert/melanoma_splice/SNAF/bam/D2.bam -r hg38 -o ../SKMEL5

bash script/ExtractHLAread.sh -s 451Lu_2 -b /projectnb/epigen/Robert/melanoma_splice/SNAF/bam/B2.bam -r hg38 -o ../451Lu
bash script/ExtractHLAread.sh -s 451Lu_1 -b /projectnb/epigen/Robert/melanoma_splice/SNAF/bam/B1.bam -r hg38 -o ../451Lu


bash script/ExtractHLAread.sh -s 1205Lu_2 -b /projectnb/epigen/Robert/melanoma_splice/SNAF/bam/A2.bam -r hg38 -o ../1205Lu
bash script/ExtractHLAread.sh -s 1205Lu_1 -b /projectnb/epigen/Robert/melanoma_splice/bam/A1.bam -r hg38 -o ../1205Lu

bash script/ExtractHLAread.sh -s WM793_2 -b /projectnb/epigen/Robert/melanoma_splice/bam/C2.bam -r hg38 -o ../WM793
bash script/ExtractHLAread.sh -s WM793_1 -b /projectnb/epigen/Robert/melanoma_splice/bam/C1.bam -r hg38 -o ../WM793

bash script/ExtractHLAread.sh -s SKMEL24_2 -b /projectnb/epigen/Robert/melanoma_splice/bam/E2.bam -r hg38 -o ../SKMEL24
bash script/ExtractHLAread.sh -s SKMEL24_1 -b /projectnb/epigen/Robert/melanoma_splice/bam/E1.bam -r hg38 -o ../SKMEL24

bash script/ExtractHLAread.sh -s SKMEL28_2 -b /projectnb/epigen/Robert/melanoma_splice/SKMEL28/bam/2.bam -r hg38 -o ../SKMEL28
bash script/ExtractHLAread.sh -s SKMEL28_1 -b /projectnb/epigen/Robert/melanoma_splice/SKMEL28/bam/1.bam -r hg38 -o ../SKMEL28

Perform SpecHLA with
bash script/whole/SpecHLA.sh -n SKMEL5_1 -1 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/SKMEL5/SKMEL5_1_extract_1.fq.gz -2 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/SKMEL5/SKMEL5_1_extract_2.fq.gz -o ../SKMEL5

bash script/whole/SpecHLA.sh -n SKMEL5_2 -1 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/SKMEL5/SKMEL5_2_extract_1.fq.gz -2 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/SKMEL5/SKMEL5_2_extract_2.fq.gz -o ../SKMEL5

bash script/whole/SpecHLA.sh -n 451Lu_2 -1 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/451Lu/451Lu_2_extract_1.fq.gz -2 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/451Lu/451Lu_2_extract_2.fq.gz -o ../451Lu

bash script/whole/SpecHLA.sh -n 451Lu_1 -1 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/451Lu/451Lu_1_extract_1.fq.gz -2 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/451Lu/451Lu_1_extract_2.fq.gz -o ../451Lu

bash script/whole/SpecHLA.sh -n 1205Lu_2 -1 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/1205Lu/1205Lu_2_extract_1.fq.gz -2 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/1205Lu/1205Lu_2_extract_2.fq.gz -o ../1205Lu

bash script/whole/SpecHLA.sh -n 1205Lu_1 -1 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/1205Lu/1205Lu_1_extract_1.fq.gz -2 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/1205Lu/1205Lu_1_extract_2.fq.gz -o ../1205Lu

bash script/whole/SpecHLA.sh -n WM793_2 -1 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/WM793/WM793_2_extract_1.fq.gz -2 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/WM793/WM793_2_extract_2.fq.gz -o ../WM793

bash script/whole/SpecHLA.sh -n WM793_1 -1 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/WM793/WM793_1_extract_1.fq.gz -2 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/WM793/WM793_1_extract_2.fq.gz -o ../WM793

bash script/whole/SpecHLA.sh -n SKMEL24_2 -1 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/SKMEL24/SKMEL24_2_extract_1.fq.gz -2 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/SKMEL24/SKMEL24_2_extract_2.fq.gz -o ../SKMEL24

#!/bin/bash -l
#$ -pe omp 8
#$ -P epigen

module load miniconda
conda activate ./spechla_env

bash script/whole/SpecHLA.sh -n SKMEL28_2 -1 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/SKMEL28/SKMEL28_2_extract_1.fq.gz -2 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/SKMEL28/SKMEL28_2_extract_2.fq.gz -o ../SKMEL28

#!/bin/bash -l
#$ -pe omp 8
#$ -P epigen

module load miniconda
conda activate ./spechla_env

bash script/whole/SpecHLA.sh -n SKMEL28_1 -1 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/SKMEL28/SKMEL28_1_extract_1.fq.gz -2 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/SKMEL28/SKMEL28_1_extract_2.fq.gz -o ../SKMEL28

#!/bin/bash -l
#$ -pe omp 8
#$ -P epigen

module load miniconda
conda activate ./spechla_env

bash script/whole/SpecHLA.sh -n SKMEL24_1 -1 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/SKMEL24/SKMEL24_1_extract_1.fq.gz -2 /projectnb/epigen/Robert/melanoma_splice/SNAF/HLA_typing/SKMEL24/SKMEL24_1_extract_2.fq.gz -o ../SKMEL24



Getting neopeptide predictions
#!/bin/bash -l
#$ -pe omp 24
#$ -P epigen

module load python3/3.7
source /projectnb/epigen/Robert/melanoma_splice/SNAF/sing/SNAF_env/bin/activate

import os,sys
import pandas as pd
import numpy as np
import anndata as ad
import snaf

# read in the splicing junction matrix
df = pd.read_csv('altanalyze_output/ExpressionInput/counts.original.pruned.txt',index_col=0,sep='\t')

# database directory (where you extract the reference tarball file) and netMHCpan folder
db_dir = '/projectnb/epigen/Robert/melanoma_splice/SNAF/sing/data'
netMHCpan_path = '/projectnb/epigen/Robert/melanoma_splice/SNAF/sing/netMHCpan-4.1/netMHCpan'

# demonstrate how to add additional control database, see below note for more
tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db}

# initiate
snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)

jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=24,add_control=add_control,outdir='result')

sample_to_hla = pd.read_csv('HLA_SNAF.txt',sep='\t',index_col=0)['hla'].to_dict()
hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]

jcmq.run(hlas=hlas,outdir='./result')

snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')



#IP/MS 
import os,sys
import pandas as pd
import numpy as np
import anndata as ad
import snaf

df = pd.read_csv('altanalyze_output/ExpressionInput/counts.original.pruned.txt',index_col=0,sep='\t')

jcmq = snaf.JunctionCountMatrixQuery.deserialize('result/after_prediction.p')
snaf.chop_normal_pep_db(fasta_path='./human_proteome_uniprot.fasta',output_path='./fasta/human_proteome_uniprot_9_10_mers_unique.fasta',mers=[9,10],allow_duplicates=False)
for sample in df.columns:  # df is the counts.original.pruned.txt file
    jcmq.show_neoantigen_as_fasta(outdir='./fasta',name='neoantigen_{}.fasta'.format(sample),stage=3,verbosity=1,contain_uid=True,sample=sample)
    snaf.remove_redundant('./fasta/neoantigen_{}.fasta'.format(sample),'./fasta/neoantigen_{}_unique.fasta'.format(sample))
    snaf.compare_two_fasta(fa1_path='./fasta/human_proteome_uniprot_9_10_mers_unique.fasta',
                        fa2_path='./fasta/neoantigen_{}_unique.fasta'.format(sample),outdir='./fasta',
                        write_unique2=True,prefix='{}_'.format(sample))


#!/bin/bash -l
#$ -pe omp 4
#$ -P epigen

module load python3/3.7
source /projectnb/epigen/Robert/melanoma_splice/SNAF/sing/SNAF_env/bin/activate

python3 neoantigen_mers.py

The above created unique fasta files for each sample. Unique means the peptides that were not repeated in the file and did not show up in the uniprot database. My plan is to run all the files through maxquant to identify neopeptides that appear in each sample. Then overlap replicates to find neopeptides found across both corin replicates and then remove those that appear in any DMSO sample. 
download bruker specific base file (https://github.com/frankligy/SNAF/tree/main/maxquant)
Frank sent the maxquant version that should be compatible ( 2.0.3.1) via dropbox since maxquant doesnt archive previous versions

#CORIN
import os,sys
import pandas as pd
import numpy as np
import anndata as ad
import snaf

dbs = ['/projectnb/epigen/Robert/melanoma_splice/SNAF/sing/fasta/J1.bed_unique2.fasta',
       '/projectnb/epigen/Robert/melanoma_splice/SNAF/sing/fasta/J2.bed_unique2.fasta',
       '/projectnb/epigen/Robert/melanoma_splice/SNAF/sing/human_proteome_uniprot.fasta']
inputs = ['/projectnb/epigen/Robert/melanoma_splice/SNAF/sing/MassSpec/Corin/2024-01-09-Corin_Slot2-45_1_2446.d']
outdir = '/projectnb/epigen/Robert/melanoma_splice/SNAF/sing/MassSpec/Corin'
snaf.proteomics.set_maxquant_configuration(base='/projectnb/epigen/Robert/melanoma_splice/SNAF/sing/mqpar.xml',dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=4,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=11)


#!/bin/bash -l
#$ -pe omp 20
#$ -P epigen

module load python3/3.7
source /projectnb/epigen/Robert/melanoma_splice/SNAF/sing/SNAF_env/bin/activate

python3 MassSpec_Corin_SKMEL5.py

#run maxquant withe the modified .xml file

#!/bin/bash -l
#$ -pe omp 20
#$ -P epigen

module load mono
module load dotnet

# run maxquant using the downloaded .exe and the modifed mqpar.xml file
mono ../../MaxQuant_2.0.3.1/bin/MaxQuantCmd.exe mqpar.xml

#DMSO
import os,sys
import pandas as pd
import numpy as np
import anndata as ad
import snaf

dbs = ['/projectnb/epigen/Robert/melanoma_splice/SNAF/sing/fasta/D1.bed_unique2.fasta',
       '/projectnb/epigen/Robert/melanoma_splice/SNAF/sing/fasta/D2.bed_unique2.fasta',
       '/projectnb/epigen/Robert/melanoma_splice/SNAF/sing/human_proteome_uniprot.fasta']
inputs = ['/projectnb/epigen/Robert/melanoma_splice/SNAF/sing/MassSpec/DMSO/2024-01-08-DMSO_Slot2-42_1_2439.d']
outdir = '/projectnb/epigen/Robert/melanoma_splice/SNAF/sing/MassSpec/DMSO'
snaf.proteomics.set_maxquant_configuration(base='/projectnb/epigen/Robert/melanoma_splice/SNAF/sing/mqpar.xml',dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=4,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=11)


#!/bin/bash -l
#$ -pe omp 20
#$ -P epigen

module load python3/3.7
source /projectnb/epigen/Robert/melanoma_splice/SNAF/sing/SNAF_env/bin/activate

python3 MassSpec_DMSO_SKMEL5.py


#!/bin/bash -l
#$ -pe omp 20
#$ -P epigen

module load mono
module load dotnet

# run maxquant using the downloaded .exe and the modifed mqpar.xml file
mono ../../MaxQuant_2.0.3.1/bin/MaxQuantCmd.exe mqpar.xml