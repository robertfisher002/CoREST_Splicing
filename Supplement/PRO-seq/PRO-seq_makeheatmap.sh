#TSS
make_heatmap -b c -l s -s s -p ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_reverse.bedGraph  -- ../TSS_anchor.txt DMSO_out_TES.matrix -500 10 275
make_heatmap -b c -l s -s s -p ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_reverse.bedGraph  -- ../TSS_anchor.txt Corin_out_TSS.matrix -500 10 275

#exon 3' exc
make_heatmap -b c -a u -l s -s s -p ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_reverse.bedGraph  -- ../exc_3_mh.txt DMSO_exc_3.matrix -500 10 100
make_heatmap -b c -a u -l s -s s -p ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_reverse.bedGraph  -- ../exc_3_mh.txt Corin_exc_3.matrix -500 10 100

#exon 3' inc
make_heatmap -b c -a u -l c -s s -p ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_reverse.bedGraph  -- ../inc_3_mh.txt DMSO_inc_3.matrix -500 10 100
make_heatmap -b c -a u -l c -s s -p ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_reverse.bedGraph  -- ../inc_3_mh.txt Corin_inc_3.matrix -500 10 100

#exon 3' background
make_heatmap -b c -a u -l c -s s -p ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_reverse.bedGraph  -- ../bkgd_3_mh.txt DMSO_bkgd_3.matrix -500 10 100
make_heatmap -b c -a u -l c -s s -p ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_reverse.bedGraph  -- ../bkgd_3_mh.txt Corin_bkgd_3.matrix -500 10 100

#exon 5' exc
make_heatmap -b c -a u -l c -s s -p ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_reverse.bedGraph  -- ../exc_5_mh.txt DMSO_exc_5.matrix -500 10 100
make_heatmap -b c -a u -l c -s s -p ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_reverse.bedGraph  -- ../exc_5_mh.txt Corin_exc_5.matrix -500 10 100

#exon 5' inc
make_heatmap -b c -a u -l c -s s -p ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_reverse.bedGraph  -- ../inc_5_mh.txt DMSO_inc_5.matrix -500 10 100
make_heatmap -b c -a u -l c -s s -p ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_reverse.bedGraph  -- ../inc_5_mh.txt Corin_inc_5.matrix -500 10 100

#exon 5' background
make_heatmap -b c -a u -l c -s s -p ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_reverse.bedGraph  -- ../bkgd_5_mh.txt DMSO_bkgd_5.matrix -500 10 100
make_heatmap -b c -a u -l c -s s -p ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_reverse.bedGraph  -- ../bkgd_5_mh.txt Corin_bkgd_5.matrix -500 10 100


#Corin Upregulated
make_heatmap -b c -l s -s s -p ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_reverse.bedGraph  -- ../SKMEL5_dom_PRO_UP.txt DMSO_UP_PRO_domTSS.matrix -500 10 275
make_heatmap -b c -l s -s s -p ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_reverse.bedGraph  -- ../SKMEL5_dom_PRO_UP.txt Corin_UP_PRO_domTSS.matrix -500 10 275

#Corin Downregulated
make_heatmap -b c -l s -s s -p ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_reverse.bedGraph  -- ../SKMEL5_dom_PRO_DOWN.txt DMSO_DOWN_PRO_domTSS.matrix -500 10 275
make_heatmap -b c -l s -s s -p ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_reverse.bedGraph  -- ../SKMEL5_dom_PRO_DOWN.txt Corin_DOWN_PRO_domTSS.matrix -500 10 275

#Corin Background
make_heatmap -b c -l s -s s -p ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1077_DMSO_1_hg38_dedup_3pr_reverse.bedGraph  -- ../SKMEL5_dom_BACKGROUND_p0.1_lf0.25.txt DMSO_BKGD_p0.1_lf0.25_PRO_domTSS.matrix -500 10 275
make_heatmap -b c -l s -s s -p ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_forward.bedGraph -m ../bedgraph/Alani001_PR1078_Corin_1_hg38_dedup_3pr_reverse.bedGraph  -- ../SKMEL5_dom_BACKGROUND_p0.1_lf0.25.txt Corin_BKGD_p0.1_lf0.25_PRO_domTSS.matrix -500 10 275
