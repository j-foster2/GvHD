## Output Directories

mkdir Figure_1
mkdir Figure_2
mkdir Figure_3
mkdir Figure_4

## Figure 1 and Extended Data Figure 1


jid1=$(sbatch --mem=10g --time=0-1 --parsable -N 1 -n 1 -o ./slurmOut/generate_geneLevel_RNA_expression_Fig1b_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif ipython generate_geneLevel_RNA_expression_Fig1b.py")
jid2=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid1 -N 1 -n 1 -o ./slurmOut/transcritToGeneConversion_mm10_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript transcritToGeneConversion_mm10.R")
jid3=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid2 -N 1 -n 24 -o ./slurmOut/bruce_RNA_laurie_foster_K4me3_data_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif ipython bruce_RNA_laurie_foster_K4me3_data.py")
jid4=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid3 -N 1 -n 1 -o ./slurmOut/Figure_1b_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_1b.R")
jid5=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid4 -N 1 -n 4 -o ./slurmOut/generate_GBA_H3K4me3_count_matrix_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript generate_GBA_H3K4me3_count_matrix.R")
jid6=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid5 -N 1 -n 1 -o ./slurmOut/GBA_H3K4me3_DESeq2_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript GBA_H3K4me3_DESeq2.R")
jid7=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid6 -N 1 -n 1 -o ./slurmOut/ExtData_Figure_1_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif ipython ExtData_Figure_1.py")
jid8=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid7 -N 1 -n 4 -o ./slurmOut/Figure_1c_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif ipython Figure_1c.py")

## Figure 2


jid9=$(sbatch --mem=100g --time=0-1 --parsable --dependency=afterok:$jid8 -N 1 -n 20 -o ./slurmOut/integrate_multiome_replicate1_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript integrate_multiome_replicate1.R ")
jid10=$(sbatch --mem=100g --time=0-1 --parsable --dependency=afterok:$jid9 -N 1 -n 1 -o ./slurmOut/integrate_multiome_replicate2.R_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript integrate_multiome_replicate2.R ")
jid11=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid10 -N 1 -n 1 -o ./slurmOut/Figure_2b_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_2b.R")
jid12=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid11 -N 1 -n 1 -o ./slurmOut/Figure_2c_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_2c.R")
jid13=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid12 -N 1 -n 1 -o ./slurmOut/Figure_2d_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_2d.R")
jid14=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid13 -N 1 -n 1 -o ./slurmOut/ExtData_Figure_2d_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript ExtData_Figure_2d.R")
jid15=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid14 -N 1 -n 1 -o ./slurmOut/post_transplant_geneSet_Enrichr_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif ipython post_transplant_geneSet_Enrichr.py")
jid16=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid15 -N 1 -n 1 -o ./slurmOut/Figure_2e_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_2e.R")
jid17=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid16 -N 1 -n 1 -o ./slurmOut/Figure_2f_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_2f.R")
jid18=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid17 -N 1 -n 1 -o ./slurmOut/Figure_2g_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_2g.R")
jid19=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid18 -N 1 -n 1 -o ./slurmOut/Figure_2h_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_2h.R")
jid20=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid19 -N 1 -n 1 -o ./slurmOut/Figure_2i_m_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_2i_m.R")
jid21=$(sbatch --mem=100g --time=0-1 --parsable --dependency=afterok:$jid20 -N 1 -n 20 -o ./slurmOut/integrate_post_transplant_replicate2_GBA_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_2h.R")
jid22=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid21 -N 1 -n 1 -o ./slurmOut/Figure_2n_o_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_2n_o.R")
jid23=$(sbatch --mem=100g --time=0-1 --parsable --dependency=afterok:$jid22 -N 1 -n 20 -o ./slurmOut/integrate_all_multiome_data_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript integrate_all_multiome_data.R")
jid24=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid23 -N 1 -n 1 -o ./slurmOut/ExtData_Figure_2e_f_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript ExtData_Figure_2e_f.R")
jid25=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid24 -N 1 -n 1 -o ./slurmOut/ExtData_Figure_2g_j_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript ExtData_Figure_2g_j.R")
jid26=$(sbatch --mem=100g --time=0-1 --parsable --dependency=afterok:$jid25 -N 1 -n 20 -o ./slurmOut/integrate_post_transplant_replicate1_GBA_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript integrate_post_transplant_replicate1_GBA.R")
jid27=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid26 -N 1 -n 1 -o ./slurmOut/ExtData_Figure_2k_l_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript ExtData_Figure_2k_l.R")


## Figure 3


# Integration of ILC2 and pcILC2 multiome data
jid28=$(sbatch --mem=100g --time=0-1 --parsable --dependency=afterok:$jid27 -N 1 -n 20 -o ./slurmOut/integrate_inVitro_ILC2_pcILC2_replicate1_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript integrate_inVitro_ILC2_pcILC2_replicate1.R")

# Figure 3d
jid29=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid28 -N 1 -n 20 -o ./slurmOut/Figure_3d_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_3d.R")

# Figure 3e
jid30=$(sbatch --mem=20g --time=0-1 --parsable -N 1 -n 1 -o ./slurmOut/mouse_ATAC_peak_union_set.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif python mouse_ATAC_peak_union_set.py")
jid31=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid30 -N 1 -n 1 -o ./slurmOut/generate_pcILC2_ILC2_ATAC_count_matrix.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript generate_pcILC2_ILC2_ATAC_count_matrix.R")
jid32=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid31 -N 1 -n 1 -o ./slurmOut/ILC2_pcILC2_ATAC_DESeq2.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript ILC2_pcILC2_ATAC_DESeq2.R")
jid33=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid32 -N 1 -n 4 -o ./slurmOut/Figure_3e.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif ipython Figure_3e.py")

# Figure 3g
jid34=$(sbatch --mem=10g --time=0-0:15 --parsable --dependency=afterok:$jid33 -N 1 -n 1 -o ./slurmOut/Figure_3g.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_3g.R")

# Figure 3h
jid35=$(sbatch --mem=30g --time=0-1 --parsable --dependency=afterok:$jid34 -N 1 -n 1 -o ./slurmOut/Figure_3h.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_3h.R")

# Figure 3i
jid36=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid35 -N 1 -n 1 -o ./slurmOut/Figure_3i.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_3i.R")

# Figure 3j
jid37=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid36 -N 1 -n 1 -o ./slurmOut/Figure_3j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_3j.R")

# Figure 3k
jid38=$(sbatch --mem=75g --time=0-2 --parsable --dependency=afterok:$jid37 -N 1 -n 4 -o ./slurmOut/Figure_3k.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_3k.R")

# Extended Data Figure 3c
jid39=$(sbatch --mem=100g --time=0-1 --parsable --dependency=afterok:$jid38 -N 1 -n 20 -o ./slurmOut/integrate_inVitro_ILC2_pcILC2_replicate2.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript integrate_inVitro_ILC2_pcILC2_replicate2.R")
jid40=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid39 -N 1 -n 1 -o ./slurmOut/ExtData_Figure_3c.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript ExtData_Figure_3c.R")

# Extended Data Figure 3d
jid41=$(sbatch --mem=35g --time=0-1 --parsable --dependency=afterok:$jid40 -N 1 -n 1 -o ./slurmOut/ExtData_Figure_3d.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript ExtData_Figure_3d.R")

# Extended Data Figure 3e-f
jid42=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid41 -N 1 -n 1 -o ./slurmOut/ExtData_Figure_3e_f.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript ExtData_Figure_3e_f.R")

# Extended Data Figure 3g
jid43=$(sbatch --mem=10g --time=0-0:15 --parsable --dependency=afterok:$jid42 -N 1 -n 1 -o ./slurmOut/ExtData_Figure_3g.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript ExtData_Figure_3g.R")

# Extended Data Figure 3h
jid44=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid43 -N 1 -n 1 -o ./slurmOut/ExtData_Figure_3h.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript ExtData_Figure_3h.R")

# Extended Data Figure 3i
jid45=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid44 -N 1 -n 1 -o ./slurmOut/ExtData_Figure_3h.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript ExtData_Figure_3i.R")

## Figure 4

# Figure 4b
jid46=$(sbatch --mem=10g --time=0-1 --parsable -N 1 -n 1 -o ./slurmOut/human_ATAC_peak_union_set_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif ipython human_ATAC_peak_union_set.py")
jid47=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid46 -N 1 -n 1 -o ./slurmOut/generate_hILC2_pc-hILC2_ATAC_count_matrix_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript generate_hILC2_pc-hILC2_ATAC_count_matrix.R")
jid48=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid47 -N 1 -n 1 -o ./slurmOut/hILC2_pc-hILC2_ATAC_DESeq2_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript hILC2_pc-hILC2_ATAC_DESeq2.R")
jid49=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid48 -N 1 -n 1 -o ./slurmOut/Figure_4b_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif ipython Figure_4b.py")

# Figure 4c-d
jid50=$(sbatch --mem=15g --time=0-3 --parsable --dependency=afterok:$jid49 -N 1 -n 1 -o ./slurmOut/hILC2_ATAC_homer_analysis_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif ipython hILC2_ATAC_homer_analysis.py")
jid51=$(sbatch --mem=15g --time=0-3 --parsable --dependency=afterok:$jid49 -N 1 -n 1 -o ./slurmOut/pc-hILC2_ATAC_homerAnalysis_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif ipython pc-hILC2_ATAC_homerAnalysis.py")
jid52=$(sbatch --mem=10g --time=0-1 --parsable --dependency=afterok:$jid50:$jid51 -N 1 -n 1 -o ./slurmOut/Figure_4c_d_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif ipython Figure_4c_d.py")

# Figure 4f-g
jid52=$(sbatch --mem=100g --time=0-2 --parsable --dependency=afterok:$jid51 -N 1 -n 20 -o ./slurmOut/integrate_human_scATAC.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript integrate_human_scATAC.R")
jid53=$(sbatch --mem=100g --time=0-3 --parsable --dependency=afterok:$jid52 -N 1 -n 1 -o ./slurmOut/addPseudoExpression_integrated_human_scATAC_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript addPseudoExpression_integrated_human_scATAC.R")
jid54=$(sbatch --mem=40g --time=0-1 --parsable --dependency=afterok:$jid53 -N 1 -n 1 -o ./slurmOut/Figure_4f_g_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure_4f_g.R")

# Figure 4h-i and Extended Data Figure 4a
jid55=$(sbatch --mem=75g --time=0-1 --parsable --dependency=afterok:$jid54 -N 1 -n 20 -o ./slurmOut/scATAC_signal_at_bulk_ILC2_pcILC2_randomILC2_peaks_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript scATAC_signal_at_bulk_ILC2_pcILC2_sharedILC2_peaks.R")
jid56=$(sbatch --mem=20g --time=0-1 --parsable --dependency=afterok:$jid55 -N 1 -n 1 -o ./slurmOut/Figure4h_i_ExtData_Fig4a_%j.txt --wrap="apptainer exec --no-mount home --bind $PWD gvhd_R4.3.1.sif Rscript Figure4h_i_ExtData_Fig4a.R")
