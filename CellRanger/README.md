# CellRanger arc (v. 2.0.0) 

## Mouse pre- and post-transplant ILC2s


```
source /nas/longleaf/home/jfoster3/.bashrc; 

module load cellranger-arc/2.0.0;

cd /proj/dllab/jfoster/serody_project/results/tenX_scMultiome/ILC1_rep1/;

ssub --mem=75g --time=1-0 -n 24 --wrap=\" \\\
cellranger-arc count \
--id=ILC1_rep1 \
--reference=/proj/dllab/jfoster/serody_project/data/mm10_data/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
--libraries=/proj/dllab/jfoster/serody_project/results/tenX_scMultiome/ILC1_rep1/ICL1_rep1_libraries.csv \
--localcores=24 \
--localmem=75\"

```
