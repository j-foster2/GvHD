ssub --mem=100g --time=0-2 -N 1 -n 20 --wrap=\"apptainer exec --no-mount home --bind $PWD ../Test_Environment/gvhd.sif Rscript integrate_multiome_replicate2.R\"
ssub --mem=10g --time=0-0:45 -N 1 -n 1 --wrap=\"apptainer exec --no-mount home --bind $PWD ../Test_Environment/gvhd.sif Rscript Figure_2b.R\"
ssub --mem=10g --time=0-0:45 -N 1 -n 1 -o ./slurmOut/ExtData_Figure_2g_j --wrap=\"apptainer exec --no-mount home --bind $PWD ../Test_Environment/gvhd.sif Rscript ExtData_Figure_2g_j.R\"
