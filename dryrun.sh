snakemake all -n --jobname "{jobid}" --jobs 100 \
		--keep-going \
		--rerun-incomplete \
		--snakefile Snakefile \
		--use-conda \
		--use-singularity \
		--singularity-args "--bind /lab" \
		--printshellcmds \
		--cluster-config config_cluster.yaml \
		--cluster "sbatch --output {cluster.output} --time {cluster.time} --mem {cluster.mem} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --parsable" \
		> snakemake.log 2>&1 &
