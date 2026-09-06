GWAS_analysis() {
    i=$1
    mkdir -p result/chr${i}
	
	plink2 --bfile UKB_gene_v3_imp_qc/UKB_gene_v3_imp_qc_chr${i}  \
		--keep data/brainimaging_british_subjs.txt \
		--extract data/snps_INFO_gt0.5.txt \
		--hwe 1e-6 \
		--maf 0.01 \
		--geno 0.05 \
		--mind 0.05 \
		--glm hide-covar \
		--pheno result/GWAS_CholinergicPathway_MNI_all.pheno \
		--covar result/GWAS_CholinergicPathway_MNI_all.cova \
		--quantile-normalize \
		--input-missing-phenotype -9 \
		--threads 20 \
		--out result/chr${i}/GWAS_chr${i}
}
toProcess=($(seq 1 1 22))
echo ${toProcess[@]}
export -f GWAS_analysis
nohup parallel -j 8 GWAS_analysis {} ::: "${toProcess[@]}" &