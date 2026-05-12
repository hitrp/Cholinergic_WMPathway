GWAS() {
    i=$1
    mkdir -p result/chr${i}
	
	plink2 --bfile UKB_gene_v3_imp_qc_chr${i}  \
		--keep british_subjs.txt \
		--extract snps_INFO_gt0.5.txt \
		--hwe 1e-6 \
		--maf 0.01 \
		--geno 0.05 \
		--mind 0.05 \
		--glm hide-covar allow-no-covars \
		--pheno phenotypes.pheno \
		--quantile-normalize \
		--input-missing-phenotype -9 \
		--threads 100 \
		--out result/chr${i}/GWAS_chr${i}
}
toProcess=($(seq 1 1 22))
echo ${toProcess[@]}
export -f GWAS
nohup parallel -j 8 GWAS {} ::: "${toProcess[@]}" &