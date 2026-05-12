#!/bin/bash
# conda activate ldsc

# Define the paths to the CSV files
resultPath=result/heritability_ldsc
secondPath=result/sumstats
outputPath=result/gc
mkdir -p $outputPath
cd $outputPath
allOutcome=(anxiety adhd BP MDD ASD SCZ PD AD MS ALS stroke_Any stroke_AIS stroke_LAS stroke_CES stroke_SVS LBD FTD)

# Loop through each metric
for metricName in $1; do
    for outcomeName in ${allmetrics[@]};do
        # Run LDSC analysis
        python ldsc/ldsc.py \
            --rg ${resultPath}/${metricName}.sumstats,${secondPath}/${outcomeName}.sumstats \
            --ref-ld-chr ldsc/data/eur_w_ld_chr/ \
            --w-ld-chr ldsc/data/eur_w_ld_chr/ \
            --out ${outputPath}/${metricName}_${outcomeName}_gc
    done
done