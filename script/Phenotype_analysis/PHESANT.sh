#!/bin/bash

baseDir="$(pwd)"

for i in $1
do
# define path
dataDir="${baseDir}/result/PhenoAsso/PheWAS"
outputDir=${dataDir}/before
# results File
CatDir="${outputDir}/${i}/results/"
mkdir -p $CatDir
codeDir="${baseDir}/software/PHESANT-master"
# expouse file and  trait of interests name
expFile="${dataDir}/${i}-all_PheWASinput.csv"
trait_name="value"
 
# ---Step1 Running a phenome scan in UK Biobank-----
## pheno data File
phenodataFile="${dataDir}/beforeImagingVisit_data_PheWASinput.csv"
confounderFile="${dataDir}/cova_PheWASinput_simple.csv"

# run PHESANT
cd $codeDir/WAS
Rscript phenomeScan.r \
--phenofile="$phenodataFile" \
--traitofinterestfile="$expFile" \
--confounderfile="$confounderFile" \
--variablelistfile="${baseDir}/result/PhenoAsso/data/outcome_info.tsv" \
--datacodingfile="${baseDir}/result/PhenoAsso/data/data-coding-ordinal-info_reordered.txt" \
--traitofinterest="$trait_name" \
--resDir="${CatDir}"  \
--userId="userId"  \
--genetic=FALSE \
--standardise=TRUE

# ---Step 2 Post-processing of results----
codeDir2=${codeDir}/resultsProcessing/
cd $codeDir2
Rscript mainCombineResults.r \
--resDir="${CatDir}"  \
--variablelistfile="${baseDir}/result/PhenoAsso/data/outcome_info.tsv"

done
