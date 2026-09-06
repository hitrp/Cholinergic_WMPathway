# logistic analysis
library(data.table)
library(speedglm)
allphenos <- as.data.frame(fread("result/PhenoAsso/PheWAS/phenoAnalysis_CholinergicPathway_MNI_all_deconf.txt"))
allphenos$IID <- NULL
colnames(allphenos)[1] <- 'eid'
prefix <- colnames(allphenos)[2:ncol(allphenos)]
covafile <- read.table("result/PhenoAsso/Disease/cova_disease.csv",sep=',',header=T)
used_covnames <- colnames(covafile)
used_covnames <- used_covnames[!(used_covnames %in% c('age0','eid'))]
used_covnames <- used_covnames[!grepl('Site|eTIV', used_covnames)]
covadata <- covafile[,c('eid','age0',used_covnames)]
targetfiles <- list.files('result/PhenoAsso/Disease/Target_data_used/',pattern = '.csv')
# in all subjs, with more covariates
list_phenofile <- c()
list_targetfile <- c()
list_numberAll <- c()
list_numberTarget <- c()
list_pvalues <- c()
list_zvalues <- c()
list_coef <- c()
list_secoef <- c()
list_FDR_separate <- c()
for(pheno in unique(prefix)){
  phenodata <- allphenos[,c('eid',pheno)]
  pvalues <- c()
  for(j in seq(1,length(targetfiles))){
    targetfile <- targetfiles[j]
    list_phenofile <- c(list_phenofile, pheno)
    list_targetfile <- c(list_targetfile, targetfile)
    # convert variables
    tempdata <- merge(phenodata, covadata, by=c('eid'))
    targetdata <- as.data.frame(fread(paste0('result/PhenoAsso/Disease/Target_data_used/',targetfile),sep=',',header=T))
    useddata <- merge(tempdata, targetdata, by=c('eid'))
    useddata$time <- useddata$BL2Target_yrs-(useddata$age2-useddata$age0)
    useddata <- useddata[useddata$time<0 | useddata$target_y==0,] # prevalent
    # logistic analysis
    cleaned_data <- useddata[,c(pheno,'target_y',used_covnames)]
    colnames(cleaned_data) <- c('pheno','status',used_covnames)
    cleaned_data <- cleaned_data[complete.cases(cleaned_data), ]
    list_numberAll <- c(list_numberAll, nrow(cleaned_data))
    list_numberTarget <- c(list_numberTarget, nrow(cleaned_data[cleaned_data$status==1,]))
    cleaned_data$pheno <- unlist(cleaned_data$pheno)
    cleaned_data$pheno <- scale(cleaned_data$pheno)

    tryCatch({
      cox_fit <- speedglm(as.formula(paste0('status~pheno+',paste(used_covnames,collapse='+'))), data=cleaned_data, family = binomial())
      result <- summary(cox_fit)$coefficients
      # Check for the possibility of infinite coefficients
      if (any(is.infinite(coefficients(cox_fit)))) {
        print("Warning: Model coefficients may be infinite.")
      }
      list_pvalues <- c(list_pvalues, as.numeric(as.character(result[2,4])))
      list_zvalues <- c(list_zvalues, result[2,3])
      list_coef <- c(list_coef, result[2,1])
      list_secoef <- c(list_secoef, result[2,2])
      pvalues <- c(pvalues, result[2,4])
    }, error = function(e) {
      list_pvalues <<- c(list_pvalues, NA)
      list_zvalues <<- c(list_zvalues, NA)
      list_coef <<- c(list_coef, NA)
      list_secoef <<- c(list_secoef, NA)
      pvalues <<- c(pvalues, NA)
    }, warning = function(w){
          list_pvalues <<- c(list_pvalues, NA)
          list_zvalues <<- c(list_zvalues, NA)
          list_coef <<- c(list_coef, NA)
          list_secoef <<- c(list_secoef, NA)
    })
  }
}
resultframe <- data.frame(pheno=list_phenofile, target=list_targetfile, numberAll=list_numberAll, numberTarget=list_numberTarget, coef=list_coef, secoef=list_secoef, zvalue=list_zvalues, pvalue=list_pvalues)
resultframe$target <- sapply(resultframe$target, function(x) gsub('.csv','',x))
resultframe$Padj_FDR_overall <- p.adjust(resultframe$pvalue,method='BH')
write.table(resultframe, 'result/PhenoAsso/Disease/allSubjs_logistic.csv',sep=',',row.names = F)