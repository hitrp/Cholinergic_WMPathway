# survial analysis
library(survival)
library(data.table)
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
# cox harzard model, in all subjs, with more covariates
list_phenofile <- c()
list_targetfile <- c()
list_numberAll <- c()
list_numberTarget <- c()
list_pvalues <- c()
list_zvalues <- c()
list_coef <- c()
list_expcoef <- c()
list_secoef <- c()
list_FDR_separate <- c()
list_up <- c()
list_down <- c()
for(pheno in unique(prefix)){
  phenodata <- allphenos[,c('eid',pheno)]
  pvalues <- c()
  for(j in seq(1,length(targetfiles))){
    targetfile <- targetfiles[j]
    list_phenofile <- c(list_phenofile, pheno)
    list_targetfile <- c(list_targetfile, targetfile)
    # convert variables
    tempdata <- merge(phenodata, covadata, by=c('eid'))
    targetdata <- read.csv(paste0('result/PhenoAsso/Disease/Target_data_used/',targetfile),sep=',',header=T)
    useddata <- merge(tempdata, targetdata, by=c('eid'))
    useddata$time <- useddata$BL2Target_yrs-(useddata$age2-useddata$age0)
    useddata <- useddata[!is.na(useddata$time) & useddata$time>0,] # Incident
    # survival analysis
    cleaned_data <- useddata[,c(pheno,'status','time',used_covnames)]
    cleaned_data <- cleaned_data[complete.cases(cleaned_data),]
    colnames(cleaned_data)[1] <- "pheno"
    list_numberAll <- c(list_numberAll, nrow(cleaned_data))
    list_numberTarget <- c(list_numberTarget, nrow(cleaned_data[cleaned_data$status==1,]))
    cleaned_data$pheno <- unlist(cleaned_data$pheno)
    cleaned_data$pheno <- scale(cleaned_data$pheno)

    if(sum(cleaned_data$status)<30){
      list_pvalues <- c(list_pvalues, NA)
      list_zvalues <- c(list_zvalues, NA)
      list_coef <- c(list_coef, NA)
      list_expcoef <- c(list_expcoef, NA)
      list_secoef <- c(list_secoef, NA)
      list_up <- c(list_up, NA)
      list_down <- c(list_down, NA)
      pvalues <- c(pvalues, NA)
      next
    }
    cox_fit <- coxph(as.formula(paste(c("Surv(time, status) ~ pheno",used_covnames),collapse='+')), data=cleaned_data)
    result <- summary(cox_fit)$coefficients
    CI_range <- summary(cox_fit)$conf.int[1,c(3,4)]
    # Check for the possibility of infinite coefficients
    if (any(is.infinite(coefficients(cox_fit)))) {
      print("Warning: Model coefficients may be infinite.")
    }
    list_pvalues <- c(list_pvalues, result[1,5])
    list_zvalues <- c(list_zvalues, result[1,4])
    list_coef <- c(list_coef, result[1,1])
    list_expcoef <- c(list_expcoef, result[1,2])
    list_secoef <- c(list_secoef, result[1,3])
    list_up <- c(list_up, CI_range[2])
    list_down <- c(list_down, CI_range[1])
    pvalues <- c(pvalues, result[1,5])
  }
}
resultframe <- data.frame(pheno=list_phenofile, target=list_targetfile, numberAll=list_numberAll, numberTarget=list_numberTarget, expcoef=list_expcoef, coef=list_coef, down=unname(list_down), up=unname(list_up), secoef=list_secoef, zvalue=list_zvalues, pvalue=list_pvalues)
resultframe$target <- sapply(resultframe$target, function(x) gsub('.csv','',x))
resultframe$Padj_FDR_overall <- p.adjust(resultframe$pvalue,method='BH')
write.table(resultframe, 'result/PhenoAsso/Disease/allSubjs_cox.csv',sep=',',row.names = F)