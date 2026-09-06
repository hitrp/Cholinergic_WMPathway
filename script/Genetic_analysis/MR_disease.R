Args<-commandArgs(TRUE)
metricname=Args[1]
MRversion=Args[2]
threshold=as.numeric(Args[3])
modality=Args[4]

library('TwoSampleMR')
library('ieugwasr')

metricPath <- 'result/plotting/'
matching_files <- list.files(metricPath, full.names = FALSE, recursive = FALSE)
metricfile <- matching_files[grepl(paste0(metricname,'_'), matching_files)]
diseasePath <- 'data/GWAS_summary/'
diseasePath_2 <- 'resources/GWAS/'
diseases <- c('anxiety','adhd','BP','MDD','ASD','SCZ','PD','AD','MS','ALS','AS','AIS','LAS','SVS','CES','FTD','LBD')
diseaseGWASs <- c('anxiety.meta.full.fs.tbl','daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta','daner_PGC_BIP32b_mds7a_0416a','daner_pgc_mdd_meta_w2_no23andMe_rmUKBB','iPSYCHPGC_ASD_Nov2017','PGC3_SCZ_wave3_clear.european.autosome.public.v3.vcf.tsv','finngen_R9_G6_PARKINSON_clear','finngen_R9_G6_ALZHEIMER_clear','finngen_R9_G6_MS_clear','Neurology/ALS/van_Rheenen2021/ALS_sumstats_EUR_only2021.txt','Neurology/Stroke/Malik2018/MEGASTROKE.1.AS.EUR.out','Neurology/Stroke/Malik2018/MEGASTROKE.2.AIS.EUR.out','Neurology/Stroke/Malik2018/MEGASTROKE.3.LAS.EUR.out','Neurology/Stroke/Malik2018/MEGASTROKE.5.SVS.EUR.out','Neurology/Stroke/Malik2018/MEGASTROKE.4.CES.EUR.out','SUM-STATS-Discovery-full_FTD_clear.txt',"Neurology/LBD/Chia2021/GCST90001390_buildGRCh38.tsv")
if(length(metricfile)==1){
outputpath <- paste0('result/metric_disease_',modality,'/',MRversion,'/')
dir.create(outputpath,recursive = TRUE)
exposurepath <- metricPath
outcomepath <- diseasePath
bfile='resources/GWAS/LD_ref_1kgenome/EUR'
outs <- diseaseGWASs
exposure_dat_0 <- read_exposure_data(
    filename = paste0(exposurepath, metricfile),
    sep = "\t",
    snp_col = "ID",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "OMITTED",
    eaf_col = "A1_FREQ",
    pval_col = "P"
)
# https://mrcieu.github.io/ieugwasr/articles/local_ld.html
exposure_dat_0$rsid <- exposure_dat_0$SNP
exposure_dat_0$trait_id <- rep(metricname,nrow(exposure_dat_0))
exposure_dat_0$exposure <- rep(metricname,nrow(exposure_dat_0))
exposure_dat_0$pval <- exposure_dat_0$pval.exposure
exposure_clumped <- ld_clump(
    exposure_dat_0,
    plink_bin = 'resources/GWAS/LD_ref_1kgenome/plink',
    bfile = bfile,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p = threshold,
    pop = "EUR"
)

for(i in seq(1,length(outs))){
    outcomefilename <- outs[i]
    outcome <- diseases[i]
   
    if(outcome %in% c('ALS','AS','AIS','LAS','SVS','CES','LBD')){
        outcomefile <- paste0(diseasePath_2, outcomefilename)
    }else{
        outcomefile <- paste0(outcomepath, outcomefilename)
    }

    snps=exposure_clumped$SNP
    if(outcome %in% c('PD','AD','MS')){
	    outcome_dat_0 <- read_outcome_data(
	        snps = snps,
	        filename = outcomefile,
	        sep = "\t",
	        snp_col = "rsids",
	        beta_col = "beta",
	        se_col = "sebeta",
	        effect_allele_col = "alt",
	        other_allele_col = "ref",
	        eaf_col = "af_alt",
	        pval_col = "pval"
	    )
    }else if(outcome %in% c('ALS')){
        outcome_dat_0 <- read_outcome_data(
            snps = snps,
            filename = outcomefile,
            sep = "\t",
            snp_col = "rsid",
            beta_col = "beta",
            se_col = "standard_error",
            effect_allele_col = "effect_allele",
            other_allele_col = "other_allele",
            eaf_col = "effect_allele_frequency",
            pval_col = "p_value"
        )
    }else if(outcome %in% c('AS','AIS','LAS','SVS','CES')){
        outcome_dat_0 <- read_outcome_data(
            snps = snps,
            filename = outcomefile,
            sep = " ",
            snp_col = "MarkerName",
            beta_col = "Effect",
            se_col = "StdErr",
            effect_allele_col = "Allele1",
            other_allele_col = "Allele2",
            eaf_col = "Freq1",
            pval_col = "P-value"
        )
    }else if(outcome %in% c('LBD')){
        outcome_dat_0 <- read_outcome_data(
            snps = snps,
            filename = outcomefile,
            sep = "\t",
            snp_col = "variant_id",
            beta_col = "beta",
            se_col = "standard_error",
            effect_allele_col = "effect_allele",
            other_allele_col = "other_allele",
            eaf_col = "effect_allele_frequency",
            pval_col = "p_value"
        )
    }else if(outcome %in% c('FTD')){
        outcome_dat_0 <- read_outcome_data(
            snps = snps,
            filename = outcomefile,
            sep = "\t",
            snp_col = "rsid",
            beta_col = "OR",
            se_col = "SE",
            effect_allele_col = "A1",
            other_allele_col = "other_allele",
            pval_col = "PVALUE"
        )        
        # the corresponding document showed that the SE are all already for log(OR)
        outcome_dat_0$beta.outcome <- log(outcome_dat_0$beta.outcome)
    }else if(outcome %in% c('anxiety')){
    	outcome_dat_0 <- read_outcome_data(
	        snps = snps,
	        filename = outcomefile,
	        sep = "\t",
	        snp_col = "SNPID",
	        beta_col = "Effect",
	        se_col = "StdErr",
	        effect_allele_col = "Allele1",
	        other_allele_col = "Allele2",
	        eaf_col = "Freq1",
	        pval_col = "P.value"
	    )
    }else if(outcome %in% c('ASD','adhd','BP','MDD')){
    	outcome_dat_0 <- read_outcome_data(
	        snps = snps,
	        filename = outcomefile,
	        sep = "\t",
	        snp_col = "SNP",
	        beta_col = "OR",
	        se_col = "SE",
	        effect_allele_col = "A1",
	        other_allele_col = "A2",
	        pval_col = "P"
	    )
        # the corresponding document showed that the SE are all already for log(OR)
	    outcome_dat_0$beta.outcome <- log(outcome_dat_0$beta.outcome)
    }else if(outcome %in% c('SCZ')){
    	outcome_dat_0 <- read_outcome_data(
	        snps = snps,
	        filename = outcomefile,
	        sep = "\t",
	        snp_col = "ID",
	        beta_col = "BETA",
	        se_col = "SE",
	        effect_allele_col = "A1",
	        other_allele_col = "A2",
	        eaf_col = "FCAS",
	        pval_col = "PVAL"
	    )
    }
    outcome_dat_0$outcome <- rep(outcome, nrow(outcome_dat_0))
    output1<-paste0(outputpath, outcome,'_',metricname,'_QC_IVW') #output data IVW
    output2<-paste0(outputpath, outcome,'_',metricname,'_QC_MR-egger') #output data MR-egger
    output3<-paste0(outputpath, outcome,'_',metricname,'_QC_outlierSNPs') #output data outlier snp
    library(RadialMR)
    dat <- harmonise_data(exposure_dat = exposure_clumped, outcome_dat = outcome_dat_0, action=2)
    dat <- dat[dat$mr_keep %in% TRUE,]
    data3<-format_radial(dat$beta.exposure,dat$beta.outcome,dat$se.exposure,dat$se.outcome,dat$SNP)
    res1<-ivw_radial(data3,0.05,1,0.0001)
    res2<-egger_radial(data3,0.05,1)
    write.table(as.data.frame(res1$data),output1,quote=F,row.names=F,col.names=T,sep="\t")
    write.table(as.data.frame(res2$data),output2,quote=F,row.names=F,col.names=T,sep="\t")
    if (class(res1$outliers)=='data.frame'){
        a1<-res1$outliers$SNP
    }else {
        a1<-c('None')
    }
    if (class(res2$outliers)=='data.frame'){
        a2<-res2$outliers$SNP
    }else{
        a2<-c('None')
    }
    a3<-unique(c(a1,a2))
    write.table(data.frame(a3),output3,quote=F,row.names=F,col.names=F,sep="\t")


    ######################################### 3. MR analysis #########################################
    output1 <- paste0(outputpath, outcome,'_',metricname,'_MR_WM') #weight median
    output2 <- paste0(outputpath, outcome,'_',metricname,'_MR_IVW') #ivw
    output3 <- paste0(outputpath, outcome,'_',metricname,'_MR_egger') #egger slope
    output4 <- paste0(outputpath, outcome,'_',metricname,'_MR_W') #weight mode
    output5 <- paste0(outputpath, outcome,'_',metricname,'_MR_robustAdjProfile') #Robust adjusted profile score
    output6 <- paste0(outputpath, outcome,'_',metricname,'_MR_waldRatio') #Wald Ratio

    library(TwoSampleMR)
    # remove outlier snps after quality control
    dat <- dat[!(dat$SNP %in% a3),]
    tsmr1<-mr(dat, method_list=c("mr_weighted_median"))
    tsmr2<-mr(dat, method_list=c("mr_ivw"))
    tsmr5<-mr(dat, method_list=c("mr_raps"))
    tsmr6<-mr(dat, method_list=c("mr_wald_ratio"))
    tsmr3<-mr(dat, method_list=c("mr_egger_regression"))
    tsmr4<-mr(dat, method_list=c("mr_weighted_mode"))
    a1<-cbind(outcome, metricname,tsmr1$nsnp,tsmr1$b,(tsmr1$b-1.96*tsmr1$se),(tsmr1$b+1.96*tsmr1$se),tsmr1$pval) #weight median
    a2<-cbind(outcome, metricname,tsmr2$nsnp,tsmr2$b,(tsmr2$b-1.96*tsmr2$se),(tsmr2$b+1.96*tsmr2$se),tsmr2$pval) #ivw
    a5<-cbind(outcome, metricname,tsmr5$nsnp,tsmr5$b,(tsmr5$b-1.96*tsmr5$se),(tsmr5$b+1.96*tsmr5$se),tsmr5$pval) #raps
    a6<-cbind(outcome, metricname,tsmr6$nsnp,tsmr6$b,(tsmr6$b-1.96*tsmr6$se),(tsmr6$b+1.96*tsmr6$se),tsmr6$pval) #WR
    CI <- 0.95
    lowerCI <- function(beta,df,SE){
    return(beta - (qt((1-CI)/2, df, lower.tail = FALSE) * SE))
    }
    upperCI <- function(beta,df,SE){
    return(beta + (qt((1-CI)/2, df, lower.tail = FALSE) * SE))
    }
    a3 <- cbind(outcome, metricname,tsmr3$nsnp,tsmr3$b, mapply(lowerCI, tsmr3$b, tsmr3$nsnp - 2, tsmr3$se), mapply(upperCI, tsmr3$b, tsmr3$nsnp - 2, tsmr3$se), tsmr3$pval) #egger slope
    a4 <- cbind(outcome, metricname,tsmr4$nsnp,tsmr4$b, mapply(lowerCI, tsmr4$b, tsmr4$nsnp - 1, tsmr4$se), mapply(upperCI, tsmr4$b, tsmr4$nsnp - 1, tsmr4$se), tsmr4$pval) #weight mode
    write.table(a1,output1,quote=F,row.names=F,col.names=F,sep="\t") #weight median
    write.table(a2,output2,quote=F,row.names=F,col.names=F,sep="\t") #ivw
    write.table(a5,output5,quote=F,row.names=F,col.names=F,sep="\t") #raps
    write.table(a6,output6,quote=F,row.names=F,col.names=F,sep="\t") #wr
    write.table(a3,output3,quote=F,row.names=F,col.names=F,sep="\t") #egger slope
    write.table(a4,output4,quote=F,row.names=F,col.names=F,sep="\t") #weight mode
    datoutput<-paste0(outputpath, outcome,'_',metricname,'_dat')
    write.table(dat, datoutput, quote=F, row.names=F,col.names=T,sep="\t")
}
}