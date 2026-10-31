# Statistical methods using fixed, de-identified analysis inputs.
# Adapted from the saved analysis script. The model formulas, estimation methods,
# interval rules, bootstrap counts and seeds below are retained.
# Separate from reproduce.py. Optional third argument: mri (MRI only) or all (default).
options(warn=1,digits=17)
suppressPackageStartupMessages({library(readr);library(dplyr)})
args<-commandArgs(trailingOnly=TRUE)
if(!length(args)%in%c(2,3))stop('Usage: Rscript fixed_input_statistics.R REPOSITORY_ROOT NEW_OUTPUT_DIRECTORY [mri|all]')
execution_scope<-if(length(args)==3)args[3] else 'all'
if(!execution_scope%in%c('mri','all'))stop('Scope must be mri or all')
if(execution_scope=='all')suppressPackageStartupMessages(library(lmerTest))
B<-normalizePath(args[1],winslash='/',mustWork=TRUE);O<-args[2]
if(dir.exists(O))stop('Choose a new output directory')
for(d in c('statistics','records','rds'))dir.create(file.path(O,d),recursive=TRUE,showWarnings=FALSE)
rd<-function(p)as.data.frame(read_csv(file.path(B,p),show_col_types=FALSE))
wr<-function(x,p)write_csv(x,file.path(O,p),na='')
mri<-rd('supplement/data/MRI_CaseData.csv')
stopifnot(nrow(mri)==18,!anyDuplicated(mri$patient_id),all(is.finite(mri$enhance_fraction)),all(mri$n_cells>0),all(mri$n_high>=0 & mri$n_high<=mri$n_cells),all(abs(mri$fraction-mri$n_high/mri$n_cells)<1e-12))
# The same confirmed rater volumes feed agreement and association analyses.
volume_cols<-c('rater1_total_volume_mm3','rater1_necrotic_volume_mm3','rater2_total_volume_mm3','rater2_necrotic_volume_mm3')
stopifnot(all(volume_cols%in%names(mri)),all(is.finite(as.matrix(mri[,volume_cols]))))
for(r in 1:2){
 total<-mri[[paste0('rater',r,'_total_volume_mm3')]];nec<-mri[[paste0('rater',r,'_necrotic_volume_mm3')]]
 stopifnot(all(total>0),all(nec>=0 & nec<total))
}
derived<-(mri$rater1_total_volume_mm3-mri$rater1_necrotic_volume_mm3)/mri$rater1_total_volume_mm3
stopifnot(all(abs(derived-mri$enhance_fraction)<1e-12))
mri$enhance_fraction<-derived
if(execution_scope=='all'){
case<-rd('main/data/Singlecell_CaseRegion.csv')
bridge<-rd('supplement/data/tissue_cell_bridge.csv')
stopifnot(nrow(case)==38,nrow(bridge)==37)
trend<-list();adj<-list();corr<-list()
for(scope in c('all38','matched37')){
 z<-if(scope=='all38')case else bridge;z$patient_id<-factor(z$patient_id)
 f<-lmer(frac_caro_high~region_ord+(1|patient_id),data=z,REML=FALSE)
 co<-coef(summary(f,ddf='Satterthwaite'))['region_ord',];ci<-co['Estimate']+c(-1,1)*qt(.975,co['df'])*co['Std. Error']
 trend[[scope]]<-data.frame(population=scope,n=nrow(z),n_patients=nlevels(z$patient_id),beta_fraction=co['Estimate'],beta_percentage_points=100*co['Estimate'],SE=co['Std. Error'],df=co['df'],CI_low=ci[1],CI_high=ci[2],P=co['Pr(>|t|)'],REML=FALSE,singular=isSingular(f),response='fraction 0-1')
 saveRDS(f,file.path(O,'rds',paste0('S3_trend_',scope,'.rds')))
}
wr(bind_rows(trend),'statistics/S3_trend_ML_Satterthwaite.csv')
for(scope in c('all37','GBM33','CWM','ED','IM','CE')){
 z<-if(scope=='all37')bridge else if(scope=='GBM33')subset(bridge,region!='CWM')else subset(bridge,region==scope)
 variable<-length(unique(z$frac_caro_high))>1
 corr[[scope]]<-data.frame(scope=scope,n=nrow(z),n_patients=length(unique(z$patient_id)),rho=if(variable)cor(z$frac_caro_high,z$bag_logit,method='spearman')else NA_real_,status=if(variable)'descriptive'else'all fractions identical')
 if(scope%in%c('all37','GBM33')){
 z$region<-factor(z$region,levels=if(scope=='all37')c('CWM','ED','IM','CE')else c('ED','IM','CE'));z$patient_id<-factor(z$patient_id)
 f<-lmer(bag_logit~fraction_per_10pp+region+(1|patient_id),data=z,REML=TRUE)
 co<-coef(summary(f,ddf='Kenward-Roger'))['fraction_per_10pp',];ci<-co['Estimate']+c(-1,1)*qt(.975,co['df'])*co['Std. Error']
 stopifnot(abs(co['Std. Error']^2-as.matrix(pbkrtest::vcovAdj(f))['fraction_per_10pp','fraction_per_10pp'])<1e-9)
 adj[[scope]]<-data.frame(population=scope,n=nrow(z),n_patients=nlevels(z$patient_id),beta=co['Estimate'],SE_KR=co['Std. Error'],df_KR=co['df'],CI_low=ci[1],CI_high=ci[2],P=co['Pr(>|t|)'],singular=isSingular(f),REML=TRUE)
 saveRDS(f,file.path(O,'rds',paste0('bridge_adjusted_',scope,'.rds')))
 wr(data.frame(scope=scope,cell_pair=paste(z$patient_id,z$region),fitted=fitted(f),residual=resid(f)) ,paste0('records/bridge_residuals_',scope,'.csv'))
 }
}
wr(bind_rows(corr),'statistics/bridge_descriptive_correlations.csv');wr(bind_rows(adj),'statistics/bridge_adjusted_REML_KR.csv')
}
q6<-function(x)as.numeric(quantile(x,c(.025,.5,.975),type=6,na.rm=TRUE))
sp<-cor.test(mri$enhance_fraction,mri$fraction,method='spearman',exact=FALSE)
set.seed(20260315);br<-replicate(5000,{i<-sample.int(18,18,replace=TRUE);suppressWarnings(cor(mri$enhance_fraction[i],mri$fraction[i],method='spearman'))});ci<-q6(br)
wr(data.frame(n=18,rho=unname(sp$estimate),P=sp$p.value,CI_low=ci[1],bootstrap_median=ci[2],CI_high=ci[3],B=5000,seed=20260315,quantile_type=6,finite_draws=sum(is.finite(br))),'statistics/MRI_Spearman.csv')
f<-glm(cbind(n_high,n_cells-n_high)~enhance_fraction,data=mri,family=quasibinomial('logit'));co<-coef(summary(f))['enhance_fraction',];be<-co[1];se<-co[2]
wr(data.frame(n=18,OR_per_10pp=exp(.1*be),CI_low=exp(.1*(be-qnorm(.975)*se)),CI_high=exp(.1*(be+qnorm(.975)*se)),P=co[4],dispersion=summary(f)$dispersion),'statistics/MRI_quasibinomial.csv')
set.seed(20260316);bo<-replicate(5000,{i<-sample.int(18,18,replace=TRUE);f0<-try(glm(cbind(n_high,n_cells-n_high)~enhance_fraction,data=mri[i,],family=quasibinomial('logit')),silent=TRUE);if(inherits(f0,'try-error'))NA_real_ else exp(.1*coef(f0)[2])});oc<-q6(bo[is.finite(bo)])
wr(data.frame(B=5000,finite_draws=sum(is.finite(bo)),seed=20260316,OR_median=oc[2],CI_low=oc[1],CI_high=oc[3],quantile_type=6),'statistics/MRI_OR_case_bootstrap.csv')
loo<-bind_rows(lapply(seq_len(18),function(i){z<-mri[-i,];tt<-cor.test(z$enhance_fraction,z$fraction,method='spearman',exact=FALSE);data.frame(omitted_patient_id=mri$patient_id[i],rho=unname(tt$estimate),P=tt$p.value)}))
wr(loo,'statistics/MRI_all18_leave_one_out.csv');wr(data.frame(rho_min=min(loo$rho),rho_max=max(loo$rho),P_min=min(loo$P),P_max=max(loo$P),all_positive=all(loo$rho>0),all_P_below_05=all(loo$P<.05)),'statistics/MRI_leave_one_out_summary.csv')
wr(data.frame(draw=1:5000,Spearman_rho=br,OR_per_10pp=bo),'statistics/MRI_bootstrap_draws.csv')

saveRDS(f,file.path(O,'rds','MRI_quasibinomial_model.rds'))
writeLines(capture.output(sessionInfo()),file.path(O,'records','sessionInfo.txt'))
writeLines(capture.output(RNGkind()),file.path(O,'records','RNGkind.txt'))
writeLines(c(paste('scope:',execution_scope),'Spearman bootstrap: 5000; seed 20260315; quantile type 6','OR bootstrap: 5000; seed 20260316; quantile type 6'),file.path(O,'records','execution_settings.txt'))
