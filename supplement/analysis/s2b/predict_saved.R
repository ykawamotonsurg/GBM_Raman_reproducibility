# Reuse fixed preprocessing and the saved native XGBoost model. No fitting.
args<-commandArgs(TRUE)
if(length(args)!=2)stop('Usage: Rscript predict_saved.R REPOSITORY_ROOT NEW_OUTPUT_DIRECTORY')
root<-normalizePath(args[1],winslash='/',mustWork=TRUE)
src<-file.path(root,'supplement/analysis/s2b');out<-args[2]
if(dir.exists(out))stop('Choose a new output directory')
stopifnot(requireNamespace('xgboost',quietly=TRUE),requireNamespace('jsonlite',quietly=TRUE))
d<-readRDS(file.path(src,'input.rds'));p<-readRDS(file.path(src,'preprocessing.rds'))
late<-d$metadata$partition=='late';x<-d$spectra[late,,drop=FALSE]
for(n in names(p$imputation_means))x[is.na(x[,n]),n]<-p$imputation_means[[n]]
x<-x[,setdiff(colnames(x),p$zero_variance_removed),drop=FALSE]
x[,names(p$normalize_means)]<-sweep(x[,names(p$normalize_means),drop=FALSE],2,p$normalize_means,'-')
x[,names(p$normalize_sds)]<-sweep(x[,names(p$normalize_sds),drop=FALSE],2,p$normalize_sds,'/')
x<-x[,p$pca_columns,drop=FALSE]
if(is.numeric(p$pca_center))x<-sweep(x,2,p$pca_center,'-')
if(is.numeric(p$pca_scale))x<-sweep(x,2,p$pca_scale,'/')
pcs<-x%*%p$pca_rotation
model<-xgboost::xgb.load(file.path(src,'model.json'))
probability<-as.numeric(predict(model,pcs));expected<-d$prediction[late]
stopifnot(length(probability)==length(expected),identical(probability,expected))
meta<-d$metadata[late,,drop=FALSE]
bag<-paste(meta$patient_id,meta$region,sep='_')
fixed<-read.csv(file.path(root,'supplement/data/s2b_temporal_late_bags.csv'),check.names=FALSE)
threshold<-read.csv(file.path(root,'supplement/data/s2b_temporal_threshold.csv'),check.names=FALSE)$threshold
comparison<-lapply(seq_len(nrow(fixed)),function(i){
 z<-probability[bag==fixed$bag_id[i]];score<-mean(sort(z,decreasing=TRUE)[seq_len(min(3,length(z)))])
 classification<-ifelse(score>=threshold,'tumor','non_tumor')
 stopifnot(length(z)==fixed$n_spectra[i],abs(score-fixed$bag_score[i])<=1e-14,
           classification==fixed$predicted_class[i])
 data.frame(bag_id=fixed$bag_id[i],bag_score=score,stored_bag_score=fixed$bag_score[i],classification=classification)
})
dir.create(out,recursive=TRUE)
write.csv(data.frame(spectrum_id=meta$spectrum_id,probability=probability),file.path(out,'late_predictions.csv'),row.names=FALSE)
write.csv(do.call(rbind,comparison),file.path(out,'bag_comparison.csv'),row.names=FALSE)
jsonlite::write_json(list(spectra=length(probability),bags=nrow(fixed),
 exact_saved_prediction_match=identical(probability,expected),maximum_spectrum_difference=max(abs(probability-expected)),
 bag_and_classification_match=TRUE,model_fitted=FALSE,threshold_selected=FALSE),
 file.path(out,'prediction_check.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
cat('PASS: saved model inference matches 502 saved probabilities and 22 bag classifications. No fit.\n')
