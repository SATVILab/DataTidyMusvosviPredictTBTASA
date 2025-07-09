rm(list=ls())
graphics.off()
library(rio)
library(pROC)
library(caret)
library(cutpointr)
library(haven)
library(DescTools)

'%!in%' = function(x,y)!('%in%'(x,y))
setwd("")#Location where the predict_tb.csv, poc_cd3_data_summary.csv, cidri_cd3_hiv_neg_data_summary.csv, hla_validation_cd3_data_summary.csv, venous_vs_fingerprick_df.csv, and Raw-Data-Harvard-Mice.csv are located
predict_master_df = import("predict_tb.csv")
POC_CD3_summary = import("poc_cd3_data_summary.csv")
CIDRI_CD3_summary_HIV_neg = import("cidri_cd3_hiv_neg_data_summary.csv")
HLA_validation_summary = import("hla_validation_cd3_data_summary.csv")

venous_vs_fingerprick_df = import("venous_vs_fingerprick_df.csv")

Mice.Normalised.Ct.Value = read.csv("Raw-Data-Harvard-Mice.csv")

#Figure 2c - Please note that the figure in the manuscript was generated using GraphPad Prism to allow for a split y-axis

par(mar=c(8,8,7,2))
par(mfrow=c(2,1))
boxplot(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV1886c"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV1886c"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq,

        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV1813c"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV1813c"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq,

        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3875"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3875"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq,

        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3619c"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3619c"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq,

        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3874"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3874"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq,

        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV0288"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV0288"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq,

        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3620c"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3620c"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq,

        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3615c"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3615c"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq,

        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
        subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq,
        ylim=c(0,0.75),ylab="Freq of IFNγ+TNF+ T cells (%)",las=2,names = c("Ag85B-HD","Ag85B-TB",
                                       "RV1813c-HD","RV1813c-TB",
                                       "ESAT-6-HD","ESAT-6-TB",
                                       "EsxV-HD","EsxV-TB",
                                       "CFP-10-HD","CFP-10-TB",
                                       "TB10.4-HD","TB10.4-TB",
                                       "EsxW-HD","EsxW-TB",
                                       "EspC-HD","EspC-TB",
                                       "MTB.LYS-HD","MTB.LYS-TB"),range = 0)

ag85B = wilcox.test(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV1886c"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV1886c"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq)

mtext(side = 3,at = 1.5, round(ag85B$p.value,digits = 4),cex = 0.7)


RV1813c = wilcox.test(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV1813c"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV1813c"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq)

mtext(side = 3,at = 3.5, round(RV1813c$p.value,digits = 4),cex = 0.7)

esat6 = wilcox.test(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3875"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3875"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq)
mtext(side = 3,at = 5.5, round(esat6$p.value,digits = 4),cex = 0.7)

esxv = wilcox.test(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3619c"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3619c"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq)
mtext(side = 3,at = 7.5, round(esxv$p.value,digits = 4),cex = 0.7)

cfp10 = wilcox.test(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3874"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3874"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq)
mtext(side = 3,at = 9.5, round(cfp10$p.value,digits = 4),cex = 0.7)

tb10.4 = wilcox.test(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV0288"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV0288"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq)
mtext(side = 3,at = 11.5, round(tb10.4$p.value,digits = 4),cex = 0.7)

esxw = wilcox.test(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3620c"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3620c"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq)
mtext(side = 3,at = 13.5, round(esxw$p.value,digits = 4),cex = 0.7)

espc = wilcox.test(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3615c"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3615c"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq)
mtext(side = 3,at = 15.5, round(espc$p.value,digits = 4),cex = 0.7)

mtblysate = wilcox.test(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="HD_Control")$bcksub_cd3_ifng_tnf_freq,
subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="TB")$bcksub_cd3_ifng_tnf_freq)
mtext(side = 3,at = 17.5, round(mtblysate$p.value,digits = 4),cex = 0.7)

mtext(side = 3,line = 5, "Figure 2c")

bck_freq_df_for_prism = as.data.frame(matrix(nrow = 39,ncol = 0))
tb_tasa_df_for_prism = as.data.frame(matrix(nrow = 39,ncol = 0))
hd_vs_tb0_stats = NULL

j=1
for(i in c("RV1886c","RV1813c","RV3875","RV3619c","RV3874","RV0288","RV3620c","RV3615c","MTBLYSATE")){
  temp_df_hd = subset(predict_master_df,predict_master_df$stimulation==i&predict_master_df$Group=="HD_Control")
  temp_df_tb_0 = subset(predict_master_df,predict_master_df$stimulation==i&predict_master_df$Group=="TB"&predict_master_df$timepoint==0)

  temp_df_hd_freq = as.data.frame(c(temp_df_hd$bcksub_cd3_ifng_tnf_freq,rep(NA,times=39-length(temp_df_hd$bcksub_cd3_ifng_tnf_freq))))
  colnames(temp_df_hd_freq) = paste("hd",i, sep = " ")

  temp_df_tb_d0_freq = as.data.frame(temp_df_tb_0$bcksub_cd3_ifng_tnf_freq)
  colnames(temp_df_tb_d0_freq) = paste("tb_0",i, sep = " ")

  temp_freq = cbind.data.frame(temp_df_hd_freq,temp_df_tb_d0_freq)
  bck_freq_df_for_prism = cbind.data.frame(bck_freq_df_for_prism,temp_freq)

  temp_freq_a_vs_b = wilcox.test(temp_df_hd$bcksub_cd3_ifng_tnf_freq,temp_df_tb_0$bcksub_cd3_ifng_tnf_freq,exact = F)
  temp_freq_auc = rbind.data.frame(subset(temp_df_hd,select=c("Group","bcksub_cd3_ifng_tnf_freq")),subset(temp_df_tb_0,select=c("Group","bcksub_cd3_ifng_tnf_freq")))
  temp_freq_auc = pROC::roc(temp_freq_auc$Group,temp_freq_auc$bcksub_cd3_ifng_tnf_freq)
  temp_stats = as.data.frame(i)
  colnames(temp_stats) = "stimulation"
  temp_stats$hd_vs_tb0_bck_cd3_freq_p_value = temp_freq_a_vs_b$p.value
  temp_stats$hd_vs_tb0_bck_cd3_freq_auc = temp_freq_auc$auc

  temp_df_hd_hladr = as.data.frame(c(temp_df_hd$tb_tasa,rep(NA,times=39-length(temp_df_hd$tb_tasa))))
  colnames(temp_df_hd_hladr) = paste("hd",i, sep = " ")

  temp_df_tb_d0_hladr = as.data.frame(temp_df_tb_0$tb_tasa)
  colnames(temp_df_tb_d0_hladr) = paste("tb_0",i, sep = " ")
  temp_meta = cbind.data.frame(temp_df_hd_hladr,temp_df_tb_d0_hladr)
  tb_tasa_df_for_prism = cbind.data.frame(tb_tasa_df_for_prism,temp_meta)

  temp_hladr_a_vs_b = wilcox.test(temp_df_hd$tb_tasa,temp_df_tb_0$tb_tasa,alternative = "less",exact = F)
  temp_hladr_auc = rbind.data.frame(subset(temp_df_hd,select=c("Group","tb_tasa")),subset(temp_df_tb_0,select=c("Group","tb_tasa")))
  temp_hladr_auc = pROC::roc(temp_hladr_auc$Group,temp_hladr_auc$tb_tasa)

  temp_stats$hd_vs_tb0_hladr_p_value = temp_hladr_a_vs_b$p.value
  temp_stats$hd_vs_tb0_auc = temp_hladr_auc$auc
  temp_stats$lower_95ci = ci.auc(temp_hladr_auc)[1]
  temp_stats$upper_95ci = ci.auc(temp_hladr_auc)[3]
  temp_stats$hd_total = length(unique(temp_df_hd$pid))
  temp_stats$hd_responders = sum(temp_df_hd$classification=="responder")
  temp_stats$hd_tasa_positive = sum(temp_df_hd$tb_tasa>10,na.rm = T)

  temp_stats$tb0_total = length(unique(temp_df_tb_0$pid))
  temp_stats$tb0_responders = sum(temp_df_tb_0$classification=="responder")
  temp_stats$tb0_tasa_positive = sum(temp_df_tb_0$tb_tasa>10,na.rm = T)

  hd_vs_tb0_stats = rbind.data.frame(hd_vs_tb0_stats,temp_stats)
}

barplot(c(sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV1886c"&predict_master_df$Group=="HD_Control")$tb_tasa)),
        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV1886c"&predict_master_df$Group=="TB")$tb_tasa)),

        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV1813c"&predict_master_df$Group=="HD_Control")$tb_tasa)),
        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV1813c"&predict_master_df$Group=="TB")$tb_tasa)),

        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3875"&predict_master_df$Group=="HD_Control")$tb_tasa)),
        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3875"&predict_master_df$Group=="TB")$tb_tasa)),

        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3619c"&predict_master_df$Group=="HD_Control")$tb_tasa)),
        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3619c"&predict_master_df$Group=="TB")$tb_tasa)),

        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3874"&predict_master_df$Group=="HD_Control")$tb_tasa)),
        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3874"&predict_master_df$Group=="TB")$tb_tasa)),

        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV0288"&predict_master_df$Group=="HD_Control")$tb_tasa)),
        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV0288"&predict_master_df$Group=="TB")$tb_tasa)),

        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3620c"&predict_master_df$Group=="HD_Control")$tb_tasa)),
        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3620c"&predict_master_df$Group=="TB")$tb_tasa)),

        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3615c"&predict_master_df$Group=="HD_Control")$tb_tasa)),
        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="RV3615c"&predict_master_df$Group=="TB")$tb_tasa)),

        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="HD_Control")$tb_tasa)),
        sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="TB")$tb_tasa)))
        ,ylim=c(0,42),ylab="Flow responders",las=2,names = c("Ag85B-HD","Ag85B-TB",
                                     "RV1813c-HD","RV1813c-TB",
                                     "ESAT-6-HD","ESAT-6-TB",
                                     "EsxV-HD","EsxV-TB",
                                     "CFP-10-HD","CFP-10-TB",
                                     "TB10.4-HD","TB10.4-TB",
                                     "EsxW-HD","EsxW-TB",
                                     "EspC-HD","EspC-TB",
                                     "MTB.LYS-HD","MTB.LYS-TB"))
mtext(side = 3,line = 5, "Figure 2d")

#Figure 3a
par(mfrow=c(2,3))
par(mar=c(8,8,7,2))
boxplot(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="HD_Control")$tb_tasa,subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="TB")$tb_tasa,xlab="",las=2,ylim=c(0,100),cex.axis=2,col="white",border=c("blue","red"),lwd=2,cex=0,range = 0)
mtext(side = 3,line = 5,"Figure 3a")
mtext(side = 3,"PredictTB")

points(rep(1,times=length(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="HD_Control")$tb_tasa)),subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="HD_Control")$tb_tasa,cex=2,col="blue")

points(rep(2,times=length(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="TB")$tb_tasa)),subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="TB")$tb_tasa,cex=2,col="red")

mtext(side = 1,line = 0.5,at = 1,"IGRA+")
mtext(side = 1,line = 0.5,at = 2,"TB")

mtext(side = 1,at = 0,line = 2, "# of samples acquired")
mtext(side = 1,at = 0,line = 3, "responders")
mtext(side = 1,at = 0,line = 4, "biomarker +")

temp_hd_vs_tb_hladr = wilcox.test(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="HD_Control")$tb_tasa,subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="TB")$tb_tasa,"less")
mtext(side = 3,line = 1,paste("p = ",round(temp_hd_vs_tb_hladr$p.value,digits = 6),sep = ""),cex=1)

mtext(side = 1,at = 1,line = 2, length(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="HD_Control")$tb_tasa))

mtext(side = 1,at = 1,line = 3, sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="HD_Control")$tb_tasa)))

mtext(side = 1,at = 1,line = 4, sum(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="HD_Control")$tb_tasa>=10,na.rm = T))

mtext(side = 1,at = 2,line = 2, length(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="TB")$tb_tasa))

mtext(side = 1,at = 2,line = 3, sum(!is.na(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="TB")$tb_tasa)))

mtext(side = 1,at = 2,line = 4, sum(subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Group=="TB")$tb_tasa>=10,na.rm = T))


#Determining TB-TASA threshold
baseline_df = subset(predict_master_df,predict_master_df$timepoint==0&predict_master_df$stimulation=="MTBLYSATE")

mtext(side = 3, line = 2.5,paste("accuracy threshold = ",round(cutpointr(baseline_df, tb_tasa, Group, metric = accuracy,na.rm=T)$optimal_cutpoint,digits = 2),sep = ""))

#Figure 3b
predict_d0_auc_cd3 = pROC::roc(baseline_df$Group,baseline_df$tb_tasa)
temp_auc = predict_d0_auc_cd3
temp_col = "darkgreen"
temp_df = baseline_df

threshold_sens =  sum(subset(temp_df,temp_df$Group=="TB")$tb_tasa>10,na.rm = T)/length(subset(temp_df,temp_df$Group=="TB"&!is.na(temp_df$tb_tasa))$tb_tasa)

threshold_spec = sum(subset(temp_df,temp_df$Group=="HD_Control")$tb_tasa<=10,na.rm = T)/length(subset(temp_df,temp_df$Group=="HD_Control"&!is.na(temp_df$tb_tasa))$tb_tasa)

set.seed(12345)
a_roc.ci.se.obj = ci.se(temp_auc, specificities = seq(0, 1, 0.01), boot.n=1000, parallel = F)
stats_a_roc = format(wilcox.test(temp_auc$cases, temp_auc$controls, paired = F, alternative = "great")$p.value, digits = 2)
transColour = alpha(temp_col,0.3)
plot(c(0,0),ylim=c(0,1),xlim=c(1,0),cex=0)
points(temp_auc$specificities,temp_auc$sensitivities,type="l",lwd=2,col=temp_col)
plot(a_roc.ci.se.obj, type="shape", col= transColour, no.roc = T, lty = 0)
points(threshold_spec,threshold_sens,pch=19,col=temp_col,cex=2)

segments(x0 = 1,x1 = 0.7,y0 = 0.9,y1=0.9,lwd = 2,lty=2)
segments(x0 = 1,x1 = 0.7,y0 = 1,y1=1,lwd = 2,lty=2)
segments(x0 = 1,x1 = 1,y0 = 0.9,y1=1,lwd = 2,lty=2)
segments(x0 = 0.7,x1 = 0.7,y0 = 0.9,y1=1,lwd = 2,lty=2)

mtext(side = 3,"Figure 3b")
text(0.2,0.15,paste("AUC = ",round(temp_auc$auc,digits = 2), "95%CI (",round(ci.auc(temp_auc)[1],digits = 2), "-",round(ci.auc(temp_auc)[3],digits = 2),sep = " "))
text(0.2,0.05,paste("sens = ",round(threshold_sens*100,digits = 2),"%",sep = ""))
text(0.2,0.0,paste("spec = ",round(threshold_spec*100,digits = 4),"%",sep = ""))


boxplot(subset(POC_CD3_summary,POC_CD3_summary$SubjectType=="LTBI")$tb_tasa,subset(POC_CD3_summary,POC_CD3_summary$SubjectType=="TB Case")$tb_tasa,xlab="",las=2,ylim=c(0,100),cex.axis=2,col="white",border=c("blue","red"),lwd=2,cex=0,range = 0)


points(rep(1,times=length(subset(POC_CD3_summary,POC_CD3_summary$SubjectType=="LTBI")$tb_tasa)),subset(POC_CD3_summary,POC_CD3_summary$SubjectType=="LTBI")$tb_tasa,cex=2,col="blue")

points(rep(2,times=length(subset(POC_CD3_summary,POC_CD3_summary$SubjectType=="TB Case")$tb_tasa)),subset(POC_CD3_summary,POC_CD3_summary$SubjectType=="TB Case")$tb_tasa,cex=2,col="red")

mtext(side = 1,line = 0.5,at = 1,"IGRA+")
mtext(side = 1,line = 0.5,at = 2,"TB")

temp_hd_vs_tb_hladr = wilcox.test(subset(POC_CD3_summary,POC_CD3_summary$SubjectType=="LTBI")$tb_tasa,subset(POC_CD3_summary,POC_CD3_summary$SubjectType=="TB Case")$tb_tasa,"less")
mtext(side = 3,line = 1,paste("p = ",round(temp_hd_vs_tb_hladr$p.value,digits = 6),sep = ""),cex=1)

mtext(side = 1,at = 0,line = 2, "# of samples acquired")
mtext(side = 1,at = 0,line = 3, "responders")
mtext(side = 1,at = 0,line = 4, "biomarker +")

mtext(side = 1,at = 1,line = 2, length(subset(POC_CD3_summary,POC_CD3_summary$SubjectType=="LTBI")$tb_tasa))

mtext(side = 1,at = 1,line = 3, sum(subset(POC_CD3_summary,POC_CD3_summary$SubjectType=="LTBI")$classification=="responders"))

mtext(side = 1,at = 1,line = 4, sum(subset(POC_CD3_summary,POC_CD3_summary$SubjectType=="LTBI")$tb_tasa>=10,na.rm = T))

mtext(side = 1,at = 2,line = 2, length(subset(POC_CD3_summary,POC_CD3_summary$SubjectType=="TB Case")$tb_tasa))

mtext(side = 1,at = 2,line = 3, sum(subset(POC_CD3_summary,POC_CD3_summary$SubjectType=="TB Case")$classification=="responders"))

mtext(side = 1,at = 2,line = 4, sum(subset(POC_CD3_summary,POC_CD3_summary$SubjectType=="TB Case")$tb_tasa>=10,na.rm = T))

abline(h=10,lwd=2,lty=2,col="red")
mtext(side = 3,line = 0,"Musvosvi et al. 2017 (Figure 3c)")


boxplot(subset(CIDRI_CD3_summary_HIV_neg,CIDRI_CD3_summary_HIV_neg$Group=="QFT")$tb_tasa,subset(CIDRI_CD3_summary_HIV_neg,CIDRI_CD3_summary_HIV_neg$Group=="TB")$tb_tasa,xlab="",las=2,ylim=c(0,100),cex.axis=2,col="white",border=c("blue","red"),lwd=2,cex=0,range = 0)
mtext(side = 1,line = 0.5,at = 1,"QFT+")
mtext(side = 1,line = 0.5,at = 2,"TB")

points(rep(1,times=length(subset(CIDRI_CD3_summary_HIV_neg,CIDRI_CD3_summary_HIV_neg$Group=="QFT")$tb_tasa)),subset(CIDRI_CD3_summary_HIV_neg,CIDRI_CD3_summary_HIV_neg$Group=="QFT")$tb_tasa,cex=2,col="blue")

points(rep(2,times=length(subset(CIDRI_CD3_summary_HIV_neg,CIDRI_CD3_summary_HIV_neg$Group=="TB")$tb_tasa)),subset(CIDRI_CD3_summary_HIV_neg,CIDRI_CD3_summary_HIV_neg$Group=="TB")$tb_tasa,cex=2,col="red")

mtext(side = 1,at = 0,line = 2, "# of samples acquired")
mtext(side = 1,at = 0,line = 3, "responders")
mtext(side = 1,at = 0,line = 4, "biomarker +")

temp_hd_vs_tb_hladr = wilcox.test(subset(CIDRI_CD3_summary_HIV_neg,CIDRI_CD3_summary_HIV_neg$Group=="QFT")$tb_tasa,subset(CIDRI_CD3_summary_HIV_neg,CIDRI_CD3_summary_HIV_neg$Group=="TB")$tb_tasa,"less")
mtext(side = 3,line = 1,paste("p = ",round(temp_hd_vs_tb_hladr$p.value,digits = 6),sep = ""),cex=1)

mtext(side = 1,at = 1,line = 2, length(subset(CIDRI_CD3_summary_HIV_neg,CIDRI_CD3_summary_HIV_neg$Group=="QFT")$tb_tasa))

mtext(side = 1,at = 1,line = 3, sum(subset(CIDRI_CD3_summary_HIV_neg,CIDRI_CD3_summary_HIV_neg$Group=="QFT")$classification=="responders"))

mtext(side = 1,at = 1,line = 4, sum(subset(CIDRI_CD3_summary_HIV_neg,CIDRI_CD3_summary_HIV_neg$Group=="QFT")$tb_tasa>=10,na.rm = T))

mtext(side = 1,at = 2,line = 2, length(subset(CIDRI_CD3_summary_HIV_neg,CIDRI_CD3_summary_HIV_neg$Group=="TB")$tb_tasa))

mtext(side = 1,at = 2,line = 3, sum(subset(CIDRI_CD3_summary_HIV_neg,CIDRI_CD3_summary_HIV_neg$Group=="TB")$classification=="responders"))

mtext(side = 1,at = 2,line = 4, sum(subset(CIDRI_CD3_summary_HIV_neg,CIDRI_CD3_summary_HIV_neg$Group=="TB")$tb_tasa>=10,na.rm = T))

abline(h=10,lwd=2,lty=2,col="red")
mtext(side = 3,"Riou et al. 2020 (Figure 3d)")

boxplot(subset(HLA_validation_summary,HLA_validation_summary$Classification=="QFT+ control")$tb_tasa,subset(HLA_validation_summary,HLA_validation_summary$Classification=="TB")$tb_tasa,xlab="",las=2,ylim=c(0,100),cex.axis=2,col="white",border=c("blue","red"),lwd=2,cex=0,range = 0)

points(rep(1,times=length(subset(HLA_validation_summary,HLA_validation_summary$Classification=="QFT+ control")$tb_tasa)),subset(HLA_validation_summary,HLA_validation_summary$Classification=="QFT+ control")$tb_tasa,cex=2,col="blue")

points(rep(2,times=length(subset(HLA_validation_summary,HLA_validation_summary$Classification=="TB")$tb_tasa)),subset(HLA_validation_summary,HLA_validation_summary$Classification=="TB")$tb_tasa,cex=2,col="red")

mtext(side = 1,line = 0.5,at = 1,"IGRA+")
mtext(side = 1,line = 0.5,at = 2,"TB")

mtext(side = 1,at = 0,line = 2, "# of samples acquired")
mtext(side = 1,at = 0,line = 3, "responders")
mtext(side = 1,at = 0,line = 4, "biomarker +")

temp_hd_vs_tb_hladr = wilcox.test(subset(HLA_validation_summary,HLA_validation_summary$Classification=="QFT+ control")$tb_tasa,subset(HLA_validation_summary,HLA_validation_summary$Classification=="TB")$tb_tasa,"less")
mtext(side = 3,line = 1,paste("p = ",round(temp_hd_vs_tb_hladr$p.value,digits = 6),sep = ""),cex=1)

mtext(side = 1,at = 1,line = 2, length(subset(HLA_validation_summary,HLA_validation_summary$Classification=="QFT+ control")$tb_tasa))
mtext(side = 1,at = 1,line = 3, sum(subset(HLA_validation_summary,HLA_validation_summary$Classification=="QFT+ control")$classification=="responders"))
mtext(side = 1,at = 1,line = 4, sum(subset(HLA_validation_summary,HLA_validation_summary$Classification=="QFT+ control")$tb_tasa>=10,na.rm = T))


mtext(side = 1,at = 2,line = 2, length(subset(HLA_validation_summary,HLA_validation_summary$Classification=="TB")$tb_tasa))
mtext(side = 1,at = 2,line = 3, sum(subset(HLA_validation_summary,HLA_validation_summary$Classification=="TB")$classification=="responders"))
mtext(side = 1,at = 2,line = 4, sum(subset(HLA_validation_summary,HLA_validation_summary$Classification=="TB")$tb_tasa>=10,na.rm = T))

abline(h=10,lwd=2,lty=2,col="red")
mtext(side = 3,"Mpande et al. 2021 (Figure 3e)")



poc_auc_cd3 = pROC::roc(POC_CD3_summary$SubjectType,POC_CD3_summary$tb_tasa)
cidri_hiv_neg_auc_cd3 = pROC::roc(CIDRI_CD3_summary_HIV_neg$Group,CIDRI_CD3_summary_HIV_neg$tb_tasa)
hla_dr_validation_auc_cd3 = subset(HLA_validation_summary,HLA_validation_summary$Group%in%c("QFT+ control","TB"))
hla_dr_validation_auc_cd3 = pROC::roc(hla_dr_validation_auc_cd3$Group,hla_dr_validation_auc_cd3$tb_tasa)

plot(c(0,0),ylim=c(0,1),xlim=c(1,0),cex=0)
j=0.8
for(i in 1:3){
  if(i==1){
    temp_title = "Musvosvi et al. 2017"
    temp_auc = poc_auc_cd3
    temp_col = "red"
    temp_df = POC_CD3_summary

    threshold_sens =  sum(subset(temp_df,temp_df$SubjectType=="TB Case")$tb_tasa>10,na.rm = T)/length(subset(temp_df,temp_df$SubjectType=="TB Case"&!is.na(temp_df$tb_tasa))$tb_tasa)

    threshold_spec = sum(subset(temp_df,temp_df$SubjectType=="LTBI")$tb_tasa<=10,na.rm = T)/length(subset(temp_df,temp_df$SubjectType=="LTBI"&!is.na(temp_df$tb_tasa))$tb_tasa)

    set.seed(12345)
    a_roc.ci.se.obj = ci.se(temp_auc, specificities = seq(0, 1, 0.01), boot.n=1000, parallel = F)
    stats_a_roc = format(wilcox.test(temp_auc$cases, temp_auc$controls, paired = F, alternative = "great")$p.value, digits = 2)
    transColour = alpha(temp_col,0.3)
    points(temp_auc$specificities,temp_auc$sensitivities,type="l",lwd=2,col=temp_col)
    plot(a_roc.ci.se.obj, type="shape", col= transColour, no.roc = T, lty = 0)
    points(threshold_spec,threshold_sens,pch=19,col=temp_col,cex=2)

  }
  if(i==2){
    temp_title = "Riou et al. 2020"
    temp_auc = cidri_hiv_neg_auc_cd3
    temp_col = "blue"
    temp_df = CIDRI_CD3_summary_HIV_neg
    threshold_sens =  sum(subset(temp_df,temp_df$Group=="TB")$tb_tasa>10,na.rm = T)/length(subset(temp_df,temp_df$Group=="TB"&!is.na(temp_df$tb_tasa))$tb_tasa)

    threshold_spec = sum(subset(temp_df,temp_df$Group=="QFT")$tb_tasa<=10,na.rm = T)/length(subset(temp_df,temp_df$Group=="QFT"&!is.na(temp_df$tb_tasa))$tb_tasa)

    set.seed(12345)
    a_roc.ci.se.obj = ci.se(temp_auc, specificities = seq(0, 1, 0.01), boot.n=1000, parallel = F)
    stats_a_roc = format(wilcox.test(temp_auc$cases, temp_auc$controls, paired = F, alternative = "great")$p.value, digits = 2)
    transColour = alpha(temp_col,0.3)
    points(temp_auc$specificities,temp_auc$sensitivities,type="l",lwd=2,col=temp_col)
    plot(a_roc.ci.se.obj, type="shape", col= transColour, no.roc = T, lty = 0)
    points(threshold_spec,threshold_sens,pch=19,col=temp_col,cex=2)

  }
  if(i==3){
    temp_title = "Mpande et al. 2021"
    temp_auc = hla_dr_validation_auc_cd3
    temp_col = "black"
    temp_df = subset(HLA_validation_summary,HLA_validation_summary$Group%in%c("QFT+ control","TB"))

    threshold_sens =  sum(subset(temp_df,temp_df$Group=="TB")$tb_tasa>10,na.rm = T)/length(subset(temp_df,temp_df$Group=="TB"&!is.na(temp_df$tb_tasa))$tb_tasa)

    threshold_spec = sum(subset(temp_df,temp_df$Group=="QFT+ control")$tb_tasa<=10,na.rm = T)/length(subset(temp_df,temp_df$Group=="QFT+ control"&!is.na(temp_df$tb_tasa))$tb_tasa)

    set.seed(12345)
    a_roc.ci.se.obj = ci.se(temp_auc, specificities = seq(0, 1, 0.01), boot.n=1000, parallel = F)
    stats_a_roc = format(wilcox.test(temp_auc$cases, temp_auc$controls, paired = F, alternative = "great")$p.value, digits = 2)
    transColour = alpha(temp_col,0.3)
    points(temp_auc$specificities,temp_auc$sensitivities,type="l",lwd=2,col=temp_col)
    plot(a_roc.ci.se.obj, type="shape", col= transColour, no.roc = T, lty = 0)
    points(threshold_spec,threshold_sens,pch=19,col=temp_col,cex=2)


  }

  text(0.2,j,paste(temp_title,"AUC = ",round(temp_auc$auc,digits = 2), "95%CI (",round(ci.auc(temp_auc)[1],digits = 2), "-",round(ci.auc(temp_auc)[3],digits = 2),sep = " "),col=temp_col)
  j=j-0.05
  text(0.2,j,paste("sens = ",round(threshold_sens*100,digits = 2),"%",sep = ""),col=temp_col)
  j=j-0.05
  text(0.2,j,paste("spec = ",round(threshold_spec*100,digits = 4),"%",sep = ""),col=temp_col)
  j=j-0.1
  }

segments(x0 = 1,x1 = 0.7,y0 = 0.9,y1=0.9,lwd = 2,lty=2)
segments(x0 = 1,x1 = 0.7,y0 = 1,y1=1,lwd = 2,lty=2)
segments(x0 = 1,x1 = 1,y0 = 0.9,y1=1,lwd = 2,lty=2)
segments(x0 = 0.7,x1 = 0.7,y0 = 0.9,y1=1,lwd = 2,lty=2)
mtext(side = 3,"Figure 3f")

#Figure 4a
baseline_all_antigens_df = subset(predict_master_df,predict_master_df$timepoint==0)
par(mar=c(7,6,5,5))
par(mfrow=c(1,1))
boxplot(subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV1886c"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV1886c"&baseline_all_antigens_df$Group=="TB")$tb_tasa,

        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV1813c"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV1813c"&baseline_all_antigens_df$Group=="TB")$tb_tasa,

        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3875"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3875"&baseline_all_antigens_df$Group=="TB")$tb_tasa,

        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3619c"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3619c"&baseline_all_antigens_df$Group=="TB")$tb_tasa,

        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3874"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3874"&baseline_all_antigens_df$Group=="TB")$tb_tasa,

        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV0288"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV0288"&baseline_all_antigens_df$Group=="TB")$tb_tasa,

        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3620c"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3620c"&baseline_all_antigens_df$Group=="TB")$tb_tasa,

        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3615c"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
        subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3615c"&baseline_all_antigens_df$Group=="TB")$tb_tasa,
        ylim=c(0,100),ylab="TB-TASA (%)",las=2,names = c("Ag85B-HD","Ag85B-TB",
                                                                           "RV1813c-HD","RV1813c-TB",
                                                                           "ESAT-6-HD","ESAT-6-TB",
                                                                           "EsxV-HD","EsxV-TB",
                                                                           "CFP-10-HD","CFP-10-TB",
                                                                           "TB10.4-HD","TB10.4-TB",
                                                                           "EsxW-HD","EsxW-TB",
                                                                           "EspC-HD","EspC-TB"
        ),range = 0)

tb_tasa_ag85B = wilcox.test(subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV1886c"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
                            subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV1886c"&baseline_all_antigens_df$Group=="TB")$tb_tasa,alternative = "less")
mtext(side = 3,at=1.5,round(tb_tasa_ag85B$p.value,digits = 4),cex=0.8)

tb_tasa_RV1813c = wilcox.test(subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV1813c"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
                              subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV1813c"&baseline_all_antigens_df$Group=="TB")$tb_tasa,alternative = "less")
mtext(side = 3,at=3.5,round(tb_tasa_RV1813c$p.value,digits = 4),cex=0.8)

tb_tasa_esat6 = wilcox.test(subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3875"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
                            subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3875"&baseline_all_antigens_df$Group=="TB")$tb_tasa,alternative = "less")
mtext(side = 3,at=5.5,round(tb_tasa_esat6$p.value,digits = 4),cex=0.8)

tb_tasa_esxv = wilcox.test(subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3619c"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
                           subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3619c"&baseline_all_antigens_df$Group=="TB")$tb_tasa,alternative = "less")
mtext(side = 3,at=7.5,round(tb_tasa_esxv$p.value,digits = 4),cex=0.8)

tb_tasa_cfp10 = wilcox.test(subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3874"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
                            subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3874"&baseline_all_antigens_df$Group=="TB")$tb_tasa,alternative = "less")
mtext(side = 3,at=9.5,round(tb_tasa_cfp10$p.value,digits = 4),cex=0.8)

tb_tasa_tb10.4 = wilcox.test(subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV0288"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
                             subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV0288"&baseline_all_antigens_df$Group=="TB")$tb_tasa,alternative = "less")
mtext(side = 3,at=11.5,round(tb_tasa_tb10.4$p.value,digits = 4),cex=0.8)

tb_tasa_esxw = wilcox.test(subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3620c"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
                           subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3620c"&baseline_all_antigens_df$Group=="TB")$tb_tasa,alternative = "less")
mtext(side = 3,at=13.5,round(tb_tasa_esxw$p.value,digits = 4),cex=0.8)

tb_tasa_espc = wilcox.test(subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3615c"&baseline_all_antigens_df$Group=="HD_Control")$tb_tasa,
                           subset(baseline_all_antigens_df,baseline_all_antigens_df$timepoint==0&baseline_all_antigens_df$stimulation=="RV3615c"&baseline_all_antigens_df$Group=="TB")$tb_tasa,alternative = "less")
mtext(side = 3,at=15.5,round(tb_tasa_espc$p.value,digits = 4),cex=0.8)
mtext(side = 3,line = 2,"Figure 4a")

#Figure 4b
par(mfrow=c(2,4))
par(mar=c(8,8,7,2))
for(i in c("RV1886c","RV1813c","RV3875","RV3619c","RV3874","RV0288","RV3620c","RV3615c")){
  plot(c(0,0),ylim=c(0,1),xlim=c(1,0),cex=0)
  if(i == "RV1886c"){
    mtext(side = 3,line = 3,"Figure 4b")
  }
  temp_df_hd = subset(predict_master_df,predict_master_df$stimulation==i&predict_master_df$timepoint==0&predict_master_df$Group=="HD_Control")
  temp_df_tb_0 = subset(predict_master_df,predict_master_df$stimulation==i&predict_master_df$timepoint==0&predict_master_df$Group=="TB")

  temp_df = rbind.data.frame(subset(temp_df_hd,select=c("Group","tb_tasa")),subset(temp_df_tb_0,select=c("Group","tb_tasa")))
  temp_auc = pROC::roc(temp_df$Group,temp_df$tb_tasa)
  colnames(temp_df) = c("Group","tb_tasa")

  threshold_sens =  sum(subset(temp_df,temp_df$Group=="TB")$tb_tasa>10,na.rm = T)/length(subset(temp_df,temp_df$Group=="TB"&!is.na(temp_df$tb_tasa))$tb_tasa)

  threshold_spec = sum(subset(temp_df,temp_df$Group=="HD_Control")$tb_tasa<=10,na.rm = T)/length(subset(temp_df,temp_df$Group=="HD_Control"&!is.na(temp_df$tb_tasa))$tb_tasa)

  set.seed(12345)
  a_roc.ci.se.obj = ci.se(temp_auc, specificities = seq(0, 1, 0.01), boot.n=1000, parallel = F)
  stats_a_roc = format(wilcox.test(temp_auc$cases, temp_auc$controls, paired = F, alternative = "great")$p.value, digits = 2)
  temp_col ="black"
  transColour = alpha(temp_col,0.3)
  points(temp_auc$specificities,temp_auc$sensitivities,type="l",lwd=2,col=temp_col)
  plot(a_roc.ci.se.obj, type="shape", col= transColour, no.roc = T, lty = 0)
  points(threshold_spec,threshold_sens,pch=19,col=temp_col,cex=2)

  segments(x0 = 1,x1 = 0.7,y0 = 0.9,y1=0.9,lwd = 2,lty=2)
  segments(x0 = 1,x1 = 0.7,y0 = 1,y1=1,lwd = 2,lty=2)
  segments(x0 = 1,x1 = 1,y0 = 0.9,y1=1,lwd = 2,lty=2)
  segments(x0 = 0.7,x1 = 0.7,y0 = 0.9,y1=1,lwd = 2,lty=2)
  mtext(side = 3,line = 0,paste(i,round(temp_auc$auc,digits = 2), "95%CI (",round(ci.auc(temp_auc)[1],digits = 2), "-",round(ci.auc(temp_auc)[3],digits = 2),")",sep = ""))
}

#Figure 5a
par(mfrow=c(2,3))
par(mar=c(6,6,8,2))
boxplot(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="IGRA+")$tb_tasa,
subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&is.na(predict_master_df$`Arm B/C`))$tb_tasa,
subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&!is.na(predict_master_df$`Arm B/C`))$tb_tasa,las=1,names = c("IGRA+","Arm A","Arm B/C"),cex=0,ylab="TB-TASA Score (%)",range = 0)
mtext(side = 3,line = 6,"Figure 5a")
abline(h=10,lwd=2,lty=2)

points(rep(1,times=length(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="IGRA+")$tb_tasa)),subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="IGRA+")$tb_tasa,cex=1.5)

points(rep(2,times=length(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&is.na(predict_master_df$`Arm B/C`))$tb_tasa)),subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&is.na(predict_master_df$`Arm B/C`))$tb_tasa,cex=1.5)

points(rep(3,times=length(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&!is.na(predict_master_df$`Arm B/C`))$tb_tasa)),subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&!is.na(predict_master_df$`Arm B/C`))$tb_tasa,cex=1.5)

mtext(side = 1,at=1,line = 2,length(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="IGRA+")$tb_tasa))

mtext(side = 1,at=1,line = 3,sum(!is.na(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="IGRA+")$tb_tasa)))

mtext(side = 1,at=1,line = 4,sum(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="IGRA+")$tb_tasa>10))

mtext(side = 1,at=2,line = 2,length(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&is.na(predict_master_df$`Arm B/C`))$tb_tasa))

mtext(side = 1,at=2,line = 3,sum(!is.na(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&is.na(predict_master_df$`Arm B/C`))$tb_tasa)))

mtext(side = 1,at=2,line = 4,sum(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&is.na(predict_master_df$`Arm B/C`))$tb_tasa>10,na.rm = T))

mtext(side = 1,at=3,line = 2,length(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&!is.na(predict_master_df$`Arm B/C`))$tb_tasa))

mtext(side = 1,at=3,line = 3,sum(!is.na(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&!is.na(predict_master_df$`Arm B/C`))$tb_tasa)))

mtext(side = 1,at=3,line = 4,sum(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&!is.na(predict_master_df$`Arm B/C`))$tb_tasa>10,na.rm = T))

mtext(side = 3,line = 0,"Baseline TB-TASA")
mtext(side = 3,line = 1,paste("IGRA+ vs Arm A p =",round(wilcox.test(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="IGRA+")$tb_tasa,subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&is.na(predict_master_df$`Arm B/C`))$tb_tasa)$p.value,digits = 4)))
mtext(side = 3,line = 2,paste("IGRA+ vs Arm B/C p =",round(wilcox.test(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="IGRA+")$tb_tasa,subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&!is.na(predict_master_df$`Arm B/C`))$tb_tasa)$p.value,digits = 4)))
mtext(side = 3,line = 3,paste("ArmA vs Arm B/C p =",round(wilcox.test(subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&is.na(predict_master_df$`Arm B/C`))$tb_tasa,subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$timepoint==0&!is.na(predict_master_df$`Arm B/C`))$tb_tasa)$p.value,digits = 4)))


#Figure 5b
arm_a_baseline = subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$`Arm A`=="Baseline")

plot(NA,NA,ylim=c(0,100),xlim=c(0,24),las=1,xlab="",xaxt = "n", yaxt = "n",ylab="TB-TASA Score (%)")
mtext(side = 3,line = 6,"Figure 5b")
axis(1, at = c(0, 8, 16, 24))
axis(2, at = c(0, 20, 40, 60, 80, 100),las=1)
mtext(side = 3,"Arm A baseline allocation")
for(i in unique(arm_a_baseline$Donor.ID)){
 temp_df = subset(arm_a_baseline,arm_a_baseline$Donor.ID==i)
 temp_df = temp_df[order(temp_df$timepoint,decreasing = F),]
 points(temp_df$timepoint,temp_df$tb_tasa,pch=ifelse(temp_df$Outcome!="success",19,1),col=ifelse(temp_df$Outcome!="success","darkred","black"),type = "o",cex=1.5)
}

j=1
for(i in c(0,8,16,24)){
if(i!=0){
  temp_df = merge(subset(arm_a_baseline,arm_a_baseline$timepoint==0,select=c("Donor.ID","tb_tasa")),subset(arm_a_baseline,arm_a_baseline$timepoint==i,select=c("Donor.ID","tb_tasa")),by="Donor.ID")
  mtext(side = 3,line = j,paste("Week 0 vs Week ",i," p value =",round(wilcox.test(temp_df$tb_tasa.x,temp_df$tb_tasa.y,paired = T)$p.value,digits = 4),sep = ""))
  j=j+1
}
  mtext(side = 1,at=i,line = 2,length(subset(arm_a_baseline,arm_a_baseline$timepoint==i)$tb_tasa))

  mtext(side = 1,at=i,line = 3,sum(!is.na(subset(arm_a_baseline,arm_a_baseline$timepoint==i)$tb_tasa)))
}

#Figure 5c
randomised_to_arm_b = subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$Arm=="B")

plot(NA,NA,ylim=c(0,100),xlim=c(0,24),las=1,xlab="",xaxt = "n", yaxt = "n",ylab="TB-TASA Score (%)")
mtext(side = 3,line = 6,"Figure 5c")
axis(1, at = c(0, 8, 16, 24))
axis(2, at = c(0, 20, 40, 60, 80, 100),las=1)
mtext(side = 3,"Randomised to Arm B")
for(i in unique(randomised_to_arm_b$Donor.ID)){
  temp_df = subset(randomised_to_arm_b,randomised_to_arm_b$Donor.ID==i)
  temp_df = temp_df[order(temp_df$timepoint,decreasing = F),]
  points(temp_df$timepoint,temp_df$tb_tasa,pch=ifelse(temp_df$Outcome!="success",19,1),col=ifelse(temp_df$Outcome!="success","darkred","black"),type = "o",cex=1.5)
}

j=1
for(i in c(0,8,16,24)){
  if(i!=0){
    temp_df = merge(subset(randomised_to_arm_b,randomised_to_arm_b$timepoint==0,select=c("Donor.ID","tb_tasa")),subset(randomised_to_arm_b,randomised_to_arm_b$timepoint==i,select=c("Donor.ID","tb_tasa")),by="Donor.ID")
  mtext(side = 3,line = j,paste("Week 0 vs Week ",i," p value =",round(wilcox.test(temp_df$tb_tasa.x,temp_df$tb_tasa.y,paired = T)$p.value,digits = 4),sep = ""))
  j=j+1
  }
  mtext(side = 1,at=i,line = 2,length(subset(randomised_to_arm_b,randomised_to_arm_b$timepoint==i)$tb_tasa))

  mtext(side = 1,at=i,line = 3,sum(!is.na(subset(randomised_to_arm_b,randomised_to_arm_b$timepoint==i)$tb_tasa)))
}

#Figure 5d
randomised_to_arm_c = subset(predict_master_df,predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Category=="TB patients"&predict_master_df$Arm=="C")

plot(NA,NA,ylim=c(0,100),xlim=c(0,24),las=1,xlab="",xaxt = "n", yaxt = "n",ylab="TB-TASA Score (%)")
mtext(side = 3,line = 6,"Figure 5d")
axis(1, at = c(0, 8, 16, 24))
axis(2, at = c(0, 20, 40, 60, 80, 100),las=1)
mtext(side = 3,"Randomised to Arm C")
for(i in unique(randomised_to_arm_c$Donor.ID)){
  temp_df = subset(randomised_to_arm_c,randomised_to_arm_c$Donor.ID==i)
  temp_df = temp_df[order(temp_df$timepoint,decreasing = F),]
  points(temp_df$timepoint,temp_df$tb_tasa,pch=ifelse(temp_df$Outcome=="success",1,ifelse(temp_df$Outcome=="TF",19,ifelse(temp_df$Outcome=="Confirmed Relapse",17,0))),col=ifelse(temp_df$Outcome!="success","darkred","black"),type = "o",cex=1.5)
}

j=1
for(i in c(0,8,16,24)){
  if(i!=0){
    temp_df = merge(subset(randomised_to_arm_c,randomised_to_arm_c$timepoint==0,select=c("Donor.ID","tb_tasa")),subset(randomised_to_arm_c,randomised_to_arm_c$timepoint==i,select=c("Donor.ID","tb_tasa")),by="Donor.ID")
  mtext(side = 3,line = j,paste("Week 0 vs Week ",i," p value =",round(wilcox.test(temp_df$tb_tasa.x,temp_df$tb_tasa.y,paired = T)$p.value,digits = 4),sep = ""))
  j=j+1
  }
  mtext(side = 1,at=i,line = 2,length(subset(randomised_to_arm_c,randomised_to_arm_c$timepoint==i)$tb_tasa))

  mtext(side = 1,at=i,line = 3,sum(!is.na(subset(randomised_to_arm_c,randomised_to_arm_c$timepoint==i)$tb_tasa)))
}


boxplot(subset(predict_master_df,predict_master_df$Arm=="A"&predict_master_df$timepoint==24&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="success")$tb_tasa,subset(predict_master_df,predict_master_df$Arm=="A"&predict_master_df$timepoint==24&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="REINFECTION")$tb_tasa,subset(predict_master_df,predict_master_df$Arm=="B"&predict_master_df$timepoint==24&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="success")$tb_tasa,subset(predict_master_df,predict_master_df$Arm=="C"&predict_master_df$timepoint==16&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="success")$tb_tasa,subset(predict_master_df,predict_master_df$Arm=="C"&predict_master_df$timepoint==16&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome!="success")$tb_tasa,cex=0,range = 0,las=1,ylab="TB-TASA Score (%)")
mtext(side = 3,line = 6,"Figure 5e")

mtext(side = 3,line = 1,paste("A suc. vs B suc = ",round(wilcox.test(subset(predict_master_df,predict_master_df$Arm=="A"&predict_master_df$timepoint==24&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="success")$tb_tasa,subset(predict_master_df,predict_master_df$Arm=="B"&predict_master_df$timepoint==24&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="success")$tb_tasa)$p.value,digits=4)))

mtext(side = 3,line = 2,paste("A suc. vs C suc = ",round(wilcox.test(subset(predict_master_df,predict_master_df$Arm=="A"&predict_master_df$timepoint==24&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="success")$tb_tasa,subset(predict_master_df,predict_master_df$Arm=="C"&predict_master_df$timepoint==16&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="success")$tb_tasa)$p.value,digits=4)))

mtext(side = 3,line = 3,paste("B suc. vs C suc = ",round(wilcox.test(subset(predict_master_df,predict_master_df$Arm=="B"&predict_master_df$timepoint==24&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="success")$tb_tasa,subset(predict_master_df,predict_master_df$Arm=="C"&predict_master_df$timepoint==16&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="success")$tb_tasa)$p.value,digits=4)))

mtext(side = 3,line = 4,paste("C suc. vs C Not Favor. = ",round(wilcox.test(subset(predict_master_df,predict_master_df$Arm=="C"&predict_master_df$timepoint==16&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="success")$tb_tasa,subset(predict_master_df,predict_master_df$Arm=="C"&predict_master_df$timepoint==16&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome!="success")$tb_tasa)$p.value,digits=4)))

for(i in 1:5){
  if(i ==1){
    temp_df = subset(predict_master_df,predict_master_df$Arm=="A"&predict_master_df$timepoint==24&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="success")$tb_tasa
    temp_arm = "Arm A"
    temp_outcome = "success"
    temp_outcome_col ="black"
    temp_outcome_pch = 1
    }

  if(i ==2){
  temp_df =subset(predict_master_df,predict_master_df$Arm=="A"&predict_master_df$timepoint==24&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="REINFECTION")$tb_tasa
  temp_arm = "Arm A"
  temp_outcome = "reinfec."
  temp_outcome_col ="darkred"
  temp_outcome_pch = 19
  }

  if(i==3){
    temp_df = subset(predict_master_df,predict_master_df$Arm=="B"&predict_master_df$timepoint==24&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="success")$tb_tasa
    temp_arm = "Arm B"
    temp_outcome = "success"
    temp_outcome_col ="black"
    temp_outcome_pch = 1
  }

  if(i==4){
    temp_df = subset(predict_master_df,predict_master_df$Arm=="C"&predict_master_df$timepoint==16&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome=="success")$tb_tasa
    temp_arm = "Arm C"
    temp_outcome = "success"
    temp_outcome_col ="black"
    temp_outcome_pch = 1
  }

  if(i==5){
    temp_df = subset(predict_master_df,predict_master_df$Arm=="C"&predict_master_df$timepoint==16&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome!="success")$tb_tasa
    temp_arm = "Arm C"
    temp_outcome = "relapse/Tx Fail"
    temp_outcome_col ="darkred"
    temp_outcome_pch = ifelse(subset(predict_master_df,predict_master_df$Arm=="C"&predict_master_df$timepoint==16&predict_master_df$stimulation=="MTBLYSATE"&predict_master_df$Outcome!="success")$Outcome=="TF",19,17)

  }

  points(rep(i,times=length(temp_df)),temp_df,pch=temp_outcome_pch,col=temp_outcome_col,cex=1.5)
  mtext(side = 1,at = i,line = 1,temp_arm,temp_df,cex=0.7)
  mtext(side = 1,at = i,line = 2,temp_outcome,temp_df,cex=0.7)
  mtext(side = 1,at = i,line = 3,length(temp_df))
  mtext(side = 1,at = i,line = 4,sum(!is.na(temp_df)))
}

temp_auc = subset(predict_master_df,predict_master_df$Arm=="C"&predict_master_df$timepoint==16&predict_master_df$stimulation=="MTBLYSATE",select=c("Outcome","tb_tasa"))
temp_auc$Outcome = ifelse(temp_auc$Outcome=="success","success","not_favourable")

temp_auc = pROC::roc(temp_auc$Outcome,temp_auc$tb_tasa)
set.seed(12345)
a_roc.ci.se.obj = ci.se(temp_auc, specificities = seq(0, 1, 0.01), boot.n=1000, parallel = F)
stats_a_roc = format(wilcox.test(temp_auc$cases, temp_auc$controls, paired = F, alternative = "great")$p.value, digits = 2)
transColour = alpha("darkred",0.3)
plot(c(0,0),ylim=c(0,1),xlim=c(1,0),cex=0)
points(temp_auc$specificities,temp_auc$sensitivities,type="l",lwd=2,col=temp_col)
plot(a_roc.ci.se.obj, type="shape", col= transColour, no.roc = T, lty = 0)
mtext(side = 3,line = 2,paste("AUC = ",round(temp_auc$auc,digits = 2), "95%CI (",round(ci.auc(temp_auc)[1],digits = 2), "-",round(ci.auc(temp_auc)[3],digits = 2),sep = " "))
mtext(side = 3,line = 1,"Arm C Relapse and Tx Failure vs Arm C Success")
mtext(side = 3,line = 6,"Figure 5f")


#Figure 6

venous_vs_fingerprick_df$bcksub_freq_for_log_plot = ifelse(venous_vs_fingerprick_df$bcksub_ifng_tnf_cd3_freq<0.001,0.001,venous_vs_fingerprick_df$raw_ifng_tnf_cd3_freq)
venous_450uL = subset(venous_vs_fingerprick_df,venous_vs_fingerprick_df$source=="450uL_venous"&venous_vs_fingerprick_df$stimulation=="ECE")
fingerprick_df = subset(venous_vs_fingerprick_df,venous_vs_fingerprick_df$source=="fingerprick"&venous_vs_fingerprick_df$stimulation=="ECE")
venous_100uL = subset(venous_vs_fingerprick_df,venous_vs_fingerprick_df$source=="100uL_venous"&venous_vs_fingerprick_df$stimulation=="ECE")


venous_450uL$bcksub_ifng_tnf_cd3_freq = ifelse(venous_450uL$bcksub_ifng_tnf_cd3_freq<0,0,venous_450uL$bcksub_ifng_tnf_cd3_freq)

fingerprick_df$bcksub_ifng_tnf_cd3_freq = ifelse(fingerprick_df$bcksub_ifng_tnf_cd3_freq<0,0,fingerprick_df$bcksub_ifng_tnf_cd3_freq)

venous_100uL$bcksub_ifng_tnf_cd3_freq = ifelse(venous_100uL$bcksub_ifng_tnf_cd3_freq<0,0,venous_100uL$bcksub_ifng_tnf_cd3_freq)


par(mfrow=c(3,2))
for(i in 1:2){
  if(i ==1){
    temp_df = venous_450uL
    temp_title = "Figure 6a"
  }

  if(i ==2){
    temp_df = fingerprick_df
    temp_title = "Figure 6b"
  }
  boxplot(subset(temp_df,temp_df$group=="Healthy QFT-")$bcksub_freq_for_log_plot,
          subset(temp_df,temp_df$group=="Healthy QFT+")$bcksub_freq_for_log_plot,
          subset(temp_df,temp_df$group=="TB")$bcksub_freq_for_log_plot,cex=0,col="white",border = c("blue","orange","red"),las=2,lwd=2,range=0,log="y",ylim=c(0.001,10))
  mtext(side = 3,line = 5,temp_title)
  y_ticks <- c(outer(1:9, 10^(-3:1), "*"))
  axis(2, at = y_ticks, labels = y_ticks, las = 1)

  points(rep(1,times=length(subset(temp_df,temp_df$group=="Healthy QFT-")$bcksub_freq_for_log_plot)),subset(temp_df,temp_df$group=="Healthy QFT-")$bcksub_freq_for_log_plot,cex=1.5)

  points(rep(2,times=length(subset(temp_df,temp_df$group=="Healthy QFT+")$bcksub_freq_for_log_plot)),subset(temp_df,temp_df$group=="Healthy QFT+")$bcksub_freq_for_log_plot,cex=1.5)

  points(rep(3,times=length(subset(temp_df,temp_df$group=="TB")$bcksub_freq_for_log_plot)),subset(temp_df,temp_df$group=="TB")$bcksub_freq_for_log_plot,cex=1.5)

  mtext(side = 1,at = 1,line = 1,"IGRA-")
  mtext(side = 1,at = 2,line = 1,"IGRA+")
  mtext(side = 1,at = 3,line = 1,"TB")

  mtext(side = 1,at = 1,line = 2,length(subset(temp_df,temp_df$group=="Healthy QFT-")$bcksub_freq_for_log_plot))
  mtext(side = 1,at = 2,line = 2,length(subset(temp_df,temp_df$group=="Healthy QFT+")$bcksub_freq_for_log_plot))
  mtext(side = 1,at = 3,line = 2,length(subset(temp_df,temp_df$group=="TB")$bcksub_freq_for_log_plot))

  mtext(side = 2,line = 2.5,"Freq. of ESAT6/CFP10/ESPC specific T cells (%)",cex=0.8)
  mtext(side = 3,unique(temp_df$source))
  abline(h=10)

  a.vs.b = wilcox.test(subset(temp_df,temp_df$group=="Healthy QFT-")$bcksub_freq_for_log_plot,
                       subset(temp_df,temp_df$group=="Healthy QFT+")$bcksub_freq_for_log_plot)

  a.vs.c = wilcox.test(subset(temp_df,temp_df$group=="Healthy QFT-")$bcksub_freq_for_log_plot,
                       subset(temp_df,temp_df$group=="TB")$bcksub_freq_for_log_plot)

  b.vs.c = wilcox.test(subset(temp_df,temp_df$group=="Healthy QFT+")$bcksub_freq_for_log_plot,
                       subset(temp_df,temp_df$group=="TB")$bcksub_freq_for_log_plot)

  mtext(side = 3, line = 1,paste("igra_neg vs igra_pos = ",round(a.vs.b$p.value,digits = 4),sep = ""))
  mtext(side = 3, line = 2,paste("igra_neg vs tb = ",round(a.vs.c$p.value,digits = 4),sep = ""))
  mtext(side = 3, line = 3,paste("igra_pos vs tb = ",round(b.vs.c$p.value,digits = 4),sep = ""))
}

venous_450uL_tb_vs_igra_pos = subset(venous_450uL,venous_450uL$group%in%c("TB","Healthy QFT+"))
fingerprick_df_tb_vs_igra_pos = subset(fingerprick_df,fingerprick_df$group%in%c("TB","Healthy QFT+"))
venous_100uL_tb_vs_igra_pos = subset(venous_100uL,venous_100uL$group%in%c("TB","Healthy QFT+"))


for(i in 1:3){
  if(i ==1){
    temp_df = venous_450uL
    temp_title = "Figure 6c"
  }

  if(i ==2){
    temp_df = fingerprick_df
    temp_title = "Figure 6d"
  }

  if(i ==3){
    temp_df = venous_100uL
    temp_title = "Figure 6e"
  }

  boxplot(subset(temp_df,temp_df$group=="Healthy QFT-")$tb_tasa,
          subset(temp_df,temp_df$group=="Healthy QFT+")$tb_tasa,
          subset(temp_df,temp_df$group=="TB")$tb_tasa,cex=0,col="white",border = c("blue","orange","red"),las=2,lwd=2,ylim=c(0,100),range = 0)
  mtext(side = 3,line = 5,temp_title)

  points(rep(1,times=length(subset(temp_df,temp_df$group=="Healthy QFT-")$tb_tasa)),subset(temp_df,temp_df$group=="Healthy QFT-")$tb_tasa,cex=1.5)

  points(rep(2,times=length(subset(temp_df,temp_df$group=="Healthy QFT+")$tb_tasa)),subset(temp_df,temp_df$group=="Healthy QFT+")$tb_tasa,cex=1.5)

  points(rep(3,times=length(subset(temp_df,temp_df$group=="TB")$tb_tasa)),subset(temp_df,temp_df$group=="TB")$tb_tasa,cex=1.5)


  mtext(side = 1,at = 1,line = 1,"IGRA-")
  mtext(side = 1,at = 2,line = 1,"IGRA+")
  mtext(side = 1,at = 3,line = 1,"TB")
  mtext(side = 1,at = 1,line = 2,length(subset(temp_df,temp_df$group=="Healthy QFT-")$tb_tasa))
  mtext(side = 1,at = 2,line = 2,length(subset(temp_df,temp_df$group=="Healthy QFT+")$tb_tasa))
  mtext(side = 1,at = 3,line = 2,length(subset(temp_df,temp_df$group=="TB")$tb_tasa))

  mtext(side = 1,at = 1,line = 3,sum(!is.na(subset(temp_df,temp_df$group=="Healthy QFT-")$tb_tasa)))
  mtext(side = 1,at = 2,line = 3,sum(!is.na(subset(temp_df,temp_df$group=="Healthy QFT+")$tb_tasa)))
  mtext(side = 1,at = 3,line = 3,sum(!is.na(subset(temp_df,temp_df$group=="TB")$tb_tasa)))

  mtext(side = 1,at = 1,line = 4,sum(subset(temp_df,temp_df$group=="Healthy QFT-")$tb_tasa>10,na.rm = T))
  mtext(side = 1,at = 2,line = 4,sum(subset(temp_df,temp_df$group=="Healthy QFT+")$tb_tasa>10,na.rm = T))
  mtext(side = 1,at = 3,line = 4,sum(subset(temp_df,temp_df$group=="TB")$tb_tasa>10,na.rm = T))


  mtext(side = 2,line = 2.5,"Proportion of HLA-DR+ ESAT6/CFP10/ESPC-specific T cells (%)",cex=0.8)
  mtext(side = 3,unique(temp_df$source))
  abline(h=10)

  a.vs.b = wilcox.test(subset(temp_df,temp_df$group=="Healthy QFT-")$tb_tasa,
                       subset(temp_df,temp_df$group=="Healthy QFT+")$tb_tasa)

  a.vs.c = wilcox.test(subset(temp_df,temp_df$group=="Healthy QFT-")$tb_tasa,
                       subset(temp_df,temp_df$group=="TB")$tb_tasa)

  b.vs.c = wilcox.test(subset(temp_df,temp_df$group=="Healthy QFT+")$tb_tasa,
                       subset(temp_df,temp_df$group=="TB")$tb_tasa)

  mtext(side = 3, line = 1,paste("igra_neg vs igra_pos = ",round(a.vs.b$p.value,digits = 4),sep = ""))
  mtext(side = 3, line = 2,paste("igra_neg vs tb = ",round(a.vs.c$p.value,digits = 4),sep = ""))
  mtext(side = 3, line = 3,paste("igra_pos vs tb = ",round(b.vs.c$p.value,digits = 4),sep = ""))
}

venous_auc_cd3_tb_vs_igra_pos_hla_dr = pROC::roc(venous_450uL_tb_vs_igra_pos$group,venous_450uL_tb_vs_igra_pos$tb_tasa)

fingerprick_auc_cd3_tb_vs_igra_pos_hla_dr = pROC::roc(fingerprick_df_tb_vs_igra_pos$group,fingerprick_df_tb_vs_igra_pos$tb_tasa)

venous100uL_auc_cd3_tb_vs_igra_pos_hla_dr = pROC::roc(venous_100uL_tb_vs_igra_pos$group,venous_100uL_tb_vs_igra_pos$tb_tasa)

plot(c(0,0),ylim=c(0,1),xlim=c(1,0),cex=0)
mtext(side = 3,"Figure 6f")
for(i in 1:3){

  if(i ==1){
    temp_auc = venous_auc_cd3_tb_vs_igra_pos_hla_dr
    temp_col ="darkgreen"
    temp_df = venous_450uL
    temp_auc_postion = 0.4
    temp_sample_type = "venous_450uL"
  }

  if(i ==2){
    temp_auc = fingerprick_auc_cd3_tb_vs_igra_pos_hla_dr
    temp_col ="purple"
    temp_df = fingerprick_df
    temp_auc_postion = 0.3
    temp_sample_type = "fingerprick"
  }

  if(i ==3){
    temp_auc = venous100uL_auc_cd3_tb_vs_igra_pos_hla_dr
    temp_col ="black"
    temp_df = venous_100uL
    temp_auc_postion = 0.2
    temp_sample_type = "venous_100uL"
  }

  set.seed(12345)
  a_roc.ci.se.obj = ci.se(temp_auc, specificities = seq(0, 1, 0.01), boot.n=1000, parallel = F)
  stats_a_roc = format(wilcox.test(temp_auc$cases, temp_auc$controls, paired = F, alternative = "great")$p.value, digits = 2)
  transColour = alpha(temp_col,0.3)
  points(temp_auc$specificities,temp_auc$sensitivities,type="l",lwd=2,col=temp_col)
  plot(a_roc.ci.se.obj, type="shape", col= transColour, no.roc = T, lty = 0)

  threshold_sens =  sum(subset(temp_df,temp_df$group=="TB")$tb_tasa>10,na.rm = T)/length(subset(temp_df,temp_df$group=="TB")$tb_tasa)

  threshold_spec = sum(sum(subset(temp_df,temp_df$group=="Healthy QFT+")$tb_tasa<=10,na.rm = T),sum(is.na(subset(temp_df,temp_df$group=="Healthy QFT+")$tb_tasa)))/length(subset(temp_df,temp_df$group=="Healthy QFT+")$tb_tasa)
  text(0.5,temp_auc_postion,paste(temp_sample_type,round(as.numeric(temp_auc$auc),digits = 2),"95%CI (",round(ci.auc(temp_auc)[1],digits = 2), " - ",round(ci.auc(temp_auc)[3],digits = 2),")","sens ",round(threshold_sens*100,digits =0),"spec ",round(threshold_spec*100,digits = 0),sep=" "),col=temp_col,cex=0.9)
  points(threshold_spec,threshold_sens,pch=19,col=temp_col,cex=2)

}

segments(x0 = 1,x1 = 0.7,y0 = 0.9,y1=0.9,lwd = 2,lty=2)
segments(x0 = 1,x1 = 0.7,y0 = 1,y1=1,lwd = 2,lty=2)
segments(x0 = 1,x1 = 1,y0 = 0.9,y1=1,lwd = 2,lty=2)
segments(x0 = 0.7,x1 = 0.7,y0 = 0.9,y1=1,lwd = 2,lty=2)


#Supp Figure 1

Mice.Normalised.Ct.Value = Mice.Normalised.Ct.Value[-c(1:3),]

Processed.DataSet = matrix(data=NA,0,0)
for(i in c(2,3,4,5,7,8,9,10,12,13,14,15,17,18,19,20,22,23,24,26,27,28,29)){
  Temp =as.matrix(Mice.Normalised.Ct.Value[,i])
  Temp.Ct = as.matrix(Temp[-c(1:3),])
  Temp.Ct = as.numeric(Temp.Ct)
  Temp.Ct =t(as.data.frame(Temp.Ct))
  Temp.1 = as.data.frame(t(Temp[1:3,]))
  Temp.data.frame= cbind(Temp.1,Temp.Ct)
  Processed.DataSet = rbind.data.frame(Processed.DataSet,Temp.data.frame)
}
rownames(Processed.DataSet)= Processed.DataSet[,1]
Processed.DataSet = Processed.DataSet[,-c(2183,2184)]

colnames(Processed.DataSet) =	Mice.Normalised.Ct.Value$X


Mean.week6_B6 = t(as.matrix(colMeans(subset(Processed.DataSet,Processed.DataSet$`Mouse Strain`=="B6"&Processed.DataSet$`Time Post Infection`=="6 weeks",select=colnames(Processed.DataSet)[4:2182]))))

Mean.week6_C3H = t(as.matrix(colMeans(subset(Processed.DataSet,Processed.DataSet$`Mouse Strain`=="C3H"&Processed.DataSet$`Time Post Infection`=="6 weeks",select=colnames(Processed.DataSet)[4:2182]))))

Antigens.of.interest = c("Rv3874","Rv3875","Rv1886c","Rv0288","Rv3620c","Rv3619c", "Rv1813c","Rv3615c","Rv0125","Rv1196","Rv2608")

Week6_Expression.DF = t(rbind.data.frame(Mean.week6_B6,Mean.week6_C3H))
colnames(Week6_Expression.DF)=c("B6_Week6","C3H_Week6")
Week6_Expression.DF = 40-Week6_Expression.DF

Week6_Expression.DF = as.data.frame(Week6_Expression.DF)
Week6_Expression.DF$Symbol.Size = ifelse(rownames(Week6_Expression.DF)%in% Antigens.of.interest,1.5,0.1)
Week6_Expression.DF$Symbol.Type = ifelse(rownames(Week6_Expression.DF)%in% Antigens.of.interest,19,1)
Week6_Expression.DF$Symbol.Color = ifelse(rownames(Week6_Expression.DF)%in% Antigens.of.interest,"red","black")

Antigens.of.interest_Week6 = subset(Week6_Expression.DF,rownames(Week6_Expression.DF)%in%Antigens.of.interest)
par(mfrow=c(1,1))
plot(Week6_Expression.DF$B6_Week6,Week6_Expression.DF$C3H_Week6,cex=Week6_Expression.DF$Symbol.Size,col=Week6_Expression.DF$Symbol.Color,pch=Week6_Expression.DF$Symbol.Type,ylim=c(5,25),xlim=c(5,25))
text(Antigens.of.interest_Week6$B6_Week6,Antigens.of.interest_Week6$C3H_Week6,rownames(Antigens.of.interest_Week6))
mtext(side = 3,"Supp Figure 1")

#Supp Figure 2
venous_vs_fingerprick_df$raw_freq_for_log_plot = ifelse(venous_vs_fingerprick_df$raw_ifng_tnf_cd3_freq<0.001,0.001,venous_vs_fingerprick_df$raw_ifng_tnf_cd3_freq)
merged_df = merge(subset(venous_vs_fingerprick_df,venous_vs_fingerprick_df$source=="450uL_venous",select=c("pid","stimulation","raw_freq_for_log_plot")),subset(venous_vs_fingerprick_df,venous_vs_fingerprick_df$source=="fingerprick",select=c("pid","stimulation","raw_freq_for_log_plot")),by=c("pid","stimulation"))

par(mfrow=c(1,1))
plot(merged_df$raw_freq_for_log_plot.x,merged_df$raw_freq_for_log_plot.y,cex=2,ylim=c(0.001,10),xlim=c(0.001,10),col=ifelse(merged_df$stimulation=="UNS","blue",ifelse(merged_df$stimulation=="ECE","darkred","black")),ylab="",xlab = "",log="xy")

mtext(side = 1, line = 3,"Freq. of ESAT-6/CFP-10/ESPC-specific T cells in 450uL venous blood (%)")
mtext(side = 2, line = 3,"Freq. of ESAT-6/CFP-10/ESPC-specific T cells in 100uL fingerprick blood (%)")

temp_ccc = CCC(merged_df$raw_freq_for_log_plot.x,merged_df$raw_freq_for_log_plot.y,na.rm = T)

mtext(side = 3,paste("CCC est. = ", round(temp_ccc$rho.c[1],digits = 2), " 95%CI (",round(temp_ccc$rho.c[2],digits = 2),"-",round(temp_ccc$rho.c[3],digits = 2),")",sep = " "))
abline(a=0,b=1,lwd=2,lty=2)
mtext(side = 3,line = 2,"Supp Figure 2")




