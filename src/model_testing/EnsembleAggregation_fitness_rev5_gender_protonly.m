function [speciesNames,y_out,y_full,y_delta,y_delta_full] = EnsembleAggregation_fitness_rev5_gender_protonly(w,w_map,w_full,exptin,exptout,inputs_idx,gender_value) 
% Generated by Netflux 0.08a on 26-Mar-2020
 
% load species parameters 
% species parameters 
speciesNames = {'AngII','AT1R','AGT','ACE','NOX','ROS','ET1','ETAR','DAG','PKC','TRPC','NE','BAR','AC','cAMP','PKA','CREB','CBP','TGFB','TGFB1R','smad3','smad7','latentTGFB','BAMBI','PDGF','PDGFR','NP','NPRA','cGMP','PKG','tension','B1int','Rho','ROCK','Ca','calcineurin','NFAT','IL6','gp130','STAT','IL1','IL1RI','TNFa','TNFaR','NFKB','PI3K','Akt','p38','TRAF','ASK1','MKK3','PP1','JNK','abl','Rac1','MEKK1','MKK4','ERK','Ras','Raf','MEK1','FAK','epac','Factin','FA','cmyc','CTGF','proliferation','SRF','EDAFN','aSMA','AP1','TIMP1','TIMP2','PAI1','proMMP14','proMMP1','proMMP2','proMMP9','fibronectin','periostin','proCI','proCIII','B3int','Src','Grb2','p130Cas','YAP','MRTF','Gactin','TNC','mTORC1','mTORC2','p70S6K','EBP1','syndecan4','proMMP3','proMMP8','proMMP12','thrombospondin4','osteopontin','contractility','RhoGEF','RhoGDI','talin','vinculin','paxillin','MLC','AT2R','TEAD2','TEAD4','Fos','Jun','RELA','CCND1','NR3C1','TCF4','SOX9','NFE2L2','AHR','ATF3','CLOCK','CTSC','CTSL','elastin','ETV4','FLI1','HMGA1','LEF1','LOXL1','MECOM','NOTCH1','P4H','POU2F1','RUNX1','RXRA','SETDB1','SMARCA2','SMARCA4','STAT5A','TBP','TP53','WT1','YY1','E2','ERB','GPR30','CyclinB1','CDK1','AMPK','ERX',}; 
tau = [1, 1.000000e-01, 10, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01, 10, 10, 1.000000e-01, 1, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01, 1, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 1, 1.000000e-01, 1.000000e-01, 10, 1.000000e-01, 10, 10, 1.000000e-01, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 10, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 10, 10, 10, 10, 10, 10, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 1, 1, 10, 10, 10, 1, 1, 1, 1, 10, 1, 1, 10, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, ]; 
ymax = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ];
 
% reaction parameters 
w_default = [1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ]; 
n = [1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.250000e+00, 1.250000e+00, 1.250000e+00, 1.250000e+00, 1.250000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, ]; 
EC50 = [6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, ]; 
% rpar = [w;n;EC50];
% params = {rpar,tau,ymax,speciesNames};
y0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ]; 
options = []; 

% load loop conditions
load(exptin,'input_scaled_protonly');
load(exptout,'output_scaled');
isFemale = contains(input_scaled_protonly.Properties.RowNames, ["TAVRA","TAVRB","TAVRC","TAVRD"]);

% match experimental inputs/outputs with model variables (w/speciesNames)
w_expand = w_default;                  % 10.28.2020: expand cluster-reduced vector to full vector
w_expand(w_full) = w(w_map);

% inputs_idx = [1 2 8 10 4 5 6];                                % inputs: w (manual)
for output = 1:length(output_scaled.Properties.VariableNames)
    idx = find(strcmpi(output_scaled.Properties.VariableNames(output),speciesNames));
    outputs_idx(output) = idx;                                  % ouputs: speciesNames
end


% run simulations: loop through all conditions
for stimulus = 1:size(input_scaled_protonly,1)
    w_stim = w_expand;         % w_stim: set to w, then specific values modified
    for species = 1:length(inputs_idx)
        w_stim(inputs_idx(species)) = input_scaled_protonly{stimulus,species};
    end
    % modify gender-specific rxns
    if isFemale(stimulus)
        w_stim(308) = gender_value * w_stim(308);
    else
        w_stim(330) = gender_value * w_stim(330);
    end
    
    rpar_stim = [w_stim;n;EC50];        % integrate w_stim into params
    params_stim = {rpar_stim,tau,ymax,speciesNames};
    tspan = [0 40];
    [~,y2] = ode23(@ODEfun,tspan,y0,options,params_stim);   % run simulation
    if stimulus == 1
        y_full = real(y2(end,:));
        y_out = real(y2(end,outputs_idx));  % y_out: end-point output levels
    else
        y_full = cat(1,y_full,real(y2(end,:)));
        y_out = cat(1,y_out,real(y2(end,outputs_idx)));
    end
end

% calculate fitness metric: MSE
% delta(Activity) = post-TAVR (rows 1:4) - pre-TAVR (rows 5:8)
delta_act = @(n) n(1:4,:) - n(5:8,:);
y_delta = delta_act(y_out);
y_delta_full = delta_act(y_full);
% expt_delta = delta_act(output_scaled{:,:});

% SE = (y_delta - expt_delta).^2;         % squared error for deltas
% MSE = mean(mean(SE));                   % average across all elements
% error = MSE;                            % assign method to fcn output
end
 
function dydt=ODEfun(t,y,params) 
% Assign names for parameters 
[rpar,tau,ymax,speciesNames]=params{:}; 
AngII = 1; 
AT1R = 2; 
AGT = 3; 
ACE = 4; 
NOX = 5; 
ROS = 6; 
ET1 = 7; 
ETAR = 8; 
DAG = 9; 
PKC = 10; 
TRPC = 11; 
NE = 12; 
BAR = 13; 
AC = 14; 
cAMP = 15; 
PKA = 16; 
CREB = 17; 
CBP = 18; 
TGFB = 19; 
TGFB1R = 20; 
smad3 = 21; 
smad7 = 22; 
latentTGFB = 23; 
BAMBI = 24; 
PDGF = 25; 
PDGFR = 26; 
NP = 27; 
NPRA = 28; 
cGMP = 29; 
PKG = 30; 
tension = 31; 
B1int = 32; 
Rho = 33; 
ROCK = 34; 
Ca = 35; 
calcineurin = 36; 
NFAT = 37; 
IL6 = 38; 
gp130 = 39; 
STAT = 40; 
IL1 = 41; 
IL1RI = 42; 
TNFa = 43; 
TNFaR = 44; 
NFKB = 45; 
PI3K = 46; 
Akt = 47; 
p38 = 48; 
TRAF = 49; 
ASK1 = 50; 
MKK3 = 51; 
PP1 = 52; 
JNK = 53; 
abl = 54; 
Rac1 = 55; 
MEKK1 = 56; 
MKK4 = 57; 
ERK = 58; 
Ras = 59; 
Raf = 60; 
MEK1 = 61; 
FAK = 62; 
epac = 63; 
Factin = 64; 
FA = 65; 
cmyc = 66; 
CTGF = 67; 
proliferation = 68; 
SRF = 69; 
EDAFN = 70; 
aSMA = 71; 
AP1 = 72; 
TIMP1 = 73; 
TIMP2 = 74; 
PAI1 = 75; 
proMMP14 = 76; 
proMMP1 = 77; 
proMMP2 = 78; 
proMMP9 = 79; 
fibronectin = 80; 
periostin = 81; 
proCI = 82; 
proCIII = 83; 
B3int = 84; 
Src = 85; 
Grb2 = 86; 
p130Cas = 87; 
YAP = 88; 
MRTF = 89; 
Gactin = 90; 
TNC = 91; 
mTORC1 = 92; 
mTORC2 = 93; 
p70S6K = 94; 
EBP1 = 95; 
syndecan4 = 96; 
proMMP3 = 97; 
proMMP8 = 98; 
proMMP12 = 99; 
thrombospondin4 = 100; 
osteopontin = 101; 
contractility = 102; 
RhoGEF = 103; 
RhoGDI = 104; 
talin = 105; 
vinculin = 106; 
paxillin = 107; 
MLC = 108; 
AT2R = 109; 
TEAD2 = 110; 
TEAD4 = 111; 
Fos = 112; 
Jun = 113; 
RELA = 114; 
CCND1 = 115; 
NR3C1 = 116; 
TCF4 = 117; 
SOX9 = 118; 
NFE2L2 = 119; 
AHR = 120; 
ATF3 = 121; 
CLOCK = 122; 
CTSC = 123; 
CTSL = 124; 
elastin = 125; 
ETV4 = 126; 
FLI1 = 127; 
HMGA1 = 128; 
LEF1 = 129; 
LOXL1 = 130; 
MECOM = 131; 
NOTCH1 = 132; 
P4H = 133; 
POU2F1 = 134; 
RUNX1 = 135; 
RXRA = 136; 
SETDB1 = 137; 
SMARCA2 = 138; 
SMARCA4 = 139; 
STAT5A = 140; 
TBP = 141; 
TP53 = 142; 
WT1 = 143; 
YY1 = 144; 
E2 = 145; 
ERB = 146; 
GPR30 = 147; 
CyclinB1 = 148; 
CDK1 = 149; 
AMPK = 150; 
ERX = 151; 
dydt = zeros(151,1); 
dydt(AngII) = (OR(rpar(1,1),AND(rpar(:,14),act(y(AGT),rpar(:,14)),act(y(ACE),rpar(:,14))))*ymax(AngII) - y(AngII))/tau(AngII); 
dydt(AT1R) = (OR(act(y(AngII),rpar(:,19)),act(y(tension),rpar(:,168)))*ymax(AT1R) - y(AT1R))/tau(AT1R); 
dydt(AGT) = (AND(rpar(:,29),inhib(y(AT1R),rpar(:,29)),act(y(p38),rpar(:,29)),inhib(y(JNK),rpar(:,29)))*ymax(AGT) - y(AGT))/tau(AGT); 
dydt(ACE) = (act(y(TGFB1R),rpar(:,50))*ymax(ACE) - y(ACE))/tau(ACE); 
dydt(NOX) = (OR(act(y(AT1R),rpar(:,20)),act(y(TGFB1R),rpar(:,99)))*ymax(NOX) - y(NOX))/tau(NOX); 
dydt(ROS) = (OR(act(y(NOX),rpar(:,21)),OR(act(y(ETAR),rpar(:,38)),inhib(y(ERB),rpar(:,312))))*ymax(ROS) - y(ROS))/tau(ROS); 
dydt(ET1) = (OR(rpar(1,9),OR(act(y(AP1),rpar(:,18)),inhib(y(NR3C1),rpar(:,247))))*ymax(ET1) - y(ET1))/tau(ET1); 
dydt(ETAR) = (act(y(ET1),rpar(:,56))*ymax(ETAR) - y(ETAR))/tau(ETAR); 
dydt(DAG) = (OR(act(y(ETAR),rpar(:,113)),act(y(AT1R),rpar(:,114)))*ymax(DAG) - y(DAG))/tau(DAG); 
dydt(PKC) = (OR(act(y(syndecan4),rpar(:,132)),AND(rpar(:,147),act(y(DAG),rpar(:,147)),act(y(mTORC2),rpar(:,147))))*ymax(PKC) - y(PKC))/tau(PKC); 
dydt(TRPC) = (OR(act(y(DAG),rpar(:,115)),act(y(tension),rpar(:,164)))*ymax(TRPC) - y(TRPC))/tau(TRPC); 
dydt(NE) = (rpar(1,7)*ymax(NE) - y(NE))/tau(NE); 
dydt(BAR) = (act(y(NE),rpar(:,55))*ymax(BAR) - y(BAR))/tau(BAR); 
dydt(AC) = (OR(act(y(BAR),rpar(:,64)),AND(rpar(:,65),act(y(AT1R),rpar(:,65)),act(y(BAR),rpar(:,65))))*ymax(AC) - y(AC))/tau(AC); 
dydt(cAMP) = (OR(act(y(AC),rpar(:,66)),act(y(ERB),rpar(:,316)))*ymax(cAMP) - y(cAMP))/tau(cAMP); 
dydt(PKA) = (OR(act(y(cAMP),rpar(:,45)),act(y(ERB),rpar(:,321)))*ymax(PKA) - y(PKA))/tau(PKA); 
dydt(CREB) = (act(y(PKA),rpar(:,53))*ymax(CREB) - y(CREB))/tau(CREB); 
dydt(CBP) = (OR(inhib(y(smad3),rpar(:,46)),inhib(y(CREB),rpar(:,47)))*ymax(CBP) - y(CBP))/tau(CBP); 
dydt(TGFB) = (OR(rpar(1,2),OR(AND(rpar(:,12),act(y(latentTGFB),rpar(:,12)),act(y(proMMP9),rpar(:,12))),OR(AND(rpar(:,13),act(y(latentTGFB),rpar(:,13)),act(y(proMMP2),rpar(:,13))),OR(inhib(y(E2),rpar(:,319)),inhib(y(AMPK),rpar(:,332))))))*ymax(TGFB) - y(TGFB))/tau(TGFB); 
dydt(TGFB1R) = (AND(rpar(:,51),act(y(TGFB),rpar(:,51)),inhib(y(BAMBI),rpar(:,51)))*ymax(TGFB1R) - y(TGFB1R))/tau(TGFB1R); 
dydt(smad3) = (OR(AND(rpar(:,30),act(y(TGFB1R),rpar(:,30)),inhib(y(smad7),rpar(:,30)),inhib(y(PKG),rpar(:,30))),OR(act(y(Akt),rpar(:,143)),OR(inhib(y(ERB),rpar(:,313)),act(y(JNK),rpar(:,317)))))*ymax(smad3) - y(smad3))/tau(smad3); 
dydt(smad7) = (OR(act(y(STAT),rpar(:,105)),AND(rpar(:,155),act(y(AP1),rpar(:,155)),inhib(y(YAP),rpar(:,155))))*ymax(smad7) - y(smad7))/tau(smad7); 
dydt(latentTGFB) = (act(y(AP1),rpar(:,68))*ymax(latentTGFB) - y(latentTGFB))/tau(latentTGFB); 
dydt(BAMBI) = (AND(rpar(:,104),act(y(TGFB),rpar(:,104)),act(y(IL1RI),rpar(:,104)))*ymax(BAMBI) - y(BAMBI))/tau(BAMBI); 
dydt(PDGF) = (rpar(1,8)*ymax(PDGF) - y(PDGF))/tau(PDGF); 
dydt(PDGFR) = (act(y(PDGF),rpar(:,63))*ymax(PDGFR) - y(PDGFR))/tau(PDGFR); 
dydt(NP) = (rpar(1,10)*ymax(NP) - y(NP))/tau(NP); 
dydt(NPRA) = (act(y(NP),rpar(:,72))*ymax(NPRA) - y(NPRA))/tau(NPRA); 
dydt(cGMP) = (act(y(NPRA),rpar(:,73))*ymax(cGMP) - y(cGMP))/tau(cGMP); 
dydt(PKG) = (act(y(cGMP),rpar(:,74))*ymax(PKG) - y(PKG))/tau(PKG); 
dydt(tension) = (OR(rpar(1,3),AND(rpar(:,162),act(y(FA),rpar(:,162)),act(y(contractility),rpar(:,162))))*ymax(tension) - y(tension))/tau(tension); 
dydt(B1int) = (OR(AND(rpar(:,44),act(y(PKC),rpar(:,44)),act(y(tension),rpar(:,44))),OR(act(y(tension),rpar(:,48)),act(y(E2),rpar(:,327))))*ymax(B1int) - y(B1int))/tau(B1int); 
dydt(Rho) = (OR(act(y(TGFB1R),rpar(:,118)),OR(AND(rpar(:,131),act(y(RhoGEF),rpar(:,131)),inhib(y(RhoGDI),rpar(:,131))),inhib(y(ERB),rpar(:,315))))*ymax(Rho) - y(Rho))/tau(Rho); 
dydt(ROCK) = (OR(act(y(Rho),rpar(:,70)),OR(act(y(TGFB),rpar(:,322)),act(y(ERB),rpar(:,325))))*ymax(ROCK) - y(ROCK))/tau(ROCK); 
dydt(Ca) = (OR(act(y(TRPC),rpar(:,116)),inhib(y(ERX),rpar(:,329)))*ymax(Ca) - y(Ca))/tau(Ca); 
dydt(calcineurin) = (act(y(Ca),rpar(:,117))*ymax(calcineurin) - y(calcineurin))/tau(calcineurin); 
dydt(NFAT) = (OR(act(y(calcineurin),rpar(:,110)),OR(act(y(cmyc),rpar(:,270)),inhib(y(cmyc),rpar(:,301))))*ymax(NFAT) - y(NFAT))/tau(NFAT); 
dydt(IL6) = (OR(rpar(1,4),OR(AND(rpar(:,15),act(y(CREB),rpar(:,15)),act(y(CBP),rpar(:,15))),OR(act(y(NFKB),rpar(:,16)),OR(act(y(AP1),rpar(:,17)),OR(act(y(YAP),rpar(:,179)),OR(act(y(AHR),rpar(:,192)),OR(act(y(HMGA1),rpar(:,199)),inhib(y(YY1),rpar(:,251)))))))))*ymax(IL6) - y(IL6))/tau(IL6); 
dydt(gp130) = (act(y(IL6),rpar(:,22))*ymax(gp130) - y(gp130))/tau(gp130); 
dydt(STAT) = (act(y(gp130),rpar(:,27))*ymax(STAT) - y(STAT))/tau(STAT); 
dydt(IL1) = (rpar(1,5)*ymax(IL1) - y(IL1))/tau(IL1); 
dydt(IL1RI) = (act(y(IL1),rpar(:,58))*ymax(IL1RI) - y(IL1RI))/tau(IL1RI); 
dydt(TNFa) = (rpar(1,6)*ymax(TNFa) - y(TNFa))/tau(TNFa); 
dydt(TNFaR) = (act(y(TNFa),rpar(:,71))*ymax(TNFaR) - y(TNFaR))/tau(TNFaR); 
dydt(NFKB) = (OR(act(y(IL1RI),rpar(:,26)),OR(act(y(ERK),rpar(:,34)),OR(act(y(p38),rpar(:,35)),OR(act(y(Akt),rpar(:,100)),OR(act(y(smad3),rpar(:,283)),inhib(y(ERB),rpar(:,328)))))))*ymax(NFKB) - y(NFKB))/tau(NFKB); 
dydt(PI3K) = (OR(act(y(TNFaR),rpar(:,28)),OR(act(y(TGFB1R),rpar(:,96)),OR(act(y(PDGFR),rpar(:,97)),act(y(FAK),rpar(:,98)))))*ymax(PI3K) - y(PI3K))/tau(PI3K); 
dydt(Akt) = (OR(AND(rpar(:,146),act(y(PI3K),rpar(:,146)),act(y(mTORC2),rpar(:,146))),act(y(ERB),rpar(:,333)))*ymax(Akt) - y(Akt))/tau(Akt); 
dydt(p38) = (OR(act(y(ROS),rpar(:,24)),OR(act(y(MKK3),rpar(:,79)),OR(act(y(Ras),rpar(:,95)),OR(AND(rpar(:,107),act(y(Rho),rpar(:,107)),inhib(y(Rac1),rpar(:,107))),inhib(y(ERB),rpar(:,314))))))*ymax(p38) - y(p38))/tau(p38); 
dydt(TRAF) = (OR(act(y(TGFB1R),rpar(:,80)),act(y(TNFaR),rpar(:,88)))*ymax(TRAF) - y(TRAF))/tau(TRAF); 
dydt(ASK1) = (OR(act(y(TRAF),rpar(:,89)),act(y(IL1RI),rpar(:,92)))*ymax(ASK1) - y(ASK1))/tau(ASK1); 
dydt(MKK3) = (act(y(ASK1),rpar(:,90))*ymax(MKK3) - y(MKK3))/tau(MKK3); 
dydt(PP1) = (act(y(p38),rpar(:,78))*ymax(PP1) - y(PP1))/tau(PP1); 
dydt(JNK) = (OR(act(y(ROS),rpar(:,25)),OR(AND(rpar(:,83),inhib(y(NFKB),rpar(:,83)),act(y(MKK4),rpar(:,83))),OR(AND(rpar(:,108),inhib(y(Rho),rpar(:,108)),act(y(MKK4),rpar(:,108))),inhib(y(ERB),rpar(:,320)))))*ymax(JNK) - y(JNK))/tau(JNK); 
dydt(abl) = (act(y(PDGFR),rpar(:,84))*ymax(abl) - y(abl))/tau(abl); 
dydt(Rac1) = (OR(act(y(abl),rpar(:,85)),AND(rpar(:,128),act(y(abl),rpar(:,128)),act(y(p130Cas),rpar(:,128))))*ymax(Rac1) - y(Rac1))/tau(Rac1); 
dydt(MEKK1) = (OR(act(y(FAK),rpar(:,67)),act(y(Rac1),rpar(:,81)))*ymax(MEKK1) - y(MEKK1))/tau(MEKK1); 
dydt(MKK4) = (OR(act(y(MEKK1),rpar(:,82)),act(y(ASK1),rpar(:,91)))*ymax(MKK4) - y(MKK4))/tau(MKK4); 
dydt(ERK) = (OR(act(y(ROS),rpar(:,23)),OR(AND(rpar(:,77),inhib(y(PP1),rpar(:,77)),act(y(MEK1),rpar(:,77))),OR(AND(rpar(:,172),act(y(ROS),rpar(:,172)),inhib(y(AT2R),rpar(:,172))),inhib(y(ERB),rpar(:,326)))))*ymax(ERK) - y(ERK))/tau(ERK); 
dydt(Ras) = (OR(act(y(AT1R),rpar(:,111)),act(y(Grb2),rpar(:,122)))*ymax(Ras) - y(Ras))/tau(Ras); 
dydt(Raf) = (act(y(Ras),rpar(:,75))*ymax(Raf) - y(Raf))/tau(Raf); 
dydt(MEK1) = (AND(rpar(:,76),inhib(y(ERK),rpar(:,76)),act(y(Raf),rpar(:,76)))*ymax(MEK1) - y(MEK1))/tau(MEK1); 
dydt(FAK) = (act(y(B1int),rpar(:,120))*ymax(FAK) - y(FAK))/tau(FAK); 
dydt(epac) = (act(y(cAMP),rpar(:,69))*ymax(epac) - y(epac))/tau(epac); 
dydt(Factin) = (AND(rpar(:,135),act(y(ROCK),rpar(:,135)),act(y(Gactin),rpar(:,135)))*ymax(Factin) - y(Factin))/tau(Factin); 
dydt(FA) = (AND(rpar(:,157),act(y(vinculin),rpar(:,157)),inhib(y(paxillin),rpar(:,157)))*ymax(FA) - y(FA))/tau(FA); 
dydt(cmyc) = (act(y(JNK),rpar(:,86))*ymax(cmyc) - y(cmyc))/tau(cmyc); 
dydt(CTGF) = (AND(rpar(:,31),act(y(CBP),rpar(:,31)),act(y(smad3),rpar(:,31)),act(y(ERK),rpar(:,31)))*ymax(CTGF) - y(CTGF))/tau(CTGF); 
dydt(proliferation) = (OR(act(y(AP1),rpar(:,52)),OR(act(y(CREB),rpar(:,54)),OR(act(y(CTGF),rpar(:,57)),OR(act(y(PKC),rpar(:,59)),OR(act(y(cmyc),rpar(:,87)),OR(AND(rpar(:,142),act(y(p70S6K),rpar(:,142)),inhib(y(EBP1),rpar(:,142))),act(y(CDK1),rpar(:,318))))))))*ymax(proliferation) - y(proliferation))/tau(proliferation); 
dydt(SRF) = (OR(act(y(MRTF),rpar(:,137)),OR(act(y(CREB),rpar(:,233)),inhib(y(CREB),rpar(:,292))))*ymax(SRF) - y(SRF))/tau(SRF); 
dydt(EDAFN) = (act(y(NFAT),rpar(:,49))*ymax(EDAFN) - y(EDAFN))/tau(EDAFN); 
dydt(aSMA) = (OR(act(y(SRF),rpar(:,112)),OR(act(y(YAP),rpar(:,170)),OR(act(y(NFKB),rpar(:,178)),OR(inhib(y(FLI1),rpar(:,235)),OR(inhib(y(LEF1),rpar(:,253)),OR(act(y(smad3),rpar(:,255)),OR(inhib(y(NFAT),rpar(:,282)),inhib(y(smad3),rpar(:,297)))))))))*ymax(aSMA) - y(aSMA))/tau(aSMA); 
dydt(AP1) = (AND(rpar(:,167),act(y(Fos),rpar(:,167)),act(y(Jun),rpar(:,167)))*ymax(AP1) - y(AP1))/tau(AP1); 
dydt(TIMP1) = (OR(act(y(AP1),rpar(:,43)),act(y(FLI1),rpar(:,225)))*ymax(TIMP1) - y(TIMP1))/tau(TIMP1); 
dydt(TIMP2) = (act(y(Jun),rpar(:,177))*ymax(TIMP2) - y(TIMP2))/tau(TIMP2); 
dydt(PAI1) = (OR(act(y(smad3),rpar(:,93)),OR(inhib(y(AHR),rpar(:,184)),OR(inhib(y(STAT),rpar(:,206)),OR(act(y(STAT5A),rpar(:,210)),OR(act(y(RUNX1),rpar(:,237)),OR(act(y(cmyc),rpar(:,250)),OR(act(y(RXRA),rpar(:,257)),OR(act(y(NOTCH1),rpar(:,265)),OR(inhib(y(cmyc),rpar(:,295)),inhib(y(NOTCH1),rpar(:,299)))))))))))*ymax(PAI1) - y(PAI1))/tau(PAI1); 
dydt(proMMP14) = (OR(act(y(AP1),rpar(:,62)),OR(act(y(NFKB),rpar(:,94)),OR(inhib(y(TEAD4),rpar(:,175)),OR(act(y(WT1),rpar(:,185)),OR(act(y(LEF1),rpar(:,194)),OR(act(y(Jun),rpar(:,197)),OR(act(y(NFAT),rpar(:,222)),OR(act(y(SMARCA2),rpar(:,232)),OR(act(y(TBP),rpar(:,234)),OR(act(y(FLI1),rpar(:,266)),inhib(y(SMARCA2),rpar(:,291))))))))))))*ymax(proMMP14) - y(proMMP14))/tau(proMMP14); 
dydt(proMMP1) = (OR(AND(rpar(:,36),inhib(y(smad3),rpar(:,36)),act(y(NFKB),rpar(:,36))),OR(AND(rpar(:,37),inhib(y(smad3),rpar(:,37)),act(y(AP1),rpar(:,37))),OR(act(y(ATF3),rpar(:,207)),OR(act(y(LEF1),rpar(:,218)),OR(inhib(y(Jun),rpar(:,254)),OR(inhib(y(ATF3),rpar(:,286)),inhib(y(LEF1),rpar(:,289))))))))*ymax(proMMP1) - y(proMMP1))/tau(proMMP1); 
dydt(proMMP2) = (OR(act(y(AP1),rpar(:,41)),act(y(NFKB),rpar(:,285)))*ymax(proMMP2) - y(proMMP2))/tau(proMMP2); 
dydt(proMMP9) = (AND(rpar(:,42),act(y(NFKB),rpar(:,42)),act(y(AP1),rpar(:,42)))*ymax(proMMP9) - y(proMMP9))/tau(proMMP9); 
dydt(fibronectin) = (OR(act(y(NFKB),rpar(:,101)),OR(act(y(TBP),rpar(:,205)),OR(act(y(NFE2L2),rpar(:,244)),OR(inhib(y(cmyc),rpar(:,274)),act(y(NFAT),rpar(:,278))))))*ymax(fibronectin) - y(fibronectin))/tau(fibronectin); 
dydt(periostin) = (OR(AND(rpar(:,32),act(y(CBP),rpar(:,32)),act(y(smad3),rpar(:,32))),OR(AND(rpar(:,33),act(y(CREB),rpar(:,33)),act(y(CBP),rpar(:,33))),OR(act(y(SOX9),rpar(:,193)),OR(inhib(y(TP53),rpar(:,245)),inhib(y(smad3),rpar(:,260))))))*ymax(periostin) - y(periostin))/tau(periostin); 
dydt(proCI) = (OR(AND(rpar(:,60),act(y(CBP),rpar(:,60)),act(y(smad3),rpar(:,60)),inhib(y(epac),rpar(:,60))),OR(act(y(SRF),rpar(:,106)),OR(act(y(CCND1),rpar(:,220)),OR(act(y(SETDB1),rpar(:,238)),act(y(ERB),rpar(:,331))))))*ymax(proCI) - y(proCI))/tau(proCI); 
dydt(proCIII) = (OR(AND(rpar(:,61),act(y(CBP),rpar(:,61)),act(y(smad3),rpar(:,61)),inhib(y(epac),rpar(:,61))),OR(act(y(SRF),rpar(:,109)),OR(act(y(ETV4),rpar(:,203)),OR(act(y(SMARCA4),rpar(:,219)),act(y(ERB),rpar(:,331))))))*ymax(proCIII) - y(proCIII))/tau(proCIII); 
dydt(B3int) = (OR(AND(rpar(:,149),act(y(tension),rpar(:,149)),inhib(y(thrombospondin4),rpar(:,149))),act(y(osteopontin),rpar(:,153)))*ymax(B3int) - y(B3int))/tau(B3int); 
dydt(Src) = (OR(act(y(B3int),rpar(:,119)),act(y(PDGFR),rpar(:,126)))*ymax(Src) - y(Src))/tau(Src); 
dydt(Grb2) = (AND(rpar(:,121),act(y(FAK),rpar(:,121)),act(y(Src),rpar(:,121)))*ymax(Grb2) - y(Grb2))/tau(Grb2); 
dydt(p130Cas) = (OR(AND(rpar(:,125),act(y(FAK),rpar(:,125)),act(y(Src),rpar(:,125))),AND(rpar(:,127),act(y(tension),rpar(:,127)),act(y(Src),rpar(:,127))))*ymax(p130Cas) - y(p130Cas))/tau(p130Cas); 
dydt(YAP) = (OR(act(y(Factin),rpar(:,129)),OR(act(y(AT1R),rpar(:,169)),inhib(y(cmyc),rpar(:,241))))*ymax(YAP) - y(YAP))/tau(YAP); 
dydt(MRTF) = (inhib(y(Gactin),rpar(:,134))*ymax(MRTF) - y(MRTF))/tau(MRTF); 
dydt(Gactin) = (inhib(y(Factin),rpar(:,136))*ymax(Gactin) - y(Gactin))/tau(Gactin); 
dydt(TNC) = (OR(act(y(NFKB),rpar(:,144)),OR(inhib(y(Jun),rpar(:,186)),OR(inhib(y(NFE2L2),rpar(:,212)),OR(inhib(y(TBP),rpar(:,230)),OR(act(y(ATF3),rpar(:,242)),inhib(y(TCF4),rpar(:,268)))))))*ymax(TNC) - y(TNC))/tau(TNC); 
dydt(mTORC1) = (act(y(Akt),rpar(:,139))*ymax(mTORC1) - y(mTORC1))/tau(mTORC1); 
dydt(mTORC2) = (inhib(y(p70S6K),rpar(:,145))*ymax(mTORC2) - y(mTORC2))/tau(mTORC2); 
dydt(p70S6K) = (act(y(mTORC1),rpar(:,140))*ymax(p70S6K) - y(p70S6K))/tau(p70S6K); 
dydt(EBP1) = (inhib(y(mTORC1),rpar(:,141))*ymax(EBP1) - y(EBP1))/tau(EBP1); 
dydt(syndecan4) = (AND(rpar(:,138),act(y(tension),rpar(:,138)),inhib(y(TNC),rpar(:,138)))*ymax(syndecan4) - y(syndecan4))/tau(syndecan4); 
dydt(proMMP3) = (AND(rpar(:,151),inhib(y(smad3),rpar(:,151)),act(y(NFKB),rpar(:,151)),act(y(AP1),rpar(:,151)))*ymax(proMMP3) - y(proMMP3))/tau(proMMP3); 
dydt(proMMP8) = (AND(rpar(:,150),inhib(y(smad3),rpar(:,150)),act(y(NFKB),rpar(:,150)),act(y(AP1),rpar(:,150)))*ymax(proMMP8) - y(proMMP8))/tau(proMMP8); 
dydt(proMMP12) = (act(y(CREB),rpar(:,154))*ymax(proMMP12) - y(proMMP12))/tau(proMMP12); 
dydt(thrombospondin4) = (act(y(smad3),rpar(:,148))*ymax(thrombospondin4) - y(thrombospondin4))/tau(thrombospondin4); 
dydt(osteopontin) = (act(y(AP1),rpar(:,152))*ymax(osteopontin) - y(osteopontin))/tau(osteopontin); 
dydt(contractility) = (AND(rpar(:,161),act(y(Factin),rpar(:,161)),act(y(vinculin),rpar(:,161)),act(y(MLC),rpar(:,161)))*ymax(contractility) - y(contractility))/tau(contractility); 
dydt(RhoGEF) = (AND(rpar(:,123),act(y(FAK),rpar(:,123)),act(y(Src),rpar(:,123)))*ymax(RhoGEF) - y(RhoGEF))/tau(RhoGEF); 
dydt(RhoGDI) = (OR(inhib(y(Src),rpar(:,124)),OR(act(y(PKA),rpar(:,130)),inhib(y(PKC),rpar(:,133))))*ymax(RhoGDI) - y(RhoGDI))/tau(RhoGDI); 
dydt(talin) = (OR(act(y(B1int),rpar(:,158)),act(y(B3int),rpar(:,159)))*ymax(talin) - y(talin))/tau(talin); 
dydt(vinculin) = (AND(rpar(:,160),act(y(tension),rpar(:,160)),act(y(talin),rpar(:,160)))*ymax(vinculin) - y(vinculin))/tau(vinculin); 
dydt(paxillin) = (AND(rpar(:,156),act(y(FAK),rpar(:,156)),act(y(Src),rpar(:,156)),act(y(MLC),rpar(:,156)))*ymax(paxillin) - y(paxillin))/tau(paxillin); 
dydt(MLC) = (act(y(ROCK),rpar(:,163))*ymax(MLC) - y(MLC))/tau(MLC); 
dydt(AT2R) = (act(y(AngII),rpar(:,171))*ymax(AT2R) - y(AT2R))/tau(AT2R); 
dydt(TEAD2) = (act(y(YAP),rpar(:,165))*ymax(TEAD2) - y(TEAD2))/tau(TEAD2); 
dydt(TEAD4) = (act(y(YAP),rpar(:,166))*ymax(TEAD4) - y(TEAD4))/tau(TEAD4); 
dydt(Fos) = (OR(act(y(ERK),rpar(:,40)),OR(act(y(JNK),rpar(:,103)),act(y(STAT),rpar(:,228))))*ymax(Fos) - y(Fos))/tau(Fos); 
dydt(Jun) = (OR(act(y(ERK),rpar(:,39)),act(y(JNK),rpar(:,102)))*ymax(Jun) - y(Jun))/tau(Jun); 
dydt(RELA) = (act(y(cmyc),rpar(:,276))*ymax(RELA) - y(RELA))/tau(RELA); 
dydt(CCND1) = (OR(inhib(y(Jun),rpar(:,213)),inhib(y(TEAD2),rpar(:,224)))*ymax(CCND1) - y(CCND1))/tau(CCND1); 
dydt(NR3C1) = (OR(act(y(TEAD2),rpar(:,198)),OR(act(y(NFKB),rpar(:,280)),inhib(y(NFKB),rpar(:,307))))*ymax(NR3C1) - y(NR3C1))/tau(NR3C1); 
dydt(TCF4) = (OR(act(y(STAT),rpar(:,200)),act(y(CREB),rpar(:,229)))*ymax(TCF4) - y(TCF4))/tau(TCF4); 
dydt(SOX9) = (OR(inhib(y(TEAD4),rpar(:,189)),OR(act(y(RELA),rpar(:,202)),act(y(NFAT),rpar(:,216))))*ymax(SOX9) - y(SOX9))/tau(SOX9); 
dydt(NFE2L2) = (OR(act(y(Jun),rpar(:,223)),act(y(RELA),rpar(:,256)))*ymax(NFE2L2) - y(NFE2L2))/tau(NFE2L2); 
dydt(AHR) = (inhib(y(smad3),rpar(:,264))*ymax(AHR) - y(AHR))/tau(AHR); 
dydt(ATF3) = (act(y(smad3),rpar(:,217))*ymax(ATF3) - y(ATF3))/tau(ATF3); 
dydt(CLOCK) = (OR(act(y(NFKB),rpar(:,273)),OR(inhib(y(cmyc),rpar(:,281)),inhib(y(NFKB),rpar(:,304))))*ymax(CLOCK) - y(CLOCK))/tau(CLOCK); 
dydt(CTSC) = (OR(act(y(FLI1),rpar(:,181)),OR(act(y(HMGA1),rpar(:,182)),OR(inhib(y(STAT5A),rpar(:,187)),OR(inhib(y(RXRA),rpar(:,195)),OR(act(y(MECOM),rpar(:,261)),inhib(y(MECOM),rpar(:,298)))))))*ymax(CTSC) - y(CTSC))/tau(CTSC); 
dydt(CTSL) = (OR(act(y(smad3),rpar(:,180)),OR(act(y(SOX9),rpar(:,188)),OR(inhib(y(RUNX1),rpar(:,236)),OR(act(y(MECOM),rpar(:,248)),inhib(y(MECOM),rpar(:,293))))))*ymax(CTSL) - y(CTSL))/tau(CTSL); 
dydt(elastin) = (OR(inhib(y(STAT),rpar(:,173)),act(y(TEAD4),rpar(:,174)))*ymax(elastin) - y(elastin))/tau(elastin); 
dydt(ETV4) = (OR(inhib(y(SRF),rpar(:,183)),inhib(y(TEAD2),rpar(:,191)))*ymax(ETV4) - y(ETV4))/tau(ETV4); 
dydt(FLI1) = (inhib(y(cmyc),rpar(:,215))*ymax(FLI1) - y(FLI1))/tau(FLI1); 
dydt(HMGA1) = (OR(act(y(TEAD4),rpar(:,269)),inhib(y(TEAD4),rpar(:,300)))*ymax(HMGA1) - y(HMGA1))/tau(HMGA1); 
dydt(LEF1) = (OR(inhib(y(SRF),rpar(:,211)),act(y(NFAT),rpar(:,231)))*ymax(LEF1) - y(LEF1))/tau(LEF1); 
dydt(LOXL1) = (OR(inhib(y(FLI1),rpar(:,190)),OR(act(y(SETDB1),rpar(:,208)),OR(inhib(y(WT1),rpar(:,227)),OR(act(y(RUNX1),rpar(:,246)),OR(act(y(CLOCK),rpar(:,249)),OR(inhib(y(SETDB1),rpar(:,287)),inhib(y(CLOCK),rpar(:,294))))))))*ymax(LOXL1) - y(LOXL1))/tau(LOXL1); 
dydt(MECOM) = (OR(act(y(Jun),rpar(:,201)),OR(act(y(STAT),rpar(:,214)),OR(act(y(CREB),rpar(:,239)),act(y(YAP),rpar(:,259)))))*ymax(MECOM) - y(MECOM))/tau(MECOM); 
dydt(NOTCH1) = (OR(inhib(y(STAT),rpar(:,258)),act(y(cmyc),rpar(:,263)))*ymax(NOTCH1) - y(NOTCH1))/tau(NOTCH1); 
dydt(P4H) = (OR(inhib(y(cmyc),rpar(:,176)),inhib(y(TP53),rpar(:,204)))*ymax(P4H) - y(P4H))/tau(P4H); 
dydt(POU2F1) = (0*ymax(POU2F1) - y(POU2F1))/tau(POU2F1); 
dydt(RUNX1) = (OR(act(y(CREB),rpar(:,243)),OR(act(y(Fos),rpar(:,271)),inhib(y(Fos),rpar(:,302))))*ymax(RUNX1) - y(RUNX1))/tau(RUNX1); 
dydt(RXRA) = (act(y(cmyc),rpar(:,240))*ymax(RXRA) - y(RXRA))/tau(RXRA); 
dydt(SETDB1) = (inhib(y(NFAT),rpar(:,284))*ymax(SETDB1) - y(SETDB1))/tau(SETDB1); 
dydt(SMARCA2) = (OR(act(y(YAP),rpar(:,221)),OR(act(y(Jun),rpar(:,226)),OR(act(y(STAT),rpar(:,262)),inhib(y(YAP),rpar(:,290)))))*ymax(SMARCA2) - y(SMARCA2))/tau(SMARCA2); 
dydt(SMARCA4) = (OR(inhib(y(TEAD2),rpar(:,196)),OR(act(y(TEAD4),rpar(:,252)),inhib(y(TEAD4),rpar(:,296))))*ymax(SMARCA4) - y(SMARCA4))/tau(SMARCA4); 
dydt(STAT5A) = (OR(act(y(RELA),rpar(:,209)),inhib(y(RELA),rpar(:,288)))*ymax(STAT5A) - y(STAT5A))/tau(STAT5A); 
dydt(TBP) = (act(y(NFAT),rpar(:,275))*ymax(TBP) - y(TBP))/tau(TBP); 
dydt(TP53) = (inhib(y(STAT),rpar(:,267))*ymax(TP53) - y(TP53))/tau(TP53); 
dydt(WT1) = (OR(act(y(cmyc),rpar(:,277)),inhib(y(cmyc),rpar(:,305)))*ymax(WT1) - y(WT1))/tau(WT1); 
dydt(YY1) = (OR(act(y(RELA),rpar(:,272)),OR(act(y(cmyc),rpar(:,279)),OR(inhib(y(RELA),rpar(:,303)),inhib(y(cmyc),rpar(:,306)))))*ymax(YY1) - y(YY1))/tau(YY1); 
dydt(E2) = (rpar(1,11)*ymax(E2) - y(E2))/tau(E2); 
dydt(ERB) = (OR(act(y(E2),rpar(:,308)),AND(rpar(:,324),act(y(cAMP),rpar(:,324)),act(y(PKA),rpar(:,324))))*ymax(ERB) - y(ERB))/tau(ERB); 
dydt(GPR30) = (act(y(E2),rpar(:,309))*ymax(GPR30) - y(GPR30))/tau(GPR30); 
dydt(CyclinB1) = (inhib(y(GPR30),rpar(:,310))*ymax(CyclinB1) - y(CyclinB1))/tau(CyclinB1); 
dydt(CDK1) = (act(y(CyclinB1),rpar(:,311))*ymax(CDK1) - y(CDK1))/tau(CDK1); 
dydt(AMPK) = (act(y(E2),rpar(:,323))*ymax(AMPK) - y(AMPK))/tau(AMPK); 
dydt(ERX) = (act(y(E2),rpar(:,330))*ymax(ERX) - y(ERX))/tau(ERX); 
end

% utility functions 
function fact = act(x,rpar) 
% hill activation function with parameters w (weight), n (Hill coeff), EC50 
    w = rpar(1); 
    n = rpar(2); 
    EC50 = rpar(3); 
    beta = (EC50.^n - 1)./(2*EC50.^n - 1); 
    K = (beta - 1).^(1./n); 
    fact = w.*(beta.*x.^n)./(K.^n + x.^n); 
    if fact>w,                 % cap fact(x)<= 1 
        fact = w; 
    end
end
 
function finhib = inhib(x,rpar) 
% inverse hill function with parameters w (weight), n (Hill coeff), EC50 
    finhib = rpar(1) - act(x,rpar);
end
 
function z = OR(x,y) 
% OR logic gate 
    z = x + y - x*y;
end
 
function z = AND(rpar,varargin) 
% AND logic gate, multiplying all of the reactants together 
    w = rpar(1); 
    if w == 0, 
        z = 0; 
    else 
        v = cell2mat(varargin); 
        z = prod(v)/w^(nargin-2);  
    end 
end
