%% Analysis for Ensemble fitting
% 10.15.2020 JR
%% *Load Data*
%% 
% * Reaction weights (x_all): from GA ensemble
% * Prediction data (y_out + y_delta): from GA aggregation script
% * Experimental data (expt_out + expt_delta): from data processing script
% * Clinical data (clin_mat): from Aguado et. al. supplemental data
%% 
% load in filepaths

homedir = "E:/Research/Aim3/ModelExpansion/";
datadir = strcat(homedir, "ModelFitting_Ensemble/1_2_rev5_gender/");
colordir = "E:/Research/Aim2/BrewerMap-master/";
exptdir = "E:/Research/Aim2/ModelFitting/ModelFitting_AnsethData/data_processed_current/";
exptdir_raw = "E:/Research/Aim2/ModelFitting/ModelFitting_AnsethData/data_raw/";
addpath(datadir, colordir, exptdir);
%% 
% load reaction weights

filedelim = "ga_finalset*";
datastr = dir(strcat(datadir, filedelim));
datanames = char(datastr.name);
for file = 1:size(datanames,1)
    filename = split(datanames(file,:));
    load(filename{1});
    if file == 1
        x_all = finalset{1};
        fval_all = finalset{2};
    else
        x_all = [x_all; finalset{1}];
        fval_all = [fval_all; finalset{2}];
    end
end
%% 
% load aggregated prediction data

% default SNM (nofit)
filename = "ga_aggregate_nofit.mat";
load(filename);

% fitted SNM
filename = "ga_aggregate_results.mat";
load(filename);

% fitted SNM - testing set
filename = "ga_aggregate_results_protonly.mat";
load(filename);

% species names (all)
filename = "ga_speciesNames.mat";
load(filename);

% aggregate fitted results across runs
y_out_mean = mean(y_out_all,3);
y_delta_mean = mean(y_delta_all,3);
SE_mean = mean(SE_all,3);

adjustPValues = @(x) reshape(mafdr(reshape(x, [], 1), "BHFDR", "true"), size(x,1), []);

%% 
% load experimental data

exptin = 'trainingdata_inputs_Anseth_11232020.mat';
exptout = 'trainingdata_outputs_Anseth_12042020.mat';
load(exptin);
load(exptout);

speciesNames_out = output_scaled.Properties.VariableNames;
expt_delta = output_scaled{1:8,:} - output_scaled{9:end,:};

% get all input/output idxs (from expt data)
inputs_other = {'tension', 'NE', 'ET1', 'E2'};
inputs_all = input_scaled.Properties.VariableNames;
inputs_all(length(inputs_all)+1:length(inputs_all)+length(inputs_other)) = inputs_other;
[~,inputs_idx] = ismember(inputs_all, speciesNames);

outputs_other = {'EDAFN', 'proMMP9', 'proMMP3', 'proMMP8', 'proMMP12', 'thrombospondin4'};
outputs_all = speciesNames_out;
outputs_all(length(outputs_all)+1:length(outputs_all)+length(outputs_other)) = outputs_other;
[~,outputs_idx] = ismember(outputs_all, speciesNames);

%% 
% load clinical data

filename = "aav3233_Data_file_S1.xlsx";
opts = detectImportOptions(strcat(exptdir_raw,filename),"Sheet","Table S1");
clin = readtable(strcat(exptdir_raw,filename),opts);

% process data: calculate delta_KCCQScore, extract numeric columns,
% calculate Pearson correlations
clin.delta_KCCQScore = clin.x1_monthPost_KCCQScore - clin.Pre_KCCQScore;
numericvars = varfun(@isnumeric, clin, "OutputFormat", "uniform");
clin_numeric = clin(:,numericvars);   % keep only training set (1-8)
clin_vars = {'Age','DaysBetweenTAVR','AoVMaxVelocity','ValveArea','MeanPressureGradient', ...
    'StrokeVolumeIndex','LVIDD','LVIDS','EjectionFraction','STSscore', ...
    'PreKCCQscore','PostKCCQscore','WBCcount','BNP','DeltaKCCQscore'};
clin_numeric.Properties.VariableNames = clin_vars;
clin_mat = clin_numeric{1:8,3:end};
clin_mat_allreps = clin_numeric{:,3:end};
clin_vars_reduced = clin_vars(3:end);

%% *Analyze pre-/post-TAVR changes in expression: 2-tailed t-tests*
% conduct ttest + calculate p-values: experimental, fitted predictions,
% non-fitted predictions
threshold = 0.05;

[~,expt_delta_p] = ttest2(output_scaled{1:8,:},output_scaled{9:end,:},"Dim",1);
[~,y_delta_p] = ttest2(y_out_mean(1:8,:),y_out_mean(9:end,:),"Dim",1);
[~,y_delta_p_nofit] = ttest2(y_out_nofit(1:8,:),y_out_nofit(9:end,:),"Dim",1);

y_TP = sum((y_delta_p <= threshold) & (expt_delta_p <= threshold));
y_FP = sum((y_delta_p <= threshold) & (expt_delta_p > threshold));
y_TN = sum((y_delta_p > threshold) & (expt_delta_p > threshold));
y_FN = sum((y_delta_p > threshold) & (expt_delta_p <= threshold));
y_acc = (y_TP + y_TN) / (y_TP + y_TN + y_FP + y_FN);
y_f1 = (2*y_TP) / ((2*y_TP) + y_FP + y_FN); 

y_TP_nofit = sum((y_delta_p_nofit <= 0.1) & (expt_delta_p <= 0.1));
y_FP_nofit = sum((y_delta_p_nofit <= 0.1) & (expt_delta_p > 0.1));
y_TN_nofit = sum((y_delta_p_nofit > 0.1) & (expt_delta_p > 0.1));
y_FN_nofit = sum((y_delta_p_nofit > 0.1) & (expt_delta_p <= 0.1));
y_acc_nofit = (y_TP_nofit + y_TN_nofit) / (y_TP_nofit + y_TN_nofit + y_FP_nofit + y_FN_nofit);
y_f1_nofit = (2*y_TP_nofit) / ((2*y_TP_nofit) + y_FP_nofit + y_FN_nofit); 

%% *Analyze output/clinical relationships: Pearson + goodness of fit*
% calculate Pearson correlations b/t model outputs + clinical variables
y_out_pre = y_out_mean(9:end,:);
y_out_pre_nofit = y_out_nofit(9:end,:);
y_out_pre_protonly = mean(y_out_all_protonly(5:end,:,:),3);
y_out_pre_allreps = [y_out_pre; y_out_pre_protonly];
expt_out_pre = output_scaled{9:end,:};
% [y_rho,y_p] = corr(y_delta_mean,clin_mat,"rows","pairwise");
% [y_rho_nofit,y_p_nofit] = corr(y_delta_nofit,clin_mat,"rows","pairwise");
% [expt_rho,expt_p] = corr(expt_delta,clin_mat,"rows","pairwise");
[y_rho, y_p] =                  corr(y_out_pre, clin_mat, "rows", "pairwise");
[y_rho_nofit, y_p_nofit] =      corr(y_out_pre_nofit, clin_mat, "rows", "pairwise");
[y_rho_allreps, y_p_allreps] =  corr(y_out_pre_allreps, clin_mat_allreps, "rows", "pairwise");
[expt_rho, expt_p] =            corr(expt_out_pre, clin_mat, "rows", "pairwise");


y_fdr = adjustPValues(y_p);
y_fdr_nofit = adjustPValues(y_p_nofit);
y_fdr_allreps = adjustPValues(y_p_allreps);
expt_fdr = adjustPValues(expt_p);

% calculate goodness of fit b/t model outputs + clinical variables
y_fit_r = zeros(size(y_out_pre, 2), size(clin_mat, 2));
y_fit_p = zeros(size(y_out_pre, 2), size(clin_mat, 2));
y_fit_r_allreps = zeros(size(y_out_pre_allreps, 2), size(clin_mat_allreps, 2));
y_fit_p_allreps = zeros(size(y_out_pre_allreps, 2), size(clin_mat_allreps, 2));
for col_clin = 1:size(clin_mat, 2)
    for col_out = 1:size(y_out_pre, 2)
        fit = fitlm(y_out_pre(:,col_out), clin_mat(:,col_clin));
        fit_allreps = fitlm(y_out_pre_allreps(:,col_out), clin_mat_allreps(:,col_clin));
        y_fit_r(col_out,col_clin) = fit.Rsquared.Ordinary;
        y_fit_p(col_out,col_clin) = fit.Coefficients.pValue(2);
        y_fit_r_allreps(col_out,col_clin) = fit_allreps.Rsquared.Ordinary;
        y_fit_p_allreps(col_out,col_clin) = fit_allreps.Coefficients.pValue(2);
    end
end

y_fit_fdr = adjustPValues(y_fit_p);
y_fit_fdr_allreps = adjustPValues(y_fit_p_allreps);

tree_rho_rows = linkage(y_rho_allreps, 'average', 'correlation');
tree_rho_cols = linkage(y_rho_allreps.', 'average', 'correlation');

%% calculate composite matrix content score (MCS)
idx_pro = [1 3 4 5 10 12 16 17 19 20 23];
idx_anti = [6 7 13 14 15];
idx_inhib = [18 21 22];

% rank training patients (sort + re-order)
rank_out = zeros(size(y_out_pre));
for col = 1:size(y_out_pre,2)
    [~,~,rank] = unique(y_out_pre(:,col));          % rank patient values (ascending)
    if any(idx_pro == col)
        rankrev = max(rank) - rank + 1;             % reverse ranks if in pro (descending)
        rank_out(:,col) = rankrev / length(rankrev);
    elseif any(idx_inhib == col)
        rankrev = max(rank) - rank + 1;             % reverse ranks if in inhib (descending)
        rank_out(:,col) = rankrev / length(rankrev);
    else
        rank_out(:,col) = rank / length(rank);      % use ascending for anti
    end
end
mcs = sum(rank_out,2);

% rank all patients
rank_out_allreps = zeros(size(y_out_pre_allreps));
for col = 1:size(y_out_pre_allreps,2)
    [~,~,rank] = unique(y_out_pre_allreps(:,col));          % rank patient values (ascending)
    if any(idx_pro == col)
        rankrev = max(rank) - rank + 1;             % reverse ranks if in pro (descending)
        rank_out_allreps(:,col) = rankrev / length(rankrev);
    elseif any(idx_inhib == col)
        rankrev = max(rank) - rank + 1;             % reverse ranks if in inhib (descending)
        rank_out_allreps(:,col) = rankrev / length(rankrev);
    else
        rank_out_allreps(:,col) = rank / length(rank);      % use ascending for anti
    end
end
mcs_allreps = sum(rank_out_allreps,2);


[mcs_rho,mcs_p] = corr(mcs, clin_mat, "rows", "pairwise", "type", "Spearman");
[mcs_rho_allreps,mcs_p_allreps] = corr(mcs_allreps, clin_mat_allreps, "rows", "pairwise", "type", "Spearman");
mcs_fdr = adjustPValues(mcs_p);
%% 
% fit linear regression model + visualize individual results (model
% predictions)

[y_lmcols,clin_lmcols] = find(y_p_allreps<0.05);

for pair = 1:length(y_lmcols)
    clin_lmcol = clin_lmcols(pair);
    y_lmcol = y_lmcols(pair);
    
    mdl = fitlm(y_out_pre_allreps(:,y_lmcol), clin_mat_allreps(:,clin_lmcol));
    
    figure("Position",[200 100 400 300]);
    mdl_plot = mdl.plot;
    xlab = split(mdl.PredictorNames,"x");
    xlabel(strcat('Pre-TAVR Activity (',speciesNames_out{y_lmcol},')'), 'Interpreter', 'tex');
    ylabel(clin_numeric.Properties.VariableNames(clin_lmcol+2));
    annot_str = "\rho = " + string(round(y_rho_allreps(y_lmcol,clin_lmcol),2)) + ...
                newline + "R^2 = " + string(round(y_fit_r_allreps(y_lmcol,clin_lmcol),2)) + ...
                newline + "p = " + string(round(y_p_allreps(y_lmcol,clin_lmcol),2,'significant'));
    annotation('textbox',[.65 .8 .2 .1], ...
                'String', annot_str, ...
                'FitBoxToText', 'on', ...
                'EdgeColor', 'none');
    legend('off');
    title(strcat('Linear Model: ',clin_numeric.Properties.VariableNames{clin_lmcol+2}), 'Interpreter', 'tex');
end

%%
% fit linear regression model + visualize individual results (expt data)

[expt_lmcols,clin_lmcols] = find(expt_p<0.05);

for pair = 1:length(expt_lmcols)
    clin_lmcol = clin_lmcols(pair);
    expt_lmcol = expt_lmcols(pair);
    
    mdl = fitlm(expt_out_pre(:,expt_lmcol), clin_mat(:,clin_lmcol));
    
    figure("Position",[200 100 400 300]);
    mdl_plot = mdl.plot;
    xlab = split(mdl.PredictorNames,"x");
    xlabel(strcat('Pre-TAVR Activity (',speciesNames_out{expt_lmcol},')'), 'Interpreter', 'tex');
    ylabel(clin_numeric.Properties.VariableNames(clin_lmcol+2));
    annot_str = "\rho = " + string(round(expt_rho(expt_lmcol,clin_lmcol),2)) + ...
                newline + "p = " + string(round(expt_p(expt_lmcol,clin_lmcol),2,'significant'));
    annotation('textbox',[.65 .8 .2 .1], ...
                'String', annot_str, ...
                'FitBoxToText', 'on', ...
                'EdgeColor', 'none');
    legend('off');
    title(strcat('Linear Model: ',clin_numeric.Properties.VariableNames{clin_lmcol+2}), 'Interpreter', 'tex');
end

%%
% Process intermediate data: filter out inputs/outputs/low \DeltaActivity
% nodes

y_delta_alloutputs = y_delta_full_all(:,outputs_idx,:);
y_delta_mean_alloutputs = mean(y_delta_alloutputs, 3);

meds_idx = 1:length(speciesNames);
meds_idx([inputs_idx outputs_idx]) = [];
y_meds = y_full_all(:,meds_idx,:);
y_meds_mean = mean(y_meds, 3);
y_delta_meds = y_delta_full_all(:,meds_idx,:);
y_delta_mean_meds = mean(y_delta_meds, 3);
% remove all-zero values
meds_zeros_idx = find(abs(mean(y_delta_mean_meds,1))<=0.01);
meds_idx(meds_zeros_idx) = [];
y_delta_mean_meds(:,meds_zeros_idx) = [];

%% *Analyze intermediate/output relationships: Pearson correlation + goodness of fit*
% calculate Pearson correlations b/t model intermediates and outputs
[y_rho_med,y_p_med] = corr(y_delta_mean_meds, y_delta_mean_alloutputs, "rows", "pairwise");

% multiple comparision correction (fwer + fdr)
y_fwer_med = y_p_med .* numel(y_p_med);
y_fdr_med = adjustPValues(y_p_med);


% calculate goodness of fit b/t model intermediates + outputs
y_fit_med_r = zeros(size(y_delta_mean_meds, 2), size(y_delta_mean_alloutputs, 2));
y_fit_med_p = zeros(size(y_delta_mean_meds, 2), size(y_delta_mean_alloutputs, 2));
for col_out = 1:size(y_delta_mean_alloutputs, 2)
    for col_med = 1:size(y_delta_mean_meds, 2)
        fit = fitlm(y_delta_mean_meds(:,col_med), y_delta_mean_alloutputs(:,col_out));
        y_fit_med_r(col_med,col_out) = fit.Rsquared.Ordinary;
        y_fit_med_p(col_med,col_out) = fit.Coefficients.pValue(2);
    end
end
y_fit_med_fdr = adjustPValues(y_fit_med_p);


% cluster correlation coeff's
tree_med_rows = linkage(y_rho_med, 'average', 'euclidean');
tree_med_cols = linkage(y_rho_med.', 'average', 'euclidean');


%% 
% reshape data + run descriptive statistics

inputs_idx = [1 2 4 5 6 8 10];
x_med = x_all;
x_med(:,inputs_idx) = [];

x_stats = mean(x_med, 1);
x_stats(2,:) = std(x_med, 1);
[~,idx] = sort(x_stats(1,:));
%% Model Fitting Visualization: parameter distributions
% visualize distributions: boxplot

figure("Position", [200 500 1000 300]);
boxplot(x_med(:,idx));
xlabel("Edge");
ylabel("Fitted Weight (a.u.)");


% visualize for clusters only: boxplot
x_inputs = x_all(:,1:4);
x_clusters = x_all(:,5:15);
[~,idx_clusters] = sort(mean(x_clusters, 1));
labels_cluster = [  "TGFB-Akt Signaling",   "AngII-YAP Signaling", ...
                    "ET1 Signaling",        "IL1 Signaling", ...
                    "NP-Calcium Signaling", "PDGF Signaling", ...
                    "JNK-ROS Signaling",    "NE Signaling", ...
                    "Integrin Signaling",   "IL6 Signaling", ...
                    "NFAT Transcription"];

figure("Position", [300 200 400 550]);
subplot(5, 1, [1,2]);
box_inputs = boxplot(fliplr(x_all(:, 1:4)), "Orientation", "horizontal");
yticklabels(flip(inputs_other));
ylabel("Input");
xlabel("Fitted Weight (a.u.)");

subplot(5, 1, [3,4,5]);
box_clusters = boxplot(x_clusters(:,idx_clusters), "Orientation", "horizontal");
yticklabels(labels_cluster(idx_clusters));
ylabel("Cluster");
xlabel("Fitted Weight (a.u.)");

%% Output Data Vlsualization: Comparison to RNA-seq Data

% Raw activity: y_out_mean + output_scaled
ydata = replace(output_scaled.Properties.RowNames,"_TAVR","-TAVR ");
figure("Position", [488 232 450 550]);
subplot(2,1,1);
fig_yout = heatmap(y_out_mean);
subplot(2,1,2);
fig_exptout = heatmap(output_scaled{:,:});
set(fig_yout,   "XDisplayLabels", nan(1, size(y_out_mean, 2)), ...
                "YData", ydata, ...
                "Title", "Fitted SNM: Predicted Output Levels");
set(fig_exptout,"XData", speciesNames_out, ...
                "YData", ydata, ...
                "Title", "Anseth Data: Output Levels");
% saveas(fig_yout, strcat(homedir,simdir,"/","y_out.fig"));
% saveas(fig_yout, strcat(homedir,simdir,"/","y_out.png"));
% saveas(fig_exptout, strcat(homedir,simdir,"/","expt_out.fig"));
% saveas(fig_exptout, strcat(homedir,simdir,"/","expt_out.png"));
%% 
% delta(activity): y_delta_mean + expt_delta

ysplit = split(ydata,"-");
ydata = unique(ysplit(:,2));
map = brewermap(21,'*RdBu');
tree = linkage(expt_delta.', 'average');

figure;
subplot(3,1,1);
[fig_dend,~,outperm] = dendrogram(tree);
set(fig_dend,   "Color", "black");
set(gca,        "Position", [0.12 0.76 0.7 0.2], ...
                "Visible", "off");
subplot(3,1,2);
fig_exptdelta = heatmap(expt_delta(:,outperm));
colorlim = max(abs(fig_exptdelta.ColorLimits));
set(fig_exptdelta,  "XDisplayLabels", nan(length(outperm),1), ...
                    "Position", [0.13 0.50 0.6826 0.25], ...
                    "YData", ydata, ...
                    "YLabel", "Experimental Data", ...
                    "Colormap", map, ...
                    "ColorLimits", [-colorlim colorlim]);
% annotation('textbox', [.82 .5 .1 .25], ...
%             'String', "\DeltaExpression_{TAVR}", ...
%             'FitBoxToText', 'on', ...
%             'EdgeColor', 'none');
        
subplot(3,1,3);
fig_ydelta = heatmap(y_delta_mean(:,outperm));
colorlim = max(abs(fig_ydelta.ColorLimits));
set(fig_ydelta, "XData", output_scaled.Properties.VariableNames(outperm), ...
                "Position", [0.13 0.21 0.6826 0.25], ...
                "YData", ydata, ...
                "YLabel", "Model Predictions", ...
                "XLabel", "Output Nodes", ...
                "Colormap",map, ...
                "ColorLimits", [-colorlim colorlim]);

% saveas(fig_ydelta, strcat(homedir,simdir,"/","y_delta.fig"));
% saveas(fig_ydelta, strcat(homedir,simdir,"/","y_delta.png"));
% saveas(fig_exptdelta, strcat(homedir,simdir,"/","expt_delta.fig"));
% saveas(fig_exptdelta, strcat(homedir,simdir,"/","expt_delta.png"));
%% 
% Squared-error: SE

figure("Position", [488 342 560 200]);
fig_se = heatmap(SE_mean(:,outperm));
set(fig_se,"XData", output_scaled.Properties.VariableNames(outperm), ...
    "YData", ydata, ...
    "Title", "Fitted SNM: Squared Errors, post-TAVR - pre-TAVR");
% saveas(fig_se, strcat(homedir,simdir,"/","SE.fig"));
% saveas(fig_se, strcat(homedir,simdir,"/","SE.png"));

for col = 1:size(SE_mean, 2)
scatter(SE_mean(:,col), SE_nofit(:,col),20,"filled");
hold on
end
axlim = max([max(get(gca, "XLim")) max(get(gca, "YLim"))]);
line([0 axlim],[0 axlim]);
xlabel("SE: Fitted Model");
ylabel("SE: Default Model");
line([0 axlim],[0 axlim],"Color",[.5 .5 .5],"LineStyle","--");
legend(output_scaled.Properties.RowNames, "Location", "southeast");
legend(ydata, "Location", "southeast");
grid on

%% Output-Clinical Visualization
% output-clinical data correlations: Pearson correlation heatmaps
fitmap = brewermap(31,'*RdYlBu');

% figure("Position", [300 200 660 550]);
figure("Position", [300 200 400 550]);

% axes("Position", [0.01 0.24 0.08 0.62]);
axes("Position", [0.74 0.24 0.08 0.705]);
[dend_row,~,outperm] = dendrogram(tree_rho_rows, 'Orientation', 'right');
set(dend_row, "Color", "black");
set(gca, "Visible", "off");

% axes("Position", [0.19 0.86 0.3 0.08]);
axes("Position", [0.20 0.94 0.55 0.05]);
[dend_col,~,clinperm] = dendrogram(tree_rho_cols);
set(dend_col, "Color", "black");
set(gca, "Visible", "off");

% axes("Position", [0.2 0.25 0.28 0.6]);
% colormap(fitmap);
% [xgrid, ygrid] = meshgrid(1:size(expt_rho,2), 1:size(expt_rho,1));
% fig_exptrho = scatter(reshape(xgrid,1,[]), ...
%                     reshape(ygrid,1,[]), ...
%                     -20*log10(reshape(expt_p(outperm,clinperm),[],1))+5, ...
%                     reshape(expt_rho(outperm,clinperm),[],1), ...
%                     "filled");
% xlim([0.5 size(expt_rho,2)+0.5]);     ylim([0.5 size(expt_rho,1)+0.5]);
% xticks(1:size(expt_rho,2));           yticks(1:size(expt_rho,1));
% xticklabels(clin_vars_reduced(clinperm));
% yticklabels(speciesNames_out(outperm));
% xtickangle(90);

% fig_exptrho = heatmap(expt_rho(flip(outperm),clinperm));
% set(fig_exptrho,    "XData", clin_vars_reduced(clinperm), ...
%                     "YData", flip(speciesNames_out(outperm)), ...
%                     "Colormap", fitmap, ...
%                     "ColorbarVisible", "off");

% axes("Position", [0.55 0.25 0.28 0.6]);
axes("Position", [0.22 0.25 0.51 0.68]);
colormap(fitmap);
[xgrid, ygrid] = meshgrid(1:size(y_rho_allreps,2), 1:size(y_rho_allreps,1));
fig_yrho = scatter(reshape(xgrid,1,[]), ...
                    reshape(ygrid,1,[]), ...
                    -40*log10(reshape(y_p_allreps(outperm,clinperm),[],1)), ...
                    reshape(y_rho_allreps(outperm,clinperm),[],1), ...
                    "filled");
xlim([0.5 size(y_rho_allreps,2)+0.5]);     ylim([0.5 size(y_rho_allreps,1)+0.5]);
xticks(1:size(y_rho_allreps,2));           yticks(1:size(y_rho_allreps,1));
xticklabels(clin_vars_reduced(clinperm));
% yticklabels(repmat("",1,length(outperm)));
yticklabels(speciesNames_out(outperm));
xtickangle(90);

% fig_yrho = heatmap(y_rho(flip(outperm),clinperm));
% set(fig_yrho,   "XData", clin_vars_reduced(clinperm), ...
%                 "YDisplayLabels", nan(length(outperm),1), ...
%                 "ColorLimits", fig_exptrho.ColorLimits, ...
%                 "Colormap", fitmap, ...
%                 "ColorbarVisible", "off");
            
% axes("Position", [0.55 0.86 0.13 0.08]);
axes("Position", [0.8 0.3 0.2 0.1]);
colormap(fitmap); cp = colorbar;
set(gca,    "Visible", "off", ...
            "CLim", [min(fig_yrho.CData) max(fig_yrho.CData)]);
set(cp,     "Position", [0.8 0.45 0.03 0.1], ...
            "AxisLocation", "in", ...
            "FontSize", 8);
cp.Label.String = "Pearson r";

pscale = [0.5 0.1 0.05 0.01];
% axes("Position", [0.72 0.88 0.1 0.02]);
axes("Position", [0.8 0.3 0.03 0.1]);
scatter(repelem(1,1,4), 1:4, flip(-40*log10(pscale)), 'black', 'filled');
% set(gca, "XAxisLocation", "top", "YColor", "none", "XTickLabelRotation", 30);
set(gca, "YAxisLocation", "right", "XColor", "none");
yticks(1:4);
yticklabels(flip(pscale));
ylabel("p-value");

% figure("Position", [300 200 660 550]);
% axes("Position", [0.2 0.25 0.28 0.6]);
% % figure("Position",[650,320,400,420]);
% fig_ydelta_p = heatmap(expt_p(flip(outperm),clinperm), ...
%     "YData",flip(speciesNames_out(outperm)), ...
%     "XData",clin_vars_reduced(clinperm), ...
%     "Colormap", pink, ...
%     "ColorLimits", [0 0.5], ...
%     "Title", "Experimental Data", ...
%     "ColorBarVisible", "off");
% axes("Position", [0.55 0.25 0.28 0.6]);
% % figure("Position",[950,320,400,420]);
% fig_exptdeltadelta_p = heatmap(y_p(flip(outperm),clinperm), ...
%     "YDisplayLabels", nan(length(outperm),1), ...
%     "XData",clin_vars_reduced(clinperm), ...
%     "Colormap", pink, ...
%     "ColorLimits", [0 0.5], ...
%     "Title", "Fitted SNM Predictions");
% 
% comparemap = brewermap(size(y_rho, 2), 'Spectral');
% figure("Position", [200 200 700 275]);
% subplot(1,2,1);
% for col = 1:size(y_rho, 2)
% scatter(expt_rho(:,col), y_rho_nofit(:,col),20,comparemap(col,:),"filled");
% hold on
% end
% line([-1 1],[-1 1],"Color",[.5 .5 .5],"LineStyle","--");
% xlabel("Pearson r (Experimental Data)");
% ylabel("Pearson r (Model Predictions)");
% title("Default Model");
% grid on; hold off;
% 
% subplot(1,2,2);
% for col = 1:size(y_rho, 2)
% scatter(expt_rho(:,col), y_rho(:,col),20,comparemap(col,:),"filled");
% hold on
% end
% line([-1 1],[-1 1],"Color",[.5 .5 .5],"LineStyle","--");
% xlabel("Pearson r (Experimental Data)");
% title("Fitted Model");
% grid on; hold off;

%%
% MCC-clinical correlations
[~, idx] = sort(mcs_p_allreps);

figure;
subplot(2,1,1);
bar(mcs_rho_allreps(idx));
set(gca,    "XTickLabel", repmat([],1,length(xticks)), ...
            "Position", [0.13 0.6 0.775 0.3]);
ylabel("Spearman \rho");

subplot(2,1,2);
bar(mcs_p_allreps(idx));
set(gca,    "XTickLabel", clin_vars_reduced(idx), ...
            "XTickLabelRotation", 45, ...
            "Position", [0.13 0.25 0.775 0.3]);
ylabel("p-value");
yline(0.05, '--', 'p=0.05', 'LabelHorizontalAlignment', 'left');
sgtitle("Spearman Correlations: Matrix Content Score");

mcs_lmcols = find(mcs_p_allreps <= 0.1);
if ~isempty(mcs_lmcols)
    for col_mcc = 1:length(mcs_lmcols)
        mdl = fitlm(mcs_allreps, clin_mat_allreps(:,mcs_lmcols(col_mcc)));

        figure("Position",[200 100 400 300]);
        mdl.plot;
        xlabel('MCS', 'Interpreter', 'tex');
        ylabel(clin_numeric.Properties.VariableNames(mcs_lmcols(col_mcc)+2));
        annotation('textbox', [.65 .8 .2 .1], ...
            'String', "\rho = " + string(round(mcs_rho_allreps(mcs_lmcols(col_mcc)), 2, "significant")) ...
            + newline + "R^2 = " + string(round(mdl.Rsquared.Ordinary, 2, "significant")) ...
            + newline + "p = " + string(round(mcs_p_allreps(mcs_lmcols(col_mcc)), 2, "significant")), ...
            'FitBoxToText', 'on', ...
            'EdgeColor', 'none');
        legend('off');
        title(strcat('Linear Model: ',clin_numeric.Properties.VariableNames{mcs_lmcols(col_mcc)+2}), 'Interpreter', 'tex');
    end
end

%% Intermediate-Output Visualization
% Intermediate-output heatmaps: Pearson correlation
fitmap = brewermap(31,'*PiYG');
figure("Position", [300 -200 550 875]);

axes("Position", [0.82 0.12 0.05 0.8175]);    % y-axis dendrogram: right
[dend_row,~,medperm] = dendrogram(tree_med_rows, 0, 'Orientation', 'right');
set(dend_row, "Color", "black");
set(gca, "Visible", "off");

axes("Position", [0.12 0.94 0.7 0.05]);    % x-axis dendrogram: top
[dend_col,~,outperm] = dendrogram(tree_med_cols, 0);
set(dend_col, "Color", "black");
set(gca, "Visible", "off");

axes("Position", [0.13 0.13 0.68 0.8]);
colormap(fitmap);
[xgrid, ygrid] = meshgrid(1:size(y_rho_med,2), 1:size(y_rho_med,1));
fig_medrho = scatter(reshape(xgrid,1,[]), ...
                    reshape(ygrid,1,[]), ...
                    -15*log10(reshape(y_fdr_med(medperm,outperm),[],1)), ...
                    reshape(y_rho_med(medperm,outperm),[],1), ...
                    "filled");
set(gca, "FontSize", 8);
xlim([0.5 size(y_rho_med,2)+0.5]);     ylim([0.5 size(y_rho_med,1)+0.5]);
xticks(1:size(y_rho_med,2));           yticks(1:size(y_rho_med,1));
xticklabels(outputs_all(outperm));
yticklabels(speciesNames(meds_idx(medperm)));
xtickangle(90);
            
axes("Position", [0.8 0.15 0.2 0.1]);
colormap(fitmap); cr = colorbar;
set(gca,    "Visible", "off", ...
            "CLim", [min(fig_medrho.CData) max(fig_medrho.CData)]);
set(cr,     "Position", [0.86 0.15 0.03 0.1], ...
            "AxisLocation", "in");
cr.Label.String = "Pearson r";


% % Intermediate-Output Heatmaps: goodness of fit
% rsquaredmap = brewermap(100,'Oranges');
% figure("Position", [300 -200 480 875]);
% 
% axes("Position", [0.81 0.1 0.08 0.82]);
% [dend_row,~,medperm] = dendrogram(tree_med_rows, 0, 'Orientation', 'right');
% set(dend_row, "Color", "black");
% set(gca, "Visible", "off");
% 
% axes("Position", [0.15 0.93 0.65 0.07]);
% [dend_col,~,outperm] = dendrogram(tree_med_cols, 0);
% set(dend_col, "Color", "black");
% set(gca, "Visible", "off");
% 
% axes("Position", [0.15 0.1 0.65 0.82]);
% fig_medrho = heatmap(y_fit_med_r(flip(medperm),outperm));
% set(fig_medrho, "XData", outputs_all(outperm), ...
%                 "YData", speciesNames(meds_idx(flip(medperm))), ...
%                 "Colormap", rsquaredmap, ...
%                 "ColorbarVisible", "off", ...
%                 "FontSize", 8);
%             
% axes("Position", [0.8 0.1 0.2 0.1]);
% colormap(rsquaredmap); cr = colorbar;
% set(gca,    "Visible", "off", ...
%             "CLim", fig_medrho.ColorLimits);
% set(cr,     "Position", [0.86 0.1 0.033 0.1], ...
%             "AxisLocation", "in");
% cr.Label.String = "Correlation R^2";
% 
% % plot fdr-adjusted p-values
% figure("Position", [300 -200 480 875]);
% axes("Position", [0.15 0.1 0.65 0.82]);
% fig_medp = heatmap(y_fit_med_fdr(flip(medperm),outperm));
% set(fig_medp,   "XData", outputs_all(outperm), ...
%                 "YData", speciesNames(meds_idx(flip(medperm))), ...
%                 "Colormap", pink, ...
%                 "ColorbarVisible", "off", ...
%                 "ColorLimits", [0 0.1], ...
%                 "FontSize", 8);
% 
% axes("Position", [0.8 0.1 0.2 0.1]);
% colormap(pink); cp = colorbar;
% set(gca,    "Visible", "off", ...
%             "CLim", fig_medp.ColorLimits);
% set(cp,     "Position", [0.86 0.1 0.033 0.1], ...
%             "AxisLocation", "in");
% cp.Label.String = "FDR-adjusted p-value";