% Model Validation: aSMA measurements
% 11.16.2020 JR

% load in filepaths
homedir = "E:/Research/Aim3/ModelExpansion/";
datadir = strcat(homedir, "ModelFitting_Ensemble/1_2_rev5_gender/");
colordir = "E:/Research/Aim2/BrewerMap-master/";
exptdir = "E:/Research/Aim2/ModelFitting/ModelFitting_AnsethData/data_processed_current/";
exptdir_raw = "E:/Research/Aim2/ModelFitting/ModelFitting_AnsethData/data_raw/";
addpath(datadir, colordir, exptdir);



% load prediction data: low-tension (lt)
% raw order: post-TAVR A-H -> pre-TAVR A-H
filename = "ga_aggregate_results.mat";
load(filename);
y_out_lt = y_out_all;
y_delta_lt = y_delta_all;
y_out_lt_mean = mean(y_out_lt,3);   % across ensemble runs
y_delta_lt_mean = mean(y_delta_lt,3);

% load prediction data: high-tension (ht)
% raw order: post-TAVR A-H -> pre-TAVR A-H
% filename = "ga_aggregate_results_tension_nocap.mat";
% load(filename);
% y_out_ht = y_out_all;
% y_delta_ht = y_delta_all;
% y_out_ht_mean = mean(y_out_ht,3);       % across ensemble runs
% y_delta_ht_mean = mean(y_delta_ht,3);

% load prediction data: p38 inhib
% raw order: post-TAVR A-H -> pre-TAVR A-H
% filename = "ga_aggregate_results_p38.mat";
% load(filename);
% y_out_p38 = y_out_all;

% load prediction data: transition
% raw order: (pre->pre pre->post) A-H
% filename = "ga_aggregate_results_transition.mat";
% load(filename);
% y_out_transition_lt = y_out_all;
% y_delta_transition_lt = y_delta_all;
% filename = "ga_aggregate_results_transition_tension.mat";
% load(filename);
% y_out_transition_ht = y_out_all;
% y_delta_transition_ht = y_delta_all;



% load experimental data
% raw order: post-TAVR A-H -> pre-TAVR A-H
exptin = 'trainingdata_inputs_Anseth_11232020.mat';
exptout = 'trainingdata_outputs_Anseth_12042020.mat';
load(exptin);
load(exptout);
speciesNames = output_scaled.Properties.VariableNames;
expt_fc = output_scaled{1:8,:} ./ output_scaled{9:end,:};

expt_1f_norm = reshape(fliplr(reshape(output_scaled{:,1},[],2)).',[],1);
expt_1g_norm = output_scaled{1:8,1} - output_scaled{9:end,1};



% load in vitro data
filename = "aav3233_Data_file_S2.xlsx";
opts = detectImportOptions(strcat(exptdir_raw,filename),"Sheet","Fig1F");
expt_1f = readtable(strcat(exptdir_raw,filename),opts);
expt_1f_mean = mean(expt_1f{:,:}, 1);

opts = detectImportOptions(strcat(exptdir_raw,filename),"Sheet","Fig1G");
expt_1g = readtable(strcat(exptdir_raw,filename),opts);
expt_1g_mean = mean(expt_1g{:,:}, 1);

% opts = detectImportOptions(strcat(exptdir_raw,filename),"Sheet","Fig4B");
% expt_4b_raw = readtable(strcat(exptdir_raw,filename),opts);
% expt_4b = reshape(expt_4b_raw{:,2:end}.',[],8);
% expt_4b_mean = mean(expt_4b, 1);

% opts = detectImportOptions(strcat(exptdir_raw,filename),"Sheet","Fig5C");
% expt_5c = readtable(strcat(exptdir_raw,filename),opts);
% expt_5c_mean = mean(expt_5c{:,:}, 1);

% opts = detectImportOptions(strcat(exptdir_raw,filename),"Sheet","Fig5D");
% expt_5d = readtable(strcat(exptdir_raw,filename),opts);
% expt_5d_mean = mean(expt_5d{:,:}, 1);

% opts = detectImportOptions(strcat(exptdir_raw,filename),"Sheet","Fig5G");
% expt_5g = readtable(strcat(exptdir_raw,filename),opts);
% expt_5g_mean = mean(expt_5g{:,:}, 1);

% opts = detectImportOptions(strcat(exptdir_raw,filename),"Sheet","Fig5H");
% expt_5h = readtable(strcat(exptdir_raw,filename),opts);
% expt_5h_mean = mean(expt_5h{:,:}, 1);



% process prediction data: reshape to match experimental
y_1f = reshape(permute(flip(reshape(y_out_lt(:,1,:),[],2,size(y_out_lt,3)),2),[2,1,3]),[],1,size(y_out_lt,3));
y_1f_mean = reshape(fliplr(reshape(y_out_lt_mean(:,1),[],2)).',[],1);

y_1g = y_delta_lt(:,1,:);
y_1g_mean = y_delta_lt_mean(:,1);

% rows_4b = ["Pre_TAVRA","Pre_TAVRB","Pre_TAVRE","Pre_TAVRF"];
% idx_4b = find(contains(output_scaled.Properties.RowNames, rows_4b));
% y_4b_raw = [y_out_lt(idx_4b,1,:); y_out_p38(idx_4b,1,1:size(y_out_lt,3))];
% y_4b = reshape(permute(flip(reshape(y_4b_raw(:,1,:),[],2,size(y_4b_raw,3)),2),[2,1,3]),[],1,size(y_4b_raw,3));
% y_4b_mean = mean(y_4b,3);

% y_5c = y_out_transition_lt(:,1,:);
% y_5c_mean = mean(y_5c,3);
% y_5d = y_delta_transition_lt(:,1,:);
% y_5d_mean = mean(y_5d,3);
% y_5g = y_out_transition_ht(:,1,:);
% y_5g_mean = mean(y_5g,3);
% y_5h = y_delta_transition_ht(:,1,:);
% y_5h_mean = mean(y_5h,3);


%%%% data visualization %%%%
% Fig. 1F (pre/post-TAVR only)
cmap_scatter = @bone;
cmap_bar = @lines;
xlabs = reshape(flipud(reshape(output_scaled.Properties.RowNames,[],2).'),[],1);
pos_1f = [250 450 1200 300];


fig_1f = makeValidationPlot(pos_1f, y_1f, y_1f_mean, ...
    expt_1f_norm, expt_1f{:,:}.', expt_1f_mean, ...
    cmap_scatter, cmap_bar, ...
    xlabs, true);

% Fig. 1G: (FC, post vs. pre-TAVR)
xsplit = split(xlabs,"_");
xlabs = unique(xsplit(:,2));
pos_1g = pos_1f;    pos_1g(3) = 1050;

fig_1g = makeValidationPlot(pos_1g, y_1g + 1, y_1g_mean + 1, ...
    expt_1g_norm + 1, expt_1g{:,:}.', expt_1g_mean, ...
    cmap_scatter, cmap_bar, ...
    xlabs, false);

% Fig. 4B: (pre-TAVR only +/- p38 inhib)
xlabs = repelem(rows_4b,2) + repmat(["","_ko"],1,4);
pos_4b = pos_1g;    pos_4b(3) = 900;

fig_4b = makeValidationPlot(pos_4b, y_4b, y_4b_mean, ...
    [], expt_4b.', expt_4b_mean, ...
    cmap_scatter, cmap_bar, ...
    xlabs, true);

% Fig. 5C: (transition pre -> pre/post-TAVR, low tension)
xsplit = split(expt_5c.Properties.VariableNames, "_");
xlabs = xsplit(:,:,2) + ": " + xsplit(:,:,1) + "->" + xsplit(:,:,3);
pos_5c = pos_4b;

fig_5c = makeValidationPlot(pos_5c, y_5c, y_5c_mean, ...
    [], expt_5c{:,:}.', expt_5c_mean, ...
    cmap_scatter, cmap_bar, ...
    xlabs, true);

% Fig. 5D: (FC, transition pre -> pre/post-TAVR, low tension)
xlabs = unique(xsplit(:,:,2));
pos_5d = pos_5c;

fig_5d = makeValidationPlot(pos_5d, y_5d, y_5d_mean, ...
    [], expt_5d{:,:}.', expt_5d_mean, ...
    cmap_scatter, cmap_bar, ...
    xlabs, false);

% Fig. 5G: (transition pre -> pre/post-TAVR, high tension)
xsplit = split(expt_5g.Properties.VariableNames, "_");
xlabs = xsplit(:,:,2) + ": " + xsplit(:,:,1) + "->" + xsplit(:,:,3);
pos_5g = pos_5d;

fig_5g = makeValidationPlot(pos_5g, y_5g, y_5g_mean, ...
    [], expt_5g{:,:}.', expt_5g_mean, ...
    cmap_scatter, cmap_bar, ...
    xlabs, true);

% Fig. 5H: (FC, transition pre -> pre/post-TAVR, high tension)
xlabs = unique(xsplit(:,:,2));
pos_5h = pos_5g;

fig_5h = makeValidationPlot(pos_5h, y_5h, y_5h_mean, ...
    [], expt_5h{:,:}.', expt_5h_mean, ...
    cmap_scatter, cmap_bar, ...
    xlabs, false);





function fig = makeValidationPlot(pos, y, y_mean, expt_norm, expt, expt_mean, cmap_scatter, cmap_bar, xlabs, isIndiv)

% create cmaps
cmap_scatter_adapted = cmap_scatter(size(y,3)+5);
if isIndiv
    cmap_bar_adapted = repelem(cmap_bar(size(y,1)/2),2,1);
else
    cmap_bar_adapted = cmap_bar(size(y,1));
end

fig = figure("Position",pos);
% plot 1 (required): simulation data
subplot(1,3,1);
for i = 1:size(y,3)
    scatter(1:size(y,1),y(:,1,i),10,cmap_scatter_adapted(i,:),'filled')
    hold on
end
y_bar = bar(y_mean(:,1),'FaceColor','flat','FaceAlpha',0.3);
y_bar.CData = cmap_bar_adapted;
if isIndiv
    ylab = "aSMA Activity (a.u.)";
else
    ylab = "\DeltaActivity (aSMA)";
end
ylabel(ylab);
title("Model Predictions");
set(gca, ...
    "XTick", 1:size(y,1), ...
    "XTickLabel", xlabs, ...
    "XTickLabelRotation", 45, ...
    "TickLabelInterpreter", "none", ...
    "FontSize", 9);
hold off

% plot 2 (optional): normalized experimental data (RNAseq)
if ~isempty(expt_norm)
    subplot(1,3,2)
    expt_bar = bar(expt_norm,'FaceColor','flat','FaceAlpha',0.3);
    expt_bar.CData = cmap_bar_adapted;
    if isIndiv
        ylab = "aSMA Expression (a.u.)";
    else
        ylab = "\DeltaExpression (aSMA)";
    end
    ylabel(ylab);
    title("Experimental Data (RNA-seq)");
    set(gca, ...
        "XTick", 1:size(expt_norm,1), ...
        "XTickLabel", xlabs, ...
        "XTickLabelRotation", 45, ...
        "TickLabelInterpreter", "none", ...
        "FontSize", 9);
    box off
end

% plot 3 (optional): raw experimental data (IF)
if ~isempty(expt) && ~isempty(expt_mean)
    if isempty(expt_norm)
        subplot(1,3,2);
    else
        subplot(1,3,3);
    end
    
    if ~isempty(expt)
        for i = 1:size(expt,2)
            scatter(1:size(expt,1),expt(:,i),10,cmap_scatter_adapted(i,:),'filled')
            hold on
        end
    end
    expt_bar = bar(expt_mean,'FaceColor','flat','FaceAlpha',0.3);
    expt_bar.CData = cmap_bar_adapted;
    if isIndiv
        ylab = "aSMA+ Cells (%)";
    else
        ylab = "Fold change (aSMA+ Cells)";
    end
    ylabel(ylab);
    title("Experimental Data (IF)");
    set(gca, ...
        "XTick", 1:size(expt,1), ...
        "XTickLabel", xlabs, ...
        "XTickLabelRotation", 45, ...
        "TickLabelInterpreter", "none", ...
        "FontSize", 9);
    hold off
end

end
