%% Sensitivity Analysis/Visualization: SNM Ensemble Results
% 12.17.2020 JR

%% *Load Data*
% load in filepaths
homedir = "E:/Research/Aim3/ModelExpansion/";
datadir = strcat(homedir, "ModelFitting_Ensemble/1_2_rev5_gender/");
colordir = "E:/Research/Aim2/BrewerMap-master/";
exptdir = "E:/Research/Aim2/ModelFitting/ModelFitting_AnsethData/data_processed_current/";
% exptdir_raw = "D:/Research/Aim2/ModelFitting/ModelFitting_AnsethData/data_raw/";
addpath(datadir, colordir, exptdir);

%%
% load d(Activity) array: ga_sensitivity_results.mat

filename = "ga_sensitivity_results.mat";
load(filename);

filename = "ga_speciesNames.mat";
load(filename);

% load experimental data
exptin = 'trainingdata_inputs_Anseth_11232020.mat';
exptout = 'trainingdata_outputs_Anseth_12042020.mat';
load(exptin);
load(exptout);
speciesNames_out = output_scaled.Properties.VariableNames;
ydata = replace(output_scaled.Properties.RowNames, "_", "-");

% get all output idxs
outputs_other = {'EDAFN', 'proMMP9', 'proMMP3', 'proMMP8', 'proMMP12', 'thrombospondin4'};
outputs_all = speciesNames_out;
outputs_all(length(outputs_all)+1:length(outputs_all)+length(outputs_other)) = outputs_other;
[~,outputs_idx] = ismember(outputs_all, speciesNames);


% load fitted SNM d(activity) data (pre-TAVR vs. post-TAVR)
filename = "ga_aggregate_results.mat";
load(filename);
% aggregate fitted results across runs
y_out_mean = mean(y_out_all,3);
y_delta_mean = mean(y_delta_all,3);
SE_mean = mean(SE_all,3);

y_out_pre = y_out_mean(9:end,:);
y_out_post = y_out_mean(1:8,:);

%%
% calculate composite matrix content score
idx_pro = [1 3 4 5 10 12 16 17 19 23 24 29];
idx_anti = [6 7 13 14 15 25 26 27 28];
idx_inhib = [18 21 22];

% rank patients (sort + re-order)
rank_out_pre = zeros(size(y_out_pre));
for col = 1:size(y_out_pre,2)
    [~,~,rank] = unique(y_out_pre(:,col));          % rank patient values (ascending)
    if any(idx_pro == col)
        rankrev = max(rank) - rank + 1;             % reverse ranks if in pro (descending)
        rank_out_pre(:,col) = rankrev / length(rankrev);
    elseif any(idx_inhib == col)
        rankrev = max(rank) - rank + 1;             % reverse ranks if in inhib (descending)
        rank_out_pre(:,col) = rankrev / length(rankrev);
    else
        rank_out_pre(:,col) = rank / length(rank);      % use ascending for anti
    end
end
mcs_pre = sum(rank_out_pre,2);

rank_out_post = zeros(size(y_out_post));
for col = 1:size(y_out_post,2)
    [~,~,rank] = unique(y_out_post(:,col));          % rank patient values (ascending)
    if any(idx_pro == col)
        rankrev = max(rank) - rank + 1;             % reverse ranks if in pro (descending)
        rank_out_post(:,col) = rankrev / length(rankrev);
    elseif any(idx_inhib == col)
        rankrev = max(rank) - rank + 1;             % reverse ranks if in inhib (descending)
        rank_out_post(:,col) = rankrev / length(rankrev);
    else
        rank_out_post(:,col) = rank / length(rank);      % use ascending for anti
    end
end
mcs_post = sum(rank_out_post,2);

mcs_delta = mcs_post - mcs_pre;

% mcc_pro = mean(y_delta_full_all(:, outputs_idx(idx_pro)), [2 3]);
% mcc_anti = mean(y_delta_full_all(:, outputs_idx(idx_anti)), [2 3]);
% mcc_inhib = mean(y_delta_full_all(:, outputs_idx(idx_inhib)), [2 3]);
% mcc = mcc_pro - mcc_anti + mcc_inhib;

adjustPValues = @(x) reshape(mafdr(reshape(x, [], 1), "BHFDR", "true"), size(x,1), []);

%% *Calculate Sensitivity Metrics*
% Calculate sensitivity/influence: distributions across runs

sensitivity_all = permute(sum(abs(y_sens_all), 1), [2 3 4 1]);   % sens x patient x run
influence_all = permute(sum(abs(y_sens_all), 2), [1 3 4 2]);     % infl x patient x run

sensitivity_mean = mean(sensitivity_all, 3);                % sens_mean x patient
influence_mean = mean(influence_all, 3);                    % infl_mean x patient

[~,sensitivity_idx] = sort(mean(sensitivity_mean, 2), 'descend');
[~,influence_idx] = sort(mean(influence_mean, 2), 'descend');

% find top nodes for each condition
[sensitivity_sort_indiv, sensitivity_idx_indiv] = sort(sensitivity_mean, 'descend');
[influence_sort_indiv, influence_idx_indiv] = sort(influence_mean, 'descend');
sensitivity_idx_top = unique(sensitivity_idx_indiv(1:10,:));
influence_idx_top = unique(influence_idx_indiv(1:10,:));

% sort top nodes for each condition (by first pre-TAVR condition)
[~,sensitivity_idx_topsorted] = sort(sensitivity_mean(sensitivity_idx_top, 9), 'descend');
[~,influence_idx_topsorted] = sort(influence_mean(influence_idx_top, 9), 'descend');

sensitivity_idx_toporder = sensitivity_idx_top(sensitivity_idx_topsorted);
influence_idx_toporder = influence_idx_top(influence_idx_topsorted);

sensitivity_topsorted_indiv = sensitivity_mean(sensitivity_idx_toporder, :);
influence_topsorted_indiv = influence_mean(influence_idx_toporder, :);

% calculate number of times each node appears in top 10
% calculate average rank of each node
influence_topranks = zeros(length(influence_idx_top), 1);
for node = 1:length(influence_idx_top)
    [noderows, ~] = find(influence_idx_indiv == influence_idx_top(node));
    meanrank = mean(noderows);
    influence_topranks(node) = meanrank;
end

sensitivity_topranks = zeros(length(sensitivity_idx_top), 1);
for node = 1:length(sensitivity_idx_top)
    [noderows, ~] = find(sensitivity_idx_indiv == sensitivity_idx_top(node));
    meanrank = mean(noderows);
    sensitivity_topranks(node) = meanrank;
end

%% *Visualize Results*
% sensitivity + influence distributions
boxmap = brewermap(2, 'Set2');
ysplit = split(ydata, "-");

% sensitivity distributions
% for patient = 1:(size(sensitivity_all, 2)/2)
for patient = 5
    figure("Position", [100 100 700 230]);
    boxplot(permute(sensitivity_all(sensitivity_idx, patient+8, :), [3 1 2]), ...
            "PlotStyle", "compact", "Colors", boxmap(1,:), "MedianStyle", "line");
    hold on
    boxplot(permute(sensitivity_all(sensitivity_idx, patient, :), [3 1 2]), ...
            "PlotStyle", "compact", "Colors", boxmap(2,:), "MedianStyle", "line");

    xticklabels(repmat([], 1, length(sensitivity_idx)));
    xlabel("Measured Node");
    ylim([0 max(sensitivity_all, [], "all")]);
    ylabel("Sensitivity");
    title(replace(ysplit(patient, 2), "TAVR", "TAVR "));
    hold off
end

% influence distributions
boxmap = brewermap(2, 'Set1');
% for patient = 1:(size(influence_all, 2)/2)
for patient = 5
    figure("Position", [100 100 700 230]);
    boxplot(permute(influence_all(influence_idx, patient+8, :), [3 1 2]), ...
            "PlotStyle", "compact", "Colors", boxmap(1,:), "MedianStyle", "line");
    hold on
    boxplot(permute(influence_all(influence_idx, patient, :), [3 1 2]), ...
            "PlotStyle", "compact", "Colors", boxmap(2,:), "MedianStyle", "line");

    xticklabels(repmat([], 1, length(influence_idx)));
    xlabel("Knock-down Node");
    ylim([0 max(influence_all, [], "all")]);
    ylabel("Influence");
    title(replace(ysplit(patient, 2), "TAVR", "TAVR "));
    hold off
end

%%
% top sensitive/influential nodes: per-patient
barmap1 = brewermap(size(influence_all, 2)/2, 'Pastel1');
barmap2 = brewermap(size(influence_all, 2)/2, 'Set1');

figure("Position", [450 80 650 400]);
figtiles = tiledlayout(1,size(influence_all, 2)/2);
for patient = 1:(size(influence_all, 2)/2)
    nexttile;
    % colormap(gca, [barmap1(patient, :); barmap2(patient, :)]);
    figbar = barh(1:size(influence_idx_top), ...
                    [flipud(influence_topsorted_indiv(:, patient)) flipud(influence_topsorted_indiv(:, patient+8))], ...
                    0.95);   % [post-TAVR pre-TAVR]
    figbar(1).FaceColor = barmap1(patient, :);  % post-TAVR (pastel)
    figbar(2).FaceColor = barmap2(patient, :);  % pre-TAVR (accent)
    yticks(1:size(influence_idx_top));
    if patient == 1
        yticklabels(speciesNames(flipud(influence_idx_top(influence_idx_topsorted))));
    else
        yticklabels(repmat([], 1, length(influence_idx_top)));
    end
    title(replace(ysplit(patient, 2), "TAVR", "TAVR "));
end
linkaxes(findall(gcf, "type", "axes"));
xlabel(figtiles, "Total Influence");
ylabel(figtiles, "Top-Ranked Node");
figtiles.TileSpacing = "compact";


figure("Position", [450 80 650 400]);
figtiles = tiledlayout(1,size(sensitivity_all, 2)/2);
for patient = 1:(size(sensitivity_all, 2)/2)
    nexttile;
    % colormap(gca, [barmap1(patient, :); barmap2(patient, :)]);
    figbar = barh(1:size(sensitivity_idx_top), ...
                    [flipud(sensitivity_topsorted_indiv(:, patient)) flipud(sensitivity_topsorted_indiv(:, patient+8))], ...
                    0.95);   % [post-TAVR pre-TAVR]
    figbar(1).FaceColor = barmap1(patient, :);      % post-TAVR (pastel)
    figbar(2).FaceColor = barmap2(patient, :);      % pre-TAVR (accent)
    yticks(1:size(sensitivity_idx_top));
    if patient == 1
        yticklabels(speciesNames(flipud(sensitivity_idx_top(sensitivity_idx_topsorted))));
    else
        yticklabels(repmat([], 1, length(sensitivity_idx_top)));
    end
    title(replace(ysplit(patient, 2), "TAVR", "TAVR "));
end
linkaxes(findall(gcf, "type", "axes"));
xlabel(figtiles, "Total Sensitivity");
ylabel(figtiles, "Top-Ranked Node");
figtiles.TileSpacing = "compact";

% CV levels for values
sensitivity_cv_pre = std(sensitivity_topsorted_indiv(:,9:end),[],2) ./ mean(sensitivity_topsorted_indiv(:,9:end),2);
sensitivity_cv_post = std(sensitivity_topsorted_indiv(:,1:8),[],2) ./ mean(sensitivity_topsorted_indiv(:,1:8),2);
influence_cv_pre = std(influence_topsorted_indiv(:,9:end),[],2) ./ mean(influence_topsorted_indiv(:,9:end),2);
influence_cv_post = std(influence_topsorted_indiv(:,1:8),[],2) ./ mean(influence_topsorted_indiv(:,1:8),2);

cvmap = brewermap(2,'Greys');
figure("Position", [450 80 200 400]);
figbar = barh(1:size(sensitivity_cv), ...
                [flip(sensitivity_cv_post) flip(sensitivity_cv_pre)], ...
                0.95);   % [post-TAVR pre-TAVR]
figbar(1).FaceColor = cvmap(1, :);      % post-TAVR (pastel)
figbar(2).FaceColor = cvmap(2, :);      % pre-TAVR (accent)
yticks(1:size(sensitivity_cv));
yticklabels(speciesNames(flipud(sensitivity_idx_top(sensitivity_idx_topsorted))));
xlabel("CV (Sensitivity)");

figure("Position", [450 80 200 400]);
figbar = barh(1:size(influence_cv), ...
                [flip(influence_cv_post) flip(influence_cv_pre)], ...
                0.95);   % [post-TAVR pre-TAVR]
figbar(1).FaceColor = cvmap(1, :);      % post-TAVR (pastel)
figbar(2).FaceColor = cvmap(2, :);      % pre-TAVR (accent)
yticks(1:size(influence_cv));
yticklabels(speciesNames(flipud(influence_idx_top(influence_idx_topsorted))));
xlabel("CV (Influence)");

%%
% Correlational analysis: sensitivity/influence vs. delta(activity)

[y_rho_infl, y_p_infl] = corr(influence_topsorted_indiv(:,9:end).', y_delta_mean);  % pre-TAVR sensitivity vs. delta(pre vs. post-TAVR)
y_fdr_infl = adjustPValues(y_p_infl);

% cluster correlation coeff's
tree_infl_rows = linkage(y_rho_infl, 'average', 'euclidean');
tree_infl_cols = linkage(y_rho_infl.', 'average', 'euclidean');

% construct clustermaps (rho/p)
fitmap = brewermap(31,'*RdYlBu');

figure("Position", [300 200 650 350]);
axes("Position", [0.81 0.22 0.08 0.67]);
[dend_row,~,inflperm] = dendrogram(tree_infl_rows, 'Orientation', 'right');
set(dend_row, "Color", "black");
set(gca, "Visible", "off");

axes("Position", [0.19 0.89 0.62 0.08]);
[dend_col,~,deltaperm] = dendrogram(tree_infl_cols);
set(dend_col, "Color", "black");
set(gca, "Visible", "off");

axes("Position", [0.2 0.23 0.6 0.65]);
colormap(fitmap);
[xgrid, ygrid] = meshgrid(1:size(y_rho_infl,2), 1:size(y_rho_infl,1));
fig_inflrho = scatter(reshape(xgrid,1,[]), ...
                    reshape(ygrid,1,[]), ...
                    -20*log10(reshape(y_p_infl(inflperm,deltaperm),[],1)), ...
                    reshape(y_rho_infl(inflperm,deltaperm),[],1), ...
                    "filled");
xlim([0.5 size(y_rho_infl,2)+0.5]);     ylim([0.5 size(y_rho_infl,1)+0.5]);
xticks(1:size(y_rho_infl,2));           yticks(1:size(y_rho_infl,1));
xticklabels(speciesNames(outputs_idx(deltaperm)));
yticklabels(speciesNames(influence_idx_toporder(inflperm)));
xtickangle(90);

axes("Position", [0.85 0.23 0.1 0.2]);
colormap(fitmap); cr = colorbar;
set(gca,    "Visible", "off", ...
            "CLim", [min(fig_inflrho.CData) max(fig_inflrho.CData)]);
set(cr,     "Position", [0.85 0.23 0.02 0.15], ...
            "AxisLocation", "in");
cr.Label.String = "Pearson r";

pscale = [0.1 0.05 0.01 0.001];
axes("Position", [0.85 0.4 0.02 0.1]);
scatter(repelem(1,1,4), 1:4, -20*log10(pscale), 'black', 'filled');
set(gca, "YAxisLocation", "right", "XColor", "none");
yticks(1:4);
yticklabels(pscale);
ylabel("p-value");


%%
% Correlational analysis: sensitivity/influence vs. MCC
[y_rho_mcs, y_p_mcs] = corr(influence_topsorted_indiv(:,9:end).', mcs_delta, ...
    "rows", "pairwise", "type", "Spearman");  % pre-TAVR sensitivity vs. mmcc
[~,idx_mcc] = sort(y_p_mcs);

figure;
subplot(2,1,1);
bar(y_rho_mcs(idx_mcc));
set(gca,    "XTick", 1:length(idx_mcc), ...
            "XTickLabel", repmat([],1,length(xticks)), ...
            "Position", [0.13 0.6 0.775 0.3]);
ylabel("Pearson r");
title("Pearson Correlations: pre-TAVR Influence vs. MCS");

subplot(2,1,2);
bar(y_p_mcs(idx_mcc));
set(gca,    "XTick", 1:length(idx_mcc), ...
            "XTickLabel", speciesNames(influence_idx_toporder(idx_mcc)), ...
            "XTickLabelRotation", 45, ...
            "Position", [0.13 0.25 0.775 0.3]);
xlabel("Top Influential Nodes");
ylabel("p-value");
yline(0.05, '--', 'p=0.05', 'LabelHorizontalAlignment', 'left');


%% *Save Data*
% Save lists of top sensitive/influential nodes

filename = "ga_sensitivity_topnodes.mat";
save(strcat(datadir, filename), "sensitivity_idx_toporder");
filename = "ga_influence_topnodes.mat";
save(strcat(datadir, filename), "influence_idx_toporder");
