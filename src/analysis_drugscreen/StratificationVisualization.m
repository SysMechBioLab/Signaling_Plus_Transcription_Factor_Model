%% Stratification Analysis/Visualization: SNM Ensemble Results
% 12.22.2020 JR

%% *Load Data*
% load in filepaths
homedir = "E:/Research/Aim3/ModelExpansion/";
datadir = strcat(homedir, "ModelFitting_Ensemble/1_2_rev5_gender/");
colordir = "E:/Research/Aim2/BrewerMap-master/";
exptdir = "E:/Research/Aim2/ModelFitting/ModelFitting_AnsethData/data_processed_current/";
% exptdir_raw = "E:/Research/Aim2/ModelFitting/ModelFitting_AnsethData/data_raw/";
addpath(datadir, colordir, exptdir);

%%
% load d(Activity) array: ga_sensitivity_results.mat

filename = "ga_stratification_results.mat";
load(filename);
y_perturb_single = y_full_all;      % perturbs x nodes x patients x runs

filename = "ga_stratification_doseresponse.mat";          % TAVR_A-H
% filename = "ga_stratification_doseresponse_protonly.mat";   % TAVR_I-L
load(filename);
if exist('y_perturb_prot', 'var')
    y_perturb_doseresponse = y_perturb_prot;
    clearvars y_perturb_prot
else
    y_perturb_doseresponse = y_perturb_all;      % perturbs x nodes x doses x patients x runs
    clearvars y_perturb_all
end

filename = "ga_speciesNames.mat";
load(filename);

filename = "ga_influence_topnodes.mat";
load(filename);

% load experimental data
exptin = 'trainingdata_inputs_Anseth_11232020.mat';       % TAVR_A-H
% exptin = 'testingdata_inputs_Anseth_12312020.mat';          % TAVR_I-L
exptout = 'trainingdata_outputs_Anseth_12042020.mat';
load(exptin);
load(exptout);
if exist('input_scaled_protonly', 'var')
    input_scaled = input_scaled_protonly;
end
speciesNames_out = output_scaled.Properties.VariableNames;
ydata = replace(input_scaled.Properties.RowNames, "_", "-");


%% *Process Data*
% Calculate sensitivity/influence: distributions across runs

mean_perturb_single = mean(y_perturb_single, [3 4]);                % perturb x node
std_perturb_single = std(mean(y_perturb_single, 4), [], 3);         % perturb x node
all_perturb_single = permute(mean(y_perturb_single, 4), [3 2 1]);   % patient x node x perturb (to match expt)

mean_perturb_doseresponse = mean(y_perturb_doseresponse, [4 5]);    % perturb x node x dose
patients_perturb_doseresponse = mean(y_perturb_doseresponse, 5);    % perturb x node x dose x patient
mean_out_single = mean(y_out_all, 4);


outputs_other = {'EDAFN', 'proMMP9', 'proMMP3', 'proMMP8', 'proMMP12', 'thrombospondin4'};
outputs_all = speciesNames_out;
outputs_all(length(outputs_all)+1:length(outputs_all)+length(outputs_other)) = outputs_other;
[~,outputs_idx] = ismember(outputs_all, speciesNames);

mean_perturb_out = mean_perturb_single(:,outputs_idx);              % outputs only
std_perturb_out = std_perturb_single(:,outputs_idx);

% Isolate perturb/output of interest (ex. STAT/periostin)
speciesPerturb = "ETAR";
speciesOutput = "proCI";
idx_perturb = find(strcmp(speciesNames(influence_idx_toporder), speciesPerturb));
idx_output = outputs_idx(strcmp(speciesNames(outputs_idx), speciesOutput));

% Calculate descriptive statistics
model_perturb_doseresponse = permute(y_perturb_doseresponse(idx_perturb,idx_output,:,:,:), [3 4 5 1 2]);
model_perturb_doseresponse_mean = mean(model_perturb_doseresponse, 3);
model_perturb_doseresponse_delta = model_perturb_doseresponse_mean - model_perturb_doseresponse_mean(1,:);

%% *Visualize Data*
% clustermap: influential nodes vs. output ndoes
tree_perturb_rows = linkage(mean_perturb_out, 'average', 'euclidean');
tree_perturb_cols = linkage(mean_perturb_out.', 'average', 'euclidean');

% construct clustermaps (rho/p)
fitmap = brewermap(31,'Reds');

figure("Position", [300 200 650 350]);
axes("Position", [0.81 0.22 0.08 0.67]);
[dend_row,~,inflperm] = dendrogram(tree_perturb_rows, 'Orientation', 'right');
set(dend_row, "Color", "black");
set(gca, "Visible", "off");

axes("Position", [0.19 0.89 0.62 0.08]);
[dend_col,~,deltaperm] = dendrogram(tree_perturb_cols);
set(dend_col, "Color", "black");
set(gca, "Visible", "off");

axes("Position", [0.2 0.23 0.6 0.65]);
fig_perturb = heatmap(std_perturb_out(flip(inflperm), deltaperm), ...
                        "Colormap", fitmap, ...
                        "XData", speciesNames(outputs_idx(deltaperm)), ...
                        "YData", speciesNames(influence_idx_toporder(flip(inflperm))), ...
                        "ColorbarVisible", "off");

axes("Position", [0.85 0.23 0.1 0.2]);
colormap(fitmap); cr = colorbar;
set(gca,    "Visible", "off", ...
            "CLim", fig_perturb.ColorLimits);
set(cr,     "Position", [0.85 0.23 0.02 0.15], ...
            "AxisLocation", "in");
cr.Label.String = "\DeltaActivity";
cr.Label.Interpreter = "tex";

%%
% Dose-response curves: isoloated perturb/output of interest
xdata = repmat(0:0.2:0.8, size(model_perturb_doseresponse, 2)/2, 1).';
ysplit = split(ydata,"-");
barmap = brewermap(size(model_perturb_doseresponse, 2)/2, 'Set1');

figure("Position", [300 300 650 275]);
figtiles = tiledlayout(1,2);
nexttile;
plot(xdata, ...
        model_perturb_doseresponse_delta(:,9:end), ...  % change for TAVR_A-H/I-L
        "Marker", "o", ...
        "LineWidth", 1.5);
set(gca, "ColorOrder", barmap);
title("pre-TAVR Conditions");
nexttile;
plot(xdata, ...
        model_perturb_doseresponse_delta(:,1:8), ...    % change for TAVR_A-H/I-L
        "Marker", "o", ...
        "LineWidth", 1.5);
set(gca, "ColorOrder", barmap);
title("post-TAVR Conditions");
legend(unique(ysplit(:,2)), "Location", "eastoutside");

linkaxes(findall(gcf, "type", "axes"));
xlabel(figtiles, strcat(speciesPerturb, " Inhibitor Dose"));
ylabel(figtiles, strcat("\Delta", speciesOutput, " Activity"));
figtiles.TileSpacing = "compact";

%% *Correlate drug response w/ exogenous inputs*
% use 'all_perturb_single' as basis
corrmap = brewermap(31, '*RdBu');
idx_output_suggested = [23 73 79 101];      % latentTGFB, TIMP1, proMMP9, osteopontin
idx_perturb_suggested = [3 4 5 6 8 11];     % gp130, NFKB, STAT, ROS, ETAR, TGFB1R

pre_perturb_single = all_perturb_single(9:end, idx_output_suggested, idx_perturb_suggested);
pre_input_scaled = input_scaled{9:end,:};

rho_perturb_single = zeros(length(idx_output_suggested), size(pre_input_scaled, 2), length(idx_perturb_suggested));
p_perturb_single = zeros(length(idx_output_suggested), size(pre_input_scaled, 2), length(idx_perturb_suggested));

for perturb = 1:length(idx_perturb_suggested)
    [rho, p] = corr(pre_perturb_single(:,:,perturb), pre_input_scaled);
    rho_perturb_single(:,:,perturb) = rho;
    p_perturb_single(:,:,perturb) = p;
end


figure;
for perturb = 1:length(idx_perturb_suggested)
    subplot(length(idx_perturb_suggested), 1, perturb);
    fig_corr = heatmap(rho_perturb_single(:,:,perturb));
    set(fig_corr,   "YData", speciesNames(idx_output_suggested), ...
                    "XData", input_scaled.Properties.VariableNames, ...
                    "Colormap", corrmap, ...
                    "ColorLimits", [-1 1]);
    title(speciesNames(influence_idx_toporder(idx_perturb_suggested(perturb))));
end
