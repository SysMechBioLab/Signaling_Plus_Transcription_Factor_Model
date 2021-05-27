% Sensitivity analysis: post-fitting model analysis
% 12.14.2020 JR

% setup variables: 
%       w_idx:      indexes of input reactions to stimulate
%       stim:       float valued 0-1 to delineate stimulus of each input
%       iterations: integer for number of k-means clustering runs
%       exptin:     string for .mat file containing (normalized) inputs

rng('shuffle');
parpool

%% *Load Data*
% load in filepaths (if executing locally)
% homedir = "D:/Research/Aim3/ModelExpansion/";
% datadir = strcat(homedir, "ModelFitting_Ensemble/1_2_rev5_gender/");
% colordir = "D:/Research/Aim2/BrewerMap-master/";
% exptdir = "D:/Research/Aim2/ModelFitting/ModelFitting_AnsethData/data_processed_current/";
% exptdir_raw = "D:/Research/Aim2/ModelFitting/ModelFitting_AnsethData/data_raw/";
% addpath(datadir, colordir, exptdir);

%%
% load reaction weights
filedelim = "ga_finalset*";
datastr = dir(filedelim);
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
% Set static parameters for simulation
clusterpath = "ClusterAnalysis_rev5.xlsx";
sheet = "reactions";
inputs_idx = [1 2 8 10 4 5 6];

exptin = 'trainingdata_inputs_Anseth_11232020.mat';
exptout = 'trainingdata_outputs_Anseth_12042020.mat';
gender_value = 0.1;
tension = 0;


[w_map, w_full] = mapClusters(clusterpath, sheet, inputs_idx);


%%
% run simulations for each set
y_default_all = zeros(151,151,16,size(x_all, 1));      % dims: node_perturbed x node_measured x patient x parameter_set
y_perturb_all = zeros(151,151,16,size(x_all, 1));      % dims: node_perturbed x node_measured x patient x parameter_set
y_sens_all = zeros(151,151,16,size(x_all, 1));      % dims: node_perturbed x node_measured x patient x parameter_set
ymax_all = zeros(151,151,16,size(x_all, 1));      % dims: node_perturbed x node_measured x patient x parameter_set

parfor run = 1:size(x_all, 1)
    w = x_all(run,:);
    
    if tension ~= 0
        if w(w_map(w_full==3)) > 1 - tension
            w(w_map(w_full==3)) = w(w_map(w_full==3)) + tension;
            fprintf("tension = %g\n", w(w_map(w_full==3)));
        else
            w(w_map(w_full==3)) = 1;
            fprintf("high tension exceeds cap; set to 1\n");
        end
    else
        fprintf("tension = %g (default)\n", w(w_map(w_full==3)));
    end
    
    [~,y_default,y_perturb,y_sens,ymax_sens] = EnsembleSensitivity_rev5_gender(w,w_map,w_full,exptin,exptout,inputs_idx,gender_value);   
    
    y_default_all(:,:,:,run) = y_default;
    y_perturb_all(:,:,:,run) = y_perturb;
    y_sens_all(:,:,:,run) = y_sens;
    ymax_all(:,:,:,run) = ymax_sens;
end


delete(gcp);

save("ga_sensitivity_results.mat", "y_default_all", "y_perturb_all", "y_sens_all", "ymax_all");