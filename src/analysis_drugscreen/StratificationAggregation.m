% Simulation of patient-specific perturbations
% 10.26.2020 JR
% Notes:    -intendended for use on Palmetto HPC
%           -should be placed in 'scripts' directory
%           -should cp ga_*.mat files from 'results' directory in sh script

% load data
rng('shuffle');
clust = parcluster();
parpool(clust, 26);
homedir = pwd;
datadir = strcat(homedir, "\\data\\training_results\\ga_ensemble\\");
clustdir = strcat(homedir, "\\data\\training_datasets_Anseth\\");
utilsdir = strcat(homedir, "\\src\\utils\\");
addpath(datadir, clustdir, utilsdir)

filedelim = strcat(datadir, "ga_finalset*");
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


% Set static parameters for simulation
clusterpath = "ClusterAnalysis_rev6.xlsx";
sheet = "reactions";
inputs_idx = [1 2 8 10 4 5 6];
[w_map, w_full] = mapClusters(clusterpath, sheet, inputs_idx);

exptin = "trainingdata_inputs_Anseth_11232020.mat";   % TAVR_A-H
% exptin = "testingdata_inputs_Anseth_12312020.mat";      % TAVR_I-L
exptout = "trainingdata_outputs_Anseth_12042020.mat";
gender_value = 0.1;

tension = 0;

inflpath = "ga_influence_topnodes.mat";
if exist(inflpath, "file")
    load(inflpath);
else
    error("MyFunction:fileNotFound", ...
        "Error: cannot find file %s. Make sure that file is located in %s directory, or run src\\analysis_sensitivity\\SensitivityVisualization.m.", ...
        inflpath, datadir)
end
perturb_lvls = flip(0.2:0.2:0.8);

% run simulations for each set
y_perturb_prot = zeros(length(influence_idx_toporder), ...   % perturbs
    151, ...                              % nodes
    length(perturb_lvls)+1, ...           % doses
    8, ...                                % sera
    size(x_all,1));                       % parameter set

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
    
    [~,y_perturb] = StratificationSimulation_doseresponse(w, ...
        w_map, ...
        w_full, ...
        exptin, ...
        exptout, ...
        inputs_idx, ...
        gender_value, ...
        influence_idx_toporder, ...
        perturb_lvls);
    
    y_perturb_prot(:,:,:,:,run) = y_perturb;
    
end

w_nofit = ones(1, size(x_all, 2));
w_nofit(w_map(w_full<=11)) = 0.1;

[~,y_perturb_nofit] = StratificationSimulation_doseresponse(w_nofit, ...
    w_map, ...
    w_full, ...
    exptin, ...
    exptout, ...
    inputs_idx, ...
    gender_value, ...
    influence_idx_toporder, ...
    perturb_lvls);

% save(strcat(homedir, "\\data\\training_results\\ga_stratification_results.mat"), "y_perturb");     % if using StratificationSimulation_singledose.m
% save(strcat(homedir, "\\data\\training_results\\ga_stratification_results_nofit.mat"), "y_perturb_nofit");
save(strcat(homedir, "\\data\\training_results\\ga_stratification_doseresponse.mat"), "y_perturb");     % if using StratificationSImulation_doseresponse.m
save(strcat(homedir, "\\data\\training_results\\ga_stratification_doseresponse_nofit.mat"), "y_perturb_nofit");

delete(gcp)

