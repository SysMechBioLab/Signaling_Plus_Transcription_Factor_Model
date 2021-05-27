% Aggregation of GA ensemble results
% 10.26.2020 JR
% Notes:    -intendended for use on Palmetto HPC
%           -should be placed in 'scripts' directory
%           -should cp ga_*.mat files from 'results' directory in sh script

% load data
% homedir = "D:/Research/Aim3/ModelExpansion/";
% datadir = strcat(homedir, "ModelFitting_Ensemble/1_2_rev5_scaling_inputs/");
% clustdir = "D:/Research/Aim2/ModelFitting/ModelFitting_AnsethData/data_processed_current/";
% addpath(clustdir)

% filedelim = strcat(datadir, "ga_finalset*");
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


% Set static parameters for simulation
clusterpath = "ClusterAnalysis_rev5.xlsx";
sheet = "reactions";
inputs_idx = [1 2 8 10 4 5 6];

exptin = 'testingdata_inputs_Anseth_11232020.mat';
exptout = 'trainingdata_outputs_Anseth_12042020.mat';
gender_value = 0.1;
% method = 'MSE';

tension = 0;
% perturb_node = 'p38';
% perturb_lvl = 0.1;

[w_map, w_full] = mapClusters(clusterpath, sheet, inputs_idx);

% remove patient-derived inputs (for rev5_scaling_inputs ONLY!!!)
% x_all_setonly = x_all;
% x_all_setonly(:, inputs_idx) = [];

% run simulations for each set
y_out_all = zeros(16,23,size(x_all,1));
y_delta_all = zeros(8,23,size(x_all,1));
SE_all = zeros(8,23,size(x_all,1));

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
    
    [~,y_out,y_full,y_delta,y_delta_full,SE] = EnsembleAggregation_fitness_rev5_gender(w,w_map,w_full,exptin,exptout,inputs_idx,gender_value);
    
    y_out_all(:,:,run) = y_out;
    y_full_all(:,:,run) = y_full;
    y_delta_all(:,:,run) = y_delta;
    y_delta_full_all(:,:,run) = y_delta_full;
    SE_all(:,:,run) = SE;
    
end

w_nofit = ones(1, size(x_all, 2));
w_nofit(w_map(w_full<=11)) = 0.1;

[speciesNames,y_out_nofit,y_full_nofit,y_delta_nofit,y_delta_full_nofit,SE_nofit] = EnsembleAggregation_fitness_rev5_gender(w_nofit,w_map,w_full,exptin,exptout,inputs_idx,gender_value);

save("ga_aggregate_results.mat", "y_out_all", "y_full_all", "y_delta_all", "y_delta_full_all", "SE_all");
save("ga_aggregate_nofit.mat", "y_out_nofit", "y_full_nofit", "y_delta_nofit", "y_delta_full_nofit", "SE_nofit");
save("ga_speciesNames.mat", "speciesNames");
