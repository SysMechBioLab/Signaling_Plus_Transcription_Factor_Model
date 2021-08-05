function ModelFittingEnsemble(iteration, iteration_max)
% SNM Model Fitting: Genetic Algorithm
% 10.29.2019 JR
% edited 03.30.2020: changed to SNM 1_2_rev2 for fitness function
%                    imports expt data from tables to retain
%                    input/output/sample metadata + parse order
% edited 02.05.2020: added input data for all measured inputs (7 total)
% edited 12.18.2019: changed scaling of input/output data
% edited 12.16.2019: using Anseth lab's input/output data for fitting (MSE)
% edited 12.05.2019: added binning method to fitness fcn, added method arg
% edited 11.04.2019: added reduced model fitness fcn, added output fcn
% edited 11.04.2019: increased population size/max generations for HPC
% edited 11.05.2019: changed upper bounds for input rxns (w(1:11))
% edited 11.17.2019: changed fitness fcn to reflect delta(Activity) errors,
%                    and changed normalization of expt_out


% Purpose: approximate reaction weights (w) for signaling network
% model (SNM) using genetic algorithm (ga)
% 
% This script sets up ga function by defining a fitness function and
% dependencies (stimulus doses, experimental data), as well as other
% settings for the solver (# of generations, upper/lower bounds for (w),
% parallelization for HPC, and population/fitness score output.
% 
% Required files:
% Fitness function: snm_1_2_rev2_fitness.m
%                   -runs SNM for each (w) across 16 perturbations
%                   -calculates mean-squared error (MSE) of predictions
%                    compared to (normalized) experimental data
% Output function: snm_1_1_rev3_gaout.m
%                   -initializes arrays to store population (history) and
%                    MSE (scores)
%                   -saves every 10th population/MSE to arrays, and exports
%                    to .mat file

clust = parcluster;
rng('shuffle')      % for multi-session job submissions
% rng('default')    % for reproducibility

% set fitness fcn inputs: stimulus doses, experimental perturbation data
filepath = 'ClusterAnalysis_rev5.xlsx';
sheet = "reactions";
inputs_idx = [1 2 8 10 4 5 6];

[w_map, w_full] = mapClusters(filepath, sheet, inputs_idx);

exptin = 'trainingdata_inputs_Anseth_11232020.mat';
exptout = 'trainingdata_outputs_Anseth_12042020.mat';
gender_value = 0.1;
% method = 'MSE';

fitnessfcn = @(w) ModelFitness(w,w_map,w_full,exptin,exptout,inputs_idx,gender_value);

% set other ga inputs: length of w, lower/upper bounds
n_vals = length(unique(w_map));
lb = zeros(1,n_vals);
ub = ones(1,n_vals);

% set ga options: parallelization, max # of generations, population size,
% output function
parpool(clust, clust.NumWorkers);
options = optimoptions('ga','UseParallel', true, ...    
    'UseVectorized', false, ...
    'MaxGenerations', 100, ...
    'PopulationSize', 500, ...
    'Display', 'diagnose');%, ...
    % 'OutputFcn', @snm_1_1_rev3_gaout, ...
    % 'PlotFcn', {'gaplotbestf', 'gaplotgenealogy'});    % for real-time visualization

% run ga function: multiple times
% iteration = 0;
% iteration_max = 1;
% x_all = zeros(iteration_max, n_vals);
% fval_all = zeros(iteration_max, 1);

while iteration < iteration_max
    fprintf('\n<<<<<<<<< GA: Iteration %d >>>>>>>>>\n', iteration)
    [x,fval,~,~,population,scores] = ga(fitnessfcn,n_vals,[],[],[],[],lb,ub,[],[],options);

    % fprintf('<<<<<<<<<<<<< GA Finished: Iteration %d >>>>>>>>>>>>>\n', iteration)
    % fprintf('The number of generations was : %d\n', output.generations);
    fprintf('Best Objective Function Value : %g\n', fval);
    
    % save final solution for each iteration
    finalset = {x,fval,population,scores};
    savepath = strcat("ga_finalset_",string(iteration),".mat");
    save(savepath, "finalset");
    
    iteration = iteration + 1;
end

delete(gcp)
end
