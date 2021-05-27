% Anseth data processing
% 12.24.2019 JR
% 
% Conversion of initial data (from GEO entry/supplement) to normalized
% values useable by model fitting scripts. 
% 
% Files needed:
% - proteomics data ("aav3233_Data_file_S1.xlsx", as log2RFU)
% - transcriptomics data ("GSE133529_ProcessedDataFile.csv", as CPM)
% - Netflux model file ("snm_1_1_rev3.xlsx", for extracting snm data)
%   - Reactions sheet: must have designated 'output' reactions in modules
%     column

%% Specify + load initial data files
datadir = "D:\Research\Aim2\ModelFitting\ModelFitting_AnsethData\data_raw";
addpath(datadir)

filepath = "GSE133529_ProcessedDataFile.csv";
opts = detectImportOptions(filepath, "Range", 2);
cpm = readtable(filepath, opts);
cpm = rmmissing(cpm,1);     % remove NaN rows

filepath = "SomaScan_Log2Transform.xlsx";
opts = detectImportOptions(filepath);
rfu = readtable(filepath, opts);

% Get model inputs/outputs from model file
modeldir = "D:\Research\Aim3\ModelExpansion";
addpath(modeldir)
filepath = "snm_1_2_rev5.xlsx";
opts = detectImportOptions(filepath,'Sheet','species');
snm_species = readtable(filepath,opts,'Sheet','species');
opts = detectImportOptions(filepath,'Sheet','reactions');
snm_rxns = readtable(filepath,opts,'Sheet','reactions');

%% Filter + Process output data (gene expression)
% ID outputs from model
output_id = strcmpi(snm_rxns.module,"output");
output_prot = split(snm_rxns.Rule(output_id),"=>");
output_prot = split(unique(output_prot(:,2))," ");
output_prot = output_prot(:,2);

% convert to gene names
% output_genes = cell.empty(length(output_prot),0);
for output = 1:length(output_prot)
    idx = find(strcmpi(snm_species.ID, output_prot(output)));
    hasMultGenes = contains(snm_species.geneName(idx),"; ");
    if hasMultGenes == true
        genes = split(snm_species.geneName(idx),"; ");
    else
        genes = snm_species.geneName(idx);
    end
    if output == 1
        output_genes = genes;
    else
        output_genes = cat(1,output_genes,genes);
    end
end
output_genes_unique = transpose(unique(output_genes));


% filter cpm data for output genes
for output = 1:length(output_genes_unique)
    idx(output) = find(strcmpi(cpm{:,1}, output_genes_unique(output)));
end
output_cpm = cpm(idx,:);
% output_cpm.Properties.RowNames = output_cpm{:,1};
% output_cpm = output_cpm(:,2:end);

% filter genes below threshold
cutoff = 0.5;
belowCuttoff = ~(mean(output_cpm{:,2:end},2)<cutoff);
output_cpm_filt = output_cpm(belowCuttoff,:);

% convert output names back to protein names
exceptions = ["FN1","TGFB1","AGT"];
for output = 1:length(output_cpm_filt{:,1})
    idx = contains(snm_species.geneName, output_cpm_filt{output,1});
    output_prot_ordered{output} = snm_species.ID(idx);
    if length(output_prot_ordered{output}) > 1
        isException = any(strcmpi(output_cpm_filt{output,1},exceptions));
        if isException == true && (output_cpm_filt{output,1} == exceptions(1))
            output_prot_ordered{output} = output_prot_ordered{output}(2);
        elseif (isException == true) && (output_cpm_filt{output,1} == exceptions(2))
            output_prot_ordered{output} = output_prot_ordered{output}(2);
        elseif (isException == true) && (output_cpm_filt{output,1} == exceptions(3))
            output_prot_ordered{output} = output_prot_ordered{output}(3);
%         elseif (isException == true) && (output_cpm_filt{output,1} == exceptions(4))
%             output_prot_ordered{output} = output_prot_ordered{output}(2);
%         elseif (isException == true) && (output_cpm_filt{output,1} == exceptions(5))
%             output_prot_ordered{output} = output_prot_ordered{output}(3);
        else
            idx = strcmpi(snm_species.geneName, output_cpm_filt{output,1});
            output_prot_ordered{output} = snm_species.ID(idx);
        end
    end
end
output_cpm_filt{:,1} = transpose(output_prot_ordered);
% remove extra CI row (COL1A2, keeping COL1A1)
output_cpm_filt = output_cpm_filt([1:3 5:end],:);


%% Filter + Process input data (serum proteomics)
% ID inputs from model
input_id = strcmpi(snm_rxns.module,"input");
input_prot = split(snm_rxns.Rule(input_id),"=>");
input_prot = split(unique(input_prot(:,2))," ");
input_prot = input_prot(:,2);

% get SomaScan protein names (manually added to model file)
input_ss = rmmissing(snm_species.inputName);
% filter SomaScan data for inputs
for input = 1:length(input_ss)
    idx = find(strcmp(rfu.Name, input_ss(input)));
    if input == 1
        idxs = idx;
    elseif isempty(idx) == false
        idxs = cat(1,idxs,idx);
    else
    end
end
input_rfu = rfu(idxs,:);

% Filter columns for (patient) samples contained in cpm dataset
idxs_in = 1;    % protein name
idxs_out = 1;
for idx = 2:length(input_rfu.Properties.VariableNames)
    var = input_rfu.Properties.VariableNames(idx);
    inCpm = any(strcmp(output_cpm_filt.Properties.VariableNames,var));
    if inCpm == true
        idxs_in = cat(1,idxs_in,idx);
    else
        if contains(var,"TAVR")
            idxs_out = cat(1,idxs_out,idx);
        end
    end
end
input_rfu_filt = input_rfu(:,idxs_in);
input_rfu_protonly = input_rfu(:,idxs_out);
input_rfu_all = join(input_rfu_filt, input_rfu_protonly);

% convert input names back to protein names
for output = 1:length(input_rfu_filt{:,1})
    idx = contains(snm_species.inputName, input_rfu_filt{output,1});
    input_prot_ordered{output} = snm_species.ID(idx);
end
input_rfu_filt{:,1} = transpose(input_prot_ordered);
input_rfu_protonly{:,1} = transpose(input_prot_ordered);
input_rfu_all{:,1} = transpose(input_prot_ordered);


%% Normalize Data
input_rfu_t = rows2vars(input_rfu_filt,'VariableNamesSource','Name');
input_rfu_t_protonly = rows2vars(input_rfu_protonly,'VariableNamesSource','Name');
% input_rfu_t_all = rows2vars(input_rfu_all,'VariableNamesSource','Name');

input_rfu_t.Properties.RowNames = input_rfu_t.OriginalVariableNames;
input_rfu_t_protonly.Properties.RowNames = input_rfu_t_protonly.OriginalVariableNames;
% input_rfu_t_all.Properties.RowNames = input_rfu_t_all.OriginalVariableNames;

input_norm_t = normalize(input_rfu_t(:,2:end), 'range', [0.1 0.6]);
% input_norm_t_protonly = normalize(input_rfu_t_protonly(:,2:end), 'range', [0.1 0.6]);
% input_norm_t_all = normalize(input_rfu_t_all(:,2:end), 'range', [0.1 0.6]);

scaleOutData = @(x, min, max) (0.1 + ((x-min)./(max-min)) .* (0.6-0.1));
input_rfu_min = varfun(@min, input_rfu_t(:,2:end));
input_rfu_max = varfun(@max, input_rfu_t(:,2:end));
input_norm_t_protonly = scaleOutData(input_rfu_t_protonly{:,2:end}, input_rfu_min{:,:}, input_rfu_max{:,:});
input_norm_t_protonly(input_norm_t_protonly > 0.95) = 0.95;
input_norm_t_protonly(input_norm_t_protonly < 0.05) = 0.05;
input_norm_t_protonly = array2table(input_norm_t_protonly, ...
                                    'VariableNames', input_rfu_t.Properties.VariableNames(2:end), ...
                                    'RowNames', input_rfu_t_protonly.Properties.RowNames);

% scaleData = @(r, c, x) (r * x{:,:}) + c;
% scaleData = @(r, sd, x) ((1 + r .* sd) .* x{:,:});
% input_sd_t = repmat(std(input_rfu_t{:,2:end}, [], 1),size(input_norm_t,1),1);
% input_scaled_t = scaleData(0.1, input_sd_t, input_norm_t);
% input_scaled_t = scaleData(0.12, 0.6, input_norm_t);
% input_scaled = array2table(input_scaled_t, ...
%     'RowNames', input_norm_t.Properties.RowNames.', ...
%     'VariableNames', input_norm_t.Properties.VariableNames.');

input_scaled = sortrows(input_norm_t, 'RowNames');
input_scaled_protonly = sortrows(input_norm_t_protonly, 'RowNames');

ouput_cpm_t = rows2vars(output_cpm_filt,'VariableNamesSource','Var1');
ouput_cpm_t.Properties.RowNames = ouput_cpm_t.OriginalVariableNames;
output_norm_t = normalize(ouput_cpm_t(:,2:end), 'range', [0.1 0.7]);

% output_scaled_t = scaleData(0.17,0.41,output_norm_t);
% output_scaled = array2table(output_scaled_t, ...
%     'RowNames', output_norm_t.Properties.RowNames.', ...
%     'VariableNames', output_norm_t.Properties.VariableNames.');
output_scaled = sortrows(output_norm_t,'RowNames');


%% Export Data
datadir = "D:\Research\Aim2\ModelFitting\ModelFitting_AnsethData\data_processed_current";
addpath(datadir)
% Specify + save data
filepath = "trainingdata_inputs_Anseth_11232020.mat";
save(strcat(datadir,"\",filepath),'input_scaled');

filepath = "trainingdata_outputs_Anseth_12042020.mat";
save(strcat(datadir,"\",filepath),'output_scaled');

filepath = "testingdata_inputs_Anseth_12312020.mat";
save(strcat(datadir,"\",filepath),'input_scaled_protonly');

