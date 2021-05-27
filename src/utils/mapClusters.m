function [w_map, w_full] = mapClusters(filepath, sheet, inputs_idx)

opts = detectImportOptions(filepath, "Sheet", sheet);
w_table = readtable(filepath, opts);

% delineate method of determining group: ("input", "cluster", "indiv")
w_table.type = repmat("indiv", size(w_table, 1), 1);    % default to "indiv"

rulesplit = split(w_table.Rule, "=>");                  % determine inputs from Rule var
w_table.type(rulesplit(:,1) == "") = "input_fit";
w_table.type(inputs_idx) = "input_set";
w_table.type(~isnan(w_table.group)) = "cluster";

% assign shorter index to map ga weights to full weight vector
w_table.index_mapped = zeros(height(w_table),1);

inInputFit = w_table.type == "input_fit";
numInputFit = numel(find(inInputFit));
w_table.index_mapped(inInputFit) = 1:numInputFit;


inCluster = w_table.type == "cluster";
w_table.index_mapped(inCluster) = w_table.group(inCluster) + numInputFit;


inIndiv = w_table.type == "indiv";
if any(inIndiv)
    groups = unique(w_table.group(inCluster));
    start = numInputFit + length(groups) + 1;    % set to after all cluster groups
    stop = numInputFit + length(groups) + length(w_table.type(inIndiv));
    w_table.index_mapped(inIndiv) = start:stop;
end

w_map = w_table.index_mapped(w_table.type ~= "input_set").';
w_full = w_table.index(w_table.type ~= "input_set").';