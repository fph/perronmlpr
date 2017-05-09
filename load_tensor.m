function R = load_tensor(name)
% loads a tensor from Gleich's set
%
% if called without arguments, returns a list of available tensor names

folder = 'mlpagerank-master/tensors/';
files = {'mtm3.mat', 'mtm4.mat', 'mtm6.mat'};
variable_names = {'R3_mats', 'R4_mats', 'R6_mats'};

if not(exist('name', 'var'))
    all_names = {};
    for i = 1:length(files)
        Rs = load(strcat(folder,files{i}));
        names = getfield(Rs, variable_names{i});
        all_names = [all_names names];
    end
    R = all_names;
    return
else
    for i = 1:length(files)
        Rs = load(strcat(folder,files{i}));
        names = getfield(Rs, variable_names{i});
        idx = find(strcmp(names, name));
        if not(isempty(idx))
            R = getfield(Rs, name);
            return
        end
    end
end
error('Tensor %s not found', name);