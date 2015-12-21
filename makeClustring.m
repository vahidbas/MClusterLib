function obj = makeClustring(varargin)

%setup parser
parser = inputParser;
algo_class_names = {'UnSpectral',...
                    'NrShiSpectral',...
                    'NrNgSpectral',...
                    'DPGMM',...
                    'ITM',...
                    'MeanShift'};
checkAlgoName = @(x) any(validatestring(x,algo_class_names));
parser.addRequired('algo_name',checkAlgoName);

% parse input
parser.parse(varargin{1});

% make object
algo_class_name = validatestring(parser.Results.algo_name,algo_class_names);
algo_class = str2func(algo_class_name);

param_class_name = [algo_class_name 'Param'];
param_class = str2func(param_class_name);

param_obj = param_class();
params = getParams(param_obj,varargin(2:end));
obj = algo_class(params);
end