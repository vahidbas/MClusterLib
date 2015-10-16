classdef DPGMM < handel
    properties (Access = private)
        hyper_params
        conc_param       % DP concenteration parameter
        cluster_params   % cluster parameters
        indicator        % cluster indicator variable
        cluster_prob     % cluster probabilities
        data_likelihood
        data
    end
    
    methods (Access = private)
        % -------------------- Constructor ---------------------------
        function obj = DPGMM(varargin)
        end
        % ------------------------------------------------------------
        % ----- Sample indicator variables for Gibbs ----------------- 
        function sampleIndicatorsGibbs(obj)
        end
        % ------------------------------------------------------------
        % ----- Sample indicator variables for Collapse Gibbs -------- 
        function sampleIndicatorsCollapseGibbs(obj)
        end
        
        function sampleClusterParameters(obj)
        end
        
        function sampleClusterProbabilities(obj)
        end
        
        function isConverged()
        end
    end
    
end