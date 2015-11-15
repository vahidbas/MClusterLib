classdef DPGMMParam < AlgoParam & handle
    %DPGMMPARAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
       
        
    end
    
    methods
        function configParser(obj)
            
            check_hyper_param = @checkHyperParam;
            check_alpha = @(x) x>0;
            
            default_plot = 'off';
            valid_plot = {'on','off'};
            checkPlot = @(x) any(validatestring(x,valid_plot));
            
            default_max_iter = 100;
            check_max_iter = @(x) x>10;
            
            defult_init_K = 2;
            check_init_K = @(x) mod(x,1) == 0 & x>=2;
            
            addRequired(obj.parser,'hyper_parameters',check_hyper_param);
            addRequired(obj.parser,'alpha',check_alpha);
            addOptional(obj.parser,'Plot',default_plot,checkPlot);
            addOptional(obj.parser,'MaxIterations',default_max_iter,check_max_iter);
            addOptional(obj.parser,'InitialK',defult_init_K,check_init_K);
  
        end
    end
    
    methods (Access = private)
        
    end
    
end

function res = checkHyperParam(hyper)
            help_message_hyper = ['The hyper parameter structure should contain ' ...
            'parameters of Normal Inverse Wishart (NIW)' ...
            'distribution:\n' ...
            '   .pseudo_observations\n' ...
            '   .expected_mean\n' ...
            '   .degree_of_freedom\n' ...
            '   .expected_covariance\n'];
        
            field_names = {'pseudo_observations', ...
                'expected_mean', ...
                'degree_of_freedom', ...
                'expected_covariance'};
            if (~isfield(hyper, field_names))
                display(help_message_hyper,class)
                res = 0;
            else
                res = 1;
            end
end
