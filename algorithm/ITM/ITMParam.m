classdef ITMParam < AlgoParam & handle
    %ITMPARAM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods (Access = protected)
        function configParser(obj)
            check_resolution = @(x) x>0;
            
            default_plot = 'off';
            valid_plot = {'on','off'};
            checkPlot = @(x) any(validatestring(x,valid_plot));
            
            addOptional(obj.parser,'Plot',default_plot,checkPlot);
            addRequired(obj.parser,'resolution_threshold',check_resolution);
        end
    end
    
end

