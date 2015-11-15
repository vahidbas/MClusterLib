classdef ITMParam < AlgoParam & handle
    %ITMPARAM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function configParser(obj)
            check_resolution = @(x) x>0;
            addRequired(obj.parser,'resolution_threshold',check_resolution);
        end
    end
    
end

