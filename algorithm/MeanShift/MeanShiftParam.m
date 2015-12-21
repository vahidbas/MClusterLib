classdef MeanShiftParam < AlgoParam & handle
    %MEANSHIFTPARAM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods (Access = protected)
        function configParser(obj)
            check_bandwidth = @(x) x>0;
            addRequired(obj.parser,'bandwidth',check_bandwidth);
        end
    end
end

