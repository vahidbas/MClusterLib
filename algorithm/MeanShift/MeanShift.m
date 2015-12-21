classdef MeanShift < ClusAlgorithm & handle
    %MEANSHIFT Summary of this class goes here
    %   Detailed explanation goes here
   properties (Access = private)
       bandwidth;
       datavec;
       dim;
   end
   methods
        % -------------------- Constructor ---------------------------
        function obj = MeanShift(params)
            obj.bandwidth = params.bandwidth; % resolution threshold
        end
        % ------------------------------------------------------------
    end
    methods (Access = protected)
        function res = clusterImp(obj)
            obj.dim = size(obj.data{1},1);
            for i = 1:length(obj.data)
                obj.datavec(:,i) = obj.data{i}; 
            end
        end
        
%         function result = getClusteringResults(obj)
%         end
        
    end
    
    methods (Access = private)
        function m = calc_meanshift_at(obj,x)
            % implementation of equation (17)
            temp = (obj.datavec - repmat(x,[1, size(obj.datavec,2)])).^2;
            
            d = sum(temp)/obj.bandwidth^2;
            
            num = sum(obj.datavec.*repmat(d,[obj.dim, 1]),2);
            denum = sum(d);
            m = num/denum-x;
        end
    end
end

