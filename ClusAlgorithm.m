classdef (Abstract) ClusAlgorithm < handle
    properties (Access = protected)
        name;
        result;
    end
    
    properties
        data;
    end
    
    methods
        function ret = getName(obj)
            ret = obj.name;
        end
        

        function res = cluster(obj,varargin)
            if numel(varargin) > 2
                error('too many input arguments')
            elseif numel(varargin) == 1
                obj.data = varargin{1};
            end
            
            obj.checkDatat(); % check the input data
                       
            obj.result = obj.clusterImp();
            
            res = obj.result;
        end
        
        function idx = getResultIdx(obj)
            numclus = max(obj.result);
            if isempty(numclus)
                error('no results is available')
            end
            
            for i=1:numclus
                idx{i} = find(obj.result == i);
            end
        end
    end
    
    methods(Access = protected)
        % check consistency of user data
        function checkDatat(obj)
            
            if ~iscell(obj.data)
                error('data must be a cell array of objects');
            end
                
            [nr, nc] = size(obj.data); % get size of input data
            
            if nr > 1 
                error('data must be an one row cell array');                
            end
            
            if nc < 2
                error('data should have at least two cells');
            end
            
        end
        
    end
    
    methods (Abstract)
        res = clusterImp(obj)
    end
end