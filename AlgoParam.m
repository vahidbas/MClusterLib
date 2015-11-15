classdef (Abstract) AlgoParam < handle
    properties (Access = protected)
        parser
        params
    end
    methods
        function p = getParams(obj,inp)
            obj.parse(inp)
            p = obj.params;
        end        
    end
    methods (Access = private)
        function parse(obj,input)
            obj.parser = inputParser;
            obj.configParser(); % make parser object
            obj.parser.parse(input{:});
            obj.params = obj.parser.Results;
        end
    end
    methods(Abstract)
        configParser(obj)
    end
end