classdef SpecClusAlgorithmParam < AlgoParam & handle
    methods 
        function configParser(obj)            
            
            check_num_cluster = @(x) mod(x,1) == 0 & x>=2;
            check_sim_func = @(x) isa(x,'function_handle');
            
            addRequired(obj.parser,'number_of_clusters',check_num_cluster);
            addRequired(obj.parser,'similarity_function',check_sim_func);
        end
    end
end