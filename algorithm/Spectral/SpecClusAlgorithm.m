classdef SpecClusAlgorithm < ClusAlgorithm & handle
    
    properties (Access = protected)
        adjacency_matrix;
        degree_matrix;
        laplacian;
        eigen_vectors;
        eigen_values;
        similarity_function;
        cluster_number;
    end
    
   %METHODS-------------------------------------------------    
    methods (Abstract)
        calcLaplacian(obj)
        calcEigen(obj)
        T = calcTransform(obj)
    end
    
    methods
        
        function obj = SpecClusAlgorithm(params)
            obj.cluster_number = params.number_of_clusters;
            obj.similarity_function = params.similarity_function;
        end
        
        function mat = getLaplacian(obj)
            mat = obj.laplacian;
        end
        function mat = getAdjMatrix(obj)
            mat = obj.adjacency_matrix;
        end
        function mat = getDegMatrix(obj)
            mat = obj.degree_matrix;
        end
        function mat = getEigVectors(obj)
            mat = obj.eigen_vectors;
        end
        function mat = getEigValues(obj)
            mat = obj.eigen_values;
        end
        
        %clustring algorithm
        function res = clusterImp(obj,varargin)                        
            obj.calcAdjMatrix();
            obj.calcDegMatrix();
            obj.calcLaplacian();
            obj.calcEigen();
            T = calcTransform(obj);
            res.indexes = kmeans(T, obj.cluster_number);            
        end
    end
    
    methods(Access = protected)
        % calculate graph adjacency matrix
        function calcAdjMatrix(obj)            
           
           n_el = length(obj.data); % get size of input data
           obj.adjacency_matrix = zeros(n_el);
           
           for i=1:n_el
               for j=(i+1):n_el
                   obj.adjacency_matrix(i,j) =...
                       obj.similarity_function(...
                       obj.data{i},...
                       obj.data{j}...
                       );
                   
                   % it is a symmetric matrix
                   obj.adjacency_matrix(j,i) = obj.adjacency_matrix(i,j);
               end
           end
        end
        
        % calculate graph adjacency matrix
        function calcDegMatrix(obj)
            obj.degree_matrix = diag(sum(obj.adjacency_matrix,2));
        end
        
        % take first k eigrn vectors (k is number of clusters)
        function mat = getKEigen(obj)
            [~, idx] = sort(obj.eigen_values);
            mat = obj.eigen_vectors(:,idx(1:obj.cluster_number));
        end
    end    
end