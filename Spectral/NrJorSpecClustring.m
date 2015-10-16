classdef NrJorSpecClustring < SpecClusAlgorithm & handle
    methods
        % constructor
        function obj = NrJorSpecClustring(varargin)
            obj = obj@SpecClusAlgorithm(varargin{:});
            obj.name = 'normalized spectral clustring (Ng, Jordan and Weiss 2002)';            
            
        end
        
        function calcLaplacian(obj)
            L = obj.degree_matrix-obj.adjacency_matrix;
            DS = sqrt(obj.degree_matrix);
            obj.laplacian = DS\L/DS;
            
        end
        
        function calcEigen(obj)
            [v,d] = eig(obj.laplacian);
            obj.eigen_vectors = v;
            obj.eigen_values = diag(d);
        end
        
        function T = calcTransform(obj)
            U = obj.getKEigen();
            s = sqrt(sum(U.^2,2));
            T = U./repmat(s,1,size(U,2));
        end
        

    end
end