classdef NrShiSpecClustring < SpecClusAlgorithm & handle
    methods
        % constructor
        function obj = NrShiSpecClustring(varargin)
            obj = obj@SpecClusAlgorithm(varargin{:});
            obj.name = 'normalized spectral clustring (Shi and Malik 2000)';            
            
        end
        
        function calcLaplacian(obj)
            obj.laplacian = obj.degree_matrix-obj.adjacency_matrix;
        end
        
        function calcEigen(obj)
            [v,d] = eig(obj.laplacian,obj.degree_matrix);
            obj.eigen_vectors = v;
            obj.eigen_values = diag(d);
        end
        
        function T = calcTransform(obj)
            T = obj.getKEigen();
        end

    end
end