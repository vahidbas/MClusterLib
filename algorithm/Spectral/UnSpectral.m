classdef UnSpectral < SpecClusAlgorithm & handle
    methods
        % constructor
        function obj = UnSpectral(params)
            obj = obj@SpecClusAlgorithm(params);
            obj.name = 'unnormalized spectral clustring';            
            
        end
        
        function calcLaplacian(obj)
            obj.laplacian = obj.degree_matrix-obj.adjacency_matrix;
        end
        
        function calcEigen(obj)
            [v,d] = eig(obj.laplacian);
            obj.eigen_vectors = v;
            obj.eigen_values = diag(d);
        end
        
        function T = calcTransform(obj)
            T = obj.getKEigen();
        end

    end
end