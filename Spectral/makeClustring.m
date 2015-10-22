function obj = makeClustring(alg_name,varargin)

switch alg_name
    case 'unnormalized spectral'
        obj = UnSpecClustring(varargin{1},varargin{2});
    case 'normalized spectral (Shi)'
        obj = NrShiSpecClustring(varargin{1},varargin{2});
    case 'normalized spectral (Ng)'
        obj = NrJorSpecClustring(varargin{1},varargin{2});
    otherwise
        error(['clustering algorithm "' alg_name '" is not recognized!'])
end
