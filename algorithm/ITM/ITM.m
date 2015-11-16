classdef ITM < ClusAlgorithm & handle
    
    properties
        last_sti; % the last recivied stimuli
    end
    properties
        n;
        s;
        num;
        rem_node;
        add_node;
        modified_node;
    end
    properties
        nodes
        adj_mat
        del_mat
        emax
    end% END PROPERTIES
    
    properties (Access = private)
        active_plot
    end
    
    methods
        %==================================================================
        % Constructor
        function obj = ITM(params)
            % modified constructor
            if nargin < 1
                error('inputs are not enough');
            end
            
            obj.emax = params.resolution_threshold; % resolution threshold
            obj.num = 0;
            
            obj.active_plot = 0;
            if isequal(params.Plot, 'on')
                obj.active_plot = 1;
            end
            
        end% Constructor
        %==================================================================
    end
    methods (Access = protected)
        function res = clusterImp(obj)
            % loop through data
            for i = 1:length(obj.data)
                step(obj, obj.data{i});
                if obj.active_plot == 1
                    if obj.num > 2
                        plot(obj, 'edge', 'on');
                        hold on
                        plot(obj.data{i}(1),obj.data{i}(2),'ob');
                        hold off
                    end                
                pause(0.1)
                end
            end
            res = getClusterResults(obj);
        end
        
    end
    methods (Access = private)
        function res = getClusterResults(obj)
            for i = 1:length(obj.data)
                obj.last_sti = obj.data{i};
                matching(obj);
                res.indexes(i) = obj.n;
                
            end
            for i=1:obj.num
                res.weight(:,i)=obj.nodes(i).w';
            end
        end
        
        function matching(obj)
            
            for  i = 1:obj.num
                % norm function is used as distance measure
                d(i) = norm(obj.nodes(i).w-obj.last_sti);
            end
            
            [~, idx] = sort(d);
            
            obj.n = idx(1);
            obj.s = idx(2);
            
        end% function matching
        %==================================================================
        function ref_vec_adapt(obj)
            %XX THIS IS NOT IN ORG ALGO
            obj.modified_node = obj.nodes(obj.n);
            obj.nodes(obj.n).w = (obj.nodes(obj.n).n * obj.nodes(obj.n).w + ...
                obj.last_sti)/(obj.nodes(obj.n).n+1);
            obj.nodes(obj.n).s = (obj.nodes(obj.n).n * obj.nodes(obj.n).s + ...
                obj.last_sti*obj.last_sti')/(obj.nodes(obj.n).n+1);
            obj.nodes(obj.n).n = obj.nodes(obj.n).n+1;
            
            N = obj.nodes(obj.n).n;
            if N>2
                %             obj.nodes(obj.n).c = obj.nodes(obj.n).s  ...
                %                 - obj.nodes(obj.n).w*obj.nodes(obj.n).w' ...
                %                 +.1*eye(size(obj.nodes(obj.n).w,1));
                obj.nodes(obj.n).c = obj.nodes(obj.n).s  ...
                    - obj.nodes(obj.n).w*obj.nodes(obj.n).w';
                
            end
            
            %XX
        end% function ref_vec_adapt
        %==================================================================
        function edg_adapt(obj)
            % (i)
            obj.adj_mat(obj.n,obj.s) = 1;
            obj.adj_mat(obj.s,obj.n) = 1;
            
            % (ii)
            m = 1;
            while m < obj.num
                
                if m == obj.n || m == obj.s || ...
                        ~obj.adj_mat(obj.n,m)
                    
                    m = m+1;
                    continue;
                end % if
                
                % if (n-w).(s-w) <= 0 w is inside of Tales sphere made by s
                % and n
                d = (obj.nodes(obj.n).w - obj.nodes(m).w)'...
                    *(obj.nodes(obj.s).w - obj.nodes(m).w);
                if d <= 0
                    % remove edge
                    obj.adj_mat(obj.n,m) = 0;
                    obj.adj_mat(m,obj.n) = 0;
                    
                    % check if there is emanating edges
                    if sum(obj.adj_mat(m,:)) == 0;
                        obj.nodes(m) = [];
                        obj.adj_mat(m,:) = [];
                        obj.adj_mat(:,m) = [];
                        
                        obj.num = obj.num - 1;
                        
                        %XX
                        obj.rem_node(end+1) = m;
                        %XX
                        
                        % adjust n and s indecies
                        if obj.n > m
                            obj.n = obj.n-1;
                        end % if
                        if obj.s > m
                            obj.s = obj.s-1;
                        end % if
                        
                    end % if
                    
                end % if
                
                m = m+1;
            end % while
            
        end% function edg_adapt
        %==================================================================
        function node_adapt(obj)
            % (i)
            d = (obj.nodes(obj.n).w - obj.last_sti)'...
                *(obj.nodes(obj.s).w - obj.last_sti);
            if d > 0 && ...
                    norm(obj.nodes(obj.n).w-obj.last_sti) > obj.emax
                % add node
                obj.num = obj.num + 1;
                obj.nodes(obj.num).w = obj.last_sti;
                
                obj.adj_mat(obj.num,obj.n) = 1;
                obj.adj_mat(obj.n,obj.num) = 1;
                
                %XX
                obj.nodes(obj.num).s = obj.last_sti*obj.last_sti';
                obj.nodes(obj.num).c = eye(size(obj.nodes(obj.num).w,1));
                obj.nodes(obj.num).n = 1;
                obj.add_node = obj.num;
                
                obj.nodes(obj.n) = obj.modified_node;
                %XX
                
                % (ii)
                if norm(obj.nodes(obj.n).w-obj.nodes(obj.s).w) < obj.emax/2
                    obj.nodes(obj.s) = [];
                    obj.adj_mat(obj.s,:) = [];
                    obj.adj_mat(:,obj.s) = [];
                    
                    obj.num = obj.num - 1;
                    
                    % adjust n and s indecies
                    if obj.n > obj.s
                        obj.n = obj.n-1;
                    end % if
                    
                    %XX
                    obj.rem_node(end+1) = obj.s;
                    obj.add_node = obj.num;
                    %XX
                    
                end
                
            end
        end% function node_adapt
        %==================================================================
        function h = plot(obj,varargin)
            label = 0;
            edge = 0;
            marker = '*';
            
            hax = 0;
            
            % input parsing
            for i =1:2:(nargin-1)
                switch varargin{i}
                    case 'label'
                        if strcmp(varargin{i+1}, 'on')
                            label = 1;
                        end
                    case 'edge'
                        if strcmp(varargin{i+1}, 'on')
                            edge = 1;
                        end
                    case 'axes'
                        hax = varargin{i+1};
                end
            end
            
            if hax == 0
                hax = gca;
            end
            
            for i = 1:obj.num
                x(i) = obj.nodes(i).w(1);
                y(i) = obj.nodes(i).w(2);
            end
            
            if obj.num > 2
                h = voronoi(hax,x,y,marker);
            else
                h = plot(hax,x,y,marker);
            end
            hold(hax, 'on')
            
            if edge == 1
                
                M = del_edges(obj);
                M = triu(M);
                M = obj.adj_mat;
                for i = 1:obj.num
                    for j = 1:obj.num
                        if M(i,j) == 1
                            h(end+1) = plot(hax,[obj.nodes(i).w(1) obj.nodes(j).w(1)],...
                                [obj.nodes(i).w(2) obj.nodes(j).w(2)],'-r');
                        end
                    end
                end
                
            end
            
            if label == 1
                for   i = 1:obj.num
                    h(end+1) = text(x(i)+.4,y(i)+.4, num2str(i),'Parent',hax);
                end
            end
            
            hold(hax, 'off')
        end% function plot
        %==================================================================
        function probs = step(obj,varargin)
            
            if nargin == 2
                sti = varargin{1};
                obj.last_sti = sti;
            elseif nargin == 3
                %XX
                sti = varargin{1};
                obj.last_sti = sti;
                sig = varargin{2};
                %XX
            end
            
            %XX
            obj.rem_node = [];
            obj.add_node = [];
            obj.del_mat = [];
            %XX
            
            if obj.num >=2
                
                matching(obj);
                ref_vec_adapt(obj);
                edg_adapt(obj);
                node_adapt(obj);
                
                %XX
            elseif obj.num == 1
                obj.n = 1;
                
                ref_vec_adapt(obj);
                if norm(obj.nodes(obj.n).w-obj.last_sti) > obj.emax
                    % add node
                    obj.num = obj.num + 1;
                    obj.nodes(obj.num).w = obj.last_sti;
                    
                    obj.adj_mat(obj.num,obj.n) = 1;
                    obj.adj_mat(obj.n,obj.num) = 1;
                    
                    obj.nodes(obj.num).s = obj.last_sti*obj.last_sti';
                    obj.nodes(obj.num).c = eye(size(obj.nodes(obj.num).w,1));
                    obj.nodes(obj.num).n = 1;
                    obj.add_node = obj.num;
                end
                
            elseif obj.num == 0
                obj.nodes(1).w = sti;
                obj.nodes(1).n = 1;
                obj.nodes(1).s = sti*sti';
                obj.nodes(1).c = eye(size(sti,1));
                
                obj.n = 1;
                obj.num = 1;
                obj.add_node = obj.num;
                
            end
            %XX
            if nargin ==3
                probs = prob_z(obj,sig);
            end
            
        end % function step
        %==================================================================
        % calc association probabilities
        function p = prob_z(obj,varargin)
            if nargin == 2
                sig = varargin{1};
                
            elseif nargin == 3
                sti = varargin{1};
                obj.last_sti = sti;
                sig = varargin{2};
            end
            
            for  i = 1:obj.num
                % norm function is used as distance measure
                d(i) = (obj.nodes(i).w-obj.last_sti)'/sig*(obj.nodes(i).w-obj.last_sti);
            end
            p = exp(-d.^2/2/obj.emax)+0.0000001;
            p = p./repmat(sum(p,2),[1 size(p,2)]);
        end
        %==================================================================
        % find full delaunay edges
        function M = del_edges(obj)
            
            if obj.num == 2
                M = [0 1; 1 0];
                return;
            end
            
            for i = 1:obj.num
                x(i) = obj.nodes(i).w(1);
                y(i) = obj.nodes(i).w(2);
            end
            
            tri = delaunay(x,y);
            
            M = zeros(size(obj.adj_mat));
            
            for i = 1:size(tri,1)
                M(tri(i,1),tri(i,2)) = 1;
                M(tri(i,2),tri(i,1)) = 1;
                
                M(tri(i,1),tri(i,3)) = 1;
                M(tri(i,3),tri(i,1)) = 1;
                
                M(tri(i,2),tri(i,3)) = 1;
                M(tri(i,3),tri(i,2)) = 1;
                
            end
            
            obj.del_mat = M;
            
        end % function del_edges
        %==================================================================
        % Make a copy of a handle object.
        function new = copy(this)
            % Instantiate new object of the same class.
            new = feval(class(this),0,0);
            
            % Copy all non-hidden properties.
            p = properties(this);
            for i = 1:length(p)
                new.(p{i}) = this.(p{i});
            end
        end % function copy
        
    end %END METHODS
end%END CLASSDEF