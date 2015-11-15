classdef DPGMM < handle
    properties (Access = private)
        hyper_params
        conc_param       % DP concenteration parameter
        cluster_params   % cluster parameters
        indicator        % cluster indicator variable
        cluster_prob     % cluster probabilities
        data_likelihood
        data
        N                % number of elements in dataset
        K           % umber of initial clusters
        p_id             % random peemutation index
    end
    
    properties (Access = private)
        help_message_hyper = ['The hyper parameter structure should contain ' ...
            'parameters of Normal Inverse Wishart (NIW)' ...
            'distribution:\n' ...
            '   .pseudo_observations\n' ...
            '   .expected_mean\n' ...
            '   .degree_of_freedom\n' ...
            '   .expected_covariance\n'];
        
        in_parser;
    end
    
    methods
        % -------------------- Constructor ---------------------------
        function obj = DPGMM(hyper,alpha)
            
            field_names = {'pseudo_observations', ...
                'expected_mean', ...
                'degree_of_freedom', ...
                'expected_covariance'};
            
            if (~isfield(hyper, field_names))
                error(obj.help_message_hyper,class(hyper))
            end
            
            obj.hyper_params = hyper;
            obj.conc_param = alpha;
            
            % configure input parser
            obj.in_parser = inputParser;
            default_plot = 'off';
            valid_plot = {'on','off'};
            checkPlot = @(x) any(validatestring(x,valid_plot));
            
            default_max_iter = 100;
            check_max_iter = @(x) x>10;
            
            defult_init_K = 2;
            check_init_K = @(x) x>1;
            
            addOptional(obj.in_parser,'Plot',default_plot,checkPlot);
            addOptional(obj.in_parser,'MaxIterations',default_max_iter,check_max_iter);
            addOptional(obj.in_parser,'InitialK',defult_init_K,check_init_K);
        end
        % ------------------------------------------------------------
        
        function clusterData(obj,data, varargin)
            obj.data = data;
            
            % parse input
            parse(obj.in_parser,varargin{:})
            % initial number of clusters
            obj.K = obj.in_parser.Results.InitialK; 
            obj.N = length(data);
            
            max_iterations = obj.in_parser.Results.MaxIterations;
            
            active_plot = 0;
            if isequal(obj.in_parser.Results.Plot, 'on')
                active_plot = 1;
            end
            
            
            % algorithm
            initialize(obj);
            iter = 0;
            while (~obj.isConverged() && iter<max_iterations)
                iterate(obj)
                
                if (active_plot)
                    plot(obj)
                end
                
                iter = iter+1;
            end
            
        end
        
        function result = getClusteringResults(obj)
            result.indexes = obj.indicator;
            result.clusters_number = obj.K;
            result.probabilities = obj.cluster_prob;
            result.parameters = obj.cluster_params;
            result.data_likelihood = obj.data_likelihood(end);
        end
        
    end
    
    methods (Access = private)
        function initialize(obj)
            K_ = obj.K;
            N_ = obj.N;
            obj.indicator = mnrand_draw(ones(1,K_)/K_,N_);
        end
        
        function iterate(obj)
            obj.p_id = randperm(obj.N); % step 1
            f_ki = calcPredictiveLik(obj); % step 2-a
            Nk = calcNumElement(obj);
            z = sampleIndicators(obj,Nk,f_ki); % step 2-b
            obj.indicator(obj.p_id) = z; % step 3, correct random permutation
            
            % increment number of clusters, if new cluster is empty it will
            % be removed in the next step
            obj.K = obj.K+1; 
            
            removeEmptyCluster(obj); % step 4
        end
        
        function f_ki = calcPredictiveLik(obj)
            K_ = obj.K;
            N_ = obj.N;
            p_i = obj.p_id;
            
            for i=1:N_
                ti = p_i(i);  % true index of element
                for k=1:K_+1
                    param = calcPredictiveParam(obj,k,i);
                    f_ki(k,i) = mvnpdf(obj.data{ti},param.m,param.C)+1e-20;
                end
            end
            
        end
        
        function param = calcPredictiveParam(obj,k,i)
            N_ = obj.N;
            p_i = obj.p_id;
            
            
            idx = find(obj.indicator == k);
            if i ~=0
                ti = p_i(i);  % true index of element
                idx(idx == ti) = []; % remove the i_th sample
            end
            L = numel(idx);
            
            xk = [];
            for j = 1:L
                xk(:,j) = obj.data{idx(j)};
            end
            
            hype = obj.hyper_params;
            kp = hype.pseudo_observations;
            tet = hype.expected_mean;
            nu = hype.degree_of_freedom;
            del = hype.expected_covariance;
            dim = length(hype.expected_mean); % dimension of data
            
            if nu <= dim+1
                error('gaussian approximation is not possible, Student-t distribution should be used.')
            end
            
            if (~isempty(xk))
                n_kp_tet_prod = kp * tet + sum(xk,2);  % (2.62)
                n_kp = kp + L; % (2.62)
                n_tet = n_kp_tet_prod / n_kp;
                
                n_nu_del_prod = nu * del + ...
                    kp * (tet * tet') - n_kp * (n_tet * n_tet') + ...
                    xk * xk';  % (2.63)
                n_nu = nu + L; % (2.63)
                n_del = n_nu_del_prod/n_nu;
                
                kp = n_kp;
                tet = n_tet;
                nu = n_nu;
                del = n_del;
            end
            
            param.m = tet;
            param.C = (kp+1)*nu/kp/(nu-dim-1)*del;
            
        end % calcPredictiveParam
        
        function Nk = calcNumElement(obj)
            alpha = obj.conc_param;
            K_ = obj.K;
            N_ = obj.N;
            p_i = obj.p_id;
            
            for i=1:N_
                ti = p_i(i);  % true index of element
                for k=1:K_
                    idx = find(obj.indicator == k);
                    idx(idx == ti) = []; % remove the i_th sample
                    Nk(k,i) = numel(idx);
                end
                Nk(obj.K+1,i) = alpha;
            end
        end % calcNumElement
        % ----- Sample indicator variables for Gibbs -----------------
        function z = sampleIndicators(obj,Nk,f_ki)
            N_ = obj.N;
            
            for i=1:N_
                p_u = Nk(:,i).*f_ki(:,i);
                p = p_u/sum(p_u);
                z(i,1) = mnrand_draw(p,1);
            end
            
        end
        % ------------------------------------------------------------
        
        function removeEmptyCluster(obj)
            K_ = obj.K;
            N_ = obj.N;
            
            new_maping = 1:K_;
            for k = 1:K_
                if sum(obj.indicator==k) == 0
                    obj.K = obj.K-1;
                    new_maping(k) = 0; 
                    for k0 = (k+1):K_
                        new_maping(k0) = new_maping(k0)-1;
                    end
                end
            end
            
            % correct order of indicators
            for k = 1:K_
               obj.indicator(obj.indicator==k) = new_maping(k);
            end
        end
        % ----- Sample indicator variables for Collapse Gibbs --------
        
        function ldl = calcDataLikelihood(obj)
            K_ = obj.K;
            N_ = obj.N;
            
            % put this latr in another function ------x
            clus_prob = zeros(1,K_);
            param = struct('m', cell(1,K_), 'C', cell(1,K_));
            for k=1:K_
                clus_prob(k) = sum(obj.indicator==k)/N_;
                param(k) = calcPredictiveParam(obj,k,0);
            end
            obj.cluster_prob = clus_prob;
            obj.cluster_params = param;
            % ----------------------------------------x
            
            for i=1:N_
                for k=1:K_
                    val(k,i) = clus_prob(k)* ...
                        mvnpdf(obj.data{i},param(k).m,param(k).C);
                end
            end
            ldl = sum(log(sum(val)));
        end
        
        function b = isConverged(obj)
            ldl = calcDataLikelihood(obj);
            obj.data_likelihood(end+1) = ldl;
            b = 0;
        end
        
        function plot(obj)
            K_ = obj.K;
            N_ = obj.N;
            

            dim = length(obj.hyper_params.expected_mean); % dimension of data
            if dim ~=2
                error('plot mode ony works for 2D data')
            end
            
            col = lines(K_);
            
            subplot(1,2,1)
            for k=1:K_
                idx = find(obj.indicator == k);
                L = numel(idx);                
                xk = [];
                for j = 1:L
                    xk(:,j) = obj.data{idx(j)};
                end
                plot(xk(1,:),xk(2,:),'.','Color',col(k,:))
                hold on
            end
            hold off
            subplot(1,2,2)
            plot(obj.data_likelihood)
            
            pause(0.1)
        end
    end
    
end