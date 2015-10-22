clear, close all
% test DPGMM class
addpath('./../DP-GMM') % add ITM class directory


hyper.pseudo_observations = .1;
hyper.expected_mean = [0; 0];
hyper.degree_of_freedom = 4;
hyper.expected_covariance = [4 .1;.1 4];

alpha = 5;

%%
% sample drown from K component mixture of 2-D Gaussian
% samples N random variable from a mixture of Gaussian with 
% normal-inverse-Wishart prior distribution on mean and covariance


K = 4;
N = 300;

col = lines(K);

for i=1:K
    [comp{i}.W, comp{i}.m] = sample_normal_iwishrnd( ...
        hyper.pseudo_observations,...
        hyper.expected_mean,...
        hyper.degree_of_freedom,...
        hyper.expected_covariance);
end

p = sample_dirichlet(alpha*ones(1,K),N);
%plot(0,0,'or')
hold on
for i=1:N
    k = mnrand_draw(p(i,:),1);
    x{i} = mvnrnd(comp{k}.m,comp{k}.W)';
    plot(x{i}(1),x{i}(2),'.','Color',col(k,:))
end
hold off
clear i K N col

%%

dp = DPGMM(hyper,alpha);
figure(2)
h = gcf;
h.Position = h.Position.*[1 1 2 1];
dp.clusterData(x,2,100)

