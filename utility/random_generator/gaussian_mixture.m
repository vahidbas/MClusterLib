function x = gaussian_mixture(hyper,alpha,K,N)
% sample drown from K component mixture Gaussian
% samples N random variable from a mixture of Gaussian with 
% normal-inverse-Wishart prior distribution on mean and covariance



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