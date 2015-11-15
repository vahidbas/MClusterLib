function theta = sample_dirichlet(alpha1, N)
% SAMPLE_DIRICHLET Sample N vectors from Dir(alpha(1), ..., alpha(k))
% theta = sample_dirichlet(alpha, N)
% theta(i,j) = i'th sample of theta_j, where theta ~ Dir

% We use the method from p. 482 of "Bayesian Data Analysis", Gelman et al.

k = length(alpha1);
theta = zeros(N, k);
scale = 1; % arbitrary?
for i=1:k
  theta(:,i) = gamrnd(alpha1(i), scale, N, 1);
  %theta(:,i) = sample_gamma(alpha(i), scale, N, 1);
end
%theta = mk_stochastic(theta);
S = sum(theta,2); 
theta = theta ./ repmat(S, 1, k);