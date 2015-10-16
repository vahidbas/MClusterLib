function [W, m] = sample_normal_iwishrnd(k,t,df,Tau)
% FROM:
% Sudderth, Erik B. (Erik Blaine), “Graphical models for visual object
% recognition and tracking.” Massachusetts Institute of Technology, 2006.
% PAGE 45-46
% k: pseudo–observations on the scale of observations
% t: expectation of mean
% df: degree of freedom in covariance matrix
% Tau: expectation of covariance matrix

W = iwishrnd(Tau,df);
m = mvnrnd(t,W/k)';