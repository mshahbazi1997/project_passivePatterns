function [G,dGdtheta] = pp1_modelpred_null(theta,Model)
% nonlinear null model
G        = ones(31)*exp(theta(1))+eye(31); 
dGdtheta = ones(31)*exp(theta(1));