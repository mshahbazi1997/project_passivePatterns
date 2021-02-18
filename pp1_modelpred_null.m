function [G,dGdtheta] = pp1_modelpred_null(theta,Model)
% nonlinear null model
numCond  = Model.numCond;
G        = ones(numCond)*exp(theta(1))+eye(numCond); 
dGdtheta = ones(numCond)*exp(theta(1));