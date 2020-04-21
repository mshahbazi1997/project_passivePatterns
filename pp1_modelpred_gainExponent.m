function [G,dGdtheta] = pp1_modelpred_gainExponent(theta,Model)
% Predicts nonlinear scaling of G = ax^b, where x = linear scaled G
% Finger params are hardcoded based on perscription (so change accordingly)
% Harvest appropriate params
gainParam = exp(theta(1));
expParam = exp(theta(2));
% make linear model out of single finger G
G_lin = Model.Gc; 
G = gainParam * real(G_lin.^expParam);

% calculate derivatives:
dGdtheta = nan(31,31,2);
% partial derivative wrt gainParam: f'(a) = x^b
dGdtheta(:,:,1) = G; % scaled deriv. by param size, so this equals: a*x^b
% partial derivative wrt expParam: f'(b) = ln(x) * ax^b
dGdtheta(:,:,2) = real(log(G_lin)).*G; % deriv will have complex vals
dGdtheta(:,:,2) = dGdtheta(:,:,2) * expParam;


% function [G,dGdtheta] = pp1_modelpred_gainExponent(theta,Model)
% % Predicts nonlinear scaling of G = ax^b, where x = linear scaled G
% % Finger params are hardcoded based on perscription (so change accordingly)
% % Harvest appropriate params
% gainParam = exp(theta(1));
% expParam = exp(theta(2));
% % make linear model out of single finger G
% G_lin = Model.Gc; 
% % x^b (accounting for negative vals in x):
% xb = sign(G_lin).*(abs(G_lin).^expParam);
% G  = gainParam * xb;  % a(x+c)^b
% 
% % calculate derivatives:
% dGdtheta = nan(31,31,2);
% % partial derivative wrt gainParam: f'(a) = x^b
% dGdtheta(:,:,1) = G; % scaled deriv. by param size, so this equals: a*x^b
% % partial derivative wrt expParam: f'(b) = ln(x) * ax^b
% dGdtheta(:,:,2) = real(log(G_lin)).*G; % deriv will have complex vals
% dGdtheta(:,:,2) = dGdtheta(:,:,2) * expParam;