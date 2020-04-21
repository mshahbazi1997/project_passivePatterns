function [G,dGdtheta] = pp1_modelpred_gainExponent_null(theta,Model)
% Predicts null G according to: G_null = a(x+c)^b, where c is constant
% background pattern
% Finger params are hardcoded based on perscription (so change accordingly)
% Harvest appropriate params
gainParam = exp(theta(1));
expParam = exp(theta(2));
addParam = exp(theta(3));
% make linear model out of single finger G
G_lin = Model.Gc; 
G  = gainParam * ((G_lin + addParam).^expParam);  % a(x+c)^b

% calculate derivatives: (these have been verified against numeric checks)
dGdtheta = nan(31,31,3);
% partial derivative wrt gainParam: f'(a) = (x+c)^b
dGdtheta(:,:,1) = G; % scale the derivative by param size
% partial derivative wrt expParam: f'(b) = ln(x+c) * a(x+c)^b
dGdtheta(:,:,2) = real(log(G_lin + addParam)).*G;
dGdtheta(:,:,2) = dGdtheta(:,:,2) * expParam;
% partial derivative wrt addParam: f'(c) = ab(x+c)^(b-1)
dGdtheta(:,:,3) = gainParam * expParam * ((G_lin + addParam).^(expParam-1));
dGdtheta(:,:,3) = dGdtheta(:,:,3) * addParam;


% function [G,dGdtheta] = pp1_modelpred_gainExponent_null(theta,Model)
% % Predicts null G according to: G_null = a(x+c)^b, where c is constant
% % background pattern
% % Finger params are hardcoded based on perscription (so change accordingly)
% % Harvest appropriate params
% gainParam = exp(theta(1));
% expParam = exp(theta(2));
% addParam = exp(theta(3));
% % make linear model out of single finger G
% G_lin = Model.Gc; 
% % now do nonlinear scaling on ALL values (including single finger params,
% % as the scaling param will adjust for this):
% G_lin = G_lin + ones(31).*addParam;
% signG = sin(G_lin);
% xb = signG.*(abs(G_lin).^expParam);
% G  = gainParam * xb;  % a(x+c)^b
% 
% % calculate derivatives: (these have been verified against numeric checks)
% dGdtheta = nan(31,31,3);
% % partial derivative wrt gainParam: f'(a) = (x+c)^b
% dGdtheta(:,:,1) = G; % scale the derivative by param size
% % partial derivative wrt expParam: f'(b) = ln(x+c) * a(x+c)^b
% dGdtheta(:,:,2) = real(log(G_lin)).*G;
% dGdtheta(:,:,2) = dGdtheta(:,:,2) * expParam;
% % partial derivative wrt addParam: f'(c) = ab(x+c)^(b-1)
% dGdtheta(:,:,3) = gainParam * expParam * (signG.*(abs(G_lin).^(expParam-1)));
% dGdtheta(:,:,3) = dGdtheta(:,:,3) * addParam;