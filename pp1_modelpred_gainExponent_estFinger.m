function [G,dGdtheta] = pp1_modelpred_gainExponent_estFinger(theta,Model)
% Predicts nonlinear scaling of G = ax^b, where x = linear scaled G.
% Note that here we also explicitly estimate finger params.
% Finger params are hardcoded based on perscription (so change accordingly)
% Harvest appropriate params
fingerParams = theta(1:14);
gainParam = exp(theta(15));
expParam = exp(theta(16));
% The first 14 paramters determine the stucture of the Finger patterns
% encoded in A. With OM = A*A'. 
indx1=double(tril(true(5),0));
indx1(indx1>0)=[-1 1:14];
indx2         = indx1';
A=zeros(5);
A(indx1==-1)=1; 
for i=1:14
    A(indx1==i | indx2==i) = fingerParams(i); 
end 
OM = A*A';
M  = Model.Ac; % finger inidicator matrix for chords
G_lin  = M*OM*M';  % linearly-scaled second moment matrix
% do nonlinearity: (Gnl = a*Gl^b)
G = gainParam * (G_lin.^expParam);

% calculate derivatives:
dGdtheta = nan(31,31,16);
% partial derivs wrt finger params:
for ii=1:14
    dA = zeros(5);
    dA(indx1==ii | indx2==ii)=1; 
    dOM = dA*A'+A*dA;
    dGdtheta(:,:,ii) = M*dOM*M';   
end
% partial derivative wrt gainParam: f'(a) = x^b
dGdtheta(:,:,ii+1) = G; % scaled deriv. by param size, so this equals: a*x^b
% partial derivative wrt expParam: f'(b) = ln(x) * ax^b
dGdtheta(:,:,ii+2) = real(log(G_lin)).*G; % deriv will have complex vals
dGdtheta(:,:,ii+2) = dGdtheta(:,:,2) * expParam;


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