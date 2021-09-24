function [G,dGdtheta] = pp1_modelpred_flexibleLinear(theta,Model)
% Predicts G-matrix from the 4 parameters of the simple
% scaling model and returns the derivative in respect to the parameters

% Harvest appropriate params
scaleParams = exp(theta(1:4));
chords     = pp1_imana('chords');
M          = chords;
numFingers = sum(M,2);
for i = 1:4
    M(numFingers==i+1,:) = M(numFingers==i+1,:).*scaleParams(i);
end

% predicted single-finger model G, given the scaling params:
G_roi = Model.Gc(:,:,1); % use group-avg pos-def G from this region
M_inv = pinv(M);
G = M_inv * G_roi * M_inv'; 

for i = 1:4 % scale param derivatives
    dM                    = zeros(size(M_inv));
    dM(:,numFingers==i+1) = M_inv(:,numFingers==i+1);
    dGdtheta(:,:,i)       = dM*G_roi*M_inv' + M_inv*G_roi*dM'; % derivative for chords with numFingers i
    dGdtheta(:,:,i)       = dGdtheta(:,:,i)*scaleParams(i); % scaled derivative 
end;

