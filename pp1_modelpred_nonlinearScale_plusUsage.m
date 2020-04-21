function [G,dGdtheta] = pp1_modelpred_nonlinearScale_plusUsage(theta,Model)
% Predicts G-matrix from the 4 parameters of the simple
% scaling model plus the chord Usage model (weighted component).  
% Returns the derivative in respect to the parameters
% Finger params are hardcoded based on perscription (so change accordingly)

% Model.Gc(:,:,1) : single-finger second moment (1:5,1:5), zeros elsewhere
% Model.Gc(:,:,2) : chord Usage G (from 'ef1_gloveData.m')

% Harvest appropriate params
scaleParams = exp(theta(1:4));
usageParam  = exp(theta(5));

G_sf = Model.Gc(1:5,1:5,1); % single finger G
G_usage = Model.Gc(:,:,2);  % chord usage G

% Build single-finger scaling G
G_sf = (G_sf+G_sf')/2;        % Symmetrize
[V,lam_G] = eig(full(G_sf));
dS    = diag(lam_G);
idx   = dS>eps; % knock out any component with negative eigenvalues to enforce matrix being semi positive definite
OM    = V(:,idx)*lam_G(idx,idx)*V(:,idx)';
% scale single finger activity by # fingers in chord
chords = pp1_imana('chords');
M      = chords;
numFingers = sum(M,2);
for i = 1:4
    M(numFingers==i+1,1:5) = M(numFingers==i+1,1:5).*scaleParams(i); %scale params
end

% Predict second moment from summated single finger G
G = M*OM*M';  % Second moment matrix
% Now add chord usage G (weighted by usageParam)
G = G + G_usage.*usageParam;

% calc derivatives of scaling params with respect to G
for i = 1:4 % scale params
    dM                      = zeros(31,5);
    dM(numFingers==i+1,1:5) = chords(numFingers==i+1,:);
    dGdtheta(:,:,i)         = dM*OM*M'+M*OM*dM'; % derivative for chords with numFingers i
    dGdtheta(:,:,i)         = dGdtheta(:,:,i)*scaleParams(i); % scaled derivative 
end
% derivative of usage G param
dGdtheta(:,:,i+1)=G_usage.*usageParam;
