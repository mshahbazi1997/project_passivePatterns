function [G,dGdtheta] = pp1_modelpred_nonlinearScale_fixedFinger(theta,Model)
% Predicts G-matrix from the 4 parameters of the simple
% scaling model and returns the derivative in respect to the parameters
% Finger params are hardcoded based on perscription (so change accordingly)
% Harvest appropriate params
scaleParams = exp(theta(1:4));
% single finger features
A = Model.Ac(:,:,1); % A = pcm_diagonalize(G_singleFinger)
G = A*A'; 
[V,lam_G] = eig(full(G));
dS    = diag(lam_G);
idx   = dS>eps;
OM    = V(:,idx)*lam_G(idx,idx)*V(:,idx)';
% activity scaling feature
chords     = pp1_imana('chords');
M          = chords;
numFingers = sum(M,2);
for i = 1:4
    M(numFingers==i+1,:) = M(numFingers==i+1,:).*scaleParams(i);
end
%OM = A*A';          
G  = M*OM*M';  % Second moment matrix

for i = 1:4 % scale params
    dM                    = zeros(size(chords));
    dM(numFingers==i+1,:) = chords(numFingers==i+1,:);
    dGdtheta(:,:,i)       = dM*OM*M'+M*OM*dM'; % derivative for chords with numFingers i
    dGdtheta(:,:,i)       = dGdtheta(:,:,i)*scaleParams(i); % scaled derivative 
end;

