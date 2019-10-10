function [G,dGdtheta] = pp1_modelpred_usageScale(theta,Model)
% Predicts distaces and G-matrix from the 18 parameters of the simple
% scaling model and returns the derivative in respect to the parameters

% single finger params
A = Model.Ac(:,:,1);
% activity scaling feature
scaleParams = exp(theta(1:4));
chords      = pp1_simulations('chords');
M           = chords;
numFingers  = sum(M,2);
for i = 1:4
    M(numFingers==i+1,:) = M(numFingers==i+1,:).*scaleParams(i);
end
OM = A*A';          
G  = M*G*M';  % Second moment matrix

% % Build dertivatives 
for i = 1:4 % scale params
    dM                    = zeros(size(chords));
    dM(numFingers==i+1,:) = chords(numFingers==i+1,:);
    dGdtheta(:,:,i)       = dM*OM*M'+M*OM*dM'; % derivative for chords with numFingers i
    dGdtheta(:,:,i)       = dGdtheta(:,:,i)*scaleParams(i); % scaled derivative 
end

