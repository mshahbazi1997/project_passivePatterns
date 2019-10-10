function [G,dGdtheta] = pp1_modelpred_nonlinearScale_fixedFinger_backgroundPattern(theta,Model)
% Predicts G-matrix from the 4 parameters of the simple
% scaling model and returns the derivative in respect to the parameters
% Finger params are hardcoded based on perscription (so change accordingly)
% Critically, we also include an additive independent background pattern (cpd2
% specific model)

% Harvest appropriate params
scaleParams = exp(theta(1:4));
addParams   = exp(theta(5:8));

G = Model.Gc(:,:,1);        % single finger features
G = blockdiag(G,1);  % additive independent pattern
G = (G+G')/2;               % Symmetrize (make pos. definite)
[V,lam_G] = eig(full(G));
dS    = diag(lam_G);
idx   = dS>eps;
OM    = V(:,idx)*lam_G(idx,idx)*V(:,idx)';

% activity scaling feature
chords     = pp1_simulations('chords');
M          = [chords,[ones(5,1);zeros(26,1)]];
numFingers = sum(chords,2);
for i = 1:4
    M(numFingers==i+1,1:5) = M(numFingers==i+1,1:5).*scaleParams(i); %scale params
    M(numFingers==i+1,6)   = addParams(i); %additive param
end
%OM = A*A';          
G  = M*OM*M';  % Second moment matrix

% calc derivatives of params with respect to G
for i = 1:4 % scale params
    dM                      = zeros(31,6);
    dM(numFingers==i+1,1:5) = chords(numFingers==i+1,:);
    dGdtheta(:,:,i)         = dM*OM*M'+M*OM*dM'; % derivative for chords with numFingers i
    dGdtheta(:,:,i)         = dGdtheta(:,:,i)*scaleParams(i); % scaled derivative 
end;
% additive params
for i = 1:4
    add                = zeros(31,1);
    add(numFingers==i+1) = 1;
    dM                 = [zeros(31,5) add];
    dGdtheta(:,:,i+4) = dM*OM*M'+M*OM*dM';
    dGdtheta(:,:,i+4) = dGdtheta(:,:,i+4)*addParams(i);
end;
