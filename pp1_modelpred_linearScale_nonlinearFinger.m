function [G,dGdtheta] = pp1_modelpred_linearScale_nonlinearFinger(theta,Model)
% Predicts G-matrix from the 14 parameters of the simple
% scaling model and returns the derivative in respect to the parameters

% Scaling is linear
% Explicitly estimates single finger params.

% Harvest appropriate params
fingerParams    = theta(1:14);
indx1           = double(tril(true(5),0));
indx1(indx1>0)  = [-1 1:14];
indx2           = indx1';
A               = zeros(5);
A(indx1==-1)    = 1; 
% finger feature
for i = 1:14
    A(indx1==i | indx2==i) = fingerParams(i); 
end
% activity scaling feature
chords     = pp1_simulations('chords');
M          = chords;
OM = A*A';          
G  = M*OM*M';  % Second moment matrix

% % Build dertivatives 
for i = 1:14 % finger params
    dA                      = zeros(5);
    dA(indx1==i | indx2==i) = 1;            % contrast for that finger
    dOM                     = dA*A'+A*dA;   % delta of omega for finger param
    dGdtheta(:,:,i)         = M*dOM*M';     % derivative of G w/ respect to omega param
end; 

