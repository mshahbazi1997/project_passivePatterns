function [G,dGdtheta] = pp1_modelpred_nonlinearScale_nonlinearFinger_rmvCM(theta,Model)
% Predicts distaces and G-matrix from the 18 parameters of the simple
% scaling model and returns the derivative in respect to the parameters

% Harvest appropriate params
fingerParams    = theta(1:14);
scaleParams     = exp(theta(15:17));
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
chords     = pp1_imana('chords');
chords     = chords(1:30,:);
M          = chords;
numFingers = sum(M,2);
for i = 1:3
    M(numFingers==i+1,:) = M(numFingers==i+1,:).*scaleParams(i);
end
OM = A*A';          
G  = M*OM*M';  % Second moment matrix

% % Build dertivatives 
for i = 1:14 % finger params
    dA                      = zeros(5);
    dA(indx1==i | indx2==i) = 1;            % contrast for that finger
    dOM                     = dA*A'+A*dA;   % delta of omega for finger param
    dGdtheta(:,:,i)         = M*dOM*M';     % derivative of G w/ respect to omega param
end; 

for i = 1:3 % scale params
    dM                      = zeros(size(chords));
    dM(numFingers==i+1,:)     = chords(numFingers==i+1,:);
    dGdtheta(:,:,14+i)      = dM*OM*M'+M*OM*dM'; % derivative for chords with numFingers i
    dGdtheta(:,:,14+i)      = dGdtheta(:,:,14+i)*scaleParams(i); % scaled derivative 
end;