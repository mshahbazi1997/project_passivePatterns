function [G,dGdtheta] = pp1_modelpred_usageScale(theta)
% Predicts distaces and G-matrix from the 18 parameters of the simple
% scaling model and returns the derivative in respect to the parameters

% finger params from usage model G
A = [0.227769349059601,-0.0347230105282332,-0.0622760019593427,-0.0692271083724356,-0.0615432281995893;
    -0.0347230105282332,0.108242506313172,0.0181421441247004,-0.0371376376887697,-0.0545240022208694;
    -0.0622760019593427,0.0181421441247004,0.0682445579701298,0.00516921636362582,-0.0292799164991133;
    -0.0692271083724356,-0.0371376376887697,0.00516921636362581,0.0633152876712225,0.0378802420263569;
    -0.0615432281995894,-0.0545240022208694,-0.0292799164991133,0.0378802420263569,0.107466904893215];
% A = [1,0.797271061785380,0.789717044525792,0.785230573262153,0.770838644824003;
%     0.797271061785380,1,0.929898611983050,0.877083465319033,0.837621292175937;
%     0.789717044525792,0.929898611983050,1,0.939389293542950,0.882864352069214;
%     0.785230573262153,0.877083465319033,0.939389293542950,1,0.952489145744138;
%     0.770838644824003,0.837621292175937,0.882864352069214,0.952489145744138,1];
% activity scaling feature
scaleParams = exp(theta(1:4));
chords      = pp1_simulations('chords');
M           = chords;
numFingers  = sum(M,2);
for i = 1:4
    M(numFingers==i+1,:) = M(numFingers==i+1,:).*scaleParams(i);
end
OM = A*A';          
G  = M*OM*M';  % Second moment matrix

% % Build dertivatives 
for i = 1:4 % scale params
    dM                    = zeros(size(chords));
    dM(numFingers==i+1,:) = chords(numFingers==i+1,:);
    dGdtheta(:,:,i)       = dM*OM*M'+M*OM*dM'; % derivative for chords with numFingers i
    dGdtheta(:,:,i)       = dGdtheta(:,:,i)*scaleParams(i); % scaled derivative 
end

