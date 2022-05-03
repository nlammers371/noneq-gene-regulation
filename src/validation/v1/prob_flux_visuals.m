% Script to generate network diagrams for Ninio-like networks
clear
close all
% initialize parameters
syms kp km ckoff kon b p1 p2 p3 p4;
% define generic transition rate matrix
R = [-kp-kon, kon ,0, kp;
     ckoff, -ckoff - kp, kp,0;
     0, km, -km-ckoff,ckoff;
     b*km,0,kon,-b*km-kon];
% pi = [p1 ; p2; p3; p4];
% Assign param values and solve for state occupancies
% find symbolic expressions for occupancies of each state
cm = jet(128);
% Rval = subs(R,[kp km ckoff kon b],[1 1 50 10 1000]);
Rval = subs(R,[kp km ckoff kon b],[1 1 5 1 50]);
[V,D] = eig(Rval');
ss_vector = V(:,diag(D)==0);
ss_vector = eval(ss_vector / sum(ss_vector));
Rflux = eval(Rval.*ss_vector);
Rflux(eye(4)==1) = 0;
% transfor to directed graph objects
s = repelem(1:4,4);
t = repmat(1:4,1,4);
fluxVec = reshape(Rflux',1,[]);
s = s(fluxVec>0);
t = t(fluxVec>0);
fluxVec = fluxVec(fluxVec>0);
G = digraph(s,t);
minNodeSize = 2;
minLineWidth = 1;
fluxWeights = fluxVec / min(fluxVec) * minLineWidth;%(1+log(fluxWeights(fluxWeights>0) / min(fluxWeights(fluxWeights>0)))) * minLineWidth;
ssWeights = ss_vector / min(ss_vector)*minNodeSize;%(1+log(ss_vector / min(ss_vector))) * minNodeSize;

plot(G,'LineWidth',fluxWeights,'MarkerSize',ssWeights,'NodeColor',cm([20 20 120 20],:))