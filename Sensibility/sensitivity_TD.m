%% Sensitivity example computation
%% Clean up
%%

%clear workspace
clear all;
%clear command window
clc;
%clear figures
clf;
%% Define symbolic quantities
%%
syms k_1 k_2 real
syms u_1 u_2 real
syms p_1 p_2 real
%the derived quantities will automatically be initialized as symbolic variables as well

%% Static equilibrium
%%
K=[k_1+k_2 -k_2; -k_2 k_2]
U=[u_1; u_2]
P=[p_1;p_2]
%% 
% Let's start by computing theses sensitivities using 
%% SYMBOLIC (MATLAB's syms)

U0=K\P
U01=simplify(U0(1))

U02=simplify(U0(2))
F=(P(1)*U0(1)+P(2)*U0(2))
gradF_k1=simplify(diff(F,k_1))
gradF_k2=simplify(diff(F,k_2))
gradU_k1=diff(U0,k_1)
gradU_k2=diff(U0,k_2)
%% 
% Numerical applications: $k_1=1$; $k_2=2$, $p_1=1$, $p_2=2$

gFk1=subs(gradF_k1,{k_1,k_2,p_1,p_2},{1,2,1,2})
gFk2=subs(gradF_k2,{k_1,k_2,p_1,p_2},{1,2,1,2})
gUk1=subs(gradU_k1,{k_1,k_2,p_1,p_2},{1,2,1,2})
gUk2=subs(gradU_k2,{k_1,k_2,p_1,p_2},{1,2,1,2})

%% 
% Now we need to compute theses sensitivities by 
%% FD (Finite Differences)
%% 

dx1=1e-6; % you can change this for evaluating the table unknwonws (see exercice)
dx2=dx1;
Kdx1=[k_1+dx1+k_2 -k_2; -k_2 k_2]
Kdx2=[k_1+k_2+dx2 -k_2-dx2; -k_2-dx2 k_2+dx2]
U1=Kdx1\P
U2=Kdx2\P
F=(P(1)*U0(1)+P(2)*U0(2))
Fdx1=(P(1)*U1(1)+P(2)*U1(2))
Fdx2=(P(1)*U2(1)+P(2)*U2(2))
gFk1_FD= (Fdx1-F)/dx1
gFk2_FD= (Fdx2-F)/dx2
gFk1_FD=eval(subs(gFk1_FD,{k_1,k_2,p_1,p_2},{1,2,1,2}))
gFk2_FD=eval(subs(gFk2_FD,{k_1,k_2,p_1,p_2},{1,2,1,2}))
%% 
% Please double click on the workspace value of  -8.999991000009000 for 
% gFk1_FD and -0.999999500000250 for gFk2_FD (WARNING DOUBLE /SINGLE PRECISION)
% 
% Now we need to compute theses sensitivities by 
%% DIRECT METHOD

invK=inv(K)
dK_k1=diff(K,k_1)
dK_k2=diff(K,k_2)
dP_k1=diff(P(1),k_1)
dP_k2=diff(P(2),k_2)
gUk1_DM=invK*(dP_k1-diff(K,k_1)*U0)
gUk2_DM=invK*(dP_k2-diff(K,k_2)*U0)
gFk1_DM=P'*gUk1_DM
gFk2_DM=P'*gUk2_DM
gFk1_DM=eval(subs(gFk1_DM,{k_1,k_2,p_1,p_2},{1,2,1,2}))
gFk2_DM=eval(subs(gFk2_DM,{k_1,k_2,p_1,p_2},{1,2,1,2}))
%% 
% Un petit test pour comparer l'exactitude des 2 méthodes

gFk1_DM==-9 
gFk1_FD==-9
gFk2_DM==-1 
gFk2_FD==-1
%% ADJOINT METHOD
% 
% 
% 
% 
% G= U'KU=U'P=U1*P1+U2*U2
% 
% 
% 
% 

G=U'*P
dG_du1=diff(G, u_1);
dG_du2=diff(G, u_2);
dG_dU=[dG_du1; dG_du2]

Lambda=inv(K')*dG_dU
Lambda=eval(subs(Lambda,{k_1,k_2,p_1,p_2},{1,2,1,2}))
U=eval(subs(U0,{k_1,k_2,p_1,p_2},{1,2,1,2}))
dG_k1=diff(G,k_1) 
dG_k2=diff(G,k_2) 
dG_k2=eval(subs(dG_k2,{k_1,k_2,p_1,p_2},{1,2,1,2}))

dG_k1=eval(subs(dG_k1,{k_1,k_2,p_1,p_2},{1,2,1,2}))
%% 
% 

dFtilde_k1=dG_k1 - Lambda'*dK_k1*U-dP_k1
dFtilde_k2=dG_k2 - Lambda'*dK_k2*U-dP_k2
gradFtilde_k1=eval(subs(dFtilde_k1,{k_1,k_2,p_1,p_2},{1,2,1,2}))
gradFtilde_k2=eval(subs(dFtilde_k2,{k_1,k_2,p_1,p_2},{1,2,1,2}))
%% 
% Un petit test pour comparer l'exactitude des 2 méthodes

gradFtilde_k1==-9 
gradFtilde_k2==-1