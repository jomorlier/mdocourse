function [c,ceq]=stress_constraints(x)
E_1=210000; %MPa
E_2=70000; %MPa
L_1=2000; %mm
L_2=1750; %mm
sl_1= 800; %MPa
sl_2= 200; %MPa
rho_1=7.8e-3; %g/mm^3
rho_2=2.7e-3; %g/mm^3
F=10000; %N
A_1=x(1); A_2=x(2);
s_1=F*L_2*E_1./(E_1*A_1*L_2+E_2*A_2*L_1); %MPa
s_2=-F*L_1*E_2./(A_1*L_2*E_1+A_2*L_1*E_2); %MPa
g_1=(abs(s_1)/sl_1-1).*(A_1>0);
g_2=(abs(s_2)/sl_2-1).*(A_2>0);
c=[g_1;g_2];ceq=[];