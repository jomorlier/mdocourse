function [c,ceq,dc,dceq]=aggregated_relaxed_stress_constraints(x,P,Epsilon)
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
barg_1rel=(abs(s_1)/sl_1-1).*(A_1/100)-Epsilon;
barg_2rel=(abs(s_2)/sl_2-1).*(A_2/100)-Epsilon;
G_ks_2=1/P*log((exp(P*(barg_1rel))+exp(P*(barg_2rel))));
 c=[G_ks_2];ceq=[];
 dgi_dAj=[-E_1^2*F*L_2^2/(x(1)*L_2*E_1+x(2)*E_2*L_1)^2/sl_1,-E_1*E_2*F*L_2*L_1/(x(1)*L_2*E_1+x(2)*E_2*L_1)^2/sl_1
     -E_1*E_2*F*L_2*L_1/(x(1)*L_2*E_1+x(2)*E_2*L_1)^2/sl_2,-E_2^2*F*L_1^2/(x(1)*L_2*E_1+x(2)*E_2*L_1)^2/sl_2];
 dc=1/(exp(P*(barg_1rel))+exp(P*(barg_2rel)))*([(abs(s_1)/sl_1-1)/100*exp(P*(barg_1rel)) (abs(s_2)/sl_2-1)/100*exp(P*(barg_2rel))]+[x(1)/100*exp(P*(barg_1rel)) x(2)/100*exp(P*(barg_2rel))]*dgi_dAj);
dc=dc';
 dceq=[];