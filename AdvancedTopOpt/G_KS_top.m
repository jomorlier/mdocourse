function [c,ceq,dc,dceq]=G_KS_top(x,nelx,nely,H,Hs,penal,KE,ft,Emin,E0,iK,jK,freedofs,F,edofMat,S,Sl,P)
global loop
loop=loop+1;
%% FE-ANALYSIS
if ft == 1
    xPhys = x;
elseif ft == 2
    xPhys(:) = (H*x(:))./Hs;
end
U = zeros(2*(nely+1)*(nelx+1),1);
sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
mse = reshape(sqrt(abs(sum((U(edofMat)*S).*U(edofMat),2))),nely,nelx); %microscopic Von Mises Stress
rcv=xPhys(:).*(mse(:)/Sl-1);
sS=reshape(S(:)*(xPhys(:)'./mse(:)'/Sl.*exp(P*(rcv(:)'-max(rcv(:)))))/sum(exp(P*(rcv(:)-max(rcv(:))))),64*nelx*nely,1);
S0=sparse(iK,jK,sS); S0=(S0+S0')/2;
dG_ksdu=S0*U;
Lambda=U;
Lambda(freedofs) = K(freedofs,freedofs)\dG_ksdu(freedofs);
dG_ksdx=(mse(:)/Sl-1).*exp(P*(rcv(:)-max(rcv(:))))/sum(exp(P*(rcv(:)-max(rcv(:)))));
% se=mse.*(Emin+reshape(xPhys,nely,nelx).^penal*(E0-Emin))/E0;
c=max(rcv(:))+1/P*log(sum(exp(P*(rcv(:)-max(rcv(:))))))-log(length(rcv(:)))/P;
ceq=[];dceq=[];
dc_du=reshape((sum((U(edofMat)*KE).*Lambda(edofMat),2)),nely,nelx);
dc_du =- penal*(E0-Emin)*xPhys(:).^(penal-1).*dc_du(:);
dc=dc_du(:)+dG_ksdx;
if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
end
figure(5)
hold on
plot(loop,(mean(xPhys(:))*100),'bo','MarkerFaceColor','b')
plot(loop,c*100,'ro','MarkerFaceColor','r')
% plot(outeriter,(1+GKSl)*VMl,'ko','MarkerFaceColor','k')
title(['Convergence volfrac = ',num2str(mean(xPhys(:))*100),', G_{KS}^l =',num2str(c*100),'%, iter = ', num2str(loop)])
grid on
legend('Volume Fraction %','G_{KSl} %')
xlabel('iter')