function [c,ceq]=G_KS_top_plot(x,nelx,nely,H,Hs,penal,KE,ft,Emin,E0,iK,jK,freedofs,F,edofMat,S,Sl,P,nodenrs)
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
mse = reshape(sqrt(sum((U(edofMat)*S).*U(edofMat),2)),nely,nelx); %microscopic Von Mises Stress
rcv=reshape(xPhys,nely,nelx).*(mse/Sl-1);
se=mse.*(Emin+reshape(xPhys,nely,nelx).^penal*(E0-Emin))/E0;
c=max(rcv(:))+1/P*log(mean(exp(P*(rcv(:)-max(rcv(:))))));
ceq=[];
  %% PLOT DENSITIES
 h = figure(1); set(h,'Color',[1 1 1]);
    [Yy,Xx]=find(nodenrs);
    Yy=nely+1-Yy;
    Xx=Xx-1;
    hold on; patchplot2 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-xPhys(:))*[1 1 1],'FaceColor','flat','EdgeColor','none'); axis equal; axis off; drawnow;hold off
    colormap(gray);
  title('density plot')
  figure(2)
  colormap(hot); hold on; patchplot3 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(se(:)),'FaceColor','flat','EdgeColor','none'); axis equal; axis off; drawnow;hold off; 
  caxis([0 max(se(:))]); axis equal; axis off;
  colorbar
  title('macroscopic stress')
  figure(3)
  colormap(hot);patchplot4 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(rcv(:)),'FaceColor','flat','EdgeColor','none'); hold off;  caxis([min(rcv(:)) max(rcv(:))]); axis equal; axis off; drawnow;
  colorbar
  title('relaxed constraint violation')
  figure(4)
  
  m=min(min(mse.*(reshape(xPhys,nely,nelx)>=0.3)));
  M=max(max(mse.*(reshape(xPhys,nely,nelx)>=0.3)));
  colormap(hot);patchplot4 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(mse(:).*(reshape(xPhys,nely*nelx,1)>=0.3)),'FaceColor','flat','EdgeColor','none'); caxis([m M]); axis equal; axis off; drawnow;
  colorbar
  title('microscopic stress in density \geq 0.3')