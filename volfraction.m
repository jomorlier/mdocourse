function [o,do]=volfraction(x,ft,H,Hs)
if ft == 1
    xPhys = x;
elseif ft == 2
    xPhys(:) = (H*x(:))./Hs;
end
o=mean(xPhys(:));
do=ones(size(xPhys))/length(xPhys(:));
if ft == 1
elseif ft == 2
    
    do(:) = H*(do(:)./Hs);
end