function [img_shift] = shift_usingkspace(img, normalized_shift)

% img = 1/2/3d image
% normalized_shift = actual shift / pixel size ( 1 means 1 pixel shift)
%                   e.g., [0.3 0.5 0] => 0.3 pixel along y(ud) / 0.5 shift along x (rl) / no shift along z (slice) 


[yn, xn, zn] = size(img);

if (mod(yn,2) ==1 )|| (mod(xn,2) ==1 )|| (mod(zn,2) ==1 )
   error('size should be even number') 
end

ys = normalized_shift(1);
xs = normalized_shift(2);
zs = normalized_shift(3);
kk = fftn(img);


yk = zeros(1,yn);
yk(1:end/2) = exp(1i*2*pi/yn*(0:yn/2-1)*ys);
yk(end/2+1) = cos(2*pi*ys/2);
yk(end/2+2:end) =exp(1i*2*pi/yn* ( (yn/2+1:yn-1)-yn)*ys);


xk = zeros(1,xn);
xk(1:end/2) = exp(1i*2*pi/xn*(0:xn/2-1)*xs);
xk(end/2+1) = cos(2*pi*xs/2);
xk(end/2+2:end) =exp(1i*2*pi/xn* ( (xn/2+1:xn-1)-xn)*xs);


zk = zeros(1,zn);
zk(1:end/2) = exp(1i*2*pi/zn*(0:zn/2-1)*zs);
zk(end/2+1) = cos(2*pi*zs/2);
zk(end/2+2:end) =exp(1i*2*pi/zn* ( (zn/2+1:zn-1)-zn)*zs);

lin_phs2 = yk'*(xk);
lin_phsz = repmat(reshape(zk,[1 1 zn]),[yn xn 1]);
lin_phs = repmat(lin_phs2,[1 1 zn]).*(lin_phsz);
kk = kk.*lin_phs;
img_shift = ifftn(kk);

end