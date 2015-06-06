clear
clc
close all
%%

nx = 30;
ny = 30;
nz = 10;

x_ = 1:nx;
y_ = 1:ny;
z_ = 1:nz;

% --- Calculate Q by Yunfei's code
%Q = prec_mat_5by5(kappa_true,alpha_true,nx,ny);

% --- Calculate Q directly from base
% --- Create 3-D base
D_base = zeros(ny,nx,nz);
D_base(1,1,1) = -6;
D_base(1,2,1) = 1;
D_base(2,1,1) = 1;
D_base(1,1,2) = 1;

D_base(1,end,1) = 1;
D_base(end,1,1) = 1;
D_base(1,1,end) = 1;

%Q_base = sqrt(nx*ny*nz)*ifft2(fft2(D_base).*fft2(D_base));
% Q_base = ifft2(fft2(D_base).*fft2(D_base));
% Q_base = real(Q_base);
% Q = circulant_block_3D(Q_base);
D = circulant_block_3D(D_base);
Q = D'*D;



% --- Generate MRF ---
L = chol(Q,'lower');
%w = randn(nx*ny*nz,1);
load rand_number.mat
v = L'\w;

w_true = reshape(v,ny,nx,nz);
for i=1:nz
subplot(2,nz/2,i)
%surf(w_true(:,:,i),'EdgeColor','none')
imagesc(w_true(:,:,i))
title(i)
end
%contourslice(w_true,[],[],[1,2,3,4],8);

% imagesc(x_,y_,reshape(w_true,ny,nx));
% colorbar;
% clim = caxis;
% figure(2)
% surf(reshape(w_true,ny,nx));



