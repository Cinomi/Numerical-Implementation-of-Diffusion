close all
clear 
clc

img = imread('camera.png');
img = im2double(img);
[height, width] = size(img);
N = height;
delta_t = 0.1;
n = 100;

% threshold
[gra_x, gra_y] = gradient(img);
lapla = sqrt(gra_x.^2 + gra_y.^2);
thre = max(lapla(:))/2;

% gamma
gamma = 1./(1+(abs(lapla)/thre).^2);
gam = spdiags(reshape(gamma,[],1),0:0,N*N,N*N);

% Dx
dx = zeros(N*N,2);
dx(:,1) = 1;
dx(:,2) = -1;
Dx = spdiags(dx,0:1,N*N,N*N);

% Dy
dy = zeros(N*N,2);
dy(:,1) = 1;
dy(:,2) = -1;
Dy = spdiags(dy,[0,N],N*N,N*N);

% PM
PM = -(Dx'*gam*Dx + Dy'*gam*Dy);

for k=1:n
    img = reshape(img,[],1);
    img = img + delta_t * PM * img;
    img = reshape(img, N, N);
    max_value = max(img(:));
    [max_x, max_y]=find(img==max_value);
    imshow(img),title({['2D Anisotropic diffusion '],['k=',num2str(k),' maximum=(',num2str(max_y),',',num2str(max_x),')']});
    hold on;
    
%     max_value = max(img(:));
%     [max_x, max_y]=find(img==max_value);
    plot(max_y, max_x, 'o','LineWidth', 2, 'MarkerEdgeColor', 'r');
    drawnow;
end