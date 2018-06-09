close all
clear 
clc

% read image
img = rgb2gray(imread('soup.jpg'));
img = im2double(img);

% get size of image
[row, col] = size(img);
N = row;

delta_t = 0.1;
[fx, fy] = gradient(img);
gradient_f = sqrt(fx^2 + fy^2);
T  = max(gradient_f(:))/2; %%%%%% can it be T = max(gradient_f(:))?

% kappa
kappa = 1./(1+(gradient_f./T).^2);
kap = spdiags(reshape(kappa,[],1),0:0,N*N,N*N);

% dx
ddx = zeros(N*N,2);
ddx(:,1) = 1;
ddx(:,2) = -1;
Dx = spdiags(ddx,[0,1],N*N,N*N);

% dy
ddy = zeros(N*N,2);
ddy(:,1)=1;
ddy(:,2)=-1;
Dy = spdiags(ddy,[0,N],N*N,N*N);

PM = - (Dx'*kap *Dx +  Dy'*kap*Dy);

output_img = img;

for k = 1:100
    output_img = reshape(output_img,[],1);
    output_img = output_img + delta_t * PM * output_img;
    
    output_img = reshape(output_img, N, N);
    imshow(output_img), title({['anisotropic diffusion 2D '],['itr=',num2str(k),' delta t=',num2str(delta_t)]});
    hold on;
    
    max_value = max(output_img(:));
    [max_row, max_col] = find(output_img == max_value);
    plot(max_col, max_row, 'o','LineWidth', 3, 'MarkerEdgeColor', 'r');
    drawnow;
end