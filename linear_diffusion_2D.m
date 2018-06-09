close all
clear 
clc

% read the image and convert to gray scale
img = imread('camera.png');
img = im2double(img);

% get size of image
[row, col] = size(img);
N = row;

% define f
fk = reshape(img, N*N, 1);

% define dd
dd = zeros(N*N, 3);
dd(:,1) = 1;
dd(:,2) = -2;
dd(:,3) = 1;

% Ax and Ay
Ax = spdiags(dd, -1:1, N*N, N*N);
Ay = spdiags(dd, [-N,0,N], N*N, N*N);

% define intentity spare martrix I
I = speye(N*N);
delta_t = 0.1;

for k = 1:100
    % semi-implicit
    f_half = (I-(delta_t/2)*Ax)\((I+(delta_t/2)*Ay)*fk);   % first half
    fk = (I-(delta_t/2)*Ay)\((I+(delta_t/2)*Ax)*f_half);   % second half
    
    % output images
    output_img = reshape(fk,N,N);
    max_value = max(output_img(:));
    [max_x, max_y] = find(output_img == max_value);
    imshow(output_img), title({['2D linear diffusion '],['k=', num2str(k), ' maximum=(',num2str(max_y),',',num2str(max_x),')']});
    hold on;
    

    plot(max_y, max_x, 'o','LineWidth', 2, 'MarkerEdgeColor', 'r');
    drawnow;
end
