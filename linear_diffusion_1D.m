close all
clear
clc

% read the image and convert to gray scale
img = rgb2gray(imread('castle.jpg'));
img = double(img);

% get the size of image
[row, col] = size(img);

% get f0, N and A
f0 = img(:,200);        % assume f0 is the 200th row of the image
N = size(f0, 1);
A = zeros(N, N);

% laplacian matrix
for i = 1:N
    for j = 1:N
        if i==j
            % first row
            if i==1
                A(i,j)=-1;
                A(i+1,j)=1;
                
            % last row
            elseif i==N
                A(i,j)=-1;
                A(i-1,j)=1;
                
            % other rows
            else
                A(i,j)=-2;
                A(i-1,j)=1;
                A(i+1,j)=1;
            end
        end
    end
end

% parameters
f_min = min(f0);
f_max = max(f0);
delta_t = 0.4;   % 0.4: max delta t to make implicit scheme stable (or 0.5?)

% std_gaussian = sqrt(2*delta_t*k);
% gaussian_filt = gaussmf(-1:1:1, [std_gaussian, 0]);

% initialize f0 
f_explicit = f0;
f_implicit = f0;
f_gaussian = f0;

for k = 1:10

    % explicit scheme
    f_explicit = explicit_scheme(f_explicit, delta_t, A);
    % figure;
    subplot(1,3,1), plot(f_explicit);
    title({['explicit scheme '],['k= ',num2str(k), ' delta t= ',num2str(delta_t)]});
    set(gca, 'XLim', [0, N]);
    set(gca, 'YLim', [f_min, f_max]);
    drawnow;
    
    % implicit scheme
    f_implicit = implicit_scheme(f_implicit, delta_t, A, N);
    subplot(1,3,2), plot(f_implicit);
    title({['implicit scheme '],['k= ',num2str(k), ' delta t= ',num2str(delta_t)]});
    set(gca, 'XLim', [0, N]);
    set(gca, 'YLim', [f_min, f_max]);
    drawnow;


end

    
    
    % convolution with gaussian filter
    sigma = sqrt(2*delta_t*10);
    gaussian_filt = fspecial('gaussian',[100,100],sigma);
    %gaussian_filt = gaussmf(-3.8:0.95:3.8, [sigma, 0]);
    %f_gaussian = conv(f_gaussian, gaussian_filt, 'same');
    f_gaussian=imfilter(f_gaussian,gaussian_filt,'replicate');
    subplot(1,3,3), plot(f_gaussian);
    title({['gaussian filter '],['k= ',num2str(10), ' delta t= ',num2str(delta_t)]});
    set(gca, 'XLim', [0, N]);
    set(gca, 'YLim', [f_min, f_max]);


%% function for explicit scheme
function [output] = explicit_scheme(fk, t, A)
    output = fk + t * A * fk;
end


% function for implicit scheme
function [output] = implicit_scheme(fk, t, A, N)
    output = (eye(N)-t*A) \ fk;
end