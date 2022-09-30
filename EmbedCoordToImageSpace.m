function [I,BW] = EmbedCoordToImageSpace(x,y,Param)

% This function takes shape coordinates and embeds them into image space 
% using a selected noise generation function.

% Obtain parameters
ImSize = Param.ImSize;
RotationAngle = Param.RotationAngle;
ScalingRatio = Param.ScalingRatio;
Persistence = Param.Persistence;
Method = Param.Method;
Motion = Param.Motion;
NbrFrames = Param.NbrFrames;

% Obtain center coordinates
[xg,yg] = meshgrid(1:ImSize,1:ImSize);
ImCenter = mean(1:ImSize);

% Perform initial rotation
rot_matrix = [cosd(RotationAngle), -sind(RotationAngle); sind(RotationAngle), cosd(RotationAngle)];
rot = rot_matrix * [x,y]';
x = reshape( rot(1,:) , [] , 1); y = reshape( rot(2,:) , [] , 1);

% Perform spline interpolation with scaling
interp_pts = fnplt(cscvn([x,y]'));
x_interp = interp_pts(1,:)';
y_interp = interp_pts(2,:)';
Max_Width = max(abs( [ max(x_interp) , min(x_interp) , ...
                       max(y_interp) , min(y_interp) ] ));
x_interp = x_interp / Max_Width * ScalingRatio * ImSize / 2;
y_interp = y_interp / Max_Width * ScalingRatio * ImSize / 2;

% Cast into pixel space (binary image == GROUND TRUTH)
BW = roipoly( zeros(ImSize) , x_interp + ImCenter , y_interp + ImCenter );

% Methods of generating cell texture
switch lower(Method)
case 'perlin'
    % Method #1: Use Perlin noise to generate cell texture
    H = fspecial('disk',3);
    P_N = perlin_noise(BW,Persistence); 
case 'fractal'
    % Method #2: Use Fractal noise to generate cell texture
    H = fspecial('disk',3);
    P_N = fractal_noise(BW,Persistence);
end
P_N = P_N - min(P_N(:)); P_N = P_N / max(P_N(:));
I = conv2( im2double(BW) , H , 'same' ) + (BW .* P_N);

% Stop here if only 1 frame is to be generated
if NbrFrames == 1, return; end;

% Create motion variables
K = double(BW); % analog for mask
L = P_N; % analog for noise
I_temp = zeros(ImSize,ImSize,NbrFrames);
BW_temp = zeros(ImSize,ImSize,NbrFrames);

% Generate motion by warping
switch lower(Motion)
case 'asympotic'
    a = @(f) 1e-2 * exp(-f^1.3);
case 'shrink-expand'
    turning_point = 0.25 + rand*0.1 + randn*0.1;
    a = @(f) ((NbrFrames-f)/NbrFrames - turning_point) * 9e-6 ;
case 'expand-shrink'
    turning_point = 0.25 + rand*0.1 + randn*0.1;
    a = @(f) -((NbrFrames-f)/NbrFrames - turning_point) * 9e-6 ;
end

% Generate translation coordinates
trans_sigma = 2;
trans_start = trans_sigma*randn(2,1);
trans_end = diag(-sign(trans_start)) * abs(trans_sigma*randn(2,1));
trans_grad = trans_end - trans_start;

for f = 1:NbrFrames
    % Apply distortion
    K = barrel_distortion(K,a(f));
    K_mask = round(K); %ceil(K-0.5);
    L = barrel_distortion(L,a(f));
    %J = barrel_distortion(J,a(f)); % Looks more natural (but is not GT)
    %J = conv2( im2double(K_mask) , H , 'same' ) + (K_mask .* P_N); 
    J = conv2( im2double(K_mask) , H , 'same' ) + (K_mask .* L); % Looks more natural

    %I_temp(:,:,f) = J;
    %BW_temp(:,:,f) = K_mask;

    % Apply translation
    t = exp(-f/(NbrFrames/2));
    tx = trans_start(1) + t*( trans_end(1)-trans_start(1) );
    ty = trans_start(2) + t*( trans_end(2)-trans_start(2) );
    I_temp(:,:,f) = imtranslate(J,[tx,ty]);
    BW_temp(:,:,f) = ceil(imtranslate(K_mask,[tx,ty]));
end
I = I_temp; BW = BW_temp;

% Show spectral plots
% [static_spectrum,static_freq] = RadialAvgPSD(P_N);
% hold on; plot(10*log10(static_freq),10*log10(static_spectrum),'b');

% Show image
% clf; imagesc(I); axis square; colormap gray;
% hold on; plot(x_interp+ImCenter,y_interp+ImCenter,'b-');
end

%% Barrel Distortion
% Reference: http://www.mathworks.com/help/images/examples/creating-a-gallery-of-transformed-images.html

function I_barrel = barrel_distortion(I,a)
% The variable 'a' controls how much barrel/pin cushion distortion there
% is. If a>0 then the image appears to shrink. If a<0 then the image
% appears to expand.

[nrows,ncols] = size(I);
[xi,yi] = meshgrid(1:ncols,1:nrows);
imid = round(size(I,2)/2);
xt = xi(:) - imid;
yt = yi(:) - imid;
[theta,r] = cart2pol(xt,yt);
s = r + a*r.^3;
[ut,vt] = pol2cart(theta,s);
u = reshape(ut,size(xi)) + imid;
v = reshape(vt,size(yi)) + imid;
tmap_B = cat(3,u,v);
resamp = makeresampler('linear','fill');
I_barrel = tformarray(I,[],resamp,[2 1],[1 2],[],tmap_B,0);

end

%% Perlin noise
% This function was modified to include persistence
% The original function was taken from:
% http://stackoverflow.com/questions/7347111/generate-procedural-perlin-noise-in-matlab
function im = perlin_noise(im,persistence)

if nargin == 1
    persistence = 0;
end

[n, m] = size(im);
i = 0;
w = sqrt(n*m);

while w > 3
    i = i + 1;
    d = interp2(randn(n, m), i-1, 'spline');
    if persistence == 0
        im = im + i * d(1:n, 1:m);
    else
        im = im + persistence^i * d(1:n, 1:m);
    end
    w = w - ceil(w/2 - 1);
end

end

%% Fractal Noise (to be precise, 1/f^p noise)
% A nice (visual) explanation can be found at http://paulbourke.net/fractals/noise/

function im = fractal_noise(im,persistence)

% If persistence wasn't defined, then pink noise (i.e. 1/f) is default
if nargin == 1
    persistence = 1;
end

% Extract size information
[N(1),N(2)] = size(im);
hr = (N(1)-1)/2;
hc = (N(2)-1)/2;

% Distance from center (normalized to 0:1)
[x,y] = meshgrid( [-hc:hc]/hc , [-hr:hr]/hr );
D = sqrt( x.^2 + y.^2 );

% Create a 1/f^p filter
H = 1 ./ D.^persistence; % filter
H = H / norm(H(:)); % normalize

im = real( ifft2( ifftshift( fft2(randn(N)) .* H )) );

end
