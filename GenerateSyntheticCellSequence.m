function [I_N,BW,I,I_DIC,B,N_static,N_dyn] = GenerateSyntheticCellSequence(param)

% Abstract:
% This program will generate synthetic cells and the respective ground
% truth data.
%
% Parameter Inputs:
%   param
%   .imsize             N variable to generate N x N image (assumed square)
%   .NbrFrames          Number of frames. If 1, then static image.
%   .epsf_M             M variable to generate DIC imaging blurring M x M kernel
%   .epsf_shear_angle   The DIC imaging blurring shear angle
%   .epsf_sigma         The DIC imaging blurring shear weighting
%   .poiss_lambda       Dynamic noise poisson distribution lambda
%   .poiss_amp          Dynamic noise poisson distribution amplitude
%   .gauss_mu           Dynamic noise gaussian distribution mean
%   .gauss_sigma        Dynamic noise gaussian distribution std dev
%   .gauss_amp          Dynamic noise gaussian distribution amplitude
%   .SNR                SNR level (in dB) of noise
%
% Outputs:
%   I_N                 Sequence of noisy DIC images (static & dynamic noise)
%   BW                  Sequence of ground truth masks
%   I                   Sequence of noiseless cell surface in the pixel-space
%   I_DIC               Sequence of noiseless DIC images
%   B                   Single image of Bias
%   N_static            Single image of static noise
%   N_dyn               Sequence of dynamic noise
%
% Note: To generate noisy DIC image with bias, use: I_N + B;
%
% The data was generated using the framework described in:
% Cell Membrane Tracking in Living Brain Tissue using Differential
% Interference Contrast Microscopy
% J.Lee, I.Kolb, C.R.Forest, C.J.Rozell
% IEEE Transactions in Image Processing

%% Load pre-learned data

% Bias data, Radially Averaged PSD data
load('BiasStaticNoiseData64.mat');

% Load PCA cell-shape data
load('Cell_PCA_data.mat');

%% Load parameters

if nargin == 1
    imsize = param.imsize;
    NbrFrames = param.NbrFrames;
    epsf_M = param.epsf_M;
    epsf_shear_angle = param.epsf_shear_angle;
    epsf_sigma = param.epsf_sigma;
    PixelspaceParam = param.PixelspaceParam;
    poiss_lambda = param.poiss_lambda;
    poiss_amp = param.poiss_amp;
    gauss_mu = param.gauss_mu;
    gauss_sigma = param.gauss_sigma;
    gauss_amp = param.gauss_amp;
    SNR = param.SNR;
else
    % Set default parameters (as in the paper)
    imsize = 64;
    NbrFrames = 100;

    % Effective Point Spread Function Parameters
    epsf_M = 5;
    epsf_shear_angle = 225;
    epsf_sigma = 0.5;

    % Cell Surface Parameters
    PixelspaceParam.ImSize = imsize;
    PixelspaceParam.RotationAngle = rand*360; % ~Uniform(0,360)
    PixelspaceParam.ScalingRatio = 0.7 + rand*0.1; % ~Uniform(0.7,0.8)
    PixelspaceParam.Method = 'fractal'; % sumofgauss perlin fractal
    PixelspaceParam.Persistence = sqrt(2); % perlin fractal
    PixelspaceParam.NbrFrames = NbrFrames;
    PixelspaceParam.Motion = 'shrink-expand'; % asympotic shrink-expand expand-shrink

    % Dynamic (Sensor) Noise Distributions Parameters
    poiss_lambda = 1e-10;
    poiss_amp = 20.6;
    gauss_mu = 0;
    gauss_sigma = 0.0181;
    gauss_amp = 0.98;
    SNR = -1;
end

%% Step #1: Generate Random Cell

[x_syn_data] = GenerateRandomCellShape(PCA_data,size(PCA_data.b,1),1);

%% Step #2: Generate Cell Surface

d = size(x_syn_data,1);

% Extract data points
x = x_syn_data(1:d/2,1);
y = x_syn_data(d/2+1:end,1);

% Obtain image using coordinates
[I,BW] = EmbedCoordToImageSpace(x,y,PixelspaceParam);

%% Step #3: Convolution against DIC Imaging Kernel

EPSF = DIC_EPSF(epsf_M,epsf_shear_angle,epsf_sigma);

I_DIC = zeros(imsize,imsize,NbrFrames);
for f = 1:NbrFrames
    I_DIC(:,:,f) = conv2(I(:,:,f),EPSF,'same');
end

%% Step #4: Extract Random Bias (assumed static)
% The bias is a component that we linear which can be added at any point in
% the simulation. We provide this as a separate component B.

rand_idx = rand*size(Bias,3);
B = Bias(:,:,ceil(rand_idx));

%% Step #5: Generate Static Noise (Cellular noise)

rand_idx = ceil(rand*size(RAPSD,2));

% Obtain Static Noise
[G] = iRadialAvgPSD(RAPSD(:,rand_idx));
N_syn_phase = exp( 1i * 2 * pi * rand(size(G)) ); % use amplitude but randomize phase
N_syn_spectrum = sqrt(2) * G .* N_syn_phase;
N_static = real(ifft2(N_syn_spectrum));

N_dyn = zeros(imsize,imsize,NbrFrames);
N_syn = zeros(imsize,imsize,NbrFrames);
I_N = zeros(imsize,imsize,NbrFrames);
signal_energy = nan(NbrFrames,1);
noise_energy = nan(NbrFrames,1);
for f = 1 : NbrFrames
    % Generate Dynamic Noise
    poiss_pd = makedist( 'Poisson' , poiss_lambda );
    gauss_pd = makedist( 'Normal' , gauss_mu , gauss_sigma );
    N_dyn = poiss_amp * random(poiss_pd,size(I)) + ...
            gauss_amp * random(gauss_pd,size(I));
    N_dyn(:,:,f) = N_dyn(1:imsize,1:imsize);

    % Noise signal
    N_syn(:,:,f) = N_static + N_dyn(:,:,f);
    %N_syn(:,:,f) = N_dyn(:,:,f); % Only dynamic noise (useful for analysis)

    % First Determine the average SNR across all of the frames
    % Note: the noise energy is the measure of the biased variance.
    % Matlab uses: sum(abs(signal(:)).^2)/length(signal(:));
    signal_energy(f) = norm(I_DIC(:,:,f),'fro')^2; % add square for power!
    noise_energy(f) = norm(N_syn(:,:,f),'fro')^2;
end

% Use unbiased estimate for the overall noise

% Scale the original signal by the average energy
I = 10^(SNR/10) * (sum(noise_energy)/sum(signal_energy)) * I;
I_DIC = 10^(SNR/10) * (sum(noise_energy)/sum(signal_energy)) * I_DIC;

% Add signal to noise
I_N = I_DIC + N_syn;

end
