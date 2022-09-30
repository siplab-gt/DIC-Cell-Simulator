function param = DefaultSyntheticCellParams

% Set default parameters (as in the paper)
param.imsize = 64;
param.NbrFrames = 100;

% Effective Point Spread Function Parameters
param.epsf_M = 5;
param.epsf_shear_angle = 225;
param.epsf_sigma = 0.5;

% Cell Surface Parameters
PixelspaceParam.ImSize = param.imsize;
PixelspaceParam.RotationAngle = rand*360; % ~Uniform(0,360)
PixelspaceParam.ScalingRatio = 0.7 + rand*0.1; % ~Uniform(0.7,0.8)
PixelspaceParam.Method = 'fractal'; % sumofgauss perlin fractal
PixelspaceParam.Persistence = sqrt(2); % perlin fractal
PixelspaceParam.NbrFrames = param.NbrFrames;
PixelspaceParam.Motion = 'shrink-expand'; % asympotic shrink-expand expand-shrink

param.PixelspaceParam = PixelspaceParam;

% Dynamic (Sensor) Noise Distributions Parameters
param.poiss_lambda = 1e-10;
param.poiss_amp = 20.6;
param.gauss_mu = 0;
param.gauss_sigma = 0.0181;
param.gauss_amp = 0.98;
param.SNR = -1;
    
end