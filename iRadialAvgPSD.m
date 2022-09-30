% This function reconstructs a synthetic FFT spectrum from a given 
% Radially Averaged PSD

function [F] = iRadialAvgPSD(PSD)

spectral_size = length(PSD(2:end)) * 2; % size of spectrum (do not count DC component)
center = [mean(1:spectral_size),mean(1:spectral_size)]; % get center coordinates

radius_space = 0 : 1 : round( 0.5 * spectral_size );
[X,Y] = meshgrid(1:spectral_size,1:spectral_size);

% Reconstructing the spectrum based on given PSD size (later, we'll resize)
F = ones(spectral_size,spectral_size) * PSD(end);

for r = length(radius_space):-1:2 % reverse order (exclude DC)
    % Set radius
    radius = radius_space(r);
    
    % Get a mask of pixels within circle of the radius
    XY = (X - center(2)).^2 + (Y - center(1)).^2 <= radius^2;
    
    % Assign to all pixels less than radius
    F(XY) = PSD(r);
end

% Get non DC-centered FFT
F = ifftshift(F);

% Assign the DC component
F(1,1) = PSD(1);

end