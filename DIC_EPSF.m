function EPSF = DIC_EPSF(M,shear_angle,epsf_sigma)

% Generates an Effective Point Spread Function (EPSF) given parameters:
%   M               : Size of blurring MxM kernel (typ M=5) 
%   shear_angle     : Shear angle in degrees
%   epsf_sigma      : Spread parameter sigma

[xg,yg] = meshgrid( [1:M] , [1:M] );
xg = xg - mean(xg(1,:));
yg = yg - mean(yg(:,1));
shear_angle_rad = shear_angle / 180 * pi;
exponent = exp( -(xg.^2 + yg.^2) / epsf_sigma^2 );
EPSF = - xg .* exponent * cos(shear_angle_rad) ...
       - yg .* exponent * sin(shear_angle_rad);

end