function [x_syn_data] = GenerateRandomCellShape(pca_data,NbrCoeffs,NbrSyntheticCells)

% This function generates synthetic cells using coefficient distribution
% Inputs:
% pca_data.b            Coefficients of training examples
% pca_data.V            Eigen-vectors
% pca_data.x_bar        Data mean
% NbrCoeffs             Number of coefficients to use
% NbrSyntheticCells     Number of Synthetic cells
%
% Outputs:
% x_syn_data            

% Extract key variables
b = pca_data.b;
V = pca_data.V;
x_bar = pca_data.x_bar;
d = size(b,1); % number of dimensions
N = size(b,2); % number of examples

dd = NbrCoeffs; % sub-dimensions of interest (select any dd < d)

x_syn_data = zeros(d,NbrSyntheticCells);

%% Covariance method (this assumes gaussian distributed coefficients)
% Unfortunately, this is an incorrect assumption.

% % Generate data covariance matrix S
% S = zeros(d);
% for i = 1:N
%     S = S + b(:,i) * b(:,i)';
% end
% S = S/(N-1);
% 
% % Obtain the cholesky factor T (for faster random generation later on)
% % Read up the function mvnrnd() for details.
% T = cholcov(S);
% 
% c = 1;
% while c <= NbrSyntheticCells
%     % Generate a random cell-shape
%     b_syn = T' * randn(size(T,1),1);
%     
%     % Resynthetize to coordinates
%     x_syn = x_bar + V * b_syn;
%     
%     % Plot
% %     plot( x_syn([d/2+1:d,d/2+1]) , x_syn([1:d/2,1]) , 'bo-' );
% %     title(['Cell #' num2str(c) ', cross-itself = ' num2str(cross_itself)]); pause;
%     if ~CheckCrossOver(x_syn), x_syn_data(:,c) = x_syn; c = c+1; end
% end

%% Kernel Density Estimation method
% This doesn't assume any gaussianity on the coefficient distributions

c = 1;
while c <= NbrSyntheticCells % Generate random cells
    
    % Generate random Eigen-values
    b_rand = zeros(d,1);

    for i = d:-1:(d-dd+1)
        clear pd;
        data = b(i,:)';
        %bw = 1.06 * min(std(data),iqr(data)/1.34) * N^-0.2;
        %pd = fitdist(data,'Kernel','Kernel','normal','Width',bw);
        pd = fitdist(data,'Kernel');
        b_rand(i) = random(pd);
    end

    % Create synthetic cell
    x_syn = x_bar + V * b_rand;

    if ~CheckCrossOver(x_syn), x_syn_data(:,c) = x_syn; c = c+1; end
end
c = c - 1;

disp(['Total cells generated: ' num2str(size(x_syn_data,2))]);
