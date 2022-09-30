clearvars;

%% Generate Synthetic Cell data

% Obtain default parameters
CellParam = DefaultSyntheticCellParams;

% Set custom parameters
CellParam.SNR = -10;
CellParam.NbrFrames = 100;

[I_N,BW,I,I_DIC,B] = GenerateSyntheticCellSequence(CellParam);

%% View Images Sequence

for f = 1:CellParam.NbrFrames
    subplot(221); imagesc(I(:,:,f)); axis square; title(['Frame = ' num2str(f)]);
    hold on; contour(BW(:,:,f),[0,0],'b');

    subplot(222); imagesc(I_DIC(:,:,f)); axis square; title('I_{DIC}');
    hold on; contour(BW(:,:,f),[0,0],'b');

    subplot(223); imagesc(I_N(:,:,f)); axis square; title('I_N');

    subplot(224); imagesc(I_N(:,:,f)+B); axis square; title('I_N + B');
    
    colormap gray;
    drawnow;
end