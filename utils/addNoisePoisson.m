function [imgNoisyRes,noiseratio] = addNoisePoisson(I, scaling)
    % 
    % [imgNoisyRes,noiseratio] = addNoisePoisson(I, scaling)
    %    I = (nxm image)
    %    sacling =  a factor (scalar)
    %    a change in the scaling by one should correspond to roughly a
    %    change of the noiseratio by 5
    % Adds Poisson Nosie to an image.
    % The scaling is complicatedly implemented in matlab.
    % Through trial and error, we found a way to circumvent to much binning
    % of the color values, by initially rescaling the intensity values
    % between 2^32 and the lower intensity value.
    % We notice, that after resacling the measured noise values seem to be
    % accruate
    % example values 
    % Scaling   17 18 19 20 21 22 23 24 25 26 27 28 29 30
    % Noise     10 12 15 19 23 32 37 42 48 57 62 68 75 80
    %
    %   TODO: The function should be rename. As Poisson noise is "applied"
    %   and not "added
    % L.Berners, 12.09.2025
    maxI = 2^16-1; %-> would corespond to uint16 (we use doulbe internally and 1e12 scaling
    IRescale = (rescale(I,2^32-2^scaling,2^32-1))/1e12;
    noisyI = imnoise(IRescale,"poisson");

    imgNoisyRes = rescale(noisyI,1,maxI);
    imgRescaleRescale = rescale(IRescale,1,maxI);

    % noiseratio = psnr(imgNoisyRes,imgRescaleRescale);
    meanSE = immse(imgNoisyRes,imgRescaleRescale);
    noiseratio = 10*log10((maxI)^2/meanSE);
end


% Use the following code, to test, how the noise is
% % % % % 
% % % % % I = imread("eight.tif");
% % % % % 
% % % % % noisyI = imnoise(I,"poisson");
% % % % % figure
% % % % % imagesc(noisyI)
% % % % % colormap gray
% % % % % noiseratio = psnr(noisyI,I);
% % % % % %
% % % % % figure
% % % % % histogram(noisyI,'DisplayName','noisy Image')
% % % % % hold on
% % % % % histogram(I,'DisplayName','image');
% % % % % title(sprintf("PSNR: %.4f",noiseratio))
% % % % % legend()
% % % % % xlim([0,255])
% % % % % %%
% % % % % maxI = 2^16-1
% % % % % powNum = 22
% % % % % mkdir("TestNoisyImgs")
% % % % % for powNum = [16:32]
% % % % %     IRescale = (rescale(I,2^32-2^powNum,2^32-1))/1e12;
% % % % %     noisyI = imnoise(IRescale,"poisson");
% % % % %     figure
% % % % %     imagesc(noisyI);
% % % % %     colormap gray
% % % % %     imgNoisyRes = rescale(noisyI,1,maxI);
% % % % %     imgRescaleRescale = rescale(IRescale,1,maxI);
% % % % %     % noiseratio = psnr(imgNoisyRes,imgRescaleRescale);
% % % % %     meanSE = immse(imgNoisyRes,imgRescaleRescale)
% % % % %     imwrite(rescale(imgNoisyRes),sprintf("TestNoisyImgs/Image_PSNR_%04.4f.png",noiseratio));
% % % % %     noiseratio = 10*log10((maxI)^2/meanSE);
% % % % %     figure
% % % % % 
% % % % %     histogram(imgNoisyRes,'DisplayName','noisy Image');
% % % % %     hold on
% % % % %     histogram(imgRescaleRescale,'DisplayName','image');
% % % % %     title(sprintf("PSNR: %.4f",noiseratio))
% % % % %     legend()
% % % % %     text(0.2,0.9,sprintf("MSE: %4e",meanSE),'Units','normalized')
% % % % %     text(0.2,0.8,sprintf("spread: 2e%i",powNum),'Units','normalized')
% % % % %     exportgraphics(gcf,sprintf("TestNoisyImgs/PSNR_%03.4f.png",noiseratio))
% % % % %     close all
% % % % % end