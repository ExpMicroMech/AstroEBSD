function setupRTMSinglePats(rtmConfigRefine,InputUser,MapInfo,MicroscopeData)


    % switch nargin

    imageFiles = dir(fullfile(InputUser.PatternFolder,'*'))
    allImgs = []
    tileNum = 1
    if isfield(InputUser,"PatternList")
        numPatterns = numel(InputUser.PatternList)
    else
        numPatterns = numel(imageFiles)
    end
    for imgNum = 1:numPatterns
        try 
            if isfield(InputUser,"PatternList")
                img = imread(fullfile(InputUser.PatternFolder,InputUser.PatternList(imgNum)));
                img = double(flipud(img));
            else
                img = double(flipud(imread(fullfile(imageFiles(imgNum).folder,imageFiles(imgNum).name))));
            end
            if size(img,3)>1
                img = mean(img,3);
            end
            allImgs = cat(3,allImgs,img);
        catch
            warning('pattern %i could not be read',imgNum)
        end

    end
    %%
    numPatterns = size(allImgs,3)
    PC = zeros(3,numPatterns);
    Phase = zeros(1,numPatterns);
    try
        PC = [MapInfo.MapData.pcx'; MapInfo.MapData.pcy';MapInfo.MapData.pcz'];
    catch
        PC = repmat(InputUser.PC',1,numPatterns);
        warning(' PC from mapdata could not be loaed. setting static PC.')
    end
    % end
    %% Lets still call them PCA Patterns to make the scritps work.
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCAPatterns", ...
                    tileNum), allImgs);
    writeDataToH5(InputUser.h5name,sprintf("/RTM/RTMPCA/%i/PCFitted", ...
                    tileNum), PC);
end
%%
% end
