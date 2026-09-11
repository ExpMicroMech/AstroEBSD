function generateSyntheticVariables(InputUser,MicroscopeData,SettingsCor,screen_int,options,noise,h5name)

predefinedOri = 0
if nargin == 5
    noise = 0
end
if nargin <= 6
    % this way we can overwrite
    h5name = InputUser.h5name % is necessary to work in parfor loops
end
if ~isfield(options,'noise')
    options.noise = 0
end
if isfield(options,'gmatrix')
    predefinedOri = 1
end

%%
disp(options)
disp(options.numOri)
if predefinedOri ==1
    gmatrix = options.gmatrix;
else
    gmatrix = rotation.rand(options.numOri,'Bunge');
end
gmatrix_conv = zeros(3,3,options.numOri);
gmatrix_in =  zeros(3,3,options.numOri);
gmatrix_out_offset =zeros(3,3, options.numOri);
gmatrix_conv_offset =zeros(3,3, options.numOri);

randomRotSmall = rotation.rand(options.numOri,'maxAngle',options.maxAngleOffset*degree,'Bunge');
% gmatrixOffset = gmatrix*randomRotSmall;
randomOffset = randn(3,options.numOri) * options.pcLimsGen;
randOffsetSmall = randn(3,options.numOri)*options.pcLims;
PC = options.pcStart' +randomOffset;
PCoffset = PC +randOffsetSmall;
Patterns = zeros(SettingsCor.screensize,SettingsCor.screensize,options.numOri);
PatternsNoisy = zeros(SettingsCor.screensize,SettingsCor.screensize,options.numOri);
PatternsOffset= zeros(SettingsCor.screensize,SettingsCor.screensize,options.numOri);
phaseId = ones(options.numOri);
angles = zeros(1,options.numOri);
pc_diff = zeros(3,options.numOri);
pc_diff_mean = zeros(3,options.numOri);
noiseRatio = ones(1,options.numOri)*1000; % 1000 means original
%%
for pNum = 1:options.numOri
    gmatrixDet = gmatrix(pNum).matrix();
    gmatrixDet= gmatrixDet'*InputUser.RTM.Rx(MicroscopeData.TotalTilt)*InputUser.EdaxConversionMatrix;
    gmatrixOffset = randomRotSmall(pNum)*gmatrix(pNum);
    gmatrixOffsetDet = gmatrixOffset.matrix();
    gmatrixOffsetDet= gmatrixOffsetDet'*InputUser.RTM.Rx(MicroscopeData.TotalTilt)*InputUser.EdaxConversionMatrix;
    % gmatrixOffset = gmatrix
    geometry = EBSP_Gnom(SettingsCor,PC(:,pNum));
    geometryOffset = EBSP_Gnom(SettingsCor,PCoffset(:,pNum));
    patOriginal = EBSP_gen(geometry,gmatrixDet,screen_int,screen_int.isHex);
    
    if noise > 0
        % Settings_CorNoise.noise = options.noise
        [patOriginalNoisy,noiseRatio(pNum)] = addNoisePoisson(patOriginal, noise);
    else
        patOriginalNoisy = patOriginal;
    end
    % patOriginal = EBSP_BGCor(patOriginal)
    patOffset = EBSP_gen(geometryOffset,gmatrixOffsetDet,screen_int,screen_int.isHex);
    gmatrix_conv(:,:,pNum) = gmatrixDet ;
    gmatrix_in(:,:,pNum) = gmatrix(pNum).matrix;
    gmatrix_conv_offset(:,:,pNum) = gmatrixOffsetDet;
    gmatrix_out_offset(:,:,pNum) = gmatrixOffset.matrix();
    Patterns(:,:,pNum) = patOriginal;
    PatternsNoisy(:,:,pNum) = patOriginalNoisy;
    PatternsOffset(:,:,pNum) = patOffset;
    angles(:,pNum ) = angle(gmatrixOffset,gmatrix(pNum))/degree;
    pc_diff(:,pNum) = PC(:,pNum) -PCoffset(:,pNum);
    pc_diff_mean(:,pNum) = sqrt(pc_diff(1,pNum)^2+ pc_diff(2,pNum)^2+pc_diff(3,pNum)^2);
    % PC(:,pNum) ==  pcStart-randomOffset'
end
% figure
% imagesc(patOriginal)
% colormap gray;
% axis equal;
%%
writeDataToH5(h5name, "/RTM/RTMPCA/1/PCAPatternsOrig",...
    Patterns);
writeDataToH5(h5name, "/RTM/RTMPCA/1/PCAPatterns",...
    PatternsNoisy);
writeDataToH5(h5name, "/RTM/RTMPCA/1/NoiseLevel",...
    noiseRatio);
writeDataToH5(h5name, "/RTM/RTMPCA/1/PatternsOffset",...
    PatternsOffset);
writeDataToH5(h5name, "/RTM/RTMPCA/1/gmatrixCommercialDetector",...
    gmatrix_conv_offset);
writeDataToH5(h5name,"/RTM/RTMPCA/1/gmatrixCommercial", ...
    gmatrix_out_offset);
writeDataToH5(h5name,"/RTM/RTMPCA/1/offset", ...
    randomRotSmall.matrix);
writeDataToH5(h5name, "/RTM/RTMPCA/1/gmatrixCommercialDetectorOrig",...
    gmatrix_conv);
writeDataToH5(h5name,"/RTM/RTMPCA/1/gmatrixCommercialOrig", ...
    gmatrix_in);
writeDataToH5(h5name,"/RTM/RTMPCA/1/PCFitted", ...
    PCoffset);
writeDataToH5(h5name,"/RTM/RTMPCA/1/pc_diff", ...
    pc_diff);
writeDataToH5(h5name,"/RTM/RTMPCA/1/pc_diff_mean", ...
    pc_diff_mean);
writeDataToH5(h5name,"/RTM/RTMPCA/1/angles", ...
    angles);



writeDataToH5(h5name,"/RTM/RTMPCA/1/PCFittedOrig", ...
    PC);
writeDataToH5(h5name,"/RTM/RTMPCA/1/Phase", ...
    PC); % the PatternRefineFFT function will epect a phase
end