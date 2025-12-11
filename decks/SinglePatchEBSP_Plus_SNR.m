%% Analyze single patterns captured at different settings, with each stored in a different submap folder.
% Suitable for TFS TruePix system
home; clear; close all;

%%
InputUser.Astro_loc='C:\Users\benja\OneDrive\Documents\GitHub\AstroEBSD';
InputUser.mtex_location='C:\Users\benja\OneDrive\Documents\MATLAB\mtex-5.11.1';

%
%start mtex if needed
try EBSD;
catch
    run(fullfile(InputUser.mtex_location,"startup.m"));
end

%start AstroEBSD if needed
try astro_loadcheck;
catch
    run(fullfile(InputUser.Astro_loc,"start_AstroEBSD.m"));
end

%% file locations

TFSProjectPath = "C:\Users\benja\Documents\EBSD\TFS_SuperRes\Res_Map";

TFSProjectLocation_ExpTime = "C:\Users\benja\Documents\EBSD\TFS_SuperRes\ExpTime";  %Exposure series Data
TFSProjectPath_SR = "C:\Users\benja\Documents\EBSD\TFS_SuperRes\Res_Map\2025-08-14T11_28_18"; %Super Resolution Data
TFSProjectPath_Spatial = "C:\Users\benja\Documents\EBSD\TFS_SuperRes\Res_Map\2025-08-14T13_02_16"; %Higher spatial magnification pattern data
TFSProjectPath_Zoom = "C:\Users\benja\Documents\EBSD\TFS_SuperRes\Res_Map\2025-08-14T11_46_24"; %zoomed in data (detector retracted)

InputUser.Phase_Input  = {'Si_20kv_HR2'}; %the phase for Radon and pattern matching etc.

ResultsDir=['C:\Users\benja\Documents\MATLAB\Superresolution\SinglePatchEBSP_Plus_SNR\results'];

%pattern file data info - exposure series
rawPatFileName = "f_00000.pak";
procPatFileName = "row_00000.bin";

if isfolder(ResultsDir) == 0
    mkdir(ResultsDir)
end

% Low level setting up for phase indexing etc. stuff - you shouldn't need to change this
RTM.Phase_Folder = fullfile(InputUser.Astro_loc,'phases'); %location of the AstroEBSD phases super-folder
RTM.Bin_loc = fullfile(RTM.Phase_Folder,'dynamic_templates'); %location of the binary files used for RTM

% initial pattern centre, big patterns
PC_In_Start=[0.4948186 0.3136417 0.4913804]; %this is for the TFSProjectPath_SR data
%% Choose some settings
% Define signal and noise radial bands (in pixels) - use for the SNR
% calculations
% these are selected assumping that a 256 pixel pattern is being used
% so half of the range is 128
signalRange = [2, 96];
noiseRange  = [96, 128];
angleToExclude = pi/12;

%Super resolution - choose the point to match the data to
% ideally this should help with projection issues at the edge of the pattern
SR_ref_ebsd_xy=[25,30]; %select a point as the reference

SR_sfac=10; %this is the scale of upsampling for the higher magnification EBSP

%% Experiment 1 - exposure series
disp("Working on Exposure Series");

%load the project info
projectJson_exp = fullfile(TFSProjectLocation_ExpTime, "project.json");
projectInfo_exp = ReadTFSJson(TFSProjectLocation_ExpTime, "project.json");

%sort out the site list
num_sites_exp=size(projectInfo_exp.Sites,1);
site_ok=true(num_sites_exp,1);
site_dir=cell(num_sites_exp,1);
for n=1:num_sites_exp
    if strcmpi(projectInfo_exp.Sites(n).Name,'LivePoint')
        site_ok(n)=0;
    else
        site_dir{n}=projectInfo_exp.Sites(n).Name;
    end
end

site_dir=site_dir(site_ok);
num_sites_exp=size(site_dir,1);

%create the containers to count the data
allExposure = zeros(num_sites_exp,1);
peakHeights    = zeros(num_sites_exp, 1);
SNRValues      = zeros(num_sites_exp, 1);
radialProfiles = cell(num_sites_exp, 1);
stepVectors    = cell(num_sites_exp, 1);


% load information from one map to help build things
TFS_DataLoc=TFSProjectLocation_ExpTime + "\" + site_dir(1);
[ebsdtemp,ebsd_header] = bIDX_to_EBSD(TFS_DataLoc);

allProcPatternsExposure = zeros(ebsd_header.PatternHeight,ebsd_header.PatternWidth,num_sites_exp);
allPatternsFFTExposure = allProcPatternsExposure;

%read the patterns
for n = 1:num_sites_exp
    TFS_DataLoc=TFSProjectLocation_ExpTime + "\" + site_dir(n);
    [ebsdtemp,ebsd_header] = bIDX_to_EBSD(TFS_DataLoc);

    [loaded_TFS_patterns] = loadPat_TFS(TFS_DataLoc,[ebsdtemp.x ebsdtemp.y],1); 
    %loadPat_TFS([file location],[list of patterns, 0 indexed, two colums x,y], [1 = procssed, 2 = raw, 3 = both])
    allProcPatternsExposure(:,:,n)=loaded_TFS_patterns.pat_proc;

    allExposure(n) = (ebsd_header.ExposureTime - 0.0003) * 1e3; % ms, removing the read out time
end

%calcuate the FFTs and PS profiles
for n=1:num_sites_exp
    thisProcPattern = allProcPatternsExposure(:,:,n);
    
    % FFT power spectrum
    thisPatternFFT = abs(fftshift(fft2(thisProcPattern))).^2;
    allPatternsFFTExposure(:,:,n) = thisPatternFFT;

    %PH
    peakHeights(n) = max(thisPatternFFT(:));

    %Radial profile
    % [r_vec, radial_prof] = radial_profile(thisPatternFFT, 2);
    [stepVectors{n}, radial_prof] = radial_profile_ExcludeAngle(thisPatternFFT, 2, angleToExclude);
    radialProfiles{n}=radial_prof;
    % Compute SNR from profile
    sig_band = radial_prof(stepVectors{n} >= signalRange(1) & stepVectors{n} <= signalRange(2));
    noi_band = radial_prof(stepVectors{n} >= noiseRange(1)  & stepVectors{n} <= noiseRange(2));
    SNRValues(n) = mean(sig_band) / mean(noi_band);

end

% save the patterns
for n=1:num_sites_exp
    exposure_dur_str=sprintf('%04d',round(allExposure(n)));
    imwrite(normalizeto16bit(flipud(allProcPatternsExposure(:,:,n))), fullfile(ResultsDir,['ExpSeries_Pat_',exposure_dur_str,'ms.tif'])); %bg corrected
    imwrite(normalizeto16bit(flipud(log10(allPatternsFFTExposure(:,:,n)))), fullfile(ResultsDir,['ExpSeries_FFT_', exposure_dur_str,'ms.tif'])); %FFT
end

disp("Exposure series done");

%% Super resolution work

disp("Working on Super Resolution");

setMTEXpref('xAxisDirection','east');      %TFS
setMTEXpref('zAxisDirection','outofplane'); %TFS

% load data
[ebsd,ebsd_header] = bIDX_to_EBSD(TFSProjectPath_SR);
ebsd.scanUnit='um'; % update scan unit
ebsd.prop.pat_number=1:ebsd_header.TotalNumberOfPoints;

pattern_loc=[ebsd.x_in,ebsd.y_in]; %note this is in 0 based indexing, whole map
[loaded_TFS_patterns] = loadPat_TFS(TFSProjectPath_SR,pattern_loc,1); %loadPat_TFS([file location],[list of patterns, 0 indexed, two colums x,y], [1 = procssed, 2 = raw, 3 = both])
allProcessedPatterns=loaded_TFS_patterns.pat_proc;
clear loaded_TFS_patterns %just to tidy up memory
ebsd=gridify(ebsd);

%select the reference point for the superresolution experiment
ref_ebsd_pt=ebsd(SR_ref_ebsd_xy(2),SR_ref_ebsd_xy(1));

% Loop the patterns and workout scale etc
t_scale=zeros(1,ebsd_header.TotalNumberOfPoints);
t_translation = zeros(2, ebsd_header.TotalNumberOfPoints);
t_rotation = zeros(1, ebsd_header.TotalNumberOfPoints);
t_peak = zeros(1, ebsd_header.TotalNumberOfPoints);

disp('Measuring the warping function and regular pattern averaging')

%measure the distortions and dewarp

imagestack=zeros(ebsd_header.PatternHeight,ebsd_header.PatternWidth,ebsd_header.TotalNumberOfPoints);
imagestack_weight=zeros(ebsd_header.PatternHeight,ebsd_header.PatternWidth,ebsd_header.TotalNumberOfPoints)+1;

OutPutSize = imref2d([ebsd_header.PatternHeight,ebsd_header.PatternWidth],1,1);
weight_single=allProcessedPatterns(:,:,1)*0+1;

SR_pattern_reference=allProcessedPatterns(:,:,ref_ebsd_pt.pat_number);

parfor pat=1:ebsd_header.TotalNumberOfPoints
    [tform(pat),t_peak(pat)] = imregcorr(allProcessedPatterns(:,:,pat),SR_pattern_reference,Method="phasecorr");
    t_scale(:,pat)=tform(pat).Scale;
    t_translation(:,pat)=tform(pat).Translation;
    t_rotation(:,pat)=tform(pat).RotationAngle;

    temp_image= imwarp(allProcessedPatterns(:,:,pat),tform(pat),'bicubic',"OutputView",OutPutSize,FillValues=0);
    temp_weight=imwarp(weight_single,tform(pat),'bicubic',"OutputView",OutPutSize,FillValues=0);

    imagestack(:,:,pat) = temp_image;
    imagestack_weight(:,:,pat) = temp_weight;
end
disp('Warping function measured and dewarping (pixel) done')


%put the warp function in the EBSD container
ebsd.prop.t_rotation=rot90(reshape(t_rotation(1,:),ebsd_header.Width,ebsd_header.Height),1);
ebsd.prop.t_translation_x=rot90(reshape(t_translation(1,:),ebsd_header.Width,ebsd_header.Height),1);
ebsd.prop.t_translation_y=rot90(reshape(t_translation(2,:),ebsd_header.Width,ebsd_header.Height),1);
ebsd.prop.t_scale=rot90(reshape(t_scale(1,:),ebsd_header.Width,60),1);
pattern_number_check=1:size(t_scale,2);
ebsd.prop.pattern_number_check=rot90(reshape(pattern_number_check,ebsd_header.Width,60),1);

%stack these back together for a single pixel wise pattern
imagestack_collapse=sum(imagestack,3);
imagestack_weight_collapse=sum(imagestack_weight,3);
superres_single_norm=imagestack_collapse./imagestack_weight_collapse;
superres_single_norm_fft = abs(fftshift(fft2(superres_single_norm))).^2;
SR_pattern_reference_fft = abs(fftshift(fft2(SR_pattern_reference))).^2;

%write the patterns to disk
imwrite(normalizeto16bit(flipud(superres_single_norm)), fullfile(ResultsDir,'SuperRes_SubPix.tif')); %bg corrected
imwrite(normalizeto16bit(flipud(log10(superres_single_norm_fft))), fullfile(ResultsDir,'SuperRes_SubPix_FFT.tif')); %FFT
imwrite(normalizeto16bit(flipud(SR_pattern_reference)), fullfile(ResultsDir,'SuperRes_Ref.tif')); %bg corrected
imwrite(normalizeto16bit(flipud(log10(SR_pattern_reference_fft))), fullfile(ResultsDir,'SuperRes_Ref_FFT.tif')); %FFT

%plot these in matlab
figure;
nexttile;
imagesc(SR_pattern_reference);
axis equal; axis xy; colormap('gray'); axis tight;
title('Reference Pattern');

nexttile;
imagesc(superres_single_norm);
axis equal; axis xy; colormap('gray'); axis tight;
title('Super Resolution Pattern');

nexttile;
imagesc(log10(SR_pattern_reference_fft));
axis equal; axis xy; colormap('gray'); axis tight;
title('Reference Pattern FFT');

nexttile;
imagesc(log10(superres_single_norm_fft));
axis equal; axis xy; colormap('gray'); axis tight;
title('Super Resolution Pattern FFT');

%% upsample
disp('Upsampling started - higher resolution')

%prepare the stack
OutPutSize = imref2d([256*SR_sfac 256*SR_sfac],1/SR_sfac,1/SR_sfac);
imagestack_sfac_new=zeros([256*SR_sfac 256*SR_sfac]);
imagestack_weight_sfac_new=imagestack_sfac_new;

%do the upsampling
parfor pat=1:ebsd_header.TotalNumberOfPoints

        tform_n_in=simtform2d(t_scale(:,pat),0,t_translation(:,pat)); %force rotation to 0
        A_tform=tform_n_in.A;
        tform_n=simtform2d(A_tform);

    temp_image= imwarp(allProcessedPatterns(:,:,pat),tform(pat),'bicubic',"OutputView",OutPutSize,FillValues=0);
    temp_weight=imwarp(weight_single,tform(pat),'bicubic',"OutputView",OutPutSize,FillValues=0);

    imagestack_sfac_new=imagestack_sfac_new+temp_image;
    imagestack_weight_sfac_new = imagestack_weight_sfac_new+temp_weight;
end

%deal with the weights, and FFT
superres_large_norm=imagestack_sfac_new./imagestack_weight_sfac_new;
superres_large_norm_fft = abs(fftshift(fft2(superres_large_norm))).^2;

%extract the central FFT to check it is the same as the other method
fft_centre=(size(superres_large_norm_fft,1)/2);
fft_centrepix_x=1:ebsd_header.PatternHeight;
fft_centrepix_y=1:ebsd_header.PatternWidth;

ycen=fft_centre+fft_centrepix_y-ebsd_header.PatternHeight/2;
xcen=fft_centre+fft_centrepix_y-ebsd_header.PatternHeight/2;
superres_large_norm_fft_centre = superres_large_norm_fft(ycen,xcen);

imwrite(normalizeto16bit(flipud(superres_large_norm)), fullfile(ResultsDir,'SuperRes_Large.tif')); %bg corrected
imwrite(normalizeto16bit(flipud(log10(superres_large_norm_fft))), fullfile(ResultsDir,'SuperRes_Large_FFT.tif')); %FFT
imwrite(normalizeto16bit(flipud(log10(superres_large_norm_fft_centre))), fullfile(ResultsDir,'SuperRes_Large_cFFT.tif')); %FFT

disp('Upsampling finished - higher resolution')

figure;
nexttile;
imagesc(superres_large_norm); colormap('gray'); axis equal; axis tight; axis xy;
title('Super resolution - up sampled')
nexttile;
imagesc(log10(superres_large_norm_fft)); colormap('gray'); axis equal; axis tight;  axis xy;
title('Super resolution - up sampled FFT')
nexttile;
imagesc(log10(superres_large_norm_fft_centre)); colormap('gray'); axis equal; axis tight;  axis xy;
title('Super resolution - up sampled FFT-centre')
nexttile;
imagesc(log10(superres_single_norm_fft));
axis equal; axis xy; colormap('gray'); axis tight;
title('Super Resolution Pattern FFT');

%% subsample a number from the original data to check how the dose/upsampling works

disp('Dose sampling started')

num_log=10; %how many patterns to sum

%allocate some variables
numPatsToSample = round(logspace(0,log10(ebsd_header.TotalNumberOfPoints),num_log));
imagestack_dose=zeros([size(imagestack,1) size(imagestack,2),num_log]);
imagestack_doseii_fft=imagestack_dose;
SNRValues_SubPixel=zeros(1,num_log);

%go through the stack and find patterns that you like
rng default; %fix the rng to be the same each time
for n=1:num_log
    if n < num_log

        %create a random list of patterns, and make sure they are unique
        selectedPatterns = randi([1,ebsd_header.TotalNumberOfPoints], numPatsToSample(n), 1);
        selectedPatterns=unique(selectedPatterns);
        while numel(selectedPatterns) < numPatsToSample(n)
            disp(['For ' int2str(n) ' iteration adding patterns to list']) %notify if the search had to grow the the random additions

            selectedPatterns2 = randi([1,ebsd_header.TotalNumberOfPoints], numPatsToSample(n), 1);
            selectedPatterns_list=unique([selectedPatterns;selectedPatterns2]);
            if numel(selectedPatterns_list) >= numPatsToSample(n)
                selectedPatterns_list = selectedPatterns_list(randperm(length(selectedPatterns_list)));
                selectedPatterns =  selectedPatterns_list(1:numPatsToSample(n));
            end
        end
    else %because we want every pattern
        selectedPatterns = 1:ebsd_header.TotalNumberOfPoints;
    end

    %stack, normalize and FFT

    imagestack_doseii=sum(imagestack(:,:,selectedPatterns),3);
    imagestack_doseii_weight=sum(imagestack_weight(:,:,selectedPatterns),3)+0.00001; %add a small weight offset to deal with the divide by 0 issue

    imagestack_dose(:,:,n)=imagestack_doseii./imagestack_doseii_weight;

    imagestack_doseii_fft(:,:,n) = abs(fftshift(fft2(imagestack_dose(:,:,n)))).^2;

    % do the radial profiles

    [r_vec, radial_prof] = radial_profile_ExcludeAngle(imagestack_doseii_fft(:,:,n), 2, angleToExclude);

    stepVectors{n}    = r_vec;
    radialProfiles{n} = radial_prof;
    % Compute SNR from profile
    sig_band = radial_prof(r_vec >= signalRange(1) & r_vec <= signalRange(2));
    noi_band = radial_prof(r_vec >= noiseRange(1)  & r_vec <= noiseRange(2));

    SNRValues_SubPixel(n) = mean(sig_band) / mean(noi_band);
end

%plot this data
figure;
for n=1:num_log
    nexttile;
    imagesc(imagestack_dose(:,:,n)); axis image; axis xy; colormap('gray');
    title(int2str(numPatsToSample(n)))
end

for n=1:num_log
    nexttile;
    imagesc(log10(imagestack_doseii_fft(:,:,n))); axis image; axis xy; colormap('gray');
    title(int2str(numPatsToSample(n)))
end

%% Sum the patterns associated with the high spatial magnification area - frame integration experiment

disp('Higher spatial magnification sampling started')
% load data
[ebsd_hm,ebsd_header_hm] = bIDX_to_EBSD(TFSProjectPath_Spatial);
ebsd_hm.scanUnit='um'; % update scan unit
ebsebsd_hmd.prop.pat_number=1:ebsd_header_hm.TotalNumberOfPoints;

pattern_loc_hm=[ebsd_hm.x_in,ebsd_hm.y_in]; %note this is in 0 based indexing, whole map loading
[loaded_TFS_patterns] = loadPat_TFS(TFSProjectPath_Spatial,pattern_loc_hm,1); %loadPat_TFS([file location],[list of patterns, 0 indexed, two colums x,y], [1 = procssed, 2 = raw, 3 = both])
allProcessedPatterns_hm=loaded_TFS_patterns.pat_proc;
clear loaded_TFS_patterns %just to tidy up memory

sumMagnifiedPatterns = mean(allProcessedPatterns_hm, 3);

figure;
imagesc(sumMagnifiedPatterns); axis image off; colormap("gray"); axis xy; axis tight;
% imwrite(normalizeto16bit(sumMagnifiedPatterns), "sumMagnifiedPatterns.png");

sumMagnifiedFFT = abs(fftshift(fft2(sumMagnifiedPatterns))).^2;
% imwrite(normalizeto16bit(log10(sumMagnifiedFFT)), "sumMagnifiedPatternsFFT.png");

numPatsToSample_hm=10;
ranges_hm = round(logspace(0,log10(ebsd_header_hm.TotalNumberOfPoints),numPatsToSample_hm));

SNRValues_Magnified=zeros(1,numPatsToSample_hm);

for n = 1:numPatsToSample_hm

     if n < numPatsToSample_hm

        %create a random list of patterns, and make sure they are unique
        selectedPatterns_hm = randi([1,ebsd_header_hm.TotalNumberOfPoints], ranges_hm(n), 1);
        selectedPatterns_hm= unique(selectedPatterns_hm);
        while numel(selectedPatterns_hm) < ranges_hm(n)
            disp(['For ' int2str(n) ' iteration adding patterns to list']) %notify if the search had to grow the the random additions

            selectedPatterns2 = randi([1,ebsd_header_hm.TotalNumberOfPoints], ranges_hm(n), 1);
            selectedPatterns_list=unique([selectedPatterns_hm;selectedPatterns2]);
            if numel(selectedPatterns_list) >= ranges_hm(n)
                selectedPatterns_list = selectedPatterns_list(randperm(length(selectedPatterns_list)));
                selectedPatterns_hm =  selectedPatterns_list(1:ranges_hm(n));
            end
        end
    else %because we want every pattern
        selectedPatterns_hm = 1:ebsd_header_hm.TotalNumberOfPoints;
    end

    meanPattern = mean(allProcessedPatterns_hm(:,:,selectedPatterns_hm), 3);
    sumMagnifiedFFT = abs(fftshift(fft2(meanPattern))).^2;

    % [r_vec, radial_prof] = radial_profile(sumMagnifiedFFT, 2);
    [r_vec, radial_prof] = radial_profile_ExcludeAngle(sumMagnifiedFFT, 2, angleToExclude);
    stepVectors{numPatsToSample_hm+2}    = r_vec;
    radialProfiles{numPatsToSample_hm+2} = radial_prof;
    % Compute SNR from profile
    sig_band = radial_prof(r_vec >= signalRange(1) & r_vec <= signalRange(2));
    noi_band = radial_prof(r_vec >= noiseRange(1)  & r_vec <= noiseRange(2));

    SNRValues_Magnified(n) = mean(sig_band) / mean(noi_band);
end

%% plot the line graph of exposure and signal to noise
fig_lineplot=figure;
plot(numPatsToSample .* 300, SNRValues_SubPixel, '.-', 'MarkerSize',30, 'LineWidth', 2);

xlabel("Exposure Time (ms)");
ylabel("SNR (FFT radial average)");
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
hold on;
plot(allExposure, SNRValues(1:num_sites_exp), '.-', 'MarkerSize',30, 'LineWidth',2);

plot(300 .* ranges_hm, SNRValues_Magnified, '.-', 'MarkerSize',30, 'LineWidth', 2);

set(gca, 'FontSize', 20);
legend(["Super Resolution" "Single Pattern" "Frame Integration"], "Location","northwest");
grid on;

% saveas(fig_lineplot,fullfile(ResultsDir,'SNR_Graph.tif'))
%% Plot the functions for the warping

%plot this data
figure;
plot(ebsd,ebsd.x); title('x')
nextAxis
plot(ebsd,ebsd.y); title('y')
nextAxis
plot(ebsd,ebsd.prop.t_translation_x); title('x translate')
nextAxis
plot(ebsd,ebsd.prop.t_translation_y); title('y translate')
nextAxis
plot(ebsd,ebsd.prop.t_scale); title('t scale')
nextAxis
plot(ebsd,ebsd.prop.t_rotation); title('t rotation')
nextAxis;
plot(ebsd,ebsd.prop.pat_number); title('Pattern Number - raw data');
nextAxis
plot(ebsd, ebsd.prop.pattern_number_check); title('Pattern Num - model build');


%% Load the zoomed in pattern
disp('Working on the camera retraction experiment')
% load map data
[ebsd_z,ebsd_header_z] = bIDX_to_EBSD(TFSProjectPath_Zoom);
ebsd_hm.scanUnit='um'; % update scan unit
ebsebsd_hmd.prop.pat_number=1:ebsd_header_z.TotalNumberOfPoints;

%load pattern data
pattern_loc_z=[ebsd_z.x_in,ebsd_z.y_in]; %note this is in 0 based indexing, whole map loading
[loaded_TFS_patterns] = loadPat_TFS(TFSProjectPath_Zoom,pattern_loc_z,1); %loadPat_TFS([file location],[list of patterns, 0 indexed, two colums x,y], [1 = procssed, 2 = raw, 3 = both])
allProcessedPatterns_z=loaded_TFS_patterns.pat_proc;
clear loaded_TFS_patterns %just to tidy up memory

%average these together
averageZoomedInPatterns = mean(allProcessedPatterns_z, 3);

% Use the high magnification pattern and fit the other patterns onto this - via the warping functions

%warp the super res to the sub region
tfun_mag=imregcorr(superres_large_norm,averageZoomedInPatterns,Method="phasecorr");
%fix the magnification
sameAsInput = affineOutputView(size(averageZoomedInPatterns),tfun_mag,"BoundsStyle","SameAsInput");

%do the warp - super resolution large vs magnified
zoom_resize_factor=2.6340; %estimate factor for the zoom function between the camera retract pattern and full pattern
mega_zoomed=imwarp(superres_large_norm,tfun_mag,'bicubic',"OutputView",sameAsInput,FillValues=0);
frame_integral_max = imresize((mean(allProcessedPatterns_hm, 3)),2.6340,Method='lanczos3');

tfun_fi=imregcorr(frame_integral_max,averageZoomedInPatterns,Method="phasecorr");
sameAsInput2 = affineOutputView(size(averageZoomedInPatterns),tfun_fi,"BoundsStyle","SameAsInput");
fi_zoomed=imwarp(frame_integral_max,tfun_fi,'bicubic',"OutputView",sameAsInput2,FillValues=0);

%do the FFTs
fft_zoom=abs(fftshift(fft2(averageZoomedInPatterns))).^2;
fft_sr_zoom=abs(fftshift(fft2(mega_zoomed))).^2;
fft_fr_zoom=abs(fftshift(fft2(fi_zoomed))).^2;

%save these
imwrite(normalizeto16bit(flipud(averageZoomedInPatterns)), fullfile(ResultsDir,'Zoomed.tif')); %bg corrected
imwrite(normalizeto16bit(flipud(log10(fft_zoom))), fullfile(ResultsDir,'Zoomed_FFT.tif')); %FFT
imwrite(normalizeto16bit(flipud(mega_zoomed)), fullfile(ResultsDir,'Zoomed_SR_Up.tif')); %bg corrected
imwrite(normalizeto16bit(flipud(log10(fft_sr_zoom))), fullfile(ResultsDir,'Zoomed_SR_Up_FFT.tif')); %FFT
imwrite(normalizeto16bit(flipud(fi_zoomed)), fullfile(ResultsDir,'Zoomed_Reg.tif')); %bg corrected
imwrite(normalizeto16bit(flipud(log10(fft_fr_zoom))), fullfile(ResultsDir,'Zoomed_Reg_FFT.tif')); %FFT


% plot these comparisons
figure;
nexttile;
imagesc(averageZoomedInPatterns); title('Zoomed in integral')
axis image; axis xy; colormap('gray');
nexttile;
imagesc(mega_zoomed); title('Super Resolution - upsample')
axis image; axis xy; colormap('gray');
nexttile;
imagesc(fi_zoomed); title('Max Frame Intensity - Auto Resize')
axis image; axis xy; colormap('gray');


nexttile;
imagesc(log10(fft_zoom)); axis image; axis xy; colormap('gray');
c_fft=clim;
nexttile;
imagesc(log10(fft_sr_zoom)); axis image; axis xy; colormap('gray');
clim(c_fft); %fix the FFT colorscale
nexttile;
imagesc(log10(fft_fr_zoom)); axis image; axis xy; colormap('gray');
clim(c_fft); %fix the FFT colorscale

%% Radon index the pattern
disp('Radon indexing beginning')
%Radon based background correction settings
Settings_CorX.gfilt=1; %use a high pass filter (do you mean high pass?)
Settings_CorX.gfilt_s=5; %low pass filter sigma
Settings_CorX.radius=1; %use a radius mask
Settings_CorX.radius_frac=0.85; %fraction of the pattern width to use as the mask
Settings_CorX.hotpixel=0; %hot pixel correction
Settings_CorX.hot_thresh=1000; %hot pixel threshold
Settings_CorX.resize=0; %resize correction
Settings_CorX.RealBG=0; %use a real BG
Settings_CorX.EBSP_bgnum=30; %number of real pattern to use for BG
Settings_CorX.SquareCrop = 0; %make square the EBSP
Settings_CorX.SplitBG=0; %deal with a split screen

%Set up the radon transform peak finder
Settings_Rad.theta_range=[-10 180 1]; %theta min, theta max, theta step - in degrees
Settings_Rad.max_peaks=8; %max number of peaks to return
Settings_Rad.num_peak=20; %number of peaks to search for - peaks will be rejected
Settings_Rad.theta_search_pix=6; %search size in theta steps
Settings_Rad.rho_search_per=0.2; %radon search in fractions
Settings_Rad.min_peak_width=0.002; %min rseperation of the peak width, in pixels

% set up the phase
[ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num ] = Phase_Builder_RTM( InputUser.Phase_Input,RTM.Phase_Folder );
 
% Background correct - mostly this is to make the Radon work, as the very
% sharp patterns causes some issues with the Radon

[EBSP_superres_Radon.PatternIn,Settings_Cor ] = EBSP_BGCor( superres_single_norm,Settings_CorX);
EBSP_superres_Radon.PC=PC_In_Start;

%Uncomment to use the GUI to check the indexing
% Astro_EBSPset(EBSP_One,Settings_CorX,Settings_Rad,Settings_PCin,InputUser);

%Perform the Radon based indexing
%Radon transform
[ EBSP_superres_Radon.Peak_Centre,EBSP_superres_Radon.Single.Peak_Set_All,EBSP_superres_Radon.Peak_Set_All,...
            EBSP_superres_Radon.R_EBSP,EBSP_superres_Radon.R_Edge,EBSP_superres_Radon.R_rho,EBSP_superres_Radon.R_theta ] ...
            = EBSP_RadHunt( EBSP_superres_Radon.PatternIn,Settings_Rad);

% Convert the bands to normal space
[ EBSP_superres_Radon.nhat_gnom] = EBSP_NormConv( EBSP_superres_Radon.Peak_Centre,size(EBSP_superres_Radon.PatternIn),EBSP_superres_Radon.PC);
   
% Index
[EBSP_superres_Radon.rotdata{1},EBSP_superres_Radon.banddata]=EBSP_Index(EBSP_superres_Radon.nhat_gnom,Crystal_LUT{1},Settings_LUT{1}.thresh_trig,Crystal_UCell{1},eye(3));

%generate the geometry for analysis
[ EBSP_superres_Radon.PatternGeometry ] = EBSP_Gnom( Settings_Cor,EBSP_superres_Radon.PC );

%plot the result
EBSP_OneFigure=Plot_SinglePattern(EBSP_superres_Radon,Crystal_UCell,Crystal_LUT,1);

saveas(EBSP_OneFigure.figure,fullfile(ResultsDir,'Super_Res_RadonIndexed.tif'))

%% Now start some pattern matching
% create a reference pattern from the Dynamics bin file
[ ~,~,~,~,~, RTM_info ] = Phase_Builder_RTM(  {InputUser.Phase_Input{1}},RTM.Phase_Folder); %#ok<CCAT1>
[screen_int] = Cube_Generate(RTM_info.bin_file,RTM_info.isHex); %use the full pattern for rendering


%%
disp('Pattern matching - started')
%setttings for RTM
RTM.screensize = 256; %size of the library patterns and the resize of the raw EBSPs
RTM.Sampling_Freq=8; %Set the SO(3) sampling freq in degrees
RTM.iterations = 12;%Set the number of iterations to do in the refinement step
RTM.LPTsize = 256; %LPT size used in pixels

%Define all rotation matrices needed in the code
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation


[ SettingsXCF] = FFT_Filter_settings( RTM.screensize, RTM.LPTsize );
SettingsXCF.filters=[2 1 64 32]; %From fiddling these seem to be the best filter settings, but there may be better ones
SettingsXCF.single=1;
[SettingsXCF.FFTfilter,SettingsXCF.hfilter] = fFilters(SettingsXCF.roisize,SettingsXCF.filters);

%resize the experiment to make it easier to work with
Settings_Cor_RTM=Settings_Cor;
Settings_Cor_RTM.radius=0;
Settings_Cor_RTM.gfilt=1; %use a high pass filter (do you mean high pass?)
Settings_Cor_RTM.gfilt_s=3; %use a high pass filter (do you mean high pass?)

%generate the geometry for analysis
[ EBSP_superres_RTM1.PatternGeometry ] = EBSP_Gnom( Settings_Cor_RTM,EBSP_superres_Radon.PC );

%Take the experiment and adjust etc.
[EBSP_superres_RTM1.PatternIn] = EBSP_BGCor( superres_single_norm,Settings_Cor_RTM);

RTM.WindowOff=1;
%Get ready for RTM
[Pat_Ref_r] = refine_prep(EBSP_superres_RTM1.PatternIn,SettingsXCF,RTM);

%refine - using the Radon index to start
[G_Refined1,regout_R1] = refine5(Pat_Ref_r,EBSP_superres_RTM1.PatternGeometry,EBSP_superres_RTM1.PatternGeometry.PC,EBSP_superres_Radon.rotdata{1}.detector,SettingsXCF,screen_int,RTM_info.isHex,RTM);

%% Use pattern matching to refine the PC

disp('Pattern matching - PC refinement started')
Refine.ss=0.2; %initial probe volume in terms of the +/- PC values to search
Refine.p=2; %order of polynomial to fit to tetrahedron
Refine.n_its=100; %number of iterations
Refine.reindex=1; %index each time
Refine.print=1; %give an output of how the search is running
Refine.PHTol=1E-7; %refinement threshold to stop, if the PH doesn't change by more than this

%fPCRefine2(PatternIn,rotmat,PatternInfo,PCRefine_Settings,SettingsXCF,screen_int,RTM_setup)
[RTM_PC_Found] = fPCRefine2(Pat_Ref_r,G_Refined1,EBSP_superres_RTM1.PatternGeometry,Refine,SettingsXCF,screen_int,RTM);
%refine
%generate the geometry for analysis
[ EBSP_superres_RTM2.PatternGeometry ] = EBSP_Gnom( Settings_Cor_RTM,RTM_PC_Found );
[G_Refined2,regout_R2] = refine5(Pat_Ref_r,EBSP_superres_RTM2.PatternGeometry,EBSP_superres_RTM2.PatternGeometry.PC,G_Refined1,SettingsXCF,screen_int,RTM_info.isHex,RTM);

%% Simulate patterns
%simulate the pattern
[ EBSP_sim_H ] = EBSP_gen( EBSP_superres_Radon.PatternGeometry,EBSP_superres_Radon.rotdata{1}.detector,screen_int); %generate the EBSP for this iteration 

[ EBSP_sim_R1 ] = EBSP_gen( EBSP_superres_RTM1.PatternGeometry,G_Refined1,screen_int ); %generate the EBSP for this iteration 
[ EBSP_sim_R2 ] = EBSP_gen( EBSP_superres_RTM2.PatternGeometry,G_Refined2,screen_int ); %generate the EBSP for this iteration 

figure;
nexttile;
pPattern(EBSP_superres_RTM1.PatternIn,EBSP_superres_RTM1.PatternGeometry); title('Experiment')
nexttile;
pPattern(EBSP_sim_H,EBSP_superres_RTM1.PatternGeometry); title('Radon')
nexttile;
pPattern(EBSP_sim_R1,EBSP_superres_RTM2.PatternGeometry); title('RTM1');
nexttile;
pPattern(EBSP_sim_R2,EBSP_superres_RTM2.PatternGeometry); title('RTM2');
nexttile;
pPattern(EBSP_sim_R2-EBSP_sim_H,EBSP_superres_RTM2.PatternGeometry); title('Radon-RTM2');
nexttile;
pPattern(EBSP_superres_RTM1.PatternIn-EBSP_sim_H,EBSP_superres_RTM2.PatternGeometry); title('Exp-Radon');
nexttile;
pPattern(EBSP_superres_RTM1.PatternIn-EBSP_sim_R2,EBSP_superres_RTM2.PatternGeometry); title('Exp-RTM1');
nexttile;
pPattern(EBSP_superres_RTM1.PatternIn-EBSP_sim_R2,EBSP_superres_RTM2.PatternGeometry); title('Exp-RTM2');

%%
fig_comp=figure;
s1=subplot(1,3,1);
Plot_EBSPAnnotated( superres_single_norm,EBSP_superres_RTM2.PatternGeometry,[],G_Refined2,Crystal_UCell{1},Crystal_LUT{1}.family_norm_list,s1);
title('Exp, indexed')
s2=subplot(1,3,2);
Plot_EBSPAnnotated( EBSP_sim_R2,EBSP_superres_RTM2.PatternGeometry,[],G_Refined2,Crystal_UCell{1},Crystal_LUT{1}.family_norm_list,s2);
title('RTM2, indexed')

s2=subplot(1,3,3);
Plot_EBSPAnnotated( EBSP_sim_R2-EBSP_superres_RTM1.PatternIn,EBSP_superres_RTM1.PatternGeometry,[],G_Refined2,Crystal_UCell{1},Crystal_LUT{1}.family_norm_list,s2);
title('RTM2-Exp, indexed')

saveas(fig_comp,fullfile(ResultsDir,'Compared_Patterns.tif'))

disp('Code finished')
