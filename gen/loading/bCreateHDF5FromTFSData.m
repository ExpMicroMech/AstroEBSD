%% Example creation of a HDF5 file
% Takes an existing file and copies it, downsizing the EBSPs
% Not all the variables need copying 
% As now bReadHDF5 gives a warning when groups of data are missing

% location_astro='C:\Users\bbrit\Documents\GitHub\AstroEBSD\';
% run(fullfile(location_astro,'start_AstroEBSD.m'));

% InputUser.HDF5_folder='D:\RKD_Data\2023-11-16T12_36_34_WSe2_orig\2023-11-16T13_36_03'; 
InputUser.HDF5_folder='D:\RKD_Data\2023-11-23T13_50_44_test_WSE2\2023-11-23T16_24_52'; 
InputUser.HDF5_file='Demo_Tianbi_EBSD_raw.h5';

LoadProcsOpts.path_patterns = fullfile(InputUser.HDF5_folder, 'patterns/');
LoadRawOpts.path_raw_patterns = fullfile(InputUser.HDF5_folder, 'D:\RKD_Data\2023-11-23T13_50_44_test_WSE2\2023-11-23T16_24_52');

%read the existing H5 file
% [ MapData,MicroscopeData,PhaseData,EBSPData ] = bReadHDF5( InputUser );

TFSJson = ReadTFSJson(InputUser.HDF5_folder, 'site.json');

% RKD_load_TFS_float32_data;
%% Write the new one

OutputUser=InputUser;

%new variable names
OutputUser.HDF5_file=[OutputUser.HDF5_file(1:end-3) '.h5'];
%we assume that this is an h5 file, hdf5 will cause issues with the (end-3)
OutputUser.DataName=InputUser.HDF5_file(1:end-3);
OutputUser.HDF5FullFile=fullfile(OutputUser.HDF5_folder,OutputUser.HDF5_file);

dtype='/EBSD/Data/'; %EBSD data location
htype='/EBSD/Header/'; %Header data location (e.g. microscope settings)

EBSPData.HDF5_loc=fullfile(OutputUser.HDF5_folder,OutputUser.HDF5_file);
EBSPData.numpats=TFSJson.Info.TotalNumberOfPoints;

EBSPData.PatternFile = strcat('/',OutputUser.DataName, dtype, "RawPatterns");

%% Map Data
%pattern centre
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'DD',TFSJson.Info.PCZ); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PCX',TFSJson.Info.PCX); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PCY',TFSJson.Info.PCY); %write the data

%indexing data
% use readpgm if you have .pgm files - which directly calls imread and is
% faster.
TFS_MEC = imread(fullfile(InputUser.HDF5_folder, 'map_fft_quality.pgm'))';
TFS_IQ = imread(fullfile(InputUser.HDF5_folder, 'map_median_electron_count.pgm'))';


h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'NIndexedBands', TFS_MEC); %write the data - replaced by counts/MEC
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'RadonQuality', TFS_IQ); %write the data

placeholder_array = zeros(size(TFS_MEC));


% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'MAD',Results.XCF); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'MAD',placeholder_array); %write the data

% don't have this - use a placeholder
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'MADPhase',placeholder_array); %write the data

%beam data
ID = linspace(0,EBSPData.numpats-1,EBSPData.numpats)';
X = mod(ID,TFSJson.Info.Width);
Y = floor(ID/TFSJson.Info.Width);

h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'X BEAM',X); %subtract 1 from the data because of the 0 indexing
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'Y BEAM',Y); %subtract 1 from the data because of the 0 indexing

h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'X SAMPLE',((X)-1)*TFSJson.Info.MapStepSize); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'Y SAMPLE',((Y)-1)*TFSJson.Info.MapStepSize); %write the data

clear ID X Y
%% Map data - orientation
% read in TFS ang file
CS = {... 
  'notIndexed',...
  crystalSymmetry('622', [3.2 3.2 12], 'X||a', 'Y||b*', 'Z||c*', 'mineral', 'Tungstenite', 'color', [0.53 0.81 0.98])};

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

% path to files
TFSAngFileName = [InputUser.HDF5_folder '\Map1.ang'];

% create an EBSD variable containing the data
ebsd = EBSD.load(TFSAngFileName,CS,'interface','ang',...
  'convertEuler2SpatialReferenceFrame');

% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PHI',ebsd.orientations.Phi / degree); %write the data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'phi2',ebsd.orientations.phi2 / degree); %write the data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'phi1',ebsd.orientations.phi1 / degree); %write the data

% clear CS TFSAngFileName ebsd

% if not, use the results from AstroEBSD Indexing.
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'PHI',Results.EulerX); %write the data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'phi2',Results.EulerZ2); %write the data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,dtype,'phi1',Results.EulerZ1); %write the data

%% Header Data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'MADMax',MicroscopeData.MADMax); %write the data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'MaxRadonBandCount',MicroscopeData.MaxRadonBandCount); %write the data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'CameraTilt',MicroscopeData.CameraTilt); %write the data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'SampleTilt',MicroscopeData.SampleTilt); %write the data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'TotalTilt',MicroscopeData.TotalTilt); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'KV',TFSJson.Info.HighVoltage / 1e3); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NCOLS',TFSJson.Info.Width); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NROWS',TFSJson.Info.Height); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'XSTEP',TFSJson.Info.MapStepSize * 1e6); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'YSTEP',TFSJson.Info.MapStepSize * 1e6); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'NPoints',TFSJson.Info.TotalNumberOfPoints); %write the data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'PatternHeight',TFSJson.Info.Geometry.PatternHeight); %write the data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'PatternWidth',TFSJson.Info.Geometry.PatternWidth); %write the data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'Magnification',MicroscopeData.Magnification); %write the data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'MapStepFactor',MicroscopeData.MapStepFactor); %write the data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'SEPixelSizeX',MicroscopeData.SEPixelSizeX); %write the data
% h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'SEPixelSizeY',MicroscopeData.SEPixelSizeY); %write the data
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'WD',TFSJson.Info.WorkingDistance); %write the data

% SEM image - source might be different.
ReadTFSRawSEMFromPgm;
h5_WritePair(OutputUser.HDF5FullFile,OutputUser.DataName,htype,'SEM Image',SEMFromPgm); %write the data

%% Pattern data

binning=2;
% dtype='uint16';

% raw patterns
RKD_load_raw_for_h5
EBSD_load_processsed_for_h5

% or processed patterns
EBSD_load_processsed_for_h5


%% Read the data back - warnings will tell you if things have failed
[ MapData_new,MicroscopeData_new,PhaseData,EBSPData ] = bReadHDF5( OutputUser );

%check that we can read this
 [ EBSDPat_back ] = bReadEBSP(EBSPData,11952);

%%

figure;
imagesc(EBSDPat_back); axis xy; axis image; axis off; 
colormap('gray');
