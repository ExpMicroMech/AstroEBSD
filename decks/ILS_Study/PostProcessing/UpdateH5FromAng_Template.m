%% 10.02.2026 L. Berners
% We want to plot the results from refining structures with different
% occupancies.
% TODO: In the future, it would be great to rework this script to work with
% arbitrary crystal structures. (e.g. load from cif and based on this plot
% everything)

% Clean up workspace
clear -regexp ^(?!.*(Phase_Folder)|eds_field_names|(Astro_loc)|rootDirOSC|MTEX_loc|rootpahtList|rootrootpath|mapFilterString|structureFilterString|additionalFilterString|rPNum|mapDirs|strucutre_dir|mapPNum|strucutre_dir|structPNum|clean_dir|h5File).
close all;
RTM.Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
RTM.Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
RTM.Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)]; %y rotation
InputUser.Phase_Folder=Phase_Folder;
InputUser.Astro_loc=Astro_loc;
warning('off','MATLAB:MKDIR:DirectoryExists');

% exportFolder = 'C:\Users\berners\sciebo\Paper6_OccupancyAndEffectsInNbNi\CompareMisalignment'
% exportFolder = "J:\Users\Berners_Lukas\clara\NbNi_septemberbatch\58Nb-42Ni_v2\Area7"

%%
InputUser.h5folder = clean_dir;
InputUser.ExportFolder =  fullfile(clean_dir,"patternExports")
mkdir(InputUser.ExportFolder)
exportFolderRoot = fullfile(InputUser.h5folder,'compareNCC');
mkdir(exportFolderRoot);
InputUser.h5refinement = h5File;
h5structure = h5info(InputUser.h5refinement,'/RTM');%,'/Patterns','min')
h5Struct = readStructFromH5(InputUser.h5refinement,"skipPatterns")
%%
InputUser.OSC_folder = strucutre_dir(structPNum).folder
phaseInfo = h5info(InputUser.h5refinement,"/UserInfo/Phase_Input"); 
Phase_Input = {}
phaseNum = 1
for phaseNum = 1:numel(phaseInfo.Datasets)
    Phase_Input{phaseNum} = h5read(InputUser.h5refinement,sprintf("/UserInfo/Phase_Input/%s",phaseInfo.Datasets(phaseNum).Name));
    phaseNum=phaseNum+1;
end
InputUser.Phase_Input = Phase_Input
InputUser.OSC_file = h5read(InputUser.h5refinement,sprintf("/UserInfo/OSC_file"));
% InputUser.OSC_file = fullfile(InputUser.OSC_file , sprintf("%s/cleaned",))
%
[ MapInfo.MapData,MicroscopeData,~,MapInfo.EBSPData,InputUser ] = loadEBSDEDAXClarity(InputUser)
%%
% this is bad this way but lets first see whether it works.
% this is bad this way but lets first see whether it works.
% maps = ["copper","bone"]
eds = struct();
fields =[]
for field = eds_field_names
    try
    %
    temp = reshape(MapInfo.MapData.props.(field),MicroscopeData.NCOLS,MicroscopeData.NROWS)';
    % temp(abs(temp-temp)>=2*std(temp(:))) = mean(temp(:)); %filter out values that are very off;
    eds.(field) =temp ;
    f = figure
    axis equal
    f = figure
    f.Position = [100 100 800 800]
    set(gcf,'windowstyle','docked')
    f.Position = [100 100 800 800]
    
    imagesc(imgaussfilt(eds.(field),1))
    axis equal
    colorbar()
    colormap('copper')
    clim([mean(temp(:))-2*std(temp(:)), mean(temp(:))+2*std(temp(:))])
    f.Position = [100 100 800 800]
    axis equal
    exportgraphics(f,fullfile(strucutre_dir.folder,sprintf("EDS_%s.png",field)))
    % close all
    fields= [fields field]
    % catch
        % fields = setdiff(fields,field)
    end
    % end
% a
end
close all
% figure
% imagesc(CPS)
% axis equal
% colormap gray
% figure
% imagesc(Ni./CPS)
% axis equal
% colormap gray
% figure
% imagesc(Nb./(Ni+Nb))
% axis equal
% colormap gray
% colorbar
% figure
% imagesc(Ni./(Ni+Nb))
% axis equal
% 
% colormap gray
% colorbar
%%
writeStructToH5(MapInfo.EBSPData.ebsd_map.prop,InputUser.h5refinement,sprintf("/MapInfo/EBSPData/props"))
tileInfo_map = h5read(InputUser.h5refinement ,sprintf('/RTM/tileInfo/MapReshaped2'));
rowfin = h5read(InputUser.h5refinement ,sprintf('/RTM/tileInfo/rowfin'));
colfin = h5read(InputUser.h5refinement ,sprintf('/RTM/tileInfo/colfin'));
rowstart = h5read(InputUser.h5refinement ,sprintf('/RTM/tileInfo/rowstart'));
colstart = h5read(InputUser.h5refinement ,sprintf('/RTM/tileInfo/colstart'));
%%
numTiles = numel(colstart);
for tileNum =1: numTiles
    compMap = tileInfo_map(rowstart(tileNum):rowfin(tileNum),colstart(tileNum):colfin(tileNum));
    PC_origIn = h5read(InputUser.h5refinement ,sprintf('/RTM/RTMPCA/%i/PCFitted',tileNum));
    numPatterns = size(PC_origIn,2);
    edsData = struct();
    % edsResults
    for field = fields
        edsValuesLoc = eds.(field);
        edsValuesLoc = edsValuesLoc(rowstart(tileNum):rowfin(tileNum),colstart(tileNum):colfin(tileNum));
        allEDSData = []
        allEDSData_std = []
        for patternNumber = 1:numPatterns
           allEDSData = [allEDSData, mean(edsValuesLoc(compMap(:)==patternNumber))];
           allEDSData_std =  [allEDSData_std, std(edsValuesLoc(compMap(:)==patternNumber))];
        end
        edsData.(field).mean = allEDSData;
        edsData.(field).std = allEDSData_std;
    end

        writeStructToH5(edsData,InputUser.h5refinement,sprintf("RTM/OriCorrectedRefinement/%i",tileNum));
    
        % end




        
        
        
        %%
        % InfoStruct = readStructFromH5(InputUser.h5refinement,'skipPatterns');
        % numTiles = numel(InfoStruct.RTM.tileInfo.colfin);
        if tileNum>1
            cotinue % not implemented for now.
        end
        Settings_Rad.theta_range=[-10 180 1]; %theta min, theta max, theta step - in degrees

        %peak hunt
        Settings_Rad.max_peaks=10; %max number of peaks to return
        Settings_Rad.num_peak=15; %number of peaks to search for - peaks will be rejected
        Settings_Rad.theta_search_pix=4; %search size in theta steps
        Settings_Rad.rho_search_per=0.1; %radon search in fractions
        Settings_Rad.min_peak_width=0.002; %seperation of the peak width, in pixels
        
        [rowVec,colVec] = meshgrid(1:size(compMap,1), 1:size(compMap,2)) ;
        SettingsCor = h5Struct.RTM.RTMRefined.Settings_Cor
        % find(tileInfo_map==1);
        iQs =  []
        iQs_exp_all = []
        % for tileNum =1: numTiles
        patOriginal = h5read(InputUser.h5refinement  ,sprintf("/RTM/RTMPCA/%i/PCAPatterns",tileNum));
            
        for pNum = 1:size(patOriginal,3)
            PCAPat = EBSP_BGCor(patOriginal(:,:,pNum),SettingsCor);
            iQs_exp = []
            patNums = find(tileInfo_map'==pNum);      %          patNums = find(tileInfo_map'==pNum);
            for i =1:10

                randomPattern = randi(numel(patNums));
                ExpPat = bReadEBSP(MapInfo.EBSPData,randomPattern);% loadCalibrationEDAX( fullfile( InputUser.OSC_folder,InputUser.OSC_file) );
                ExpPat = EBSP_BGCor(ExpPat,SettingsCor);
                [ ~,Peak_Quality,~,~,~,~,~ ]  = EBSP_RadHunt( ExpPat,Settings_Rad );
                iQs_exp = [iQs_exp Peak_Quality(1)];
                imwrite(flipud(rescale(ExpPat)),fullfile(InputUser.ExportFolder ,sprintf("PCAAndexpComparison_comp_%ipat_%i_%.3f.tif",pNum,randomPattern,Peak_Quality(1))))
            end
            iQs_exp_all = cat(2, iQs_exp_all, iQs_exp(:));
            % iQs_exp_all = cat(2, iQs_exp_all, iQs_exp(:));
            [ ~,Peak_Quality,~,~,~,~,~ ]  = EBSP_RadHunt( PCAPat,Settings_Rad );
            disp(Peak_Quality)
            iQs = [iQs Peak_Quality(1)]
            imwrite(flipud(rescale(ExpPat)),fullfile(InputUser.ExportFolder ,sprintf("PCAAndexpComparison_comp_%i.tif",pNum)))
        end
    
        %%
        IQDAta.IQs = iQs
        IQData.iQs_exp_all=iQs_exp_all
        IQData.Settings_Rad = Settings_Rad
        % figure
        % plot(iQs)
        writeStructToH5(IQDAta,InputUser.h5refinement,sprintf("RTM/OriCorrectedRefinement/%i",tileNum));
         % writeStructToH5(IQDAta,InputUser.h5refinement,sprintf("RTM/OriCorrectedRefinement/%i",tileNum));
    
        % end
end