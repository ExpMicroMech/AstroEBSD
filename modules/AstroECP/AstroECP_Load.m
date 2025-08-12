function [ECP_Pat] = AstroECP_Load(Input_Data)
%ASTROECP_LOAD Load the experimental ECP
%
% Input_Data.image_folder = folder with the image
% Input_Data.image_frame = file name
%
% Set up to load the Tescan data - OR - a regular RGB image

try
    ExpImage_Filename=fullfile(Input_Data.image_folder,Input_Data.image_name);

    if isfield(Input_Data,'image_frame')
        ECP_Pat.frame=Input_Data.image_frame;
        [ExpImage_Image,data1,td] = TescanFrame_Load(ExpImage_Filename,Input_Data.image_frame);
        ECP_Pat.pattern=flipud(double(ExpImage_Image));
        ECP_Pat.filename=ECP_Pat;
        ECP_Pat.size = size(ExpImage_Image); %size of the library patterns and the resize of the raw EBSP
    else
        try %try to load it as a TFS frame
            ExpImage_Image=imread(ExpImage_Filename);
            info1 = imfinfo(ExpImage_Filename);
            i1_text=info1.UnknownTags.Value;
            i1_text_pairs=regexp(i1_text, '[\f\n\r]', 'split');

            if size(i1_text_pairs,2) > 2 %assume that this is a TFS file if this loader works
                nonEmptyCells = ~cellfun('isempty', i1_text_pairs);
                i1_text_pairs_full=i1_text_pairs(nonEmptyCells);

                % [ResolutionX] =
                % str2double(fReadTFSHeaderPair('ResolutionX',i1_text_pairs_full));
                % %not needed
                [ResolutionY] = str2double(fReadTFSHeaderPair('ResolutionY',i1_text_pairs_full));
                ExpImage_Image=ExpImage_Image(1:ResolutionY,:); %crop the data bar

                ECP_Pat.pattern=double(flipud(ExpImage_Image));
                ECP_Pat.size =size(ECP_Pat.pattern);
                ECP_Pat.frame=0;
            end

        catch

            %load this frame - normal image loader
            warning('An image is loaded with limited other information')
            ExpImage_Image=imread(ExpImage_Filename);
            if size(ExpImage_Image,3) == 3
                ExpImage_Image=rgb2gray(Input_Data.ECP_Pat);
            end
            ECP_Pat.pattern=double(flipud(ExpImage_Image));
            ECP_Pat.size =size(ECP_Pat.pattern);
            ECP_Pat.frame=0;
        end
    end

catch
    warning('Failed to load the experimental ECP, creating a blank frame');
    ECP_Pat.pattern=zeros([100,100]);
    ECP_Pat.size =size(ECP_Pat.pattern);
end

% set some default filters
ECP_Pat.Settings_Cor.gfilt=0; %use a low pass filter
ECP_Pat.Settings_Cor.gfilt_s=7; %low pass filter sigma
ECP_Pat.Settings_Cor.radius=0;
ECP_Pat.Settings_Cor.radius_frac=0.6;  %smaller radius to crop the black portion and to avoid abberation effects on sides 

ECP_Pat.Settings_Cor.max_var_pc_x=0;
ECP_Pat.Settings_Cor.max_var_pc_y=0;
ECP_Pat.Settings_Cor.max_var_pc_z=1;
ECP_Pat.Settings_Cor.max_var_ori=3*degree;
ECP_Pat.Settings_Cor.pattern_crop_factor=4;  % 1 uses orginal resolution, 2 reduces it to half with faster processing

[ECP_Pat.ECP_Pat_BG,Settings_Cor_out] = EBSP_BGCor( ECP_Pat.pattern,ECP_Pat.Settings_Cor );

ECP_Pat.Settings_Cor_refine=Settings_Cor_out;
ECP_Pat.Settings_Cor_refine.radius=1;
ECP_Pat.Settings_Cor_refine.radius_frac=0.6;
ECP_Pat.Settings_Cor_refine.resize=1;

% ECP_Pat.size=size(Input_Data.ECP_Pat)/pattern_crop_factor;
% ECP_Pat.size=uint32(pattern_info.size);
% ECP_Pat.Settings_Cor_refine.size=pattern_info.size;
