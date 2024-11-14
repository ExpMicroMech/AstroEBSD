function [EBSDpat,EBSDpat_unprocessed]=loadEBSP_T(fname,EBSP_size,pattern_number)

%% 
% EBSP_PW=double(dataset_header..Pattern_Width);
% EBSP_PH=double(dataset_header..Pattern_Height);
% loc=header.fileloc;
% fname=header.fname;
% gnum=size(h5info(fname).Groups,1);
pname=h5info(fname).Groups(1).Name;

% if gnum>else
% 
Pat = [pname '/EBSD/Data/Processed Patterns'];
EBSDpat=double((rot90(shiftdim(h5read(fname,Pat,[1 1 pattern_number],[EBSP_size(1) EBSP_size(2) 1])))));
if nargout == 2
    try
        UPat = [h5info(fname).Groups(1).Name '/EBSD/Data/Unprocessed Patterns'];
        EBSDpat_unprocessed=double((rot90(shiftdim(h5read(fname,UPat,[1 1 pattern_number],[EBSP_size(1) EBSP_size(2) 1])))));
    catch
        EBSDpat_unprocessed=[];
    end
end
%have to install the filter, as per the answer here
% https://www.mathworks.com/matlabcentral/answers/393621-how-can-i-import-a-lzf-compressed-hdf5-dataset-in-matlab







