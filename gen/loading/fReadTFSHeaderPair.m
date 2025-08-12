function [text_retrieved] = fReadTFSHeaderPair(text_to_find,tif_info_string)
%FREADTFSHEADERPAIR Read a TFS extended information header pair of data
% INPUTS
% tif_info_string the string of the file info
% text_to_find - string of the text that you want to find, do not add the '=' symbol
% 
% OUTPUTS
% text_retrieved - text after the equals sign, in a string

arguments (Input)
    text_to_find
    tif_info_string
end

arguments (Output)
    text_retrieved 
end

text_to_find_eq=[text_to_find '='];

contains_text = cellfun(@(x)~isempty(strfind(x,text_to_find_eq)), tif_info_string, 'UniformOutput', false);
contains_text=cell2mat(contains_text);
text_found=tif_info_string{contains_text};
text_retrieved_sloc=strfind(text_found,text_to_find_eq);
text_retrieved=text_found(text_retrieved_sloc+numel(text_to_find_eq):end);

end