function [pat_raw] = rEBSD_TFS_Raw_Single(data_all, pat_num)
% Created by T. Ben Britton and edited by Tianbi Zhang, November 2024

% reads from TFS data - will give you back a pattern

s_size=256*256;
n = pat_num;

pat_raw=zeros(549,549);

range_min=1+(n-1)*512*512;
range_max=512*512+(n-1)*512*512;
data=data_all(range_min:range_max);

%seperate the chips
data_0=data(0*s_size+1:1*s_size); data_0=reshape(data_0,[256 256]);
data_1=data(1*s_size+1:2*s_size); data_1=reshape(data_1,[256 256]);
data_2=data(2*s_size+1:3*s_size); data_2=reshape(data_2,[256 256]);
data_3=data(3*s_size+1:4*s_size); data_3=reshape(data_3,[256 256]);
%arrange these in the array
pat_raw([1:256]+35,1:256)=rot90(data_2,3);
pat_raw(1:256,257:512)= rot90(data_1,2);
pat_raw([257:512]+35,[1:256]+35)= rot90(data_3,4);
pat_raw(257:512,[257:512]+35)=rot90(data_0, 2);


%cut out the boundary around the edge
pat_raw=pat_raw(36:512,36:512);
end