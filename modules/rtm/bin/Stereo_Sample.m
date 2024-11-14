function [ebsp_i] = Stereo_Sample(xjs,yjs,zjs,screen_int_s)
%STEREO_SAMPLE Summary of this function goes here
%   Detailed explanation goes here

r_s=[xjs(:),yjs(:),zjs(:)];
r_s_norm=r_s./repmat(sqrt(dot(r_s,r_s,2)),1,3);

rn_x=r_s_norm(:,1);
rn_y=r_s_norm(:,2);
rn_z=r_s_norm(:,3);

s_inds_p=find(rn_z >= 0);
s_inds_m=find(rn_z  < 0);

pn_xp= rn_x(s_inds_p)./(1+rn_z(s_inds_p));
pn_yp= rn_y(s_inds_p)./(1+rn_z(s_inds_p));

pn_xm=-rn_x(s_inds_m)./(1-rn_z(s_inds_m));
pn_ym=-rn_y(s_inds_m)./(1-rn_z(s_inds_m));

ebsp_p=screen_int_s.pp(pn_xp,pn_yp);
ebsp_m=screen_int_s.pm(pn_xm,pn_ym);

ebsp_i=zeros(size(rn_x));
ebsp_i(s_inds_p)=ebsp_p;
ebsp_i(s_inds_m)=ebsp_m;


end

