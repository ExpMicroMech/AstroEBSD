function [ ebsd,GrainData, PhaseData,cs_all ] = fGConstruct( Data,PhaseData,GrainID_Setup )
%FGCONSTRUCT Extract the reference locations and the grain data
%Utilises MTEX codes - working with v5


%read the phase based data
GrainData.num_phase=size(PhaseData,2);
GrainData.cifs=cell(GrainData.num_phase,1);
%             Data.colors={'r','g','b','y'};

%CMM Edit: Call a translator to correct names from Astro to MTEX cif.

[ GrainData,cs_all, PhaseData] = fGTranslate(GrainData, PhaseData);

% OLD CODE BELOW:
% for phasen = 1:GrainData.num_phase
%     %hard code here - should have a translation tool properly
%     if strcmpi(PhaseData(phasen).Name,'Titanium') == 1
%         GrainData.cifs{phasen}='Ti-Titanium-alpha.cif';
%     elseif strcmpi(PhaseData(phasen).Name,'Titanium-beta') == 1
%         GrainData.cifs{phasen}='Ti-Titanium-beta.cif';
%     elseif strcmpi(PhaseData(phasen).Name,'Ferrite, bcc (New)') == 1
%         GrainData.cifs{phasen}='Fe-Iron-alpha.cif';
%     else
%         error(['Phase = ' PhaseData(phasen).Name ' not in data base - please correct']);
%     end
%     cs=loadCIF(GrainData.cifs{phasen});
%                     cs.color=Data.colors{phasen};
%     if phasen == 1
%         cs_all=cs;
%     else
%         cs_all={cs_all,cs};
%     end
% end

%build an EBSD container for MTEX
prop.x = double(Data.XSample);
prop.y = double(Data.YSample);
prop.xi =double(Data.XBeam_Map);
prop.yi =double(Data.YBeam_Map);
prop.PMap=double(Data.PMap);

Phi1=double(Data.Phi1);
Phi=double(Data.Phi);
Phi2=double(Data.Phi2);
MADPhase=double(Data.MADPhase);
MADPhase(MADPhase==0)=1; %force the phase to be at least 1

prop.quality=double(Data.RadonQuality);

ori = rotation('Euler',Phi1,Phi,Phi2);

ebsd = EBSD(ori,MADPhase,[cs_all],'options',prop);
ebsd2 = ebsd('indexed');
%create the grains
[grains,ebsd2.grainId,ebsd2.mis2mean]=calcGrains(ebsd2,'angle',GrainID_Setup.mis_tol*degree);

%clean up the grain size - using the minimum threshold
grains2=grains(grains.grainSize > GrainID_Setup.GrainID_minsize);
num_grains=size(grains2,1);
%%  calculate the grains, reference points and dilate the areas
for n=1:num_grains
    %subdivide
    ebsd_grain=ebsd2(ebsd2.grainId==grains2(n).id);
    grain_inds=sub2ind(size(ebsd.prop.x),ebsd_grain.prop.yi, ebsd_grain.prop.xi);
    
    %subset the qualitydata
    quality=ebsd_grain.prop.quality;
    quality_min=quantile(quality,GrainID_Setup.QualTrim);
    quality_ok=find(quality>quality_min);
    grain_inds_qual_ok=grain_inds(quality_ok);
    
    %calculate the euclidian grain boundary distance
    grain_image=ones(size(ebsd.prop.x));
    grain_image(grain_inds) =0;
    grain_edis=bwdist(grain_image);
    grain_dis_set=grain_edis(grain_inds_qual_ok);
    
    grain_dis_min=quantile(grain_dis_set,GrainID_Setup.P_Edge_Dist);
    grain_dis_ok=find(grain_dis_set>=grain_dis_min);%CMM EDIT - might break
    grain_inds_qual_plusdis_ok=grain_inds_qual_ok(grain_dis_ok);
 
    %extract the orientations for this set
    ebsd_sub=ebsd2(grain_inds_qual_plusdis_ok);
    
    %put them into one "grain"
%     [grain_cloud,ebsd_sub.grainId]=calcGrains(ebsd_sub,180*degree);
%     
    mean_ori=mean(ebsd_sub.orientations);
    
    %find a point with an orientation close to this
    angle_set=angle(ebsd_sub.orientations,mean_ori)*180/pi;
    [~,point_ok]=min(angle_set);
    
    
    graindata_x=ebsd_grain.xi;
    graindata_y=ebsd_grain.yi;
    
    %dilate the grain data
    
    %dilate in X
    graindata_xp_x=repmat(graindata_x,1,GrainID_Setup.Edge_Dil+1)+repmat(0:GrainID_Setup.Edge_Dil,size(graindata_x,1),1);
    graindata_xm_x=repmat(graindata_x,1,GrainID_Setup.Edge_Dil+1)+repmat([-GrainID_Setup.Edge_Dil:0],size(graindata_x,1),1);
    
    %join the arrays
    graindata_xdil_y=repmat(graindata_y,1,2*(GrainID_Setup.Edge_Dil+1));
    graindata_xdil_x=[graindata_xm_x graindata_xp_x];
    
    %column
    graindata_xdil_y(:)=graindata_xdil_y;
    graindata_xdil_x(:)=graindata_xdil_x;
    
    %clean the overspan
    graindata_xdil_y=graindata_xdil_y(graindata_xdil_x > 0);
    graindata_xdil_x=graindata_xdil_x(graindata_xdil_x > 0);
    
    graindata_xdil_y=graindata_xdil_y(graindata_xdil_x < Data.NumPoints(1)+1);
    graindata_xdil_x=graindata_xdil_x(graindata_xdil_x < Data.NumPoints(1)+1);
    
    graindata_xdil_i=sub2ind(size(Data.PMap),graindata_xdil_y(:),graindata_xdil_x(:));
    %             graindata_xdil_xi=Data_InputMap.XBeam_Map(sub2ind(size(Data_InputMap.PMap),graindata_xdil_y(:),graindata_xdil_x(:)));
    
    graindata_xdil_iunique=unique(graindata_xdil_i);
    
    graindata_x=Data.XBeam_Map(graindata_xdil_iunique);
    graindata_y=Data.YBeam_Map(graindata_xdil_iunique);
    %             graindata_bx=Data_InputMap.XSample(graindata_xdil_iunique);
    %             graindata_by=Data_InputMap.YSample(graindata_xdil_iunique);
    %             graindata_p=Data_InputMap.PMap(graindata_xdil_iunique);
    %             graindata{n}=[graindata_p graindata_x graindata_y graindata_bx graindata_by];
    %
    
    %repeat for Y
    graindata_yp_y=repmat(graindata_y,1,GrainID_Setup.Edge_Dil+1)+repmat(0:GrainID_Setup.Edge_Dil,size(graindata_y,1),1);
    graindata_ym_y=repmat(graindata_y,1,GrainID_Setup.Edge_Dil+1)+repmat([-GrainID_Setup.Edge_Dil:0],size(graindata_y,1),1);
    
    %join the arrays
    graindata_ydil_x=repmat(graindata_x,1,2*(GrainID_Setup.Edge_Dil+1));
    graindata_ydil_y=[graindata_ym_y graindata_yp_y];
    
    %column
    graindata_ydil_x(:)=graindata_ydil_x;
    graindata_ydil_y(:)=graindata_ydil_y;
    
    %clean the overspan
    graindata_ydil_x=graindata_ydil_x(graindata_ydil_y > 0);
    graindata_ydil_y=graindata_ydil_y(graindata_ydil_y > 0);
    
    graindata_ydil_x=graindata_ydil_x(graindata_ydil_y < Data.NumPoints(2)+1);
    graindata_ydil_y=graindata_ydil_y(graindata_ydil_y < Data.NumPoints(2)+1);
    
    graindata_ydil_i=sub2ind(size(Data.PMap),graindata_ydil_y(:),graindata_ydil_x(:));
    graindata_ydil_iunique=unique(graindata_ydil_i);
    
    
    graindata_x=Data.XBeam_Map(graindata_ydil_iunique);
    graindata_y=Data.YBeam_Map(graindata_ydil_iunique);
    graindata_bx=Data.XSample(graindata_ydil_iunique);
    graindata_by=Data.YSample(graindata_ydil_iunique);
    graindata_p=Data.PMap(graindata_ydil_iunique);
    
    GrainData.RefPoint(n)=ebsd_sub(point_ok);

    %populate the arrays
    GrainData.GrainPts{n}=[graindata_p graindata_x graindata_y graindata_bx graindata_by graindata_ydil_iunique];
    
end
end

