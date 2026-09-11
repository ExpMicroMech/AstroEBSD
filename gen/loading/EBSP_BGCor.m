function [ EBSP2,Settings_Cor_out] = EBSP_BGCor( EBSP,Settings_Cor )
%EBSP_BGCOR Background correct the EBSP
%Use as [ EBSP2,Settings_Cor ] = EBSP_BGCor( EBSP,Settings_Cor )
%Inputs
%EBSP - array that contains the EBSP
%Settings_Cor - structure that contains information
%               on how backgroudn correction should be performed
%
%               As an example set:
%
% %background correction
% Settings_Cor.gfilt=1; %use a low pass filter
% Settings_Cor.gfilt_s=5; %low pass filter sigma
% 
% %radius mask
% Settings_Cor.radius=1; %use a radius mask
% Settings_Cor.radius_frac=0.9; %fraction of the pattern width to use as the mask
% 
% %hold pixel
% Settings_Cor.hotpixel=1; %hot pixel correction
% Settings_Cor.hot_thresh=1000; %hot pixel threshold
% 
% %resize
% Settings_Cor.resize=1; %resize correction
% Settings_Cor.size=150; %image width
% 
% Settings_Cor.RealBG=0; %use a real BG
% Settings_Cor.EBSP_bgnum=30; %number of real pattern to use for BG
% 
% Settings_Cor.SquareCrop=1; %crop to a square
%
%Settings_Cor
%   EBSP2 = corrected ESBP array (as double)
%   Settings_Cor = information from corrections as a structure
%

%% Versioning
%v1 - TBB 14/04/2017

%% Start Code
if isempty(Settings_Cor)
    Settings_Cor=struct;
end

EBSP2=EBSP;

%check fields exist & create if needed - this is ordered in the order of
%operations to aid with debugging & adding new correction routines 
%as needed

if ~isfield(Settings_Cor,'NaNClean')
   Settings_Cor.NaNClean=0;
end
if ~isfield(Settings_Cor,'InvSquareCrop')
   Settings_Cor.InvSquareCrop=0;
end
if ~isfield(Settings_Cor,'MeanCentre')
    Settings_Cor.MeanCentre=0;
end

if ~isfield(Settings_Cor,'SquareCrop')
    Settings_Cor.SquareCrop=0;
end

if ~isfield(Settings_Cor,'hotpixel')
    Settings_Cor.hotpixel=0;
end

if ~isfield(Settings_Cor,'resize')
    Settings_Cor.resize=0;
end

if ~isfield(Settings_Cor,'gaussfit')
    Settings_Cor.gaussfit=0;
end

if ~isfield(Settings_Cor,'hot_thresh')
    Settings_Cor.hot_thresh=0;
end


if ~isfield(Settings_Cor,'gfilt_s')
    Settings_Cor.gfilt_s=0;
end

if ~isfield(Settings_Cor,'blur')
    Settings_Cor.blur=0;
end

if ~isfield(Settings_Cor,'RealBG')
    Settings_Cor.RealBG=0;
end

if ~isfield(Settings_Cor,'radius_frac')
    Settings_Cor.radius_frac=1;
end

if ~isfield(Settings_Cor,'radius')
    Settings_Cor.radius=0;
end

if ~isfield(Settings_Cor,'gfilt')
    Settings_Cor.gfilt=0;
end

if ~isfield(Settings_Cor,'gfiltQuad')
    Settings_Cor.gfiltQuad= 0;
end
if ~isfield(Settings_Cor,'SplitBG')
    Settings_Cor.SplitBG=0;
end

if ~isfield(Settings_Cor,'Square')
    Settings_Cor.Square=0;
end

if ~isfield(Settings_Cor,'SatCor')
    Settings_Cor.SatCor=0;
end

if ~isfield(Settings_Cor,'TKDBlur2')
    Settings_Cor.TKDBlur2=0;
end

if ~isfield(Settings_Cor,'LineError')
   Settings_Cor.LineError=0; 
end

if ~isfield(Settings_Cor,'MeanCentre')
   Settings_Cor.LineError=0; 
end
if ~isfield(Settings_Cor,'cross')
   Settings_Cor.cross=0; 
end
if ~isfield(Settings_Cor,'gaussBlur')
    Settings_Cor.gaussBlur=0;
end

if ~isfield(Settings_Cor,'gaussian_range')
    Settings_Cor.gaussian_range=0;
end
if ~isfield(Settings_Cor,'quadrants')
   Settings_Cor.quadrants=0; 
end

if ~isfield(Settings_Cor,'rkd')
    Settings_Cor.rkd=0;
end



%% Start the corrections
% if Settings_Cor
if Settings_Cor.SquareCrop == 1 %crop the image to a square
    EBSP_size=min(size(EBSP2));
    crop.centre = [size(EBSP2, 2)/2 size(EBSP2, 1)/2];
    s=1:EBSP_size; s=s-mean(s,"omitnan");

    sy=s+crop.centre(2)+0.5; 
    sx=s+crop.centre(1)+0.5;
    
    %helps deal with issue with when the patterns are not even in size to start
    sx=floor(sx);
    sy=floor(sy);
    EBSP2=EBSP2(sy,sx);
end
%%
if Settings_Cor.InvSquareCrop == 1 %crop the image to a square
    [EBSP_size,xvsy]=max(size(EBSP));
    crop.centre = [size(EBSP, 2)/2 size(EBSP, 1)/2];
    s=1:EBSP_size; s=s-mean(s,"omitnan");
    s2=(max(size(EBSP))-min(size(EBSP)))/2+1:min(size(EBSP))+(max(size(EBSP))-min(size(EBSP)))/2;
    % s2=s2-mean(s2,"omitnan");
    if xvsy==2
        sy=s2;
        sx=s+crop.centre(1)+0.5;
    else
        sy=crop.centre(2)+0.5; 
        sx=s2;
    end
    
    %helps deal with issue with when the patterns are not even in size to start
    sx=floor(sx);
    sy=floor(sy);
    %%
    EBSP2=ones(EBSP_size,EBSP_size);
    %%
    % EBSP2=EBSP2*mean(EBSP)
    EBSP2(sy,sx)=EBSP;
end
%%
%%
if Settings_Cor.NaNClean == 1 %NaN correction
    [y,x]=find(isnan(EBSP2));
    
    [threexthreey,threexthreex]=meshgrid(-1:1,-1:1);
    
    if numel(y) > 0
        for s=1:numel(y)
            y_v=y(s)+threexthreey;
            x_v=x(s)+threexthreex;
            
            y_v(y_v<0)=[];
            x_v(x_v<0)=[];
            
            y_v(y_v<0)=[];
            x_v(x_v<0)=[];
            
            y_v(y_v>size(EBSP2,1))=[];
            x_v(x_v>size(EBSP2,2))=[];
                        
            sv=sub2ind(size(EBSP2),y_v(:),x_v(:));
            EBSP2(y,x)=median(EBSP2(sv),'omitnan');
        end
    end
    
end
%%
if Settings_Cor.SatCor == 1 %saturation correction
    I_noise= imnoise(zeros(size(EBSP2)),'gaussian');
    I_noise=mean(EBSP2(:),"omitnan")+std(EBSP2(:),"omitnan")*I_noise;
    mvals=find(EBSP2 == max(EBSP2(:)));
    EBSP2(mvals)=I_noise(mvals);
end
%%
if Settings_Cor.RealBG == 1
    % EBSP2=fix_mean(EBSP2);
    EBSP2=EBSP2./Settings_Cor.EBSP_bg;
end

if Settings_Cor.SplitBG == 1
    EBSPw=size(EBSP2,2);
    
    EBSP2a=EBSP2(:,1:EBSPw/2);
    EBSP2b=EBSP2(:,EBSPw/2+1:end);
    
    gf=Settings_Cor.gfilt_s*size(EBSP2a,1)/100;
    EBSP2a_g = imgaussfilt(EBSP2a,gf);
    EBSP2b_g = imgaussfilt(EBSP2b,gf);
    
    EBSP2=[EBSP2a./EBSP2a_g EBSP2b./EBSP2b_g];
   
end

%cor the pattern for hot pixels
if Settings_Cor.hotpixel == 1
    [EBSP2,Settings_Cor.hotpixl_num]=cor_hotpix(EBSP2,Settings_Cor.hot_thresh);
end

%%
if Settings_Cor.LineError==1
    
    meanval=mean(EBSP2(:),"omitnan");
    patstd=std(EBSP2(:),"omitnan");
    
    loc1=floor(0.9667*size(EBSP,1));
    loc2=ceil(0.98*size(EBSP,1));
    
    for i=loc1:loc2
        line=rand([1,EBSP_size]);
        EBSP2(i,:)=meanval+line.*(3*patstd)-3*patstd./2;
    end
end
%%
%fix the mean and nanstd
EBSP2=fix_mean(EBSP2);
%%
%resize the image
if Settings_Cor.resize == 1
    if isscalar(Settings_Cor.size)
        Settings_Cor.size(2)=Settings_Cor.size(1);
    end

    cs=floor([Settings_Cor.size(1) Settings_Cor.size(2)*size(EBSP2,2)/size(EBSP2,1)]);
    EBSP2 = imresize(EBSP2,cs(1:2));
    cs=size(EBSP2);
else
    cs=size(EBSP2);
end
Settings_Cor.size=cs;
%%
if Settings_Cor.Square == 1
    %square crop on minimum dimension
    ms=min(Settings_Cor.size);
    cv=floor(Settings_Cor.size/2);
    stepv=(1:ms)-floor(ms/2);
    EBSP2=EBSP2(cv(1)+stepv,cv(2)+stepv);
    cs=size(EBSP2);
    Settings_Cor.size=cs;
end
%%

if Settings_Cor.gfilt == 1
    EBSP2=EBSP2-min(EBSP2(:))+1;
    gf=Settings_Cor.gfilt_s*size(EBSP2,1)/100;
    EBSP2B = imgaussfilt(EBSP2,gf);
    EBSP2 = EBSP2./EBSP2B;
    EBSP2=fix_mean(EBSP2);
end

%%
if Settings_Cor.gaussfit == 1
    
    EBSPData.PW=size(EBSP2,2);
    EBSPData.PH=size(EBSP2,1);
    [bg2,Settings_Cor.gaussparams]=bg_fit(EBSP2,EBSPData);
    EBSP2=EBSP2./bg2;
    EBSP2=fix_mean(EBSP2);
end
%%
if Settings_Cor.blur == 1
    ix=size(EBSP2,2);
    Iblur = imgaussfilt(EBSP2, Settings_Cor.blurf(1),'filtersize',Settings_Cor.blurf(2));
    EBSP2=EBSP2-Iblur;
end
if Settings_Cor.gaussBlur ~=0
    % ix=size(EBSP2,2);
    EBSP2 = imgaussfilt(EBSP2, Settings_Cor.gaussBlur);
    % EBSP2=EBSP2-Iblur;
end



if Settings_Cor.radius == 1
    %%
    % EBSP2=fix_mean(EBSP2);
    r_thresh=Settings_Cor.radius_frac*4/3*cs(1)/2;
    
    [xgrid,ygrid]=meshgrid(1:cs(2),1:cs(1));
    xgrid=double(xgrid);
    ygrid=double(ygrid);
    r_grid=sqrt(double(xgrid-size(EBSP2,2)/2).^2+double(ygrid-size(EBSP2,1)/2).^2);
    %%
    EBSP2(r_grid>=r_thresh) = 0 ;
    % EBSP2(r_grid>=r_thresh) = mean(EBSP2(:)) ;
    % 
    %
    % EBSP2(r_grid<r_thresh)=  (EBSP2(r_grid<r_thresh)-mean(EBSP2(r_grid<r_thresh),"omitnan"))/std(EBSP2(r_grid<r_thresh),"omitnan");
    EBSP2(r_grid<r_thresh)=  (EBSP2(r_grid<r_thresh)-min(EBSP2(r_grid<r_thresh)))/max(EBSP2(r_grid<r_thresh));
    % EBSP2(r_grid>=r_thresh) = 0 ;
    % EBSP2(r_grid>=r_thresh) = 0.5 ;
    %%
    % normalize
else
    EBSP2=fix_mean(EBSP2);
end
%%
if Settings_Cor.TKDBlur2 == 1
    gf=Settings_Cor.gfilt_s*size(EBSP2,1)/200;
    EBSP2B = imgaussfilt(EBSP2,gf);
    EBSP2=EBSP2./EBSP2B;
    for n=1:3
        EBSP_med = medfilt2(EBSP2);
        hpix=find(abs(EBSP_med-EBSP2) > 0.8);
        EBSP2(hpix)=EBSP_med(hpix);
    end
%     EBSP2=fix_mean(EBSP2);
end
%%
if Settings_Cor.MeanCentre==1
    %taken from fix fixmean function, without the make positive part
    EBSP2=(EBSP2-mean(EBSP2(:),"omitnan"))/std(EBSP2(:),"omitnan");
end

%%
if Settings_Cor.quadrants==1
    %%
    cross_size=Settings_Cor.cross;
    center=uint16(size(EBSP2)/2);
    mask=zeros(size(EBSP2));
    mask2 = mask;
    mask2(ceil(center(1)/3):center(1)-cross_size+1 , ceil(center(2)/3):center(2)-cross_size+1)=1;
    mask(1:center(1)-cross_size+1 , 1:center(2)-cross_size+1)=1;
    quadDims = [center(1)-cross_size+1,center(2)-cross_size+1];
    EBSP2 = rescale(EBSP2);
    for i=1:4
        % disp(i)
        
       
        % figure
        % imagesc(mask)
        
        if Settings_Cor.gfiltQuad == 1
            gf=Settings_Cor.gfilt_s*size(EBSP2,1)/50;
            filtImage = reshape(EBSP2(~mask==0),quadDims(1),quadDims(2));
            EBSP2B = imgaussfilt(filtImage,gf,'Padding','replicate');
            EBSP2B=reshape(EBSP2B,quadDims(1),quadDims(2));
            EBSP2(~mask==0) = filtImage./EBSP2B;
        end
      EBSP2(~mask==0)=(EBSP2(~mask==0)-mean(EBSP2(~mask2==0),"omitnan"))/std(EBSP2(~mask2==0),"omitnan");
         mask=rot90(mask);
        mask2 = rot90(mask2);
         % figure
    % imagesc(EBSP2)
    % clim([-3 3])
    % axis xy
    % axis equal
    % colormap gray
    end
    EBSP2=(EBSP2-mean(EBSP2(:),"omitnan"))/std(EBSP2(:),"omitnan");
    
    % if Settings_Cor.gfilt == 1
    % gf=Settings_Cor.gfilt_s*size(EBSP2,1)/100;
    % EBSP2B = imgaussfilt(EBSP2,gf);
    % EBSP2 = EBSP2./EBSP2B;
    % end
    % if Settings_Cor.MeanCentre==1
    % %taken from fix fixmean function, without the make positive part
    % EBSP2=(EBSP2-mean(EBSP2(:),"omitnan"))/std(EBSP2(:),"omitnan");

    % end
end




if Settings_Cor.gaussian_range>0

    gaussian_range=Settings_Cor.gaussian_range*max(size(EBSP2))/2;
    mask =  gauss2d(size(EBSP2),gaussian_range);
    mask = rescale(mask);
    % figure
    % imagesc(mask)
    % size(mask)
    % colorbar();
    EBSP2 = mask.*EBSP2;
end

if Settings_Cor.cross>=1
    cross_size=Settings_Cor.cross;
    center=uint16(size(EBSP2)/2);
    % slice=[center(1)-cross_size:center(1)+cross_size,center(2)-cross_size,center(2)+cross_size]
    % EBSP2(slice)=0;
    % EBSP2=(EBSP2-mean(EBSP2(:),"omitnan"))/std(EBSP2(:),"omitnan");
    mask=ones(size(EBSP2));
    mask(center(1)-cross_size+1:center(1)+cross_size,:)=0;
    mask(:,center(2)-cross_size+1:center(2)+cross_size)=0;
    % EBSP2(slice)=0;
    
end

if Settings_Cor.MeanCentre==1
    %taken from fix fixmean function, without the make positive part
    EBSP2=(EBSP2-mean(EBSP2(:),"omitnan"))/std(EBSP2(:),"omitnan");
end

Settings_Cor_out=Settings_Cor;
end


function [EBSP2,num_hot]=cor_hotpix(EBSP,hthresh)

%Image correction - modified from JH code
%median filter

EBSP_med = medfilt2(EBSP);

%hot pixel correct
h_pix=find(abs(EBSP-EBSP_med)>hthresh); %1000 is chosen as an arbitary number
EBSP2=EBSP;
EBSP2(h_pix)=EBSP_med(h_pix);
num_hot=size(h_pix);
%%
end




function EBSP2=fix_mean(EBSP)
%zero mean & fix stdev
% EBSP2=(EBSP-nanmean(EBSP(:)))/nanstd(EBSP(:));
EBSP2=(EBSP-trimmean(EBSP(:),10))/nanstd(EBSP(:));
%make positive
% EBSP2=EBSP2-min(EBSP2(:))+1;
end

function [bg2,params1]=bg_fit(EBSP2,EBSPData)
%build an image grid for bg fitting
ygv=1:1:EBSPData.PH;
xgv=1:1:EBSPData.PW;
[xgr,ygr] = meshgrid(xgv,ygv);

xsize=EBSPData.PW;

[xmax,ixv]=max(EBSP2);
[mv,iy]=max(xmax);
ix=ixv(iy);
range_e=mv-min(EBSP2(:));

params=[      ix            EBSPData.PH/4       iy            EBSPData.PH/4     range_e 1      range_e  1];

bg_dif=@(params,xgr,ygr,xsize,EBSP2)(abs(x_bg( params,xgr,ygr,xsize)./EBSP2-1));
singleval=@(x)(sum(x(:)));

fun_bg_solve=@(params)(singleval(bg_dif( params,xgr,ygr,xsize,EBSP2)));

[params1,fval] = fminsearch(fun_bg_solve,params);

bg2=x_bg( params1,xgr,ygr,xsize);
end

function [ bg ] = x_bg( params,xgr,ygr,xsize)
%X_BG Summary of this function goes here
%   Detailed explanation goes here

bg= exp( - ((((xgr-params(1)).^2)./(2*(params(2).^2))) + (((ygr-params(3)).^2)./(2*(params(4).^2)))));

bg(:,1:xsize/2)=params(5)*bg(:,1:xsize/2)+params(6);
bg(:,xsize/2+1:end)=params(7)*bg(:,xsize/2+1:end)+params(8);

end

  function G = gauss2d(size2,sigma)

  %% 2d "classic" gauss filter
        [X, Y] = meshgrid(1:size2(1), 1:size2(2));
        A=1/sqrt(2*pi*sigma^2);
        G = 1*exp(-1/(sigma^2)*((Y-size2(1)/2).^2 + (X-size2(2)/2).^2));
  end
  
  function G=gauss2_disc(size2,sigma,disc_size)

  %% circular gauss filter
     G = zeros(size2(1),size2(2));
     r_thresh=disc_size*min(size2)/2;
     [xgrid,ygrid]=meshgrid(1:size2(1),1:size2(2));
     xgrid=double(xgrid);
     ygrid=double(ygrid);
     r_grid=sqrt(double(xgrid-size2(1)/2).^2+double(ygrid-size2(2)/2).^2);
     G(r_grid<=r_thresh) = 1;
     G=imgaussfilt(G,sigma);
  end


  function G=gauss2_square(size2,sigma,disc_size)

  %% square gauss filter
     G = zeros(size2(1),size2(2));
     G( uint16(size2(1)*(1-disc_size)/2) : uint16(size2(1)*(1-(1-disc_size)/2)) ,uint16( size2(2)*(1-disc_size)/2) : uint16( size2(2)*(1-(1-disc_size)/2)) ) = 1;
     
%%
     G=imgaussfilt(G,sigma);
%      Y(size(2)*(disc_size/2):size(2)*(1-disc_size/2)=1
  end
