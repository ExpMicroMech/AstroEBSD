function [EBSP,XCF_data_fill] = refine_prep(EBSP_cor,SettingsXCF,RTM_setup)
%REFINE_PREP Prepare a pattern to be refined

screensize=RTM_setup.screensize;
if isfield(RTM_setup,'tkd_onaxis')
    %change the settings for doing TKD analysis
    
    rmin=RTM_setup.tkd_onaxis(4);
    rmax=RTM_setup.tkd_onaxis(5);
    
    EBSP_cor(~RTM_setup.grid_ok)=0;
    
else
    rmin=10;
    rmax=screensize/sqrt(2);
end

if isfield(RTM_setup,'WindowOff')
else
    RTM_setup.WindowOff=0;
end

if RTM_setup.WindowOff == 1
    SettingsXCF.hfilter=SettingsXCF.hfilter*0+1;
end

% LPTsize=RTM_setup.LPTsize;
[EBSP.FFT,XCF_data_fill]  =fROIEx2(EBSP_cor,SettingsXCF);

EBSP.logp = logsample(EBSP_cor, rmin, rmax/sqrt(2), screensize/2, screensize/2, RTM_setup.LPTsize, RTM_setup.LPTsize); %Transform the reference image into LPT space, logsample is in logsample
EBSP.pat_in=EBSP_cor;

if isfield(SettingsXCF,'single')
else
    SettingsXCF.single=0;
end

if SettingsXCF.single == 1
    EBSP.logp=single(EBSP.logp);
    EBSP.FFT=single(EBSP.FFT);
end

end

