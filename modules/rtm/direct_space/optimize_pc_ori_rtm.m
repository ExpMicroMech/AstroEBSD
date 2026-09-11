function [pc_out,eangs,peakhight,lattice_strain,barrelDistortion,optiparams] =optimize_pc_ori_rtm(pc_in,eangs,PatternIn,Settings_Cor,screen_int,options,varargin)
% OPTIMIZE_PC_ORI
%
%Optimizes a orientation and patterncenter based on the normalised
%crosscorrelation normxcorr2
%[pc_out,eangs,peakhight]=
%OPTIMIZE_PC_ORI(pc_in,eangs,PatternIn,Settings_Cor,screen_int,varargin)
% - flag 'Mode', with either 'fmincon' or 'fminsearch' to select for
% optimisation method
%- Settings_Cor: define either max_var_pc or max_var_pc_x,max_var_pc_y and
%max_var_pc_z as field names to set the constraints for the refienement.
%- pc_in: 3x1 double, containing pcx, pcy, DD
%- eangs: 3x1 double, contianing Euleranglex in ZXZ convention in radiants
%- screen_int, contains interpolants of the masterpattern
%- Settings_Cor, struct, with cropping and correction infos, as well as
%refinement constraints.
%
%
%Outputs
%pc_out: refined pattern_center
%eangs: refined orientation
%peakhigh: normalised cross product between simulation and experiment after
%refinement
%
%

% Settings_Cor.size=256;
% Settings_Cor.SquareCrop=0;
% eangs=extractions.(tile_name).oris(p_cur,:);
% eangs=conv_G_to_EA(gmatrix_pseudo);

% TODO: implement.
debugWarnings = 0;

if debugWarnings==0
    warning('off','all')
end


mode='fmincon';
refine_lattice_parameter=0;
refine_patternCentre=0;
refine_barrel = 0;
backgroundVariation = 0;
gradients= 0;
doPlot = 0;
doCorrection= 0;
% pcScaling = 100;
% eulerScaling = 15;




for i=1:numel(varargin(:))
    if(strcmp(varargin{i},'Mode'))
        mode=varargin{i+1};
    end
    if(strcmp(varargin{i},'gradients'))
        gradients = 1;
    end
    if strcmp(varargin{i},'LPRatio')
        refine_lattice_parameter=1;

    end
    if strcmp(varargin{i},'refinePC')
        refine_patternCentre = 1;
    end
    if strcmp(varargin{i},'refineBarrel')
       refine_barrel = 1;
    end
    if strcmp(varargin{i},'BackgroundVariation')
        backgroundVariation = 1;
    end
    if strcmp(varargin{i},'doCorrection')
        doCorrection = 1;
    end
end
if ~isfield(options,'metric')
    metric = struct();
    metric.metric = 'NCC';
    options.metric = metric;
end
if strcmp(class(options.metric),'string')
    metric = struct();
    metric.metric = options.metric;
    options.metric = metric;
end
if ~isfield(options,'minFunVal')
    options.metric.minVal=-1;
else
    options.metric.minVal = options.minFunVal;
end
options.gradients = gradients;
%%





%% scale variables
try
    scalingFactor = options.scalingFactor;
catch
    scalingFactor = 1;
end
try
    pcScaling = options.pcScaling;
catch
    options.pcScaling = scalingFactor*250;
end
try
    eulerScaling = options.oriScaling;
catch
    options.oriScaling = scalingFactor*250;
end
% TODO: lattice strain.




pc_in = pc_in* options.pcScaling;
eangs = eangs* options.oriScaling;


lattice_strain = [1 1 1];
pcx_in=double(pc_in(1));
pcy_in=double(pc_in(2));
pcz_in=double(pc_in(3));

try
    barrel = Settings_Cor.BarrelDistortion;
catch
    barrel = 0;
end

fields={'max_var_pc_x','max_var_pc_y','max_var_pc_y'};

if any(isfield(Settings_Cor,fields))
    try
    var_pcx=Settings_Cor.max_var_pc_x*options.pcScaling;
    var_pcy=Settings_Cor.max_var_pc_y*options.pcScaling;
    var_pcz=Settings_Cor.max_var_pc_z*options.pcScaling;
    catch
        error('you have spacefied a constrained for one pc variblae but not all of them')
    end
else
    var_pcx=Settings_Cor.max_var_pc*options.pcScaling;
    var_pcy=Settings_Cor.max_var_pc*options.pcScaling;
    var_pcz=Settings_Cor.max_var_pc*options.pcScaling;
end

phi1_in=double(eangs(1));
PHI_in=double(eangs(2));
phi2_in=double(eangs(3));
var_phi1=Settings_Cor.max_var_ori*options.oriScaling;
var_PHI=Settings_Cor.max_var_ori*options.oriScaling;
var_phi2=Settings_Cor.max_var_ori*options.oriScaling;
barrelTemp = 0;

if isfield(Settings_Cor,'BarrelDistortion')
    barrelTemp= Settings_Cor.BarrelDistortion;
    Settings_Cor.BarrelDistortion = 0;
end
%%
if doCorrection == 1
    [EBSP_single.ebsp_cor,EBSP_single.PatternInfo]=EBSP_BGCor(PatternIn,Settings_Cor );
else
    [~,EBSP_single.PatternInfo]=EBSP_BGCor(PatternIn,Settings_Cor );
    EBSP_single.ebsp_cor = PatternIn;
end
%%
% if backgroundVariation ==1
%     initialSimulationGeometry = EBSP_Gnom(EBSP_single.PatternInfo,pc_in); %you can change PC_in if you want
%     Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
%     Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
%     Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)];
%     gmatrix=Rz(phi2_in)*Rx(PHI_in)*Rz(phi1_in);
%     output_pattern= EBSP_gen(initialSimulationGeometry,gmatrix,screen_int);
%     EBSP_single.ebsp_cor = estimateBackgroundDifference(output_pattern,EBSP_single.ebsp_cor,Settings_Cor)
%     % EBSP_ref =EBSP_BGCor(EBSP_ref, Settings_Cor);
%     % diffPat =  EBSP_single.ebsp_cor - output_pattern;
%     % diffPatFiltered = imgaussfilt(diffPat,20);
%     % ebsp_cor = EBSP_single.ebsp_cor -diffPatFiltered;
%     % [ebsp_cor,~]=EBSP_BGCor(ebsp_cor,Settings_Cor );
%     % imagesc(diffPat-diffPatFiltered);
% end

%%
Settings_Cor.BarrelDistortion = barrelTemp;
EBSP_single.PatternInfo.size=size(EBSP_single.ebsp_cor);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                    % Function calls for solvres depending on
                    % desired refinement
                    % TODO: this needs to be streamlined.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if refine_patternCentre==1 &  refine_lattice_parameter==1

    x0=[pcx_in,pcy_in,pcz_in,phi1_in,PHI_in,phi2_in,lattice_strain(1),lattice_strain(2),lattice_strain(3)];
    f=@(x)(optimization_with_lattice_and_patternCentre(x,screen_int,Settings_Cor,EBSP_single.PatternInfo,EBSP_single.ebsp_cor,options));


elseif refine_lattice_parameter==1

    x0=[phi1_in,PHI_in,phi2_in,lattice_strain(1),lattice_strain(2),lattice_strain(3)];
    f=@(x)(optimization_with_lattice(x,screen_int,Settings_Cor,EBSP_single.PatternInfo,EBSP_single.ebsp_cor,pcx_in,pcy_in,pcz_in,options));


elseif refine_patternCentre==1

    x0=[pcx_in,pcy_in,pcz_in,phi1_in,PHI_in,phi2_in];
    f=@(x)(optimization_with_patternCentre(x,screen_int,Settings_Cor,EBSP_single.PatternInfo,EBSP_single.ebsp_cor,options));
    l=[pcx_in-var_pcx;pcy_in-var_pcy;pcz_in-var_pcz;phi1_in-var_phi1;PHI_in-var_PHI;phi2_in-var_phi2];
    u=[pcx_in+var_pcx;pcy_in+var_pcy;pcz_in+var_pcz;phi1_in+var_phi1;PHI_in+var_PHI;phi2_in+var_phi2];



elseif refine_barrel == 1 
    
    x0 = [barrel];
    f = @(x)(optimization_Barrel(x,screen_int,Settings_Cor,EBSP_single.PatternInfo,EBSP_single.ebsp_cor,pcx_in,pcy_in,pcz_in,phi1_in,PHI_in,phi2_in,options));
    l = [-50];
    u = [50];


else
    
    x0=[phi1_in,PHI_in,phi2_in];
    f=@(x)(optimization(x,screen_int,Settings_Cor,EBSP_single.PatternInfo,EBSP_single.ebsp_cor,pcx_in,pcy_in,pcz_in,options));
    l=[phi1_in-var_phi1;PHI_in-var_PHI;phi2_in-var_phi2];
    u=[phi1_in+var_phi1;PHI_in+var_PHI;phi2_in+var_phi2];
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                    % Options for different algorithms


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






opts = optimoptions('fmincon');
try
    opts.FunValCheck = options.FunValCheck;
catch
    warning('No objective Limit specified')
end
% try
    % opts.FunVa

if strcmp(mode,'fminsearch')
    opts = optimset('fminsearch');
end


if strcmp(mode,'genetic')
    % mutationRate = 0.7
    % if isfield(options,'mutationRate')
    %     mutationRate = options.MutationRate
    % end
    % crossOverFraction = 0.7
    if doPlot == 1
        opts = optimoptions('ga','PlotFcn', @gaplotbestf);
    else
        opts = optimoptions('ga')
    end
end

if strcmp(mode,'particleswarm')
    % mutationRate = 0.7
    % if isfield(options,'mutationRate')
    %     mutationRate = options.MutationRate
    % end
    % crossOverFraction = 0.7
    if doPlot == 1
        opts = optimoptions('particleswarm','PlotFcn', @pswplotbestf);
    else
        opts = optimoptions('particleswarm')
    end
end

if strcmp(mode,'patternsearch')
    % mutationRate = 0.7
    % if isfield(options,'mutationRate')
    %     mutationRate = options.MutationRate
    % end
    % crossOverFraction = 0.7
    if doPlot == 1
        opts = optimoptions('patternsearch','PlotFcn', @psplotbestf);
    else
        opts = optimoptions('patternsearch');
    end
end


%%%%%% Particle swarm
try
    opts.SelfAdjustmentWeight = options.SelfAdjustmentWeight;
catch
    try
    opts.SelfAdjustmentWeight  = 1.49;
    catch
        warning('SelfAdjustmentWeight not implemented for %s  algorithm',mode)
    end
end

try
    opts.SocialAdjustmentWeight = options.SocialAdjustmentWeight;
catch
    try
    opts.SocialAdjustmentWeight  = 1.49;
    catch
        warning('SocialAdjustmentWeight not implemented for  %s  algorithm',mode)
    end
end





%%% GENETIC
try
    opts.CrossoverFraction = options.CrossoverFraction;
catch
    try
    opts.CrossoverFraction  = 0.7;
    catch
        warning('CrossoverFraction not implemented for %s  algorithm',mode)
end

try
    opts.MaxGenerations = options.MaxGenerations;
catch
    try
    opts.MaxGenerations  = 50;
    catch
        warning('MaxGenerations not implemented for  %s  algorithm',mode)
    end
end


%%% FMINCON
try
    opts.StepTolerance = options.StepTolerance;
catch
    try
    opts.StepTolerance = 1e-6;
    catch
        warning('StepTolerance not implemented for  %s  algorithm',mode)
    end
end

try
    opts.OptimalityTolerance = options.OptimalityTolerance;
catch
    try
    opts.OptimalityTolerance = 1e-6;
    catch
        warning('not implemented for  %s  algorithm',mode)
    end
end


%%% Fminsearch
try
    opts.TolFun = options.TolFun;
catch
    try
    opts.StepTolerance = 1e-6;
    catch
        warning('Tol Fun not implemented for %s  algorithm',mode)
    end
end


try
    opts.TolX = options.TolX;
catch
    try
    opts.TolX = 1e-6;
    catch
        warning('TolX not implemented for %s  algorithm',mode)
    end
end



%%% ALL Global
try
    opts.MaxIter = options.MaxIter;
catch
    try
    opts.MaxIter = 300;
    catch
        warning(' MaxIter not implemented for genetic %s  algorithm',mode)
    end
end
try
    opts.MaxIterations = options.MaxIter;
catch
    try
    opts.MaxIterations = 300;
    catch
        warning(' MaxIterations not implemented for genetic %s  algorithm',mode)
    end
end

%%% all parallelizable
try
    opts.UseParallel = options.UseParallel;
catch
    try
        opts.UseParallel = 0;
          catch
            warning(' UseParallel not implemented for %s  algorithm',mode)
    end
    end
end

%%%%% display options
% opts.Display = 'final-detailed';
try
    opts.Display = options.Display;
catch
    opts.Display = 'final';
end

optiparams = struct();











%%%% start the optimisers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                    % Call the different solvers


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if strcmp(mode,'fmincon')
    try
        [x,fval,exitflag,output] = fmincon(f,x0,[],[],[],[],l,u,[],opts);
    catch
        warning('fmincon did not run, likely your function value is too bad')
        x = x0;
        exitflag = '';
        output = '';
        fval = 10; %peakhight = 1-fVal, so fVal neads to be a high positive number % lets make it not to high to not cut off numbers from the other values.
    end
    optiparams.exitflag = exitflag;
    optiparams.output = output;
elseif strcmp(mode,'fminsearch')
    try
        [x,fval,optiparams.exitflag,optiparams.output] = fminsearch(f,x0,opts);%,l,u)%,[],options)
    catch
        warning('fmincon did not run, likely your function value is too bad')
        x = x0;
        exitflag = '';
        output = '';
        fval = 10;
    end
elseif strcmp(mode,'fminunc')
    [x,fval] = fminunc(f,x0);%,l,u)%,[],options)

elseif strcmp(mode,'genetic')
    [x,fval,exitflag,output,population,scores] = ga(f,numel(x0),[],[],[],[],l,u,[],[],opts);

    optiparams.exitflag = exitflag;
    optiparams.output = output;
    optiparams.population = population;
    optiparams.scores = scores;

elseif strcmp(mode,'particleSwarm')
    [x,fval,exitflag,output,points] = particleswarm(f,numel(x0),l,u,opts);
    optiparams.exitflag = exitflag;
    optiparams.output = output;
    % optiparams.population = population;
    optiparams.points = points;

elseif strcmp(mode, 'patternsearch')
    [x,fval,exitflag,output] = patternsearch(f,x0,[],[],[],[],l,u,[],opts);
    optiparams.exitflag = exitflag;
    optiparams.output = output;

elseif strcmp(mode,'snobfit')
    [x,fval,exitflag,output,population,scores] = snobfit(f,numel(x0),[],[],[],[],l,u,[],[],opts);

else
    error('you have specified a none supported mode')

end

%%
if refine_patternCentre==1 
    pc_out=[x(1)/options.pcScaling,x(2)/options.pcScaling,x(3)/options.pcScaling];
    eangs=[x(4)/options.oriScaling,x(5)/options.oriScaling,x(6)/options.oriScaling];

elseif ~refine_barrel == 1
    pc_out = pc_in/options.pcScaling;
    eangs=[x(1)/options.oriScaling,x(2)/options.oriScaling,x(3)/options.oriScaling];

end


if refine_lattice_parameter==1

    if refine_patternCentre ==1
        lattice_strain=[x(7), x(8), x(9)];

    else
        lattice_strain=[x(4), x(5), x(6)];
    
    end
end


%% needs to be updated, if we refine this alongside not only pc and ori. 
if refine_barrel == 1 && ~refine_patternCentre==1 
    barrelDistortion = x(1);
    pc_out=pc_in;
    eangs=eangs;
else
    barrelDistortion = Settings_Cor.BarrelDistortion;
end

peakhight=1-fval;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    % Function Definitions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







function outMetric=optimization_with_patternCentre(x,screen_int,Settings_Cor,PatternInfo,Pattern_in,options)
    % opTIMIZATION
    % internal function, calculates the normalised cross correlation
    % between two input images
    % 
    %%
    pcx=x(1)/options.pcScaling;
    pcy=x(2)/options.pcScaling;
    pcz=x(3)/options.pcScaling;
    phi1=x(4)/options.oriScaling;
    PHI=x(5)/options.oriScaling;
    phi2=x(6)/options.oriScaling;
        % disp(pcScaling);
        % disp([pcx pcy pcz])
        % disp([phi1 PHI phi2])
    % disp(eulerScaling);
    % disp(x);
    compare_pat=generate_pattern(screen_int,Settings_Cor,PatternInfo,pcx,pcy,pcz,phi1,PHI,phi2);
    % imwrite(compare_pat,'before.png')
    compare_pat=EBSP_BGCor(compare_pat,Settings_Cor);
    % imwrite(compare_pat,'after.png')
     Pattern_in_BGsub = Pattern_in;
    
    if isfield(options, 'BGVar')
        if options.BGVar == 1
        tempBarrel = Settings_Cor.BarrelDistortion;
        % tempBarrel = 0;
        Settings_Cor.BarrelDistortion = 0;
        % gFilt = str2num(get(h_BackgroundVariationValue,'String'))
        Pattern_in_BGsub = estimateBackgroundDifference(compare_pat,Pattern_in,Settings_Cor,options.BGVarVal);
        Settings_Cor.BarrelDistortion = tempBarrel;
        % Pattern_Sim_diff = Pattern_EXP - Pattern_Sim;
        end
    end
   
    outMetric = 1- getMatchMetric(Pattern_in_BGsub,compare_pat,options.metric);
% error('end')
%%
    end %


function outMetric = optimization_with_patternCentreAndBarrelDistortion(x,screen_int,Settings_Cor,PatternInfo,Pattern_in,options)
    % opTIMIZATION
    % internal function, calculates the normalised cross correlation
    % between two input images
    % 
    %%
    pcx=x(1)/options.pcScaling;
    pcy=x(2)/options.pcScaling;
    pcz=x(3)/options.pcScaling;
    phi1=x(4)/options.oriScaling;
    PHI=x(5)/options.oriScaling;
    phi2=x(6)/options.oriScaling;

    barrelCorr = x(7);
    Settings_Cor.BarrelDistortion = barrelCorr;
    compare_pat=generate_pattern(screen_int,Settings_Cor,PatternInfo,pcx,pcy,pcz,phi1,PHI,phi2);
    compare_pat=EBSP_BGCor(compare_pat,Settings_Cor);
    %%
    Pattern_in_BGsub = Pattern_in;
    if isfield(options, 'BGVar')
        if options.BGVar ==1
        tempBarrel = Settings_Cor.BarrelDistortion;
        % tempBarrel = 0;
        Settings_Cor.BarrelDistortion = 0;
        % gFilt = str2num(get(h_BackgroundVariationValue,'String'))
        Pattern_in_BGsub = estimateBackgroundDifference(compare_pat,Pattern_in,Settings_Cor,options.BGVarVal);
        Settings_Cor.BarrelDistortion = tempBarrel;
        end
        % Pattern_Sim_diff = Pattern_EXP - Pattern_Sim;
    end
    % normcross=normxcorr2(Pattern_in_BGsub,compare_pat);
    outMetric = 1-getMatchMetric(Pattern_in_BGsub,compare_pat,options.metric);
end



function outMetric=optimization_with_lattice_and_patternCentre(x,screen_int,Settings_Cor,PatternInfo,Pattern_in,options)
    % opTIMIZATION
    % internal function, calculates the normalised cross correlation
    % between two input images
    % 
    %
    pcx=x(1)/options.pcScaling;
    pcy=x(2)/options.pcScaling;
    pcz=x(3)/options.pcScaling;
    phi1=x(4)/options.oriScaling;
    PHI=x(5)/options.oriScaling;
    phi2=x(6)/options.oriScaling;
    l_x=x(7);
    l_y=x(8);
    l_z=x(9);
    compare_pat=generate_pattern_lattice(screen_int,Settings_Cor,PatternInfo,pcx,pcy,pcz,phi1,PHI,phi2);
    compare_pat=EBSP_BGCor(compare_pat,Settings_Cor);
    %%
    
    % normcross=normxcorr2(Pattern_in,compare_pat);
    outMetric = 1-getMatchMetric(Pattern_in_BGsub,compare_pat,options.metric);
end



function outMetric=optimization(x,screen_int,Settings_Cor,PatternInfo,Pattern_in,pcx,pcy,pcz,options)
    % opTIMIZATION
    % internal function, calculates the normalised cross correlation
    % between two input images
    % 
    %%
    % pcx=x(1);
    % pcy=x(2);
    % pcz=x(3);
    phi1=x(1)/options.oriScaling;
    PHI=x(2)/options.oriScaling;
    phi2=x(3)/options.oriScaling;
    pcx = pcx/options.oriScaling;
    pcy = pcy/options.oriScaling;
    pcz = pcz/options.oriScaling;
    compare_pat=generate_pattern(screen_int,Settings_Cor,PatternInfo,pcx,pcy,pcz,phi1,PHI,phi2);
    % imwrite(compare_pat,'before.png')
    compare_pat=EBSP_BGCor(compare_pat,Settings_Cor);
    % imwrite(compare_pat,'after.png')
     Pattern_in_BGsub = Pattern_in;
    if isfield(options, 'BGVar')
        if options.BGVar ==1
        tempBarrel = Settings_Cor.BarrelDistortion;
        % tempBarrel = 0;
        Settings_Cor.BarrelDistortion = 0;
        % gFilt = str2num(get(h_BackgroundVariationValue,'String'))
        Pattern_in_BGsub = estimateBackgroundDifference(compare_pat,Pattern_in,Settings_Cor,options.BGVarVal);
        Settings_Cor.BarrelDistortion = tempBarrel;
        end
        % Pattern_Sim_diff = Pattern_EXP - Pattern_Sim;
    end
    % barrelCorr = x(1);
    % Settings_Cor.BarrelDistortion = barrelCorr;
    
    % compare_pat=EBSP_BGCor(Pattern_in_BGsub,Settings_Cor);
    %%
    
    % normcross=normxcorr2(Pattern_in,compare_pat);
    outMetric = 1-getMatchMetric(Pattern_in_BGsub,compare_pat,options.metric);
end



function outMetric=optimization_Barrel(x,screen_int,Settings_Cor,PatternInfo,Pattern_in,pcx,pcy,pcz,phi1,PHI,phi2,options)
    % opTIMIZATION
    % internal function, calculates the normalised cross correlation
    % between two input images
    
    barrelCorr = x(1);
    Settings_Cor.BarrelDistortion = barrelCorr;
    compare_pat=generate_pattern(screen_int,Settings_Cor,PatternInfo,pcx,pcy,pcz,phi1,PHI,phi2);
    compare_pat=EBSP_BGCor(compare_pat,Settings_Cor);
    %%
    %%
    Pattern_in_BGsub = Pattern_in;
    if isfield(options, 'BGVar')
        if options.BGVar ==1
        tempBarrel = Settings_Cor.BarrelDistortion;
        % tempBarrel = 0;
        Settings_Cor.BarrelDistortion = 0;
        % gFilt = str2num(get(h_BackgroundVariationValue,'String'))
        Pattern_in_BGsub = estimateBackgroundDifference(compare_pat,Pattern_in,Settings_Cor,options.BGVarVal);
        Settings_Cor.BarrelDistortion = tempBarrel;
        end
        % Pattern_Sim_diff = Pattern_EXP - Pattern_Sim;
    end
    % normcross=normxcorr2(Pattern_in_BGsub,compare_pat);
    
    % normcross=normxcorr2(Pattern_in,compare_pat);
    outMetric = 1-getMatchMetric(Pattern_in_BGsub,compare_pat,options.metric);
  
end




function outMetric=optimization_with_lattice(x,screen_int,Settings_Cor,PatternInfo,Pattern_in,pcx,pcy,pcz,options)
    % opTIMIZATION
    % internal function, calculates the normalised cross correlation
    % between two input images
    % 
    %
    % pcx=x(1);
    % pcy=x(2);
    % pcz=x(3);
    phi1=x(1)/options.oriScaling;
    PHI=x(2)/options.oriScaling;
    phi2=x(3)/options.oriScaling;
    l_x=x(4);
    l_y=x(5);
    l_z=x(6);
    compare_pat=generate_pattern_lattice(screen_int,Settings_Cor,PatternInfo,pcx,pcy,pcz,phi1,PHI,phi2);
    compare_pat=EBSP_BGCor(compare_pat,Settings_Cor);
    %%
    
    outMetric = 1-getMatchMetric(Pattern_in_BGsub,compare_pat,options.metric);
   
end

%%


function output_pattern=generate_pattern(screen_int,Settings_Cor,PatternInfo,pcx,pcy,pcz,phi1,PHI,phi2)
    % GENERATE_PATTERN
    % calculates the gnomonic projection based on dynamical simulation and
    % pattern geometry information
    %
    %

    Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
    Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
    Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)];
    gmatrix=Rz(phi2)*Rx(PHI)*Rz(phi1);
    pc=[pcx,pcy,pcz];
    % output_pattern=
    %%
    geometry=EBSP_Gnom(PatternInfo,pc); %you can change PC_in if you want
    output_pattern= EBSP_gen(geometry,gmatrix,screen_int);
end



function output_pattern=generate_pattern_lattice(screen_int,Settings_Cor,PatternInfo,pcx,pcy,pcz,phi1,PHI,phi2,a,b,c)
    % GENERATE_PATTERN
    % calculates the gnomonic projection based on dynamical simulation and
    % pattern geometry information
    %
    %

    Rz=@(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1]; %z rotation
    Rx=@(theta)[1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)]; %x rotation
    Ry=@(theta)[cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)];
    gmatrix=Rz(phi2)*Rx(PHI)*Rz(phi1);
    pc=[pcx,pcy,pcz];
    % output_pattern=
    %%
    geometry=EBSP_Gnom(PatternInfo,pc); %you can change PC_in if you want
    output_pattern= EBSP_gen(geometry,gmatrix,screen_int,1,'lattice_scaling',[a b c]);
end

    % function normcross_quadrant_wise()

end
% function stop = stopFunction(x,optimValues,state) 
% stop = false;
% if x(1) > 11 
% stop = true; 
% disp('Stopping, x > 11') 
% end
% end