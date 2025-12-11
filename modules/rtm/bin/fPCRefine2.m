function [PC_best] = fPCRefine2(PatternIn,rotmat,PatternInfo,PCRefine_Settings,SettingsXCF,screen_int,RTM_setup)
%% fPCRefine2 - McAuliffe PC Refinement using gradient descent
%
% Inputs -
% PC_start = [PCx, PCy, PCz]
%
% Refine = structure
% Refine.ss=0.08; %initial probe volume
% Refine.p=2; %order of polynomial to fit to tetrahedron
% Refine.n_its=10; %number of interations
%
% EBSD_geom = geometry of the detector
%
% update 2025 - TBB - reducing the inputs

%rotation matricies

if isfield(PCRefine_Settings,'PHTol')
else
    PCRefine_Settings.PHTol = 1E-3;
end


phbest=0;

%create the descent kernel
xs=[1,3,2,2,2,2,2];
ys=[2,2,1,3,2,2,2];
zs=[2,2,2,2,1,3,2];

%start the PC values
PC_start=PatternInfo.PC;
PC_best=PC_start;

%initalise a counter
skipgradcalc=0;

n=1;

%count for a number of iterations
while n < PCRefine_Settings.n_its
    phlast=phbest;

    phs=zeros(3,3,3);

    for i=1:7
        x=xs(i);
        y=ys(i);
        z=zs(i);

        x_itval=(x-2).*PCRefine_Settings.ss;
        y_itval=(y-2).*PCRefine_Settings.ss;
        z_itval=(z-2).*PCRefine_Settings.ss;

        PC_refined=[PC_best(1)+x_itval,PC_best(2)+y_itval,PC_best(3)+z_itval];
        [ EBSD_geom ] = EBSP_Gnom( PatternInfo,PC_refined );

        %Optimise the fit, allowing for small orientation updates
        [~,regout] = refine5(PatternIn,EBSD_geom,PC_start,rotmat,SettingsXCF,screen_int,screen_int.isHex,RTM_setup);

        %extract the PH from the refinement
        phs(x,y,z)=regout(4);
    end

    %fit to iterpolants and find gradients
    xcol=reshape(phs(:,2,2),1,3);
    ycol=reshape(phs(2,:,2),1,3);
    zcol=reshape(phs(2,2,:),1,3);
    sv=[-PCRefine_Settings.ss,0,PCRefine_Settings.ss];

    %sum peakheights across the axes and fit a polynomial
    [fitx]=polyfit(sv,xcol,PCRefine_Settings.p); %xaxis
    [fity]=polyfit(sv,ycol,PCRefine_Settings.p); %yaxis
    [fitz]=polyfit(sv,zcol,PCRefine_Settings.p); %zaxis

    %evaluate the gradient at x=0
    grad_x=fitx(PCRefine_Settings.p);
    grad_y=fity(PCRefine_Settings.p);
    grad_z=fitz(PCRefine_Settings.p);

    %descend
    grad=0.1*[grad_x.*PCRefine_Settings.ss,grad_y.*PCRefine_Settings.ss,grad_z.*PCRefine_Settings.ss];

    PCRefine_Settings.PC_it(n,:)=PC_start;
    PCRefine_Settings.PH_it(n)=phs(2,2,2);

    phtrial=phs(2,2,2);

    % decide whether to continue
    if n<PCRefine_Settings.n_its && phtrial>phbest

        %update the current best PC and peak height
        %move the trial PC (for the next iteration)
        PC_best=PC_start;
        phbest=phtrial;



        if skipgradcalc==0
            %set the next PC to follow the gradient
            PC_start=[PC_best(1)+grad(1),PC_best(2)+grad(2),PC_best(3)+grad(3)];
        else
        end

        skipgradcalc=0;
        counter=0;

    elseif n<PCRefine_Settings.n_its %else reduce the step size
        PCRefine_Settings.ss=0.9*PCRefine_Settings.ss;
        %and move slightly randomly
        reduction=(5-mod(counter,5))/5;

        grad=reduction.*0.01.*[(rand-0.5),(rand-0.5),(rand-0.5)].*PC_start; % -0.25 to 0.25 percent of PC_start
        PC_best=[PC_best(1)+grad(1),PC_best(2)+grad(2),PC_best(3)+grad(3)];
        skipgradcalc=1;

        counter=counter+1;



    end

    if PCRefine_Settings.print == 1
        fprintf('PH = %0.5f, I = %3.0f, PC = [ %0.5f %0.5f %0.5f ] \n',phbest,n,PC_best(1),PC_best(2),PC_best(3));
    end

    %increment the loop
    p=n;
    n=n+1;
    
    %break the loop if the 
    if phbest-phlast < PCRefine_Settings.PHTol %break the search if the PH is not getting better
        n=PCRefine_Settings.n_its;
        
    end

end
fprintf('PH = %0.5f, I = %3.0f, PC = [ %0.5f %0.5f %0.5f ] END \n',phbest,p,PC_best(1),PC_best(2),PC_best(3));

end


