function [varargs,varargNum,optimisationOptions]= checkRealSpaceOptions(optimisationOptions)


if ~isfield(optimisationOptions,'fminsearch')
    optimisationOptions.fminsearch=0;
end

if ~isfield(optimisationOptions,'refinePC')
    optimisationOptions.refinePC=1;
    % this is actually supposed to be one as default
end
if ~isfield(optimisationOptions,'genetic')
    optimisationOptions.genetic=0;
end
if ~isfield(optimisationOptions,'particleSwarm')
    optimisationOptions.particleSwarm=0;
end
if ~isfield(optimisationOptions,'patternsearch')
    optimisationOptions.patternsearch=0;
end
if ~isfield(optimisationOptions,'snobfit')
    optimisationOptions.snobfit=0;
end
if ~isfield(optimisationOptions,'gradients')
    optimisationOptions.gradients = 0;
end
if ~isfield(optimisationOptions,['doPlot' ...
        '' ...
        ''])
    optimisationOptions.doPlot = 0;
end
if ~isfield(optimisationOptions, 'parFor')
    optimisationOptions.parFor = 1; % parallel default;
    optimisationOptions.UseParallel = 0; %do no tuse parallel in inner loops.               
end
% optimisationOptions.refinePC ==1



% set up variable argumetns
varargNum = 1;
% varargs = {""}
if optimisationOptions.fminsearch ==1
    varargs{varargNum} = 'Mode';
    varargs{varargNum+1} = 'fminsearch';
    % varargNum = varargNum+
end
% this way the first option will be overwritten, by the later one.
if optimisationOptions.genetic ==1
    varargs{varargNum} = 'Mode';
    varargs{varargNum+1} = 'genetic';
end
if optimisationOptions.snobfit ==1
    varargs{varargNum} = 'Mode';
    varargs{varargNum+1} = 'snobfit';
end
if optimisationOptions.particleSwarm ==1
    varargs{varargNum} = 'Mode';
    varargs{varargNum+1} = 'particleSwarm';
end
if optimisationOptions.patternsearch ==1
    varargs{varargNum} = 'Mode';
    varargs{varargNum+1} = 'patternsearch';
end
varargNum = varargNum+2;


if optimisationOptions.gradients ==1
    varargs{varargNum} = 'gradients';
    varargNum = varargNum+1;
end
if optimisationOptions.refinePC ==1
    varargs{varargNum} = 'refinePC';
    varargNum = varargNum+1;
end
if optimisationOptions.doPlot ==1
    varargs{varargNum} ='doPlot'
    varargNum = varargNum+1;
end
if ~exist("varargs")
    varargs = "";
end

disp(varargs);

end