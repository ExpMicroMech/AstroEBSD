function metricOUT = getMatchMetric(patSimIn,patIn,metricIn)
    gradients = 0;
    minVal = -1000;
    if strcmp(class(metricIn),'struct')
        try 
            metric = metricIn.metric;
        catch
            metric = "NCC";
        end
        if isfield(metricIn,'gradients')
            gradients = metricIn.gradients;
        end
        if isfield(metricIn,'minVal')
            minVal = metricIn.minVal;
        end
       
    else
        metric = metricIn;
    end
    % if nargin ==3
    if gradients == 1

       [GxSim,GySim] = imgradientxy((patSimIn));
       [GxPat,GyPat] = imgradientxy((patIn));
       metricX =  patternMetric(GxSim,GxPat,metric);
       metricY =  patternMetric(GySim,GyPat,metric);
       metricOUT = (metricX+metricY)/2;
    %%
    else
        metricOUT= patternMetric(patSimIn,patIn,metric);
    end
    metricOUT = double(metricOUT); % solvers might bother, but to save memory we sometimes use single precisoin interally.
    if metricOUT<minVal
        metricOUT = NaN;
    end
end
function metricOUT= patternMetric(patSim,pat,metric)
    if strcmp(metric,'NCC')
        metricOUT = nxCorrelation(patSim,pat);
%         metricOUT = max(ncc(:));
    % end
    elseif strcmp(metric,'MI')
        metricOUT = 1-mutualInfo(patSim,pat);
    % end
    elseif strcmp(metric,'SSIM')
        metricOUT = ssim(patSim,pat);
    elseif strcmp(metric,'ndp')
        metricOUT = normdotprod(patSim,pat);
    elseif strcmp(metric,'RMSE')
        % elseif strcmp(metric,'ndp')
        metricOUT = rmse(patSim(:),pat(:));
    else
        metricOut = double(-1);
    end
end