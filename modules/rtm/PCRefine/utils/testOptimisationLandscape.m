function metric = testOptimisationLandscape(PCGrid,rotGrid,screen_int,pat,SettingsCor,metricName,gradients)
    options.metric = metricName;
%     disp(nargin)
    if nargin==7
        options.gradients = gradients;
    end
    if ~isfield(SettingsCor,'useParallelOptiLandscape')
        SettingsCor.useParallelOptiLandscape = 0;
    end
    %%
        % hw = waitbar(0,'Running Test of optimisation landscape....');
    metric = zeros(size(PCGrid,1),size(PCGrid,2));
    PCGrid = reshape(PCGrid,numel(metric),3);
    PCGrid = PCGrid();
    centre = ceil(size(metric,1)/2);

    if ~isfield(SettingsCor, "Debug")
        SettingsCor.Debug = 0;
    end
    numMetric = numel(metric)
    metricSize = size(metric)
    %%
    if SettingsCor.useParallelOptiLandscape == 1
        parfor gridNum = 1:numMetric
    %         gridNum = 8
            [row,col] = ind2sub(metricSize,gridNum);
            
            gmatrix = rotGrid(gridNum).matrix;
            pc = PCGrid(gridNum,:);
            geometry = EBSP_Gnom(SettingsCor,pc);
            
            patSim = EBSP_gen(geometry,gmatrix,screen_int,screen_int.isHex);
            patSim = EBSP_BGCor(patSim,SettingsCor);
            %%
            
            metric(gridNum) = getMatchMetric(patSim,pat,options);
        end
       
    else
        for gridNum = 1:numMetric
    %         gridNum = 8
    %%
            [row,col] = ind2sub(metricSize,gridNum);
            
            gmatrix = rotGrid(gridNum).matrix;
            pc = PCGrid(gridNum,:);
            geometry = EBSP_Gnom(SettingsCor,pc);
            
            patSim = EBSP_gen(geometry,gmatrix,screen_int,screen_int.isHex);
            patSim = EBSP_BGCor(patSim,SettingsCor);
            %%
            
            metric(gridNum) = getMatchMetric(patSim,pat,options);
          
        end
    end