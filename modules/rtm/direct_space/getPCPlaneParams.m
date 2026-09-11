function [PC_x_mod,PC_y_mod,PC_z_mod] = getPCPlaneParams(InputUser,PatternData,std_devsForCleaning)

    x_unfiltered = h5read(InputUser.AnalysisPCHDF5File,'/UserInfo/XGrid');
    y_unfiltered = h5read(InputUser.AnalysisPCHDF5File,'/UserInfo/YGrid');
    if isfield(PatternData,'MapData')
        ind = sub2ind([max(PatternData.MapData.XBeam),max(PatternData.MapData.YBeam)],x_unfiltered,y_unfiltered);
        x_unfiltered = (PatternData.MapData.XSampleCor(ind))
        y_unfiltered = (PatternData.MapData.YSampleCor(ind))
    end
    % RefinementNumber = 9
    % PatternData.DataName
    NCC = PatternData.NCC;
    % %------------ REmove NCC outliers -----------
    % NCC = h5read(InputUser.AnalysisHDF5File,'/UserInfo/YGrid');
    outliers = find(NCC > mean(NCC - std_devsForCleaning*std(NCC)));% ..
    PCX = PatternData.PC_out_(1,:);
    PCY = PatternData.PC_out_(2,:);
    PCZ = PatternData.PC_out_(3,:);
    
    PCX = PCX(:,outliers);
    PCY = PCY(:,outliers);
    PCZ = PCZ(:,outliers);
    x = x_unfiltered(outliers);
    y = y_unfiltered(outliers);
    
    filler = zeros(size(PCX'));
    RHS =  [x(:),y(:),filler+1];
    [PC_x_mod,~]=robustfit(RHS,PCX','bisquare',4.685,0);
    [PC_y_mod,~]=robustfit(RHS,PCY','bisquare',4.685,0);
    [PC_z_mod,~]=robustfit(RHS,PCZ','bisquare',4.685,0);
    
    
    PCXGrid_unfiltered = x_unfiltered*PC_x_mod(1)+y_unfiltered*PC_x_mod(2)+PC_x_mod(3);
    PCYGrid_unfiltered = x_unfiltered*PC_y_mod(1)+y_unfiltered*PC_y_mod(2)+PC_y_mod(3);
    PCZGrid_unfiltered = x_unfiltered*PC_z_mod(1)+y_unfiltered*PC_z_mod(2)+PC_z_mod(3);
     
    PCXGrid = PCXGrid_unfiltered(outliers);
    PCYGrid = PCYGrid_unfiltered(outliers);
    PCZGrid = PCZGrid_unfiltered(outliers);
    size(x_unfiltered)
    dim1 = size(x,1); % to get not co parallel vectors
    dim2 = size(x,2); %
    x0 = [PCXGrid(1), PCYGrid(1),PCZGrid(1)];
    x1 = [PCXGrid(2), PCYGrid(2), PCZGrid(2)];
    x2 = [PCXGrid(numel(x)), PCYGrid(numel(x)), PCZGrid(numel(x))];
    
    n = cross([x0-x1],[x0-x2]);
    deltaPC = dot(([PCX',PCY',PCZ']-x1)',repmat(n',size(PCX)))./norm(n);
    deltaPCFromPoints = vecnorm([PCX'-PCXGrid(:),PCY'-PCYGrid(:),PCZ'-PCZGrid(:)],2,2);
    %% ------------- Remove outliers -------------------- 
    
    outliers2 = deltaPCFromPoints < mean(deltaPCFromPoints - std_devsForCleaning*std(deltaPCFromPoints)) ...
                | deltaPCFromPoints > mean(deltaPCFromPoints + std_devsForCleaning*std(deltaPCFromPoints));
    
    PCX = PatternData.PC_out_(1,:);
    PCY = PatternData.PC_out_(2,:);
    PCZ = PatternData.PC_out_(3,:);
    x = h5read(InputUser.AnalysisPCHDF5File,'/UserInfo/XGrid');
    y = h5read(InputUser.AnalysisPCHDF5File,'/UserInfo/YGrid');
    if isfield(PatternData,'MapData')
        ind = sub2ind([max(PatternData.MapData.XBeam),max(PatternData.MapData.YBeam)],x,y);
        x = (PatternData.MapData.XSampleCor(ind))
        y = (PatternData.MapData.YSampleCor(ind))
    end
    PCX = PCX(:,outliers(~outliers2));
    PCY = PCY(:,outliers(~outliers2));
    PCZ = PCZ(:,outliers(~outliers2));
    x = x(outliers(~outliers2));
    y = y(outliers(~outliers2));
    filler = zeros(size(PCX'));
    RHS =  [x(:),y(:),filler+1];
    [PC_x_mod,~]=robustfit(RHS,PCX','bisquare',4.685,0);
    [PC_y_mod,~]=robustfit(RHS,PCY','bisquare',4.685,0);
    [PC_z_mod,~]=robustfit(RHS,PCZ','bisquare',4.685,0);

    deltax = acos(n(1)/norm(n))/degree-90;
deltay = acos(n(2)/norm(n))/degree-90;
deltaz = acos(n(3)/norm(n))/degree-90;
end