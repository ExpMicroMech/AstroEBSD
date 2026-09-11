function structConfig = defaultStructureFactorConfig(options)
     % Return default configuration parameters for the structure factor,
     % such as Voltage, wavelenght and number of bands to plot
     %
     %


    switch nargin
         case 0
             options = struct();
     end
     
     if ~isfield(options,'Voltage')
         options.Voltage = 20E3;
     end
    
    structConfig.V_in = options.Voltage;
    structConfig.lattice_param_scale=1E-10;
    % lambda
    structConfig.lambda = calc_lambda(structConfig.V_in);
    structConfig.lambda_p=structConfig.lambda/structConfig.lattice_param_scale;
    structConfig.maxInd = 30;
    structConfig.cutoff = 0.05;
    structConfig.minD = 0.2;
    % structConfig
    cset=[1,0,0; % Colors of the Bands in rgb (1 =255)
    0,1,0;
    0,0,1;
    1,0,1;
    1,1,1]; %yellow does not work !
    cset=[cset;cset*0.5];
    structConfig.cset = cset;
    structConfig.AlignEdax = 1;
    structConfig.LineStyle = '-';
    structConfig.LineWidth = 1.2;
    structConfig.alpha = 0.4;
end