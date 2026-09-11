function color = rwthMaiGruen(intensity)
    switch intensity 
        case 100 
            color = [189, 205, 0];
        case 75 
            color = [208, 217, 92];
        case 50 
            color = [224, 230, 154];
        case 25  
            color = [240, 243, 208];
        case 10  
            color = [249, 250, 237];
        otherwise  
            warning('Color not defined for this intensity, returning default')
            color = [189, 205, 0]; % Default to full intensity (100)
    end

    color = color/256;
end