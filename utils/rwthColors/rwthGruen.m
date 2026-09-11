function color = rwthGruen(intensity)
    switch intensity
        case 100 
            color = [87, 171, 39];
        case 75 
            color = [141, 192, 96];
        case 50 
            color = [184, 214, 152];
        case 25 
            color = [221, 235, 206];
        case 10 
            color = [242, 247, 236];
        otherwise  
            warning('Color not defined for this intensity. Returning default')
            color = [87, 171, 39]; % Default to full intensity (100)
    end

    color = color/256;
end