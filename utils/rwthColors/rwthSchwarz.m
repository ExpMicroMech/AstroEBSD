function color = rwthSchwarz(intensity)
    switch intensity
        case 100 
            color = [0, 0, 0];
        case 75
            color = [100, 101, 103];
        case 50
            color = [156, 158, 159];
        case 25
            color = [207, 209, 210];
        case 10
            color = [236, 237, 237];
        otherwise
            warning('Color not defined for this intensity, returning full intensity')
            color = [0, 0 ,0]; % Default to full intensity (100)
    end
        color = color/256;
end