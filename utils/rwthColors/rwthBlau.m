function color = rwthBlau(intensity)
    switch intensity
        case 100 
            color = [0, 84, 159];
        case 75
            color = [64, 127, 183];
        case 50
            color = [142, 186, 229];
        case 25
            color = [199, 221, 242];
        case 10
            color = [232, 241, 250];
        otherwise
            warning('Color not defined for this intensity. Returning default')
            color = [0, 84, 159]; % Default to full intensity (100)


    end
    color = color/256;
end