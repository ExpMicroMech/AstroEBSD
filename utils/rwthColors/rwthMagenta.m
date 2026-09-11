function color = rwthMagenta(intensity)
    switch intensity
        case 100 
            color = [227, 0 ,102];
        case 75
            color = [233 ,96 ,136];
        case 50
            color = [241 ,158 ,177];
        case 25
            color = [249 ,210 ,218];
        case 10
            color = [253 ,238 ,240];
        otherwise
            warning('Color not defined for this intensity, returning full intensity')
            color = [227 ,0 ,102]; % Default to full intensity (100)
    end

    color = color/256;
end
