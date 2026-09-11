function color = rwthLila(intensity)
    switch intensity
        case 100
            color = [122, 111, 172];
        case 75
            color = [155, 145, 193];
        case 50
            color = [188, 181, 215];
        case 25
            color = [222, 218, 235];
        case 10
            color = [242, 240, 247];
        otherwise
            warning('Color not defined for this intensity. Returning default')
            color = [122, 111, 172]; % Default to full intensity (100)
    end

    color = color/256;
end