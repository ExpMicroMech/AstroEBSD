function color = rwthOrange(intensity)
    switch intensity      
        case 100      
            color = [246, 168, 0];
        case 75      
            color = [250, 190, 80];
        case 50      
            color = [253, 212, 143];
        case 25      
            color = [254, 234, 201];
        case 10      
            color = [255, 247, 234];       
        otherwise      
            warning('Color not defined for this intensity. Returning default')       
            color = [246 ,168 ,0]; % Default to full intensity (100)        
    end      

    color = color/256;
end   