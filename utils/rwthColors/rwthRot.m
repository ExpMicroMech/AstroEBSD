function color = rwthRot(intensity)
    switch intensity      
        case 100     
            color = [204 ,7 ,30];
        case 75      
            color = [216 ,92 ,65];
        case 50      
            color = [230 ,150 ,121];
        case 25      
            color = [243 ,205 ,187];
        case 10       
            color = [250 ,235 ,227 ];
        otherwise   
             warning('Color not defined for this intensity. Returning default')       
             color=	[204;7;30]; % Default to full intensity (100)        
    end
    color = color/256;
end