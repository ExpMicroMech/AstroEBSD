function color = rwthGelb(intensity)
    switch intensity
        case 100 
            color = [255 ,237 ,0 ];
        case 75 
            color = [255 ,240 ,85 ];
        case 50 
            color = [255 ,245 ,155 ];
        case 25  
              color=	[255	250	209	];	
        case	10	
              color=	[255	253	238	];	
        otherwise	
              warning('Color not defined for this intensity. Returning full intensity.')
              color=	[255	  237	  0]; % Default to full brightness (100)	
    end	
    color = color/256;

end	
