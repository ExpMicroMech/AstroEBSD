function color = rwthBordeaux(intensity)
    switch intensity         
         case 100          
             color=	[161,16,53];         
         case 75         
                color=	[182,82,86];         
            case 50          
                color=	[205,139,135];         
            case 25
                color=	[229,197,192];         
            case 10
                color=	[245,232,229];         
            otherwise         
               warning('Color not defined for this intensity. Returning default')         
                  color=	[161 ,16 ,53 ]; % Default to full brightness (100)         
    end         
        color = color/256
end      