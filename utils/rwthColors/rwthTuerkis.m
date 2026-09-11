
function	color=rwthTuerkis(intensity)			
         switch intensity			
             case 100			
                    color=[0,152,161];	
             case 75				
                   color=[0, 177, 183];			
             case 50			
                  color=[137,204,207];				
             case 25
                  color=[202,231,231];
             case 10
                  color=[235,246,246];
             otherwise
                  warning('color not defined for this intensitiy. returnig default') ,
            color=[0, 152, 161];
         end
         color = color/256;
end