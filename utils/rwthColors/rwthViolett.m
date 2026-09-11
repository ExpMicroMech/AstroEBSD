function color = rwthViolett(intensity)
    switch intensity         
         case 100           
                 color =[97 ,33 ,88];       
             case 75          
                 color =[131 ,78 ,117];     
            case 50           
                 color =[168 ,133 ,158];     
             case 25 
                 color =[210,192,205];     
             case 10 
                 color =[237,229,234];        
        otherwise           
                warning('Color not defined for this intesity.returning default')    
                color= [97,33,88];%Default to maximum brightness(100);    
        end
% end      
    color = color/256;
end 