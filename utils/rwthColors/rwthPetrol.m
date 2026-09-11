function color =	rwthPetrol(intensity)	
       switch	intensity	
          case 100	
              color= [0	 , 97,	  101];	
          case 75	
              color= [45	  ,127,	  131];	
          case 50	
              color= [125,	  164	,  167];	
          case 25		
             color=[191  	,208,  	209 ];		
          case 10		
             color=[230 , 	236  ,	236 ];		
          otherwise		
             warning('Color not defined for this intensity. Returning full brightness.')		
              color= [0	 , 97,	  101];	% Default to full brightness (100)		
        end		
      color = color/256;
end		
