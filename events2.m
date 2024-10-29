function [value,isterminal,direction] = events2(t,y,C)
   
    value=y-C;
    isterminal=0; 
    %is terminal =0 non si ferma all'evento 
    % isterminal =1 si ferma all'evento
    direction=0; % non gli interessa se cresce o decresce
    
end

