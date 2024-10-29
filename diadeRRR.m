function [p3] = diadeRRR(p1,p2,d13,d23,config)

% function che risolve il problema della diade RRR 
% chiede in input 2 punti le distanze dal 3 punto e la configurazione
    r1=d13;
    r2=d23;
    d12=sqrt( (p2(1) - p1(1) )^2 + (p2(2) - p1(2) )^2 );
    lambda = ( 1/2 ) * ( 1 + ( r1^2 - r2^2 ) / (d12^2) );

    if config == 1
        mu= sqrt( (r1^2) / (d12^2) - lambda^2 );
    else
        mu= - sqrt( (r1^2) / (d12^2) - lambda^2 );
    end
    
    %sistema locale
    p13x=lambda*d12;
    p13y=mu*d12;
    p13loc=[p13x,p13y];
    %calcolo rotazione
    m= ( p1(2) - p2(2) )/ ( p1(1) - p2(1) );
    x = - atan(m);
    %ruoto il SdR 
    matricedirotazione = [ cos(x) -sin(x); sin(x) cos(x)]; 
    p13=p13loc*matricedirotazione;
    %sommo i vettori 01 13
    p3=p1+p13;
    
end

