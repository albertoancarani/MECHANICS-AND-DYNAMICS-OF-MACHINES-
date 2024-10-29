function varargout = fun04c(t,y,flag,A,B,C)
%C= 95% del regime
%varargout gli argomenti in out sono dinamici
    switch flag %è una sorta di if
        case '' %'' se il flag è la stringa vuota mi da varargout{1}=..
            varargout{1}=f(t,y,A,B);
        case 'events' % nel caso il flag sia questa fa..
            [varargout{1:3}]=events(t,y,C);
        otherwise
            error(['Error message']);
    end
    


end



