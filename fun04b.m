function varargout = fun04b(t,y,flag,A1,B1,A2,B2)
%varargout gli argomenti in out sono dinamici
    switch flag %è una sorta di if
        case '' %'' se il flag è la stringa vuota mi da varargout{1}=..
            varargout{1}=f(t,y,A1,B1,A2,B2);
        case 'events' % nel caso il flag sia questa fa..
            [varargout{1:3}]=events(t,y);
        otherwise
            error(['Error message']);
    end
    


end

