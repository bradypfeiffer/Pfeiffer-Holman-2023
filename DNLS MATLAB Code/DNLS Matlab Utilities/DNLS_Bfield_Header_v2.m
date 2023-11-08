function DNLS_Bfield_Header_v2(PStr, WorA, k, time)
%% Writes Trial header to B-field *.txt file
%

    %  WorA =   'w' to overwrite or create new file ... for header
    %           'a' to append to existing file ... for post run diagnostics
    function [Line] = PrintLine(DummyCon)
        stars = '*************************************\r\n';
        us = '_____________________________________\r\n';
        nl = '\r\n';
        
        %Dummy Identifiers
        DummyTag = [DummyCon('title'),':\r\n'];
        keys = DummyCon.keys;
        values = DummyCon.values;
        sz = size(keys);
        NumLines = sz(2) - 1;
        KeySz = sz(2);
        Line = [stars,DummyTag,us];
        Tspace = zeros(1,NumLines);
        for ii = 1:KeySz
            if ~strcmp(keys(ii),'title')
                Tspace(ii) = length(keys{ii});
            end
        end
        Tmax = max(Tspace);
        TS = cell(1,sz(2) - 1);
        for ii = 1:NumLines
            TS{ii} = repmat(' ',1,(Tmax - Tspace(ii))); 
            TS(ii) = {[TS(ii),'  =  ']};
        end 
        
        lineindex = 0;
        for ii = 1:KeySz
            ii;
            key = keys{ii};
            valuec = values{ii};
            if ~strcmp(key,'title')
                lineindex = lineindex + 1;
%                 if class(valuec) == cell
%                     valuec = valuec{:};
%                 end
                Line = [Line,num2str(key), num2str(TS{lineindex}{1}), num2str(TS{lineindex}{2}), num2str(valuec),nl];
            end
        end
        Line = [Line,nl,nl];
    end
    Trial = PStr.TCon;
    Prfl  = PStr.PrflCon;
    Physics = PStr.PCon;
    Numerics = PStr.NCon;
    Folder = PStr.FCon;
    Diagnostics = containers.Map('null', 'blah');
    DConBool = 0;
    if isfield(PStr,'DCon')
        DConBool = 1;
        Diagnostics = PStr.DCon;
    end
    file = Folder('Bfieldfile');
    elapsed = num2str(cputime - time); t = cputime;
    
    if DConBool == 1
        [DiagnosticsLine] = PrintLine(Diagnostics);
        
        Bfile = fopen(file,WorA);
        fprintf(Bfile,DiagnosticsLine);
        fclose(Bfile);  
    else        
        [TrialLine] = PrintLine(Trial);
        [PrflLine] = PrintLine(Prfl);
        [PhysicsLine] = PrintLine(Physics);
        [NumericsLine] = PrintLine(Numerics);

        Bfile = fopen(file,WorA);
        fprintf(Bfile,TrialLine);
        fprintf(Bfile,PrflLine);
        fprintf(Bfile,PhysicsLine);
        fprintf(Bfile,NumericsLine);
        fclose(Bfile); 
    end
    elapsed = num2str(cputime - t); t = cputime;
end