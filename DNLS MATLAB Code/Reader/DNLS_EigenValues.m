function DNLS_EigenValues(PStr, BoolCon, EvalCon, ArchiveBase, TrialFolderRel, ...
            TrialNumber, RunNumber, LineNumber, ...
            Energy, Signs1p, NReEvs, ReEvs, LREval, CEVList, debugDNLS_Eigenvalues,...
            msgHalley, Xh, Bh, dispersion, alpha)
    %% Functions
    function [Enrg] = Energy_Calc(Lh,Xh,Bh,Bo)

        DX = Xh(2) - Xh(1);
        stps = length(Bh);

        Enrg = 0;
        for id = 1:stps
            Enrg = Enrg + 0.5*(abs(Bh(id))^2 - abs(Bo(1))^2)*DX;
        end
        Enrg;
        
    end
    
    function [HelicityTotal] = Helicity_Calc(Xh,Bh)
        
        DX = Xh(2) - Xh(1);
        stps = length(Bh);
        
        Enrg = 0;
        HR = 0;
        QR = 0;
        HelicityTotal = 0;
        for idOut = 1:stps
            HR = HR + conj(Bh(idOut))*DX;
            QR = QR + Bh(idOut)*DX;              
            HelicityTotal = HelicityTotal + 1i/2*(Bh(idOut)*HR - conj(Bh(idOut))*QR)*DX;
        end
        
    end

    if debugDNLS_Eigenvalues ==1
        ['Within call to write eigenvalues to file: Start']
    end
    EVFolderRel = [TrialFolderRel,'\EV']; %e.g. '\Trial 11\EV'
    EVFolderAbs = [ArchiveBase,EVFolderRel]; %e.g. 'F:\ ...\Trial 11\'
    EVfilename = [EVFolderAbs, '\EV_',TrialNumber,'_',RunNumber,'.txt'];
    if debugDNLS_Eigenvalues ==1
        ['Within call to write eigenvalues to file: About to create/validate EV folder']
    end
    
    utilities(EVFolderAbs, EVFolderRel); %Creates EV folder if needed
    if debugDNLS_Eigenvalues ==1
        ['Within call to write eigenvalues to file: Declared folders.']
        ['Within call to write eigenvalues to file: Write eigenvalues to screen.']
    end
    
    %% Extract Parameter Containers by category 
        Physics     = PStr.PCon;
        Prfl        = PStr.PrflCon;
        Numerics    = PStr.NCon;
        Graphics    = PStr.GCon;
        Trial       = PStr.TCon;
        Diagnostics = PStr.DCon;
                 
        N = PStr.NCon('N');
        L = PStr.NCon('L');
        dt = PStr.NCon('dt');
        dx = PStr.NCon('dx');
        TimestepAttempted = PStr.NCon('TimestepsAttempted');
        timestepsPerFrame = PStr.NCon('timestepsPerFrame');
        
        TimePerLine = dt*timestepsPerFrame;
        
        zmax = Graphics('zmax');
        zmin = Graphics('zmin');
        
        
        exact = PStr.GCon('DisplayExactSolution');

        TrialID = PStr.TCon('TrialID');
        TrialName = PStr.TCon('TrialName');
        TimeStamp = PStr.TCon('TimeStamp');   
        
        boolBgraph = BoolCon('boolBgraph');
        boolHodogram = BoolCon('boolHodogram');
        boolVideo = BoolCon('boolVideo');
        boolREV = BoolCon('boolREV');
        EVZoom = BoolCon('EVZoom');
        boolEnergy = BoolCon('boolEnergy');
        boolRoot = BoolCon('boolRoot');
        boolCEV = BoolCon('boolCEV');
        boolCEVgraph = BoolCon('boolCEVgraph');
        boolRAD = BoolCon('boolRAD');
        boolRADgraph = BoolCon('boolRADgraph');
        boolHalley = BoolCon('boolHalley');
        SFlag = BoolCon('REV_ErrFlag');
    %% Create Header for Eigenvalue Txt File
    nl = '\r\n\r\n';
    if ~exist(EVfilename)
        EVfile = fopen(EVfilename,'a');
        out = {TrialNumber,RunNumber};
        out = out.';
        fprintf(EVfile,'Results for Trial # %s, Run # %s \n',out{:});
        fclose(EVfile);
    end
    %% Prints real eigenvalues to command window
%     clc;
%     if  boolRoot || boolEnergy

        Xrange = [Xh(1), Xh(end)];
        Bconditions = [Bh(1), Bh(end)];
        fprintf(['*******************************************************\n']);
        fprintf(['Results for:']);
        fprintf(' Trial #%s, ',TrialNumber);
        fprintf('Run #%s, ', RunNumber);
        fprintf(' Line #%s\n', num2str(LineNumber));
        fprintf(['Evaluated on ',datestr(datetime('now')),'\n']);
        fprintf(['Xmin = ',num2str(Xrange(1)),', Xmax = ',num2str(Xrange(2)),' --> DeltaX = ',num2str(Xrange(2) - Xrange(1)),'\n']);
        fprintf(['B_left = ',num2str(Bconditions(1)),', B_right = ',num2str(Bconditions(2)),'\n']);
        fprintf(['*******************************************************\n\n']);
%     end    
       

        if NReEvs > 0

            fprintf(['Number of Real Eigenvalues = ',num2str(NReEvs),'\n']);
            for ev = 1:(length(ReEvs) -1)
                fprintf(' %2.5f, ',ReEvs(ev));
            end
            fprintf(' %2.5f  \n',ReEvs(length(ReEvs)));

    %         fprintf(['Signs of Real Eigenvalues = ',num2str(NReEvs),'\n']);
            for ev = 1:(length(ReEvs) -1)
                if Signs1p(ev) >0
                    fprintf(' bright, ');
                elseif Signs1p(ev) >0
                    fprintf(' dark, ');
                else
                    fprintf(' ????, ');
                end
            end

            if Signs1p(length(Signs1p)) >0
                fprintf(' bright \n');
            elseif Signs1p(length(Signs1p)) <0
                fprintf(' dark \n');
            else
                fprintf(' ???? \n');
            end

            BCratio = abs(Bconditions(1)/Bconditions(2));
            if boolEnergy && ( BCratio < 1.03) && ( BCratio > 0.97 )
                Lh = [Xh(1), Xh(end)];
                Eintegrate = 2*Energy_Calc(Lh,Xh,Bh,Bh(1));
                E_solitons = zeros(1,NReEvs);
                Helicity_solitons = zeros(1,NReEvs);
                for iREidx = 1: NReEvs
                        E_solitons(iREidx) = dispersion/alpha*(2*pi()*(1 + sign(Signs1p(iREidx)))-4*acos(ReEvs(iREidx)/abs(Bh(1))));
                        Helicity_solitons(iREidx) = -2*(dispersion/alpha)^2*sqrt(abs(Bh(end))^2 -ReEvs(iREidx)^2 )/(abs(Bh(end))^2*ReEvs(iREidx));
                end
    %             ratio_energy = E_solitons/Eintegrate;
                Helicity = Helicity_Calc(Xh,Bh);

                fprintf(['Energy = ',num2str(Eintegrate),'\n\n']);
                % Display Energy Analysis
                fprintf('Energy Analysis: \n');
                for iREidx = 1: NReEvs
                    fprintf(['                 ',char(955)','(',num2str(iREidx),') = ',num2str(ReEvs(iREidx)),', sign = ',num2str(sign(Signs1p(iREidx))),', E(',num2str(ReEvs(iREidx)),') = ',num2str(E_solitons(iREidx)),' ... ',num2str(E_solitons(iREidx)/Eintegrate*100),' %% of total energy of ',num2str(Eintegrate),'\n']);
                end
                ratio_energy = sum(E_solitons)/Eintegrate;
                if ( ratio_energy < 1.1 ) && ( ratio_energy > 0.9 )
                    fprintf(['Total soliton energy = ',num2str(sum(E_solitons)),' is ',num2str(sum(E_solitons)/Eintegrate*100),' %% of total energy of ',num2str(Eintegrate),'\n\n']);
                else
                    fprintf(2,['Total soliton energy = ',num2str(sum(E_solitons)),' is ',num2str(sum(E_solitons)/Eintegrate*100),' %% of total energy of ',num2str(Eintegrate),'\n\n']);
                end

                % Display Helicity Analysis
                fprintf('\n');
                fprintf('Helicity Analysis:  \n');
                for iREidx = 1: NReEvs
                    fprintf(['                 ',char(955)','(',num2str(iREidx),') = ',num2str(ReEvs(iREidx)),', helicity = ',num2str(Helicity_solitons(iREidx)),' ... ',num2str(Helicity_solitons(iREidx)/Helicity*100),' %% of total helicity of ',num2str(Helicity),'\n']);
                end
                ratio_helicity = sum(Helicity_solitons)/Helicity;
                if ( ratio_helicity < 1.1 ) && ( ratio_helicity > 0.9 )
                    fprintf(['Total soliton helicity = ',num2str(sum(Helicity_solitons)),' is ',num2str(ratio_helicity*100),' %% of total helicity of ',num2str(Helicity),'\n\n']);
                else
                    fprintf(2,['Total soliton helicity = ',num2str(sum(Helicity_solitons)),' is ',num2str(ratio_helicity*100),' %% of total helicity of ',num2str(Helicity),'\n\n']);
                end

            end   


        else
            fprintf(['No Real Eigenvalues found between ',char(955),'= ',num2str(LREval(1)),' and ',num2str(LREval(2)),' with boolREV = ',num2str(boolREV),'\n']);
        end
    
    if SFlag == 0
        fprintf(['There were no errors above tolerance.\n']);
    elseif SFlag > 0
        fprintf(['There were  ',num2str(SFlag),' ',char(955),' evaluations with errors above tolerance.\n']);
    end
    
    if boolRoot
%         fprintf(EVfile,' \n');
        fprintf(' \n');
%         fprintf(EVfile,'********************************************* \n');
        fprintf('boolRoot = 1 Complex Eigenvalue ... \n');
%         fprintf(EVfile,'********************************************* \n');
        
        fprintf(' %s ',num2str(CEVList{1}));
        if iscellstr(CEVList(1))
            if ~contains(CEVList{1},'<') && ( prod(abs(imag(CEVList{1})) < 0.001 ) ) && prod( real(CEVList{1}) > 0 )
                fprintf('<-- Probable Real Eigenvalue??? \n');
            end
        elseif abs(imag(CEVList{1})) < 0.001 && real(CEVList{1}) > 0
            fprintf('<-- Probable Real Eigenvalue??? \n');
        end
%         if ~contains(CEVList{1},'<') && abs(imag(CEVList{1})) < 0.001 && real(CEVList{1}) > 0
%             fprintf('<-- Probable Real Eigenvalue??? \n');
%         end
        
        fprintf(' \n');
%         fprintf(' \n');
    if boolHalley == 1
        fprintf('________________________________________ \n');
        fprintf('Complex root diagnostics from Halley \n');
        fprintf('________________________________________ \n');
        for nh = 1:length(msgHalley)
            fprintf(msgHalley{nh});
            fprintf('\n');
        end
    end
    else
        fprintf('boolRoot = 0 no search for Complex Eigenvalues\n');
%         fprintf(' \n');
    end
    fprintf(['_______________________________________________________\n']);
    fprintf(['_______________________________________________________\n\n']);
    %% Archives eigenvalues results to Txt file
    if debugDNLS_Eigenvalues ==1
        ['Within call to write eigenvalues to file: Declaring']
    end
%     nl = '\r\n\r\n';
    EVfile = fopen(EVfilename,'a');
    S1ptype = strings(1,NReEvs);
    for ii = 1:length(Signs1p)
        if real(Signs1p(ii)) > 0
            S1ptype(ii) = ' bright ';
        elseif real(Signs1p(ii)) < 0
            S1ptype(ii) = '  dark  ';
        else
            S1ptype(ii) = '  ????  ';
        end
    end
    fprintf(EVfile,' \n');
    fprintf(EVfile,' \n');
    fprintf(EVfile,'**************************************** \n');
    fprintf(EVfile,['Results for Line # ',num2str(LineNumber),' evaluated on ',datestr(datetime('now')),'\n']);
    fprintf(EVfile,['Xmin = ',num2str(Xrange(1)),', Xmax = ',num2str(Xrange(2)),' --> DeltaX = ',num2str(Xrange(2) - Xrange(1)),'\n']);
    fprintf(EVfile,['B_left = ',num2str(Bconditions(1)),', B_right = ',num2str(Bconditions(2)),'\n']);
    fprintf(EVfile,'**************************************** \n\n');

    if NReEvs > 0     
        fprintf(EVfile,'Number of real eigenvalues found between lambda = %s and %s --> N = %s \n',num2str(LREval(1)), num2str(LREval(2)),num2str(NReEvs) );
%         fprintf(EVfile,[' and ',num2str(LREval(2))]);
%         fprintf(EVfile,[' N = %s \n',num2str(NReEvs)]);
%         fprintf(EVfile,' %s',num2str(ReEvs(evn)));
        for evn = 1:NReEvs
            if NReEvs == 1
                fprintf(EVfile,'%s',num2str(ReEvs(evn)));
            end
            if NReEvs > 1 && evn < NReEvs
                fprintf(EVfile,'%s,',num2str(ReEvs(evn)));
            end
            if NReEvs > 1 && evn == NReEvs
                fprintf(EVfile,'%s',num2str(ReEvs(evn)));
            end        
        end
        fprintf(EVfile,'\n');
        
        for evn = 1:NReEvs
            if NReEvs == 1
                fprintf(EVfile,' %s ',S1ptype(evn));
            end
            if NReEvs > 1 && evn < NReEvs
                fprintf(EVfile,' %s,',S1ptype(evn));
            end
            if NReEvs > 1 && evn == NReEvs
                fprintf(EVfile,' %s ',S1ptype(evn));
            end  
            
        end
%         fprintf(EVfile,' %s  \n\n',S1ptype(length(Signs1p)));
%         fprintf(EVfile,' \n');
        fprintf(EVfile,' \n');
    else
        fprintf(EVfile,['No real eigenvalues found between lambda =',num2str(LREval(1))]);
        fprintf(EVfile,[' and ',num2str(LREval(2))]);
        fprintf(EVfile,[', boolREV = %s \n',num2str(boolREV)]);
        fprintf(EVfile,' \n');
    end
    
    BCratio = abs(Bconditions(1)/Bconditions(2));
    if boolEnergy && ( BCratio < 1.03) && ( BCratio > 0.97 )
%         fprintf(EVfile,'Energy =  %s \n\n',num2str(Energy));
        Lh = [Xh(1), Xh(end)];
        Eintegrate = 2*Energy_Calc(Lh,Xh,Bh,Bh(1));
        E_solitons = zeros(1,NReEvs);
        Helicity_solitons = zeros(1,NReEvs);
        for iREidx = 1: NReEvs
                E_solitons(iREidx) = dispersion/alpha*(2*pi()*(1 + sign(Signs1p(iREidx)))-4*acos(ReEvs(iREidx)/abs(Bh(1))));
                Helicity_solitons(iREidx) = -2*(dispersion/alpha)^2*sqrt(abs(Bh(end))^2 -ReEvs(iREidx)^2 )/(abs(Bh(end))^2*ReEvs(iREidx));
        end
%             ratio_energy = E_solitons/Eintegrate;
        Helicity = Helicity_Calc(Xh,Bh);

        fprintf(EVfile,['Energy = ',num2str(Eintegrate),'\n\n']);
        % Display Energy Analysis
        fprintf(EVfile,'Energy Analysis: \n');
        for iREidx = 1: NReEvs
            fprintf(EVfile,['lambda(',num2str(iREidx),') = ',num2str(ReEvs(iREidx)),', sign = ',num2str(sign(Signs1p(iREidx))),', E(',num2str(ReEvs(iREidx)),') = ',num2str(E_solitons(iREidx)),' ... ',num2str(E_solitons(iREidx)/Eintegrate*100),' %% of total energy of ',num2str(Eintegrate),'\n']);
        end
        ratio_energy = sum(E_solitons)/Eintegrate;
        if ( ratio_energy < 1.1 ) && ( ratio_energy > 0.9 )
            fprintf(EVfile,['Total soliton energy = ',num2str(sum(E_solitons)),' is ',num2str(sum(E_solitons)/Eintegrate*100),' %% of total energy of ',num2str(Eintegrate),'\n']);
        else
            fprintf(EVfile,['Total soliton energy = ',num2str(sum(E_solitons)),' is ',num2str(sum(E_solitons)/Eintegrate*100),' %% of total energy of ',num2str(Eintegrate),'\n']);
        end

        % Display Helicity Analysis
        fprintf(EVfile,'\n');
        fprintf(EVfile,'Helicity Analysis:  \n');
        for iREidx = 1: NReEvs
            fprintf(EVfile,['lambda(',num2str(iREidx),') = ',num2str(ReEvs(iREidx)),', helicity = ',num2str(Helicity_solitons(iREidx)),' ... ',num2str(Helicity_solitons(iREidx)/Helicity*100),' %% of total helicity of ',num2str(Helicity),'\n']);
        end
        ratio_helicity = sum(Helicity_solitons)/Helicity;
        if ( ratio_helicity < 1.1 ) && ( ratio_helicity > 0.9 )
            fprintf(EVfile,['Total soliton helicity = ',num2str(sum(Helicity_solitons)),' is ',num2str(ratio_helicity*100),' %% of total helicity of ',num2str(Helicity),'\n\n']);
        else
            fprintf(EVfile,['Total soliton helicity = ',num2str(sum(Helicity_solitons)),' is ',num2str(ratio_helicity*100),' %% of total helicity of ',num2str(Helicity),'\n\n']);
        end
        
%         if NReEvs > 1
%             idx = linspace(1,NReEvs,NReEvs);
%             C = cell(1,NReEvs);
%             for idxN = 1:(NReEvs)
%                 C{idxN} = nchoosek(idx,idxN);
%             end
%             C
% 
%         end
    end
    
    if SFlag == 0
        fprintf(EVfile,'There were no errors above tolerance.\n');
    elseif SFlag > 0
        fprintf(EVfile,['There were %s  lambda-evaluations with errors above tolerance.\n',num2str(SFlag),char(955)]);
    end
    if boolCEV
        CEVTemp = '';
        fprintf(EVfile,[nl,'CEVList',nl]);
        for ci = 1: length(CEVList)
            if iscell(CEVList)
                ceentry = CEVList{1,ci};
            end
            if isnumeric(CEVList)
                ceentry = CEVList(1,ci);
            end
            
            if ~isempty(ceentry)
                if isnumeric(ceentry)
                    fprintf(EVfile,' \n');
                    dlmwrite(EVfilename, ceentry, 'precision', '%2.5f','-append', 'delimiter', ' ');
                    fprintf(EVfile,' \n');
%                     CEVTemp = [CEVTemp, num2str(ceentry)];
                end
                if ischar(ceentry)
%                     dlmwrite(EVfilename,[nl, ceentry, nl]);
                    fprintf(EVfile,[nl,ceentry,nl]);
%                     CEVTemp = [CEVTemp, ceentry];
                end
            end
        end
%         dlmwrite(EVfilename,CEVTemp,'precision','%.4f','-append','delimiter',' ');
    end
    if boolRoot
%         fprintf(EVfile,' \n');
%         fprintf(EVfile,' \n');
%         fprintf(EVfile,'********************************************* \n');
        fprintf(EVfile,'Complex Eigenvalue with boolRoot = 1 ... \n');
%         fprintf(EVfile,'********************************************* \n');
        
        fprintf(EVfile,' %s ',num2str(CEVList{1}));
        
        if iscellstr(CEVList{1})
            if ~contains(CEVList{1},'<') && abs(imag(CEVList{1})) < 0.001 && real(CEVList{1}) > 0
                fprintf(EVfile,'<-- Probable Real Eigenvalue??? \n');
            end
        elseif abs(imag(CEVList{1})) < 0.001 && real(CEVList{1}) > 0
            fprintf(EVfile,'<-- Probable Real Eigenvalue??? \n');
        end
        
%         if ~contains(CEVList{1},'<') && abs(imag(CEVList{1})) < 0.001 && real(CEVList{1}) > 0
%             fprintf(EVfile,'<-- Probable Real Eigenvalue??? \n');
%         end
        
        fprintf(EVfile,' \n');
%         fprintf(EVfile,' \n');
        if boolHalley == 1
            fprintf(EVfile,'________________________________________ \n');
            fprintf(EVfile,'Complex root diagnostics from Halley \n');
            fprintf(EVfile,'________________________________________ \n');
            for nh = 1:length(msgHalley)
                fprintf(EVfile,msgHalley{nh});
                fprintf(EVfile,'\n');
            end
        end
    else
        fprintf(EVfile,'boolRoot = 0, No search for Complex Eigenvalues\n');
%         fprintf(EVfile,' \n');
    end
    fprintf(EVfile,'________________________________________ \n');
    status = fclose('all');
    status;


end

