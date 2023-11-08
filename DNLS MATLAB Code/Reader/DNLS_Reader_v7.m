function  [RealEV_List, EVLeft_return, EVRight_return, PStr] = DNLS_Reader_v7(BoolCon, EvalCon, ArchiveBase, TrialNumber, RunNumber, LineNumber, LNstart, LNend, XLeft, XRight )
%% Primary change: Implementation of Newton's Method

% see Comments: 7/10/2017   
%% Comments: 7/10/2017
% v7 attempts to extend the allowed boundary conditions to the fully
% asymmetric case.
%
% v6 added ability to read PStr parameters as well as individual B-field
% lines from a '*.mat' mat file format. This file format is native to
% Matlab and: is an OS independent archiving format; streamlines and
% optimizes read/write time and memory usage.
%% Comments: 6/1/2017
% 1) It was observed that N = 4096, L = 90, Courant = 0.03 for the
% 2p-soliton with lambda = 0.5 + i*0.5 (intended to reproduce Fig. 2 of the
% Magnetic Hole Formation paper) created sufficiently high numerical error
% that the complex eigenvalue dynamics did not match those published.
%
% However, by increasing N 4x, the numerical error in the time development
% of the archived B-field was found to have the complex eigenvalue at least
% reasonably close to the published results.
%
% It seems likely the error was created on the 'Writer' side as the
% 'Reader' uses an error-checking integration method which has an interpolation 
% feature.
%
% --> Add a feature to the Reader eigenvalue method which allows it to vary
% the resolution of the B-field. For example if B(n) = B(xn) with xn = -L +
% n*Dx, then allow B(mn) = B(xmn)  with xmn = -L + (m*n)*Dx where m could
% take on a single value of m = 1, 2, 4, 8, 16 etc. to which n would range
% from n = 1 up to (N/m). Currently the method has a default setting of m =
% 4 as the Writer method oversamples by 4x allowing 3/4 of the
% corresponding high k index components of the FFT to be zero'd out to
% prevent aliasing errors. Allowing m to be 8 or 16 would GREATLY speed up
% the process of evaluating eigenvalues. It would also allow a
% determination to be made if this results in an acceptable numerical error
% on the Reader side.
%
% 2) It is observed that the algorythms used in determining  both the real
% and the complex eigenvalues are quite slow. The real eigenvalue method
% uses the bisection method which is not terribly inefficient, but perhaps
% could be improved even so.
%
% The complex eigenvalue algorythm is extremely slow. It currently
% evaluates the scattering coefficient, s11(lambda), at every point on a
% grid in some rectangular region of the complex plane. The current default
% is for the grid to have Nc*Nc evaluations where Nc might be 25, for
% example. For N, mentioned in (1) above, with N = 4096/4, this has been
% observed to take as long as 30 minutes or more to evaluate this rather
% course grid. If the grid size is '1' on a side the grid spacing is 0.04
% which produces an eigenvalue result with an error that might be as much
% as 10% for a middling target value of, say, 0.5 + i*0.5. To reduce this
% error by a factor of 2 requires a factor of 2 increase of Nc which would
% increase the computation time by 4-fold. A more efficient algorythm is
% required.
%
% 2a)   Kayla Toll has suggested an intuitive 'boxing' method which
% evaluates the signs of the real and imaginary parts of s11 just at points 
% on the boundary of a region. There are a handful of cases where one or
% more zeroes are inside the 'box' or just outside it. Each of these cases
% creates conditions, or patterns, of the signs of Re(s11) and Im(s11) ...
% or rather how they change ... as one moves continuously around the
% boundary of the box. If this check indicates no eigenvalue is enclosed,
% the algorythm is terminated and required less than 4*Nc evaluations of
% s11.
%
% If there is one (or more) eigenvalue(s) enclosed, divide the box in half
% with (Nc-2) additional evalutions and check the left hand side then, if
% necessary, the right hand side. Whichever side has an eigenvalue, divide
% that box in half with an additional (Nc/2 - 2) evaluations of s11 and
% repeat. A sketch of this process indicates the potential to produce a
% result comparable to the current method with about 1/4th the calcuations,
% so perhaps in a time about 1/4th the time required in the current method.
% Depending on the precision required the time required could perhaps be
% either more or less than this 1/4th estimate. To be sure, 10 minutes is
% much better than 40 minutes! This is mostly likely true if there is a
% single complex eigenvalue. If there are multiple complex eigenvalues the
% speed-up will likely be reduced.
%
% 2b) There are algorythms for finding zeroes of complex functions such as
% Newton's Method or Halley's Method (and many others). Each of these
% methods requires evaluation of the derivative of s11 at each point in the
% iteration. For a method requiring just the first derivate, each iteration
% requires 3 evaluations of s11. A second derivate also requires 3
% evaluations. A third derivative requires 5 evaluations of s11. Newton's
% method (needs f') has 3 s11 evaluations per iteration. Halley's method
% (needs f''') has 5 evaluations of s11 per iteration. It is claimed that
% Newton's method has quadratic convergence and Halley's method has cubic
% convergence. It will require testing to see which would be more efficient
% or effective and the effectiveness surely depends on how close the
% initial guess is for the first iterant, z0.
%
% It is possible that a hybrid approach with (2a) providing a narrowed
% region for the zero, that is an initial estimate for z0, followed by (2b) 
% to produce a convergent estimate.
%
% 2c) An alternate approach is to use the 'Argument principle':
% https://en.wikipedia.org/wiki/Argument_principle
% It would need to be tested to determine if this provides any efficiency
% beyond Kayla's method (2a). It is essentially a mathematically formal
% approach to the same 'boxing' method as Kayla's. Kayla's method only
% requires knowledge of the signs of Re(s11) and Im(s11) around the
% boundary while the Argument principle requires numerical values with
% some, unknown precision, and with some, unknown, sampling distance, dz.
% There are two potential advantages:
% --> The result of the Argument principle, if executed with sufficient
% accuracy, will yield the number of eigenvalues in a region.
% --> The integral required of the Argument principle can be evaluated
% using built in matlab routines which monitor precision and where the
% precision can be specified by the user.
%
% For the immediate task, we have an expectation of what the initial
% complex eigenvalue should be so ...
%
%****************************************
%****************************************
% Implement Newton's Method taking an initial z0 specified by the user for 
% the initial profile (aka 'line 1). Use the result of Newton's Method for 
% line 1 as the initial z0 for the next time step to be evaluated.
%****************************************
%**************************************** 
%% Comments: 6/13/2016
    %6/13/2016
    %Primary goal: allow integration over a specified range:
    %           x: EVL --> EVR
    % the range where the B-field goes asymptotically from
    % B(x = EVL) ~ Bo*Exp[I*theta_M] at the 'far left' ... to ...
    % B(x = EVR) ~ Bo*Exp[I*theta_P] at the 'far right' 
    % At the moment EVL and EVR are specified by the 'Writer' program and
    % stored in the Excel Archive, but these values could be determined
    % within 'Reader' if useful.
    %  --> The values will be passed as a single argument, Lh = [EVL, EVR]
    %      keeping the argument structure of the previous functions
    %      unchanged.
    %  --> The one place this range is required is in the function Psiminus 
    %      which passes them to ode113 through the argument 'range',
    %      previously range = [-L, +L], but now simply assign range = Lh
    %  --> The slightly more detailed task of modifying 'B4 and X4' to
    %  correspond to the range, Lh, also needs to be addressed. It seems
    %  sensible to determine the indices of X4 corresponding to Lh and
    %  define Bh and Xh as corresponding sections of B4 and X4.    
    %1) Began work to handle the various trial parameters in a more
    %systematic way. Created a structure array, 'PStr' keyed on the parameter
    %containers: Diagnostics, Folder, Graphics, Numerics, Physics etc.
    %The named containers are assigned values in a call to
    %'DNLS_Excel_read_v2.m' and, in turn, are assigned as keys to the
    %single structure, PStr, which is the single thing returned.
    % --> It may be worthwhile to create a mirror of this structure in the
    %Writer program.
    % --> It will be necessary to archive the results found from the Reader
    % program along with the various parameters (e.g. range of complex
    % plane searched, grid size, values for 'opts' used in ODE45, range
    % along the real and imaginary branch for which the radiative
    % scattering data was evaluated, archived trial # etc. being analyzed
    % ...) Presumably these 'Reader' parameters will be written to an Excel
    % file analogous to the 'Writer' version and the results stored in a
    % 'Results' subfolder within the same Trial folder being studied
    %
    %Future Work:
    % 0) Validate method for determining eigenvalues for RABC using
    % 'helical' initial profile.
    % 1) Integrate all user declared variables to be set from 'Call_DNLS_Reader'
    % 2) Implement Zach's EV evaluation loop with method to write results
    % to file for archiving
    % 3) Develop labeling features for graphs: Titles, axes, parameters,
    % updates, location of complex eigenvalues
    % 4) Develop method to evaluate contribution of  'Radiative' Scattering Data
    % to 'energy' and net B-field rotation.
    % 5) Develop method for evaluating solitons positions in x ... and
    % recording these to an archive file            
%% FUNCTION DEFINITIONS  
% *************************************************************

    function [xf, bf, bfCalc] = MaxMinLinearRegression(x, fidx, f)
        xf = x(fidx);
        Xf = [ones(length(xf),1), xf'];
        Yf = f(fidx);        
        bf = Xf\Yf';
        bfCalc = (Xf*bf)';
    end
    function [fabsMinIdxs, fabsMaxIdxs] = FindMaxMinIndices(f, Np)
        Nmm = length(f);
        % if length(f) is odd, reduce its length by 1
        % ... by dropping it's last data point
        if mod(Nmm,2) == 1
            f = f(1: end -1);
            Nmm = Nmm - mod(Nmm,2);
        end
        Kmmx = [0:Nmm/2-1,-Nmm/2:-1] * pi/L; % Nmm defined above to be even
        fabs = abs(f);
        fabsPrimed = real( ifft( 1i*Kmmx.*fft( fabs ) ) );
        fabsDoublePrimed = real( ifft( -Kmmx.*Kmmx.*fft( fabs ) ) );
        fabsPrimedShift = circshift(fabsPrimed,1);
        fabsPrimeSignChanges = fabsPrimed.*fabsPrimedShift;
        fabsMaxMinIdxs = find(fabsPrimeSignChanges < 0);
        
        
        fabsMaxValues = -1*ones(1,Np);
        fabsMaxIdxs = zeros(1,Np);
        idxMValues = 0;
        for idxcktmp = 1:length(fabsMaxMinIdxs) 
            idxck = fabsMaxMinIdxs(idxcktmp); %idxCL is the actual index of PLb that may be a max or min
            maxtmp = fabs(idxck);
            imaxtmp = idxck;
            if idxck>2 && idxck < (length(fabs) - 2) && fabsDoublePrimed(idxck) < 0 % check that idxCK is not on the boundary and that it is a maximum

                if idxMValues == Np
                    idxcktmp = length(fabsMaxMinIdxs)+2;
                    continue;
                end
                              
                for idxtmp = (idxck -2):(idxck + 2)
                    if fabs(idxtmp) > maxtmp
                        maxtmp = fabs(idxtmp);
                        imaxtmp = idxtmp;
                    end                  
                end
                
                idxMValues = idxMValues + 1;
                fabsMaxValues(idxMValues) = maxtmp;
                fabsMaxIdxs(idxMValues) = imaxtmp;                               
            end          
        end
        
        fabsMaxValues = fabsMaxValues( fabsMaxValues > 0 );
        fabsMaxIdxs   = fabsMaxIdxs( fabsMaxIdxs > 0);
        
        fabsMinValues = -1*ones(1,Np);
        fabsMinIdxs = zeros(1,Np);
        idxMValues = 0;
        for idxcktmp = 1:length(fabsMaxMinIdxs)
            idxck = fabsMaxMinIdxs(idxcktmp); %idxCL is the actual index of PLb that may be a max or min
            mintmp = fabs(idxck);
            imintmp = idxck;
            if idxck>2 && idxck < (length(fabs) - 2) && fabsDoublePrimed(idxck) > 0 % check that idxCK is not on the boundary and that it is a minimum
                
                if idxMValues == Np
                    idxcktmp = length(fabsMaxMinIdxs)+2;
                    continue;
                end
                
                for idxtmp = (idxck -2):(idxck + 2)
                    if fabs(idxtmp) < mintmp
                        mintmp = fabs(idxtmp);
                        imintmp = idxtmp;
                    end
                end
                
                idxMValues = idxMValues + 1;
                fabsMinValues(idxMValues) = mintmp;
                fabsMinIdxs(idxMValues)   = imintmp;
            end
        end
        
        fabsMinValues = fabsMinValues( fabsMinValues > 0 );
        fabsMinIdxs   = fabsMinIdxs( fabsMinIdxs > 0);       
        
    end
    function [fabsMaxMinIdxs, fabsDoublePrimed] = FindMaxMinIndicesOriginal(f)
        Nmm = length(f);
        % if length(f) is odd, reduce its length by 1
        % ... by dropping it's last data point
        if mod(Nmm,2) == 1
            f = f(1: end -1);
            Nmm = Nmm - mod(Nmm,2);
        end
        Kmmx = [0:Nmm/2-1,-Nmm/2:-1] * pi/L; % Nmm defined above to be even
        fabs = abs(f);
        fabsPrimed = real( ifft( 1i*Kmmx.*fft( fabs ) ) );
        fabsDoublePrimed = real( ifft( -Kmmx.*Kmmx.*fft( fabs ) ) );
        fabsPrimedShift = circshift(fabsPrimed,1);
        fabsPrimeSignChanges = fabsPrimed.*fabsPrimedShift;
        fabsMaxMinIdxs = find(fabsPrimeSignChanges < 0);        
    end

    function [B4fun, BFlinesfun] = FindBfieldLines(BfieldMatfile, BFieldFile, LineNumber, N)
        % Modify to include LineNumberStart, LineNumberEnd and
        % LineIncrement to return B4All ... all F-field lines
        % LNStart:LIncrement:LNEnd
        
        ['From: ',thisFile,' Begin Read B-field Line ...  ']
        if exist(BfieldMatfile,'file') && 0
            ['From: ',thisFile,' Begin Read B-field Line. Trying Matfile Archive  ']
            BFMat = matfile(BfieldMatfile);
            BfieldNames = fieldnames(BFMat);
            idxRun = strfind(BfieldNames,'Run');
            if str2num(RunNumber) <= length(BFMat.PStr(:,1))
                runN = ['Run',RunNumber];
            end
            if sum(~cellfun('isempty',strfind(BfieldNames,runN))) > 0
                if LineNumber <= length(BFMat.(runN)(:,1))              
                    B4 = BFMat.(runN)(LineNumber,:);
                    B4All = BFMat.(runN);
                    ;
                else
                    ['Line # not in Matfile']
                    length(BFMat.(runN)(:,1))
                    return;
                end
            else
                ['Run # not in file/nr']
                ['Field Names:']
                fieldnames(BFMat)
                return;
            end
        else 
            %Read all lines from Bfield file for given TrialNumber and RunNumber
            ['From: ',thisFile,' Begin Read B-field Line. Trying text file archive']
            BFh = fopen(BFieldFile);
            if BFh < 1
                [BFieldFile, ' does not exist.']
                return
            end
            BFline = fgetl(BFh); % each call to fget1(BFh) assigns the next line in the text file to BFline
            BFlinesfun = cell(0,1);
            NumOfLinesfun = 0;
            while ischar(BFline) % This loop parses through until it finds a B-field line
                                 % and appends it to BFlinesfun
                BFline = fgetl(BFh);
                if strfind(BFline,'LineStart')>0
                    NumOfLinesfun = NumOfLinesfun + 1;
                    BFlinesfun{end+1,1} = fgetl(BFh);
                    ['From: ',thisFile,' verifying B-field Line # ',num2str(NumOfLinesfun)]
                end
            end
            fclose(BFh);
        
            %Find Bfield profile for Line n  
            LineStartNfun = ['LineStart',num2str(LineNumber)];
            numlines = length(BFlinesfun); %returns number of lines in BFh
            IndexBFfun = strfind(BFlinesfun, LineStartNfun); %returns line # of flag before Bfield Line n
            lineN = find(not(cellfun('isempty',IndexBFfun))) + 1; %gives line # of nth Bfield Line
            if LineNumber < (NumOfLinesfun + 1)
                ['From: ',thisFile,' loading B-field Line # ',num2str(LineNumber)]
                temp = BFlinesfun(LineNumber);
                C = textscan(temp{1},'%f');
                B4fun = C{1}.'; %B is an N/4 complex valued, double precision array
            else
                ['Line # ', num2str(LineNumber),' is greater than the ', num2str(NumOfLinesfun), ' lines in the file.' ]
                return;
            end
        end
        if length(B4fun)> (N/4)
            B4temp = B4fun(1:round(length(B4fun)/(N/4)):end);
            B4fun = B4temp;
        end

    end

    function [LReBvals, RhoR, AbsRhoR, LImBvals, RhoI, AbsRhoI] = Rho(NReBranch, NImBranch, LReBranchEval, LImBranchEval, RBeps,IBeps, Lh, Bh, Xh,B0,theta_M,theta_P,opts)
        %Evaluate RhoI = s21/s11 to right of imaginary branch
        LReBvals = linspace(LReBranchEval(1),LReBranchEval(2),NReBranch);        
        RhoR = zeros(1,NReBranch);
        t = cputime;
        star1 = '';
        for kk = 1:NReBranch/1
            star1 = [star1,'*'];
        end
        WhatsHappening = 'Evaluating Scattering Data above real branch cut:';
        star1 = [star1,'\r'];
        star2 = '';
        TimeLeft = 'Remaining Time ... ';
        WittyComment = '';
        Taverage = 0;
        for ii = 1:NReBranch
                l =LReBvals(ii) + 1i*RBeps;
                S = S11_S21(l,Lh,Xh,Bh,BC,opts);
                s11 = S*[1,0]';
                s21 = S*[0,1]';
                RhoR(ii) = s21/s11;
            clc
            elapsed = cputime - t; t = cputime;
            increase = elapsed/(Taverage + 1.01);
            Taverage = (Taverage*(ii-1) + elapsed)/ii;
            RemainingTime = (NReBranch - ii)*Taverage;
            WittyComment = Snappy(RemainingTime,increase);
            TimeTag = [WhatsHappening,'\r',WittyComment,TimeLeft, num2str(RemainingTime),' s'];
            if mod(ii,1) == 0
                star2 = [star2,'*'];
            end
            fprintf([TimeTag, '\r',star1,star2,'\r',star1]);
        end  
        AbsRhoR = abs(RhoR);
        
         %Evaluate RhoI = s21/s11 to right of imaginary branch
        LImBvals = linspace(LImBranchEval(1),LImBranchEval(2),NImBranch);        
        RhoI = zeros(1,NImBranch);
        t = cputime;
        star1 = '';
        for kk = 1:NImBranch/1
            star1 = [star1,'*'];
        end
        WhatsHappening = 'Evaluating Scattering Data along imaginary branch cut:';
        star1 = [star1,'\r'];
        star2 = '';
        TimeLeft = 'Remaining Time ... ';
        WittyComment = '';
        Taverage = 0;
        for ii = 1:NImBranch
                l =1i*LImBvals(ii) + IBeps;
                S = S11_S21(l,Lh,Xh,Bh,BC,opts);
                s11 = S*[1,0]';
                s21 = S*[0,1]';
                RhoI(ii) = s21/s11;
            clc
            elapsed = cputime - t; t = cputime;
            increase = elapsed/(Taverage + 1.01);
            Taverage = (Taverage*(ii-1) + elapsed)/ii;
            RemainingTime = (NImBranch - ii)*Taverage;
            WittyComment = Snappy(RemainingTime,increase);
            TimeTag = [WhatsHappening,'\r',WittyComment,TimeLeft, num2str(RemainingTime),' s'];
            if mod(ii,1) == 0
                star2 = [star2,'*'];
            end
            fprintf([TimeTag, '\r',star1,star2,'\r',star1]);
        end  
        AbsRhoI = abs(RhoI);    
    end
        
    function [CEVList, NCE] = ComplexEVFinder(LComplexEval,NCEval,CevState)
        %Sort through CevState to find the (ii,jj) indices at the center of
        %any 3x3 grid containing all four states.
        %Store resulting (ii,jj) indices in a 2xNCEval matrix = IndexListTemp
        IndexListTemp = zeros(2,NCEval);
        EVtempcount = 0;
         for ii = 2:(NCEval - 1)
             for jj = 2:(NCEval - 1)
                 Factor = CevState(ii,jj)*CevState(ii,jj+1)*CevState(ii,jj-1)*...
                     CevState(ii-1,jj)*CevState(ii-1,jj+1)*CevState(ii-1,jj-1)*...
                     CevState(ii+1,jj)*CevState(ii-1,jj+1)*CevState(ii+1,jj-1);
                 if mod(Factor,210) == 0
                     EVtempcount = EVtempcount + 1;
                     IndexListTemp(:,EVtempcount) = [ii,jj];
                 end
             end
         end
%       IndexListTemp  
%       EVtempcount
        if EVtempcount == 0
            CEVList(1,1) = 0;NCE = 0;
            return
        end
        if EVtempcount == 1
            II = IndexListTemp(1,1);
            JJ = IndexListTemp(2,1);
            l =ComplexEV(LComplexEval, NCEval, II,JJ);
            CEVList(1,1) = l;NCE = 1;
            return
        end
        
%      IndexListTemp
%      EVtempcount
         %Create IndexList by removing blank elements from IndexListTemp
         
         IndexList = zeros(2,EVtempcount);
         count = 0;
         for ii = 1:length(IndexListTemp(1,:))
             if IndexListTemp(1,ii) ~=0
                 count = count + 1;
                 IndexList(:,count) = IndexListTemp(:,ii);
             end
         end
%          IndexList
%          EVtempcount

         %Remove multiple indice listings of single eigenvalues
         %Sort through list and pass over any candidate index that is adjacent to
         %another entry on the list
         %If multiple entries exist, the resulting array will have zeroes
         SortedIndexListTemp = zeros(2,EVtempcount);
         EVCount = 0;
         for ii = 1:length(IndexList)
             EVi = IndexList(:,ii);
             keep = 1;
             for jj = (ii + 1): length(IndexList)
                 EVj = IndexList(:,jj);
                 if ((EVi(1) - EVj(1))^2 + (EVi(2) - EVj(2))^2) < 1.5^2
                     keep = keep - 1;
                 end
             end
             if keep == 1          
                 EVCount = EVCount + 1;
                 SortedIndexListTemp(:,EVCount) = EVi;
             end
         end
%           SortedIndexListTemp
%           EVCount

         %Remove zero listings from SortedIndexListTemp and create a single row
         %array whose entries are the complex valued eigenvalues
         %The final result will be CEVList
         CEVList = zeros(1,EVCount);
         count = 0;
         for ii = 1:length(SortedIndexListTemp(1,:))
             if SortedIndexListTemp(1,ii) ~=0
                 count = count + 1;
                 II = SortedIndexListTemp(1,ii);
                 JJ = SortedIndexListTemp(2,ii);
                 l =ComplexEV(LComplexEval, NCEval, II,JJ);
                 CEVList(1,count) = l;
             end
         end   
%          CEVList
    end

    function [Enrg] = Energy_Calc(Lh,Xh,Bh,Bo)
        
        DX = Xh(2) - Xh(1);
        stps = length(Bh);
        
        Enrg = 0;
        for id = 1:stps
            Enrg = Enrg + 0.5*(abs(Bh(id))^2 - abs(Bo(1))^2)*DX;
        end
        Enrg;  
    end

    function [FELan] = FELandau_Calc(Lh,Xh,Bh,Bo)
        
        stps = length(Bh);
        if mod(stps,2) ==1
            Bh = Bh(2:end);
            stps = length(Bh);
            Xh = Xh(2:end);
        end
        DX = Xh(2) - Xh(1);
        LH = (  Xh(end) - Xh(1)  )/2;
        stps = length(Bh);
        Kh = [0:Nh/2-1,-Nh/2:-1] * pi/LH;

        Fbsq = fft(abs(Bh.^2));
        FHilbertbsq = -1i*sign(Kh).*Fbsq;
        Hilbertbsq = ifft(FHilbertbsq);
        Fbsqx = -1i*Kh.*Fbsq;
        bsqx = ifft(Fbsqx);
        
        FELan = 0;
        for id = 1:stps
            FELan = FELan - 0.5*(Hilbertbsq(id).*bsqx(id))*DX;
            % FELan = 42;
        end
        FELan;  
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

    function [CEVList] = Newton(z,zev,Lh,Xh,Bh,BC,opts,Nnewt)
        
        CEVList = cell(1,3);
        TempEV = zeros(1,3);
        %zev = 0.5 + 1i*0.5;
        z0 = z;% the initial guess for the root
        ['z',num2str(0),' = ', num2str(real(z0)), ' + i(',num2str(imag(z0)),')']
        del = 0.001*abs(z0);
        %del = (zev - z0)/10;
        RelErrFactor = 1;
        RelErrOld = 0;
        evn = 0; % number of times iteration results had decreasing relative errors for three consecutive iterations
        for ni = 1:Nnewt
            %del = (1i)^ni*0.00001;
            
            %fo = Newton(z0,Lh,Xh,Bh,B0,theta_M,theta_P,opts);
            fo = S11_S21(z0,Lh,Xh,Bh,BC,opts)*[1,0]';
            %fr = Newton(z0 + del,Lh,Xh,Bh,B0,theta_M,theta_P,opts);
            fr = S11_S21(z0 + del,Lh,Xh,Bh,BC,opts)*[1,0]';
            %fl = Newton(z0 - del,Lh,Xh,Bh,B0,theta_M,theta_P,opts);
            fp = (fr - fo)/(del);
            z1 = z0 - fo/fp;
            
            RelErrNew = abs((z1-z0)/z0);
            if RelErrNew < RelErrOld
                RelErrFactor = RelErrFactor*2;
            else
                RelErrFactor = 1;
            end
            RelErrOld = RelErrNew;
            if RelErrFactor >= 8 && evn <= 7
                evn = evn + 1;
                TempEV(1,evn) = z1;
                RelErrFactor = 1;
                
            elseif evn >= 8
                break;
            end
            
            if real(z1) >=0
                ['z',num2str(ni),' = ', num2str(real(z1)), ' + i(',num2str(imag(z1)),') with percent relative error = ', num2str(abs((z1-z0)/z0)*100),' and percent error = ', num2str(abs((z1 - zev)/zev)*100)]
            else
                ['Error:', 'z',num2str(ni),' = ', num2str(real(z1)), ' + i(',num2str(imag(z1)),'), found with Re(z) < 0 --> Not a valid eigenvalue']
                z1 = -1;
                break;
            end
            
%                ['z',num2str(ni),' = ', num2str(real(z1)), ' + i(',num2str(imag(z1)),') with percent relative error = ', num2str(abs((z1-z0)/z0)*100),' and percent error = ', num2str(abs((z1 - zev)/zev)*100)]
         
% %             del = abs(z1 - z0);
            z0 = z1;
            
            
            
        end
        CEVList(1,1) = {z0};
        if real(z0) < 0
            CEVList(1,1) = {'Iteration failed: Re(z) < 0'};
        elseif evn == 0
            CEVList(1,1) = {'Iteration failed: Did not converge.'};
        elseif evn == 1
            CEVList(1,1) = {TempEV(1,1)};
            NewtMsg = '';
        elseif evn == 2
            CEVList(1,1) = {(TempEV(1,1) + TempEV(1,2))/2};
            errlist = abs((TempEV(1,1) - TempEV(1,2))/(TempEV(1,1) + TempEV(1,2)))*100;
            CEVList(1,2) = {['Iteration had two convergent sequences. Average reported. Percent relative error = ',num2str(errlist)]};

        elseif evn == 3
            errlist = [2*(TempEV(1,2)-TempEV(1,1))/(TempEV(1,2)+TempEV(1,1)),2*(TempEV(1,3)-TempEV(1,1))/(TempEV(1,3)+TempEV(1,1)),2*(TempEV(1,3)-TempEV(1,2))/(TempEV(1,3)+TempEV(1,2))]
            CEVList(1,1) = {(TempEV(1,1) + TempEV(1,2) + TempEV(1,3))/3};
            CEVList(1,2) = {['Iteration had three convergent sequences. Average reported. Percent relative error = ',num2str(max(abs(errlist)))]};

        end
        ['End of Newton'];
        
    end

    function [CEVList, HalleyResults] = Halley(z,zev,Lh,Xh,Bh,B0,theta_M,theta_P,opts, Nhalley, LEval)
        
        CEVList = cell(1,3);
%         HalleyResultsTmp = cell(1,Nhalley);
%         CEVList(1,1) = 0;
        z0 = z;% the initial guess for the root
        ['z',num2str(0),' = ', num2str(real(z0)), ' + i(',num2str(imag(z0)),')']
        zold = 0.9*z0+0.01*z0*(1-i*1)/sqrt(2);
%         del = (zev - z0)/10;
%         del = 0.00001;
        for ni = 1:Nhalley
            del = (z0 - zold)/150*(1 - i*1)/sqrt(2);
%             fo = S11_S21(z0,Lh,Xh,Bh,B0,theta_M,theta_P,opts)*[1,0]';
            fo = S11_S21(z0,Lh,Xh,Bh,BC,opts)*[1,0]';
%             fr = S11_S21(z0 + del,Lh,Xh,Bh,B0,theta_M,theta_P,opts)*[1,0]';
            fr = S11_S21(z0 + del,Lh,Xh,Bh,BC,opts)*[1,0]';
%             fl = S11_S21(z0 - del,Lh,Xh,Bh,B0,theta_M,theta_P,opts)*[1,0]';
            fl = S11_S21(z0 - del,Lh,Xh,Bh,BC,opts)*[1,0]';
            fp = (fr - fl)/(2*del);
            
            fpp = (fr -2*fo + fl)/(del^2);
            den = (fp - fpp*fo/(2*fp));
            z1 = z0 - fo/(fp - fpp*fo/(2*fp));
            msgHalley = 'Wow, this failed before it even started!';
            if ( real(z1) > LEval(1,1) && real(z1) <  LEval(1,2) ) && (   imag(z1) > LEval(2,1) && imag(z1) <  LEval(2,2)  )
                
                msgHalley = ['z',num2str(ni),' = ', num2str(real(z1)), ' + i(',num2str(imag(z1)),') with percent relative error = ', num2str(abs((z1-z0)/z0)*100),' and percent error = ', num2str(abs((z1 - zev)/zev)*100)];
                msgHalley
                HalleyResults{ni,1} = msgHalley;
            elseif isnan(z1)
                msgHalley = ['Error:', 'z evaluated to NaN. This could be due to a floating point error if the relative error is extremely small.'];
                msgHalley
                HalleyResults{ni,1} = msgHalley;
                z1 = -1;
                break;
            else
                msgHalley = ['Error:', 'z',num2str(ni),' = ', num2str(real(z1)), ' + i(',num2str(imag(z1)),'), is outside of the evaluation window.'];
                msgHalley
                HalleyResults{ni,1} = msgHalley;
                z1 = -1;
                break;
            end
            
            
            zold = z0;
            z0 = z1;         
        end 
        CEVList(1,1) = {z0};
    end

    function [CEVList] = Halley24(z,zev,Lh,Xh,Bh,B0,theta_M,theta_P,opts,Nhalley24)
        
%         CEVList(1,1) = 0;
        CEVList = cell(1,3);
        z0 = z;% the initial guess for the root
        ['z',num2str(0),' = ', num2str(real(z0)), ' + i(',num2str(imag(z0)),')']
%         del = (zev - z0)/10;
        zold = 0.9*z0;
        for ni = 1:Nhalley24
            del = (z0 - zold)/10;
%             foo = S11_S21(z0,Lh,Xh,Bh,B0,theta_M,theta_P,opts)*[1,0]';
%             fro = S11_S21(z0 + del,Lh,Xh,Bh,B0,theta_M,theta_P,opts)*[1,0]';
            foo = S11_S21(z0,Lh,Xh,Bh,BC,opts)*[1,0]';
            fro = S11_S21(z0,Lh,Xh,Bh,BC,opts)*[1,0]';
            fpo = (fro - foo)/(del);
            
            y1 = z0 - foo/fpo;
            
            %if ni == 1
                del = (y1 - z0)/10;
            %end
            
%             fo1 = S11_S21(y1,Lh,Xh,Bh,B0,theta_M,theta_P,opts)*[1,0]';
%             fr1 = S11_S21(y1 + del,Lh,Xh,Bh,B0,theta_M,theta_P,opts)*[1,0]';
            fo1 = S11_S21(z0,Lh,Xh,Bh,BC,opts)*[1,0]';
            fr1 = S11_S21(z0,Lh,Xh,Bh,BC,opts)*[1,0]';
            fp1 = (fr1 - fo1)/(del);
            
            den = 2*foo*fp1^2 - fpo^2*fo1 + fpo*fp1*fo1;
            num = 2*foo*fo1*fp1;

            z1 = y1 - num/den;
            ['z',num2str(ni),' = ', num2str(real(z1)), ' + i(',num2str(imag(z1)),') with percent relative error = ', num2str(abs((z1-z0)/z0)*100),' and percent absolute error = ', num2str(abs((z1 - zev)/zev)*100)]
            zold = z0;
            z0 = z1;         
        end 
        CEVList(1,1) = {z0};
    end

    function [ Sval, SAtime, S11Err, SAErrFlag, TotalSteps, FieldPrepped, Lh, Bh, Xh] = ScatteringDataGrid(NCEval, FieldPrepped, LComplexEval,Lh, Bh, Xh, BC, SMatrixMethod, NFracSS, eps_AS, Rfactor, opts) 
        % ToDo:
        % 1) Check Newton's Method for 'Root Find'. Results do not seem to
        % be consistent with graph of S11 vs. lambda
        %
        % 2) Be sure graph GUI terminates smoothly and passes any complex
        % eigenvalues found to be archived
        %
        % 3) Look into method and requirements to estimate # of eigenvalues
        % within region by numerical evaluation of:
        %
        %    1/(2*pi*i)PathIntegral(f'(z)/f(z)*dz)
        %
        % where the path is closed and counterclockwise around a region
        % with Im(lambda*zeta) > 0.
        % Also look into extension to allow parts of path to cross:
        % --> real axis in region where real eigenvalues can exist and also regions
        % --> real axis beyond point where real eigenvalues can exist
        % --> imaginary axis
        %
        % 4)For efficiency, look into methods and changes needed to allow
        % grid to have varying step sizes. The idea being that this would
        % allow for faster 'integration' through regions where the B-field
        % is, or is expected to be, constant. Perhaps have the user select
        % regions where B-field is constant. 
        %
        % 5) ToDo for the Writer methods: It appears that in splices there
        % may be a very slight discontinuity which ought to be
        % insignificant. It may be that the method used to check that the
        % B-field is resolved finely enough so that the first derivatives
        % are continuous may pick up on this slight splicing mismatch and
        % drive the step size, and thus the number of grid points to
        % unhelpful and numerially inefficient values. So,
        % --> check if this is the case by creating a profile with a very
        % slight splicing discontinuity.
        % --> If this does happen, look into identifying and smoothing such
        % splices during the creation of a spliced profile.
        %
        % 6) [Done]Be sure that xL and xR archived in profile analyses either for
        % real or complex eigenvalues correspond to the original, archived
        % profile and not an adapted grid used in the analysis.
        %
        % [ReEvFinal, Signs1p, NReEvs, Xrange, Bconditions,Xh, Bh, SFlag]
        % S11Err = percent error between S_BigStep & S_1/2_BigStep
        eps_AS_original = eps_AS;
        eps_AS_tmp = eps_AS;
        Sval = zeros(NCEval);
        t = cputime; 
        star1 = '';
        for kk = 1:NCEval
            star1 = [star1,'*'];
        end
        WhatsHappening = 'Evaluating Scattering Data:';
        star1 = [star1,'\r'];
        star2 = '';
        TimeLeft = 'Remaining Time ... ';
        WittyComment = '';
        Taverage = 0;
        
        SFlag = 1;

        lamdas     = zeros(NCEval);
%         S11s       = zeros(NCEval);
        SAtime     = zeros(NCEval);
        S11Err     = zeros(NCEval);
        SAErrFlag  = zeros(NCEval);
        alpha_ratio  = zeros(NCEval);
        TotalSteps = zeros(NCEval);
        Eps_Max = zeros(1,NCEval);
        eps_record = cell(1,NCEval);
        err_record = cell(1,NCEval);
        lambda_Complex = zeros(NCEval);
        
%         for ii = 1:NCEval
        ii = 0;
        Max_eps = 10;
        [SvalTemp, SFlag, Xh0, Bh0,  RelPctErr, ErrStep, StepsTaken, SA_Time, msg_S_terse, msg_S_verbose,SA_ERR] = S_Adaptive_Complete(Xh, Bh, 1, NFracSS, 0.4+1i*0.4, 0.1, Rfactor);

        while ii< NCEval
            ii = ii + 1;
            jj = 0;
            while jj < ( NCEval )
                jj = jj + 1;
%             for jj = 1:NCEval
                l =ComplexEV(LComplexEval, NCEval, ii,jj);
                lambda_Complex(ii,jj) = l;
                
                % For each new eigenvalue it is necessary to re-choose XL and XR:
                [xL_is_alphaXR] = alpha_XR(real(l),imag(l),abs(Bh(1)),abs(Bh(end)));        
                XR = ( Xh(end) - Xh(1))/(1 + xL_is_alphaXR);
                XL = - xL_is_alphaXR*XR;
                if alpha_ratio(ii,jj) == 0
                    alpha_ratio(ii,jj) = xL_is_alphaXR;
                end
                Xh = linspace(XL, XR, length(Bh));
                NumberSteps = 2^floor(log2(length(Bh))/2 + 1);
                
                if SMatrixMethod == 1
                    if ( ii == 1 ) && ( jj == 1 ) && ~FieldPrepped
                        BF_Prep = 1;
                        FieldPrepped = 1;
                    else
                        BF_Prep = 0;
                    end
                    
                    ErrGood = 0;
                    XhTmp = Xh0; % Default to original field resolution
                    BhTmp = Bh0;
                    BfieldResolve = 0;
                    IJVerbose = 0;
                    while ~ErrGood
                        [SvalTemp, SFlag, XhTmp, BhTmp,  RelPctErr, ErrStep, StepsTaken, SA_Time, msg_S_terse, msg_S_verbose,SA_ERR] = S_Adaptive_Complete(XhTmp, BhTmp, BfieldResolve, NFracSS, l, eps_AS_tmp, Rfactor);
                        %[S_Adapt1121, ErrFlag, Xuniform, b2nsymm,  RelPctErr, ErrStep, StepsTaken, S_Adaptive_Time, msg_S_terse, msg_S_verbose, SA_ERR] = S_Adaptive_Complete(Xas, Bas, FieldPrep, NFractionSTEPS, lambda_as, eps_as)
                        Sval(ii,jj) = transpose(SvalTemp)*[1,0]';
                        NumberSteps = 2^floor(log2(length(BhTmp))/2 + 1);
                        Err_ratio = abs(RelPctErr)/(100*abs(eps_AS_original));
                        lamdas(ii,jj) = l;
                        SAtime(ii,jj) = SA_Time;
                        S11Err(ii,jj) = real(RelPctErr);
                        SAErrFlag(ii,jj) = SA_ERR;
                        TotalSteps(ii,jj) = sum(NumberSteps./StepsTaken);
                        ErrString = sprintf('%4.6f',Err_ratio);
                        s11string = sprintf('%4.6e + i(%4.6e)',real(Sval(ii,jj)), imag(Sval(ii,jj)));
                         if Err_ratio > 1
                             if IJVerbose
                                 fprintf(['i = (',num2str(ii),'/',num2str(NCEval),') j = (',num2str(jj),'/',num2str(NCEval),'), S11 = ',s11string,', Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(lamdas(ii,jj)),', Avg # steps =  ', num2str(TotalSteps(ii,jj)),', Number of grid points = ',num2str(length(BhTmp)),'\r']);
                             end
                                                        ErrGood = 0;

                            eps_AS_tmp = eps_AS_tmp/min((Err_ratio),5)*0.25; % Lowers eps by a minimum of 1/4 and a maximum of 1/20
                            
                            if eps_AS_tmp < 10*(-2)*eps_AS_original % This increase B-field resolution if decreasing tolerance is not sufficient
                                eps_AS_tmp = eps_AS_original;
                                Rfactor = Rfactor*0.1;
                                BfieldResolve = 1;
                            end


                        elseif ( Err_ratio > 0.7 ) && ( Err_ratio <= 1.0 )
                            if IJVerbose
                                fprintf(['i = (',num2str(ii),'/',num2str(NCEval),') j = (',num2str(jj),'/',num2str(NCEval),'), S11 = ',s11string,', Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(lamdas(ii,jj)),', Avg # steps =  ', num2str(TotalSteps(ii,jj)),', Number of grid points = ',num2str(length(BhTmp)),'\r']);
                            end
                                eps_AS_tmp = eps_AS_tmp*0.5;
                            ErrGood = 1;
                            BfieldResolve = 0;
                        elseif ( Err_ratio >= 0.5 ) && ( Err_ratio <= 0.7 )
                            if IJVerbose
                                fprintf(['i = (',num2str(ii),'/',num2str(NCEval),') j = (',num2str(jj),'/',num2str(NCEval),'), S11 = ',s11string,', Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(lamdas(ii,jj)),', Avg # steps =  ', num2str(TotalSteps(ii,jj)),', Number of grid points = ',num2str(length(BhTmp)),'\r']);
                            end
                                eps_AS_tmp = eps_AS_tmp*0.75;
                            ErrGood = 1;
                            BfieldResolve = 0;
                        elseif Err_ratio < 0.5
                            if IJVerbose
                                fprintf(['i = (',num2str(ii),'/',num2str(NCEval),') j = (',num2str(jj),'/',num2str(NCEval),'), S11 = ',s11string,', Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(lamdas(ii,jj)),', Avg # steps =  ', num2str(TotalSteps(ii,jj)),', Number of grid points = ',num2str(length(BhTmp)),'\r']);
                            end
                            eps_AS_tmp = eps_AS_tmp*min(2/Err_ratio,3);
                            ErrGood = 1;
                            BfieldResolve = 0;
                        end     
                        
                        
                        
                    end
                    Sval(ii,jj) = transpose(SvalTemp)*[1,0]';
  
                    if BF_Prep == 1
                        NumberSteps = 2^floor(log2(length(Bh))/2 + 1);
                    end
%                     lamdas(ii,jj) = l;
%                     SAtime(ii,jj) = SA_Time;
%                     S11Err(ii,jj) = real(RelPctErr);
%                     SAErrFlag(ii,jj) = SA_ERR;
%                     TotalSteps(ii,jj) = sum(NumberSteps./StepsTaken);
%                     Err_ratio = max(abs(S11Err(ii,jj)))/(100*abs(eps_AS_original));
%                     ErrString = sprintf('%4.6f',Err_ratio);
%                     fprintf(['i = (',num2str(ii), '/',num2str(NCEval),'), j = (',num2str(jj), '/',num2str(NCEval),'), error ratio = ',ErrString,', # steps taken = ',num2str(TotalSteps(ii,jj)),'\r']);
                else
                    [SvalTemp, SFlag] = S11_S21(l,Lh,Xh,Bh,BC, opts);
                    % SvalTemp = [1.1298 - 0.5676i,  -0.0003 - 0.0006i];
                    % SvalTemp*[1,0]' = [1.1298 - 0.5676i,  -0.0003 -0.0006i]*(1) = 1.1298 - 0.5676i
                    %                                                         (0)
                    Sval(ii,jj) = SvalTemp*[1,0]';
                end
                
%                 Sval(ii,jj) = S11_S21(l,Lh,Xh,Bh,BC,opts)*[1,0]';
            end
            clc
            elapsed = cputime - t; t = cputime;
            increase = elapsed/(Taverage + 1.01);
            Taverage = (Taverage*(ii-1) + elapsed)/ii;
            RemainingTime = (NCEval - ii)*Taverage;
            WittyComment = Snappy(RemainingTime,increase);
            TimeTag = [WhatsHappening,'\r',WittyComment,TimeLeft, num2str(RemainingTime),' s'];
            
            Err_ratio = max(abs(S11Err(ii,:)))/(100*abs(eps_AS_original));
            
            eps_list = eps_record{ii};
            eps_list = [eps_list,eps_AS_tmp];
            eps_record{ii} = eps_list;
            
            err_list = err_record{ii};
            err_list = [err_list, Err_ratio];
            err_record{ii} = err_list;  
            
            
            if Err_ratio > 1
%                 star2 = [star2,'*'];
                fprintf([TimeTag, '\r',star1,star2,'\r',star1]);
                fprintf(['For line, ',num2str(ii),'/',num2str(NCEval),' --> Max. err = ',num2str(max(abs(S11Err(ii,:)))),' (Max. err)/tolerance = ',num2str(Err_ratio),' > 1 and this line will be repeated with a smaller tolerance\r']);
                
                ii = ii - 1;
%                 if eps_AS_tmp < Max_eps
%                     Max_eps = .3*Max_eps + 0.7*eps_AS_tmp;
%                 end
%                 Max_eps = eps_AS_tmp;
                if eps_AS_tmp < eps_AS_original
                    eps_AS_tmp = eps_AS_tmp/(Err_ratio)*.5;
                else
                    eps_AS_tmp = eps_AS_original/(Err_ratio);
                end

            elseif ( Err_ratio > 0.7 ) && ( Err_ratio <= 1.0 )
                star2 = [star2,'*'];
                fprintf([TimeTag, '\r',star1,star2,'\r',star1]);
                fprintf(['For line, ',num2str(ii),'/',num2str(NCEval),' --> Max. err = ',num2str(max(abs(S11Err(ii,:)))),' (Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(ComplexEV(LComplexEval, NCEval, ii,1)),' --> ',num2str(ComplexEV(LComplexEval, NCEval, ii,jj)),', Avg # steps =  ', num2str(mean(TotalSteps(ii,:))),'\r']);
%                 if eps_AS_tmp < Max_eps
%                     Max_eps = eps_AS_tmp;
%                 end
                eps_AS_tmp = eps_AS_tmp*0.5;
            elseif ( Err_ratio >= 0.5 ) && ( Err_ratio <= 0.7 )
                star2 = [star2,'*'];
                fprintf([TimeTag, '\r',star1,star2,'\r',star1]);
                fprintf(['For line, ',num2str(ii),'/',num2str(NCEval),' --> Max. err = ',num2str(max(abs(S11Err(ii,:)))),' (Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(ComplexEV(LComplexEval, NCEval, ii,1)),' --> ',num2str(ComplexEV(LComplexEval, NCEval, ii,jj)),', Avg # steps =  ', num2str(mean(TotalSteps(ii,:))),'\r']);                
                eps_AS_tmp = eps_AS_tmp*0.75;
            elseif Err_ratio < 0.5
                star2 = [star2,'*'];
                fprintf([TimeTag, '\r',star1,star2,'\r',star1]);
                fprintf(['For line, ',num2str(ii),'/',num2str(NCEval),' --> Max. err = ',num2str(max(abs(S11Err(ii,:)))),' (Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(ComplexEV(LComplexEval, NCEval, ii,1)),' --> ',num2str(ComplexEV(LComplexEval, NCEval, ii,jj)),', Avg # steps =  ', num2str(mean(TotalSteps(ii,:))),'\r']);                
                eps_AS_tmp = eps_AS_tmp/Err_ratio/2.0;
%                 if eps_AS_tmp > Max_eps
%                     eps_AS_tmp =  Max_eps;
%                 end
                    
%                 if eps_AS_tmp > Max_eps
%                     eps_AS_tmp = Max_eps;
%                 end
                    
            end
            if ii > 0
                Eps_Max(ii) = Max_eps;
            end
            fprintf(['Number of grid points = ',num2str(length(Bh0)),'\r']);
%             Err_ratio = max(abs(S11Err(ii,:)))/(100*eps_AS_original);
%             eps_AS_tmp = eps_AS_original/Err_ratio^2;
%             max(abs(S11Err(ii,:)))  ;    
        end   
        
        eps_one = zeros(1,length(eps_record));
        for epsIdx =1:length(eps_record)
            epsTmp = eps_record{epsIdx};
            eps_one(epsIdx) = epsTmp(1);
        end
        
        err_one = zeros(1,length(err_record));
        MaxNumErr = -1;
        for errIdx =1:length(eps_record)
            errTmp = err_record{errIdx};
            MaxNumtmp = length(errTmp);
            if MaxNumtmp > MaxNumErr
                MaxNumErr = MaxNumtmp;
            end
            err_one(errIdx) = errTmp(1);
        end
        
        EpsFig = figure('units','normalized','outerposition',[(NMonitors -1 +.00) 0.05 0.3 0.5]);%,'title','Error as fraction of Tolerance'); % create figure to display color map               
        plot(eps_one);
        t_EpsFig = title('eps');

        ErrFig = figure('units','normalized','outerposition',[(NMonitors -1 +.7) 0.05 0.3 0.5]);%,'title','Total Number of Steps for Evaluation'); % create figure to display color map               
        plot(err_one);
        t_ErrFig = title('Err Ratio');
        
        err_record
        err_record
    end
                              
    function [CevState,PrntStrng,Xr, Yr, Xg, Yg, Xb, Yb, Xw, Yw] = CevStateMatrix(LComplexEval,NCEval,Sval)
        %Search through values of s11 and determine the 'state' of each value
        %on the grid:
        %   if Re(s11) > 0 and Im(s11) > 0 --> 'PP' and color = red and factor = 2
        %   if Re(s11) > 0 and Im(s11) < 0 --> 'PN' and color = green and factor = 3
        %   if Re(s11) < 0 and Im(s11) > 0 --> 'NP' and color = blue and factor = 5
        %   if Re(s11) < 0 and Im(s11) < 0 --> 'NN' and color = white and factor = 7
        % Output is matrix with 'state' of each eigenvalue
            NCEvs = 0;
            Xr = []; Xg = []; Xb = []; Xw = [];
            Yr = []; Yg = []; Yb = []; Yw = [];
            CevState = zeros(NCEval,NCEval); %'state' represented by 2, 3, 5, 7 as described above
            PP = 2; PN = 3; NP = 5; NN = 7;
            for ii = 1:NCEval
                for jj = 1:NCEval
                    l = ComplexEV(LComplexEval, NCEval, ii,jj);
                    Rl = real(l);Il = imag(l);
                    s11 = Sval(ii,jj);
                    if sign(real(s11)) > 0 && sign(imag(s11)) > 0
                        CevState(ii,jj) = 2;
                        Xr = horzcat(Xr,Rl);
                        Yr = horzcat(Yr,Il);
                    end

                    if sign(real(s11)) > 0 && sign(imag(s11)) < 0
                        CevState(ii,jj) = 3;
                        Xg = horzcat(Xg,Rl);
                        Yg = horzcat(Yg,Il);
                    end

                    if sign(real(s11)) < 0 && sign(imag(s11)) > 0
                        CevState(ii,jj) = 5;
                        Xb = horzcat(Xb,Rl);
                        Yb = horzcat(Yb,Il);
                    end

                    if sign(real(s11)) < 0 && sign(imag(s11)) < 0
                        CevState(ii,jj) = 7;
                        Xw = horzcat(Xw,Rl);
                        Yw = horzcat(Yw,Il);
                    end 
                end
            end 
            % plot(Xr,Yr,'+r',Xg,Yg,'og',Xb,Yb,'xb',Xw,Yw,'^k')
        PrntStrng = 'plot(';
        boolr = ~isempty(Xr) && ~isempty(Yr);% true if there are 'r' points
        boolg = ~isempty(Xg) && ~isempty(Yg);% etc.
        boolb = ~isempty(Xb) && ~isempty(Yb);
        boolw = ~isempty(Xw) && ~isempty(Yw);
        if boolr % if'r' points, add their plot to print string
            PrntStrng = [PrntStrng,['Xr, Yr, '],'''',['+r'],''''];
        end
        if boolg %if 'g' points, add them too
            if boolr %if 'r' plot already there, separate with ','
                PrntStrng = [PrntStrng,','];
            end
            PrntStrng = [PrntStrng,['Xg, Yg, '],'''',['og'],''''];
        end
        if boolb % if 'b' points, add them also
            if boolr || boolg % if either 'r' or 'g' points already there, add ','
                PrntStrng = [PrntStrng,','];
            end
            PrntStrng = [PrntStrng,['Xb, Yb, '],'''',['xb'],''''];
        end
        if boolw % if 'w' points, add them as well
            if boolr || boolg || boolb % if either 'r' or 'g' or 'b' points already there, add ','
                PrntStrng = [PrntStrng,','];
            end
            PrntStrng = [PrntStrng,['Xw, Yw, '],'''',['^k'],''''];
        end
        PrntStrng = [PrntStrng,')']; %close with ')'
        
    end

    function out = Snappy(in, increase)
        out = '\r';
        if in > 60*60
            out = 'To see a world in a grain of sand and heaven in a wild flower Hold infinity in the palms of your hand and eternity in an HOUR!';
        end
        
        if in > 30*60 && in < 60*30
            out = 'Beam me up, Scotty!\r';
        end
        if in < 30*60 && in > 20*60
            out = 'Well, might as well go get some coffee!\r';
        end
        if in < 20*60 && in > 10*60
            out ='Thank you for waiting, an operator will be with you shortly.\r';
        end
        
        if in < 30
            out = 'Wait for it ...\r';
        end
        if in < 10
            out = 'Bazinga!\r';
        end
        if increase > 1
            out = ['I','''','m with Star Fleet. We don','''','t lie!\r'];
        end
    end
       
    function l = ComplexEV(LComplexEval, NCEval, ii,jj)
        l = LComplexEval(1,1) + (LComplexEval(1,2) - LComplexEval(1,1))/(NCEval - 1)*(ii - 1) + ...
                    1i*(LComplexEval(2,1) + (LComplexEval(2,2) - LComplexEval(2,1))/(NCEval - 1)*(jj - 1));
    end

    function [ReEvFinal, Signs1p, NReEvs, Xrange, Bconditions,Xh, Bh, SFlag, S11s, IS21s] = RealEVs(LREval, NEval, NBSIter, Xh, Bh, NFracSS, eps_AS, SMatrixMethod, opts)
        %(LREval, NREval, NBSIter, Lh, Xh, Bh, opts)
        %
        % SMatrixMethod --> 1 for S_Adaptive_Complete
        %               --> 2 for S11_S21
        % NFracSS       --> # of fractional steps to take. Used in call to S_Adaptive_Complete
        % eps_AS        --> as in relative error = 1 + eps_AS, maximum relative
        %                   error for Scattering matrix output of S_Adaptive_Complete
        %Return list of real eigenvalues in range within LREval
        %NEval specifies the resolution
        eps_AS_original = eps_AS;
        eps_AS_tmp = eps_AS;
        Xrange = [Xh(1),Xh(end)];
        Bconditions = [Bh(1), Bh(end)];
        Lh = [-1, 1]*(  Xh(end) - Xh(1)  )/2;
        BC = [Bh(1), Bh(end)];
        NBrkts = 0; %counter for # of pairs of bracket points
        LEvalGrid = linspace(LREval(1),LREval(2),NEval);
        reDel = 0.1*(LREval(2) - LREval(1))/NEval;
        %*****************************
        %S11(l) evalutated on the real axis for 0 < l < Bo must be real valued
        %We will use the Bisection method which relies on knowing the
        %locations, ln, such that s11(ln)*s11(ln+1) < 0.
        %the first loop below sweeps through a fairly fine grid of points
        %and records one value of ln every time the sign is observed to change
        %See: http://www.math.niu.edu/~dattab/MATH435.2013/ROOT_FINDING.pdf
        ReBrkts = zeros(2,NEval);%store bracket points in pairs (a, b)     
        ReBrktsLeftIdx = zeros(1,NEval);%store bracket points in pairs (a, b)
        SFlag = 1;

        lamdas = zeros(1,length(LEvalGrid));
        S11s = zeros(1,length(LEvalGrid));
        IS21s = zeros(1,length(LEvalGrid));
        SAtime = zeros(1,length(LEvalGrid));
        S11Err = zeros(1,length(LEvalGrid));
        SAErrFlag = zeros(1,length(LEvalGrid));
        TotalSteps = zeros(1,length(LEvalGrid));
        
        AdjustXRXLforLambda = 0;
        if AdjustXRXLforLambda
            [xL_is_alphaXR] = alpha_XR(real(LEvalGrid(:,1)),imag(LEvalGrid(:,1)),abs(Bh(1)),abs(Bh(end)));        
            XR = ( Xh(end) - Xh(1))/(1 + xL_is_alphaXR);
            XL = - xL_is_alphaXR*XR;          
        else
            XR = ( Xh(end) - Xh(1))/(1 + 1);
            XL = - 1*XR;
        end
        Nh = length(Bh);
        Xh = linspace(XL, XR, Nh);
        
        if SMatrixMethod == 1
            PrecisionMet = 0;
            IsFieldPrepped = 0;
            Rfactor = 2*10^(-4);
            while ~PrecisionMet
                [S, SFlag, Xh, Bh,  RelPctErr, ErrStep, StepsTaken, SA_Time, msg_S_terse, msg_S_verbose,SA_ERR] = S_Adaptive_Complete(Xh, Bh, ~IsFieldPrepped, NFracSS, LEvalGrid(:,1), eps_AS_tmp, Rfactor);
                IsFieldPrepped = 1;
                Err_ratio = abs(RelPctErr)/(100*abs(eps_AS_original));
                eps_AS_tmp = eps_AS_tmp/Err_ratio^1;
                if Err_ratio < 1
                    PrecisionMet = 1;
                end
            end
            NumberSteps = 2^floor(log2(length(Bh))/2 + 1);
            lamdas(1) = LEvalGrid(:,1);
            S11s(1) = real(S(1));
            IS21s(1) = real(1i*S(2));
            SAtime(1) = SA_Time;
            S11Err(1) = real(RelPctErr);
            SAErrFlag(1) = SA_ERR;
            TotalSteps(1) = sum(NumberSteps./StepsTaken);
        else
            [S, SFlag] = S11_S21(LEvalGrid(:,1),Lh,Xh,Bh,BC, opts);
            S11s(1) = real(S(1));
            IS21s(1) = real(1i*S(2));
        end
        
        if SFlag < 1
            ['Error in evaluation of Scattering Data: Psiminus']
            return
        end
        s11 = S(1);
        signold = sign(real(s11));
%         signold
        Sdebug = zeros(2,200);
%         ReEvFinal = [];
%         Signs1p = [];
        % ReEvFinal(1) = {'No real eigenvalues found'};%default case valid if NBrkts = 0
        ReEvFinal = [];%default case valid if NBrkts = 0
        Signs1p = 0;
        NReEvs = 0;
        
%         for ii = 2:NEval
        ii = 1;
        while ii < NEval
            ii = ii + 1;
%             [S, SFlag] = S11_S21(LEvalGrid(:,ii),Lh,Xh, Bh, BC,opts);
            if AdjustXRXLforLambda
                [xL_is_alphaXR] = alpha_XR(real(LEvalGrid(:,1)),imag(LEvalGrid(:,1)),abs(Bh(1)),abs(Bh(end)));        
                XR = ( Xh(end) - Xh(1))/(1 + xL_is_alphaXR);
                XL = - xL_is_alphaXR*XR;          
            else
                XR = ( Xh(end) - Xh(1))/(1 + 1);
                XL = - 1*XR;
            end
            Nh = length(Bh);
            Xh = linspace(XL, XR, Nh);
            if SMatrixMethod == 1
                [S           , SFlag , ~       , ~      ,  RelPctErr, ErrStep, StepsTaken, SA_Time        , msg_S_terse, msg_S_verbose,SA_ERR] = S_Adaptive_Complete(Xh, Bh, 0, NFracSS, LEvalGrid(:,ii), eps_AS_tmp, Rfactor);
  %             [S_Adapt1121, ErrFlag, Xuniform, b2nsymm,  RelPctErr, ErrStep, StepsTaken, S_Adaptive_Time, msg_S_terse, msg_S_verbose, SA_ERR]

                lamdas(ii) = LEvalGrid(:,ii);
                S11s(ii) = real(S(1));
                IS21s(ii) = real(1i*S(2));
                SAtime(ii) = SA_Time;
                S11Err(ii) = abs(RelPctErr);
                SAErrFlag(ii) = SA_ERR;
                TotalSteps(ii) = sum(NumberSteps./StepsTaken);
            else
                lamdas(ii) = LEvalGrid(:,ii);
                [S, SFlag] = S11_S21(LEvalGrid(:,ii),Lh,Xh, Bh, BC,opts);
                S11s(ii) = real(S(1));
                IS21s(ii) = real(1i*S(2));
            end
            if SFlag < 1
                ['Error in evaluation of Scattering Data: Psiminus']
                % ReEvFinal(1) = {'Error in evaluation of Scattering Data: Psiminus.'};%default case valid if NBrkts = 0
                ReEvFinal(1) = [];%default case valid if NBrkts = 0
                Signs1p = 0;
                NReEvs = 0;
                return
            end
            s11 = S(1);
            Sdebug(:,ii) = transpose(S);
            signnew = sign(real(s11));
            s11string = sprintf('%4.4E',real(S(1)));
%             signnew
            if SMatrixMethod == 1
                Err_ratio = abs(RelPctErr)/(100*abs(eps_AS_original));
            
                if Err_ratio > 1
                    fprintf(['For line, ',num2str(ii),'/',num2str(NEval),', Re(S11) = ',s11string,' with (Max. err)/tolerance = ',num2str(abs(RelPctErr)),'/',num2str((100*abs(eps_AS_original))),' --> (Max. err)/tolerance = ',num2str(Err_ratio),' > 1 and this line will be repeated with a smaller tolerance\r']);                
                    ii = ii - 1;

                    eps_AS_tmp = eps_AS_tmp/(Err_ratio)*.25;
    %                 if eps_AS_original > eps_AS_tmp
    %                     eps_AS_tmp = eps_AS_tmp/(Err_ratio)*.5;
    %                 else
    %                     eps_AS_tmp = eps_AS_original/(Err_ratio)*.92653589793;
    %                 end

                elseif ( Err_ratio > 0.7 ) && ( Err_ratio <= 1.0 )
                    fprintf(['For line, ',num2str(ii),'/',num2str(NEval),', Re(S11) = ',s11string,' with (Max. err)/tolerance = ',num2str(abs(RelPctErr)),'/',num2str((100*abs(eps_AS_original))),' --> (Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(LEvalGrid(:,ii)),', Avg # steps =  ', num2str(mean(TotalSteps(ii))),', Number of grid points = ',num2str(length(Bh)),'\r']);
                    eps_AS_tmp = eps_AS_tmp*0.5;
                elseif ( Err_ratio >= 0.5 ) && ( Err_ratio <= 0.7 )
                    fprintf(['For line, ',num2str(ii),'/',num2str(NEval),', Re(S11) = ',s11string,' with (Max. err)/tolerance = ',num2str(abs(RelPctErr)),'/',num2str((100*abs(eps_AS_original))),' --> (Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(LEvalGrid(:,ii)),', Avg # steps =  ', num2str(mean(TotalSteps(ii))),', Number of grid points = ',num2str(length(Bh)),'\r']);
                    eps_AS_tmp = eps_AS_tmp*0.75;
                elseif Err_ratio < 0.5
                    fprintf(['For line, ',num2str(ii),'/',num2str(NEval),', Re(S11) = ',s11string,' with (Max. err)/tolerance = ',num2str(abs(RelPctErr)),'/',num2str((100*abs(eps_AS_original))),' --> (Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(LEvalGrid(:,ii)),', Avg # steps =  ', num2str(mean(TotalSteps(ii))),', Number of grid points = ',num2str(length(Bh)),'\r']);
                    eps_AS_tmp = eps_AS_tmp*min(2/Err_ratio,3);

                end
            
                if Err_ratio < 1
                    if signnew ~= signold
                        NBrkts = NBrkts + 1;
                        ReBrkts(1,NBrkts) = LEvalGrid(:,ii-1);%record la of this bracket pair   
                        ReBrkts(2,NBrkts) = LEvalGrid(:,ii);%record lb of this bracket pair 
                        ReBrktsLeftIdx(NBrkts) = ii; % Record the index of LEvalGrid just to the left of a sign change in S11
                        signold = signnew; 
                    end
                end
            else
                if signnew ~= signold
                    NBrkts = NBrkts + 1;
                    ReBrkts(1,NBrkts) = LEvalGrid(:,ii-1);%record la of this bracket pair   
                    ReBrkts(2,NBrkts) = LEvalGrid(:,ii);%record lb of this bracket pair 
                    ReBrktsLeftIdx(NBrkts) = ii; % Record the index of LEvalGrid just to the left of a sign change in S11
                    signold = signnew; 
                end
                
            end
                      
%             fprintf(['Number of grid points = ',num2str(length(Bh)),'\r']);

            
        end
        SFlag = sum(SAErrFlag);
%          ReBrkts
%          NBrkts
%   
        % If there are real eigenvalues, NBrkts > 0, then proceed to
        % determine their values and the signs of the corresponding solitons
        if NBrkts > 0
            %% First Attempt the Adaptive Step Size Method to refine eigenvalues and determine signs
            AdaptiveStepSizeMethod = 0; % Set AdaptiveStepSizeMethod to 1 to default to this method
            OldBracketMethod = 1;       % Upon execution of AdaptiveStepSizeMethod, if it fails, it will turn on OldBracketMethod
            NewBracketMethod = 0;       % If sign( real(iS21/S11') ) is ambiguous, this method will be used
            if AdaptiveStepSizeMethod
                t = cputime;
                lambdaGrid = linspace(LEvalGrid(:,1), LEvalGrid(:,end),NEval);
                lambdaGridFine = linspace(LEvalGrid(:,1), LEvalGrid(:,end),50*NEval);
                S11Fine = interp1(lambdaGrid,S11s,lambdaGridFine,'spline');
                IS21Fine = interp1(lambdaGrid,IS21s,lambdaGridFine,'spline');
                S11FineShiftOne = circshift(S11Fine,-1);
                IS21FineShiftOne = circshift(IS21Fine,-1);
                S11Zeros = find((S11Fine(2:end-1).*S11FineShiftOne(2:end-1))<0);
                IS21Zeros = find((IS21Fine(2:end-1).*IS21FineShiftOne(2:end-1))<0);
                NReEvs = 0;
                UpResolution = 0; % Flag set to 1 if a zero of S11
                if ~isempty(S11Zeros) 
                    S11Zeros = S11Zeros + 1;   % indices of zeros of S11
                    IS21Zeros = IS21Zeros + 1; % indices of zeros of IS21
                    NReEvs = length(S11Zeros);
                end
                
%                 if NReEvs > 0
%                     for ResIdx = 1:NReEvs
%                         if min(S11Zeros(ResIdx) - IS21Zeros) < 2
%                         end
%                     end
%                 end
                
                S11Slope = (S11FineShiftOne - S11Fine)/(lambdaGridFine(2) - lambdaGridFine(1));
%                 if NBrkts && NBrkts ==
%                     NReEvs = length(S11Zeros);
%                 else
%                     NReEvs = 0;
%                 end
%                 ReEvFinal = [];
%                 Signs1p = [];
                WholeLineCurveFitTime = cputime - t;
                WholeLineCurveFit = 0;
                if WholeLineCurveFit
                    % If the finer resolution, S11Fine, is a faithful
                    % representation of S11s, then the number and location
                    % of places it changes sign should coincide with the
                    % number and location of the sign changes found using
                    % the bracket method.
                    ReEvFinal = lambdaGridFine(S11Zeros);
                    Signs1p = IS21Fine(S11Zeros)./S11Slope(S11Zeros);
                else
                    % If NBrkts ~= NReEvs, it sometimes happens that the
                    % spline interpolation method forces rapid
                    % oscillations, and false zeros, in sections where S11s
                    % is nearly flat and zero though perhaps without actually 
                    % crossing zero. This has been observed to happen in the LHS 
                    % of the spectrum, creating spurious sign changes if S11s 
                    % has a large amplitude oscillation on the RHS of the spectrum.
                    %
                    % In this case, evaluate S11Fine only across regions
                    % nearby the bounds recorded in ReBrkts(1,ev) and
                    % ReBrkts(2, ev) for ev = 1, 2, 3, ... NBrkts.
                    t = cputime;
                    ReEvFinal = zeros(1,NBrkts);
                    Signs1p   = zeros(1,NBrkts);
                    DelLambda2111 = (-1)*ones(1,NBrkts);
                    NReEvs = 0;
%                         ev = 1;
                    for ev = 1:NBrkts    

                        % Establish values for lambdaIdxLeft and lambdaIdxRight
                        % to use for finer resolution curve fit for S11 and S21
                        % in vicinity of a bracket point (that is, a location where
                        % the sign of S11 changes).
                        IdxStepLeftTry  = 2;
                        IdxStepRightTry = 2;
                        LeftIdxTest = ReBrktsLeftIdx(ev) - IdxStepLeftTry;
                        RightIdxTest = (ReBrktsLeftIdx(ev) + IdxStepRightTry);
                        if ( LeftIdxTest > 0 )
                            lambdaIdxLeft = LeftIdxTest;
                        else
                            lambdaIdxLeft = 1;
                        end
                        if ( RightIdxTest < (NEval + 1) )
                            lambdaIdxRight = RightIdxTest;
                        else
                            lambdaIdxRight = NEval;
                        end

                        % Define section of lambda spectrum adjacent to a sign
                        % change recorded in ReBrktsLeftIdx
                        lambdaGridBrkt = linspace(lambdaGrid(lambdaIdxLeft), lambdaGrid(lambdaIdxRight),( lambdaIdxRight - lambdaIdxLeft + 1 ));

                        % Define a corresponding section of the S11s values
                        S11sBrkt = S11s(lambdaIdxLeft:lambdaIdxRight);

                        % Define a corresponding section of the IS21s values
                        IS21sBrkt = IS21s(lambdaIdxLeft:lambdaIdxRight);

                        % Define a finer resolution spectrum for the same section
                        lambdaGridFineBrkt = linspace(lambdaGrid(lambdaIdxLeft), lambdaGrid(lambdaIdxRight),50*( lambdaIdxRight - lambdaIdxLeft + 1 ));

                        % Create a fine scale interpolation, S11FineBrkt, of
                        % S11s for section surrounding a sign change
                        S11FineBrkt = interp1(lambdaGridBrkt,S11sBrkt,lambdaGridFineBrkt,'spline');

                        % Create a fine scale interpolation, IS21FineBrkt, of
                        % IS21s for section surrounding a sign change
                        IS21FineBrkt = interp1(lambdaGridBrkt,IS21sBrkt,lambdaGridFineBrkt,'spline');
                        IS21FineBrktShiftOne = circshift(IS21FineBrkt,-1);
                        IS21BrktZeros = find((IS21FineBrkt(2:end-1).*IS21FineBrktShiftOne(2:end-1))<0) + 1;


                        S11FineBrktShiftOne = circshift(S11FineBrkt,-1);
                        S11BrktZeros = find((S11FineBrkt(2:end-1).*S11FineBrktShiftOne(2:end-1))<0) + 1;
                        deltalambda = (lambdaGridFineBrkt(2) - lambdaGridFineBrkt(1));
                        S11BrktSlope = (S11FineBrktShiftOne - S11FineBrkt)/deltalambda;
                        
                        if ~isempty(IS21BrktZeros) && ~isempty(S11BrktZeros)
                            DelLambda2111(ev) = min(abs(IS21BrktZeros - S11BrktZeros));
                        elseif isempty(IS21BrktZeros) && ~isempty(S11BrktZeros)
                            DelLambda2111(ev) = min(abs(S11BrktZeros));
                        end
                        
                        
                        if ~isempty(S11BrktZeros)
%                             S11BrktSlope = (S11FineBrktShiftOne - S11FineBrkt)/(lambdaGridFineBrkt(2) - lambdaGridFineBrkt(1));
                            ReEvFinal(ev) = lambdaGridFineBrkt(S11BrktZeros(1));
                            if DelLambda2111(ev) > 1
                                Signs1p(ev) = IS21FineBrkt(S11BrktZeros(1))./S11BrktSlope(S11BrktZeros(1));
                            else
                                Signs1p(ev) = 0;
                            end
                            NReEvs = NReEvs + 1;
                        else
                            OldBracketMethod = 1; % Go to old bracket method if this modified curve fitting approach fails
                        end
                    end
                    BrktCurveFitTime = cputime - t;
                    
                end
                
                CurveFitCWOutput = 1;
                if CurveFitCWOutput
                    figure;
                    plot(lamdas,S11s/max(abs(S11s)),'kx');
                    hold 'on';
                    plot(lambdaGridFine,S11Fine/max(abs(S11Fine)),'color',[ 0 0 0]);
                    hold 'on';
                    plot(lambdaGridFine,IS21Fine/max(abs(IS21Fine)),'color',[ 0 0 0]);
                    hold 'on';
                    plot(lamdas,S11Err/max(abs(S11Err)),'color',[ 1 0 0]);
                    hold 'on';
                    plot(lamdas,SAtime/max(SAtime),'color',[ 0 0 1]);




                    [MaxErr, MaxErrIdx] = max(abs(S11Err));

                    Eintegrate = 2*Energy_Calc(Lh,Xh,Bh,Bh(1));
                    E_solitons = zeros(1,NReEvs);
                    Helicity_solitons = zeros(1,NReEvs);
                    for iREidx = 1: NReEvs
                            E_solitons(iREidx) = dispersion/alpha*(2*pi()*(1 + sign(Signs1p(iREidx)))-4*acos(lambdaGridFine(S11Zeros(iREidx))/abs(Bh(1))));
                            Helicity_solitons(iREidx) = -2*(dispersion/alpha)^2*sqrt(abs(Bh(end))^2 -ReEvFinal(iREidx)^2 )/(abs(Bh(end))^2*ReEvFinal(iREidx));
                    end
        %             ratio_energy = E_solitons/Eintegrate;
                    Helicity = Helicity_Calc(Xh,Bh);
                    fprintf('**************************************************************************************\n*\n');
                    fprintf(['*   From RealEVs: After For Loop evaluating s11/s21 from ',char(955),'_min = ',num2str(lamdas(1)),' to ',char(955),'_max = ',num2str(lamdas(end)),'\n*\n']);
                    fprintf(['*   S11 max error = ',num2str(MaxErr(1)),' %% at ',char(955),' = ', num2str(lamdas(MaxErrIdx(1))),'\n*\n']);

                    if NReEvs > 0
                        fprintf('*   Eigenvalues: \n');

                        for iREidx = 1: NReEvs
                            evString = sprintf('%2.5f',ReEvFinal(iREidx));
                            if DelLambda2111(iREidx) < 1
                                fprintf(['*                 ',char(955)','(',num2str(iREidx),') = ',evString,', sign is indeterminant at resolution of delta lambda = ',num2str(deltalambda),'\n']);
                            else
                                fprintf(['*                 ',char(955)','(',num2str(iREidx),') = ',evString,', sign = ',num2str(sign(Signs1p(iREidx))),'\n']);
                            end
                                
                        end


                        if ( abs(Bh(1)/Bh(end)) < 1.03 ) && ( abs(Bh(1)/Bh(end)) > .97 ) % if 'energy can be evaluated

                            % Display Energy Analysis
                            fprintf('*\n');
                            fprintf('*   Energy Analysis: \n');
                            for iREidx = 1: NReEvs
                                fprintf(['*                 ',char(955)','(',num2str(iREidx),') = ',num2str(ReEvFinal(iREidx)),', sign = ',num2str(sign(Signs1p(iREidx))),', E(',num2str(ReEvFinal(iREidx)),') = ',num2str(E_solitons(iREidx)),' ... ',num2str(E_solitons(iREidx)/Eintegrate*100),' %% of total energy of ',num2str(Eintegrate),'\n']);
                            end
                            ratio_energy = sum(E_solitons)/Eintegrate;
                            if ( ratio_energy < 1.1 ) && ( ratio_energy > 0.9 )
                                fprintf(['*   Total soliton energy = ',num2str(sum(E_solitons)),' is ',num2str(sum(E_solitons)/Eintegrate*100),' %% of total energy of ',num2str(Eintegrate),'\n*\n']);
                            else
                                fprintf(2,['*   Total soliton energy = ',num2str(sum(E_solitons)),' is ',num2str(sum(E_solitons)/Eintegrate*100),' %% of total energy of ',num2str(Eintegrate),'\n*\n']);
                            end

                            % Display Helicity Analysis
                            fprintf('*\n');
                            fprintf('*   Helicity Analysis:  \n');
                            for iREidx = 1: NReEvs
                                fprintf(['*                 ',char(955)','(',num2str(iREidx),') = ',num2str(ReEvFinal(iREidx)),', helicity = ',num2str(Helicity_solitons(iREidx)),' ... ',num2str(Helicity_solitons(iREidx)/Helicity*100),' %% of total helicity of ',num2str(Helicity),'\n']);
                            end
                            ratio_helicity = sum(Helicity_solitons)/Helicity;
                            if ( ratio_helicity < 1.1 ) && ( ratio_helicity > 0.9 )
                                fprintf(['*   Total soliton helicity = ',num2str(sum(Helicity_solitons)),' is ',num2str(ratio_helicity*100),' %% of total helicity of ',num2str(Helicity),'\n*\n']);
                            else
                                fprintf(2,['*   Total soliton helicity = ',num2str(sum(Helicity_solitons)),' is ',num2str(ratio_helicity*100),' %% of total helicity of ',num2str(Helicity),'\n*\n']);

                            end

                        end


                    else
                        fprintf('*   No real eigenvalues found. \n');
                    end
                    fprintf('**************************************************************************************\n\n');
                end
            end


                if NewBracketMethod
                    %% Finalize the eigenvalue brackets found above
                    if NBrkts > 0 %Pair down the array to match the actual # of points
                        ReBrktsfinal = zeros(2,NBrkts);
                        for ii = 1:NBrkts
                            ReBrktsfinal(:,ii) = ReBrkts(:,ii);
                        end   
                    end

%                     ReEvFinal(1) = {'No real eigenvalues found.'};%default case valid if NBrkts = 0


                    %% Background on Bisection Method used to refine real eigenvalues
                    % Apply the Bisection method to the bracketting l-values recorded in 'ReEvlist'
                    %to improve precision. 
                    % The Bisection Method:
                    % Suppose la and lb are listed sequentially in ReEvlist, then if
                    % fa = sign(real(S11(la)) and fb = sign(real(S11(lb)), it must be
                    % that fa*fb < 0
                    %Take c = (a+b)/2, then there are two cases ...
                    % 1) if fa*fc < 0, then the zero must be between a and c, so narrow 
                    %the bracket with b = c and repeat.
                    % 2) else fa and fc have the same sign and the zero must be between
                    % c and b, then update a = c
                    % Iterate a fixed number of times or until targeted precision is
                    % reached. The method will work if there is just a single zero
                    % between each pair of bracket points. So it is possible it could
                    % be blind to multiple zeroes ... and so overlook zeros.

                    %% Refine Real Eigenvalues, if any, and evaluate signs of corresponding solitons
                    NReEvs = NBrkts;
%                     Signs1p = 0;
                    if NBrkts > 0 % if at least one pair of bracketting points is found
                                  % then refine eigenvalues and evaluate their signs

                        ReEvFinal = zeros(1,NBrkts); % the # real EV's = # of pairs of bracketting points  
%                         Signs1p = zeros(1,NBrkts);

                        % Evaluate final real eigenvalues using the Bisection Method
%                         Rfactor = 2*10^(-4);
                        RfactorOriginal = 2*10^(-4);                        
                        for ii = 1:NBrkts
                            if Signs1p(ii) == 0
                                a = ReBrktsfinal(1,ii);
                                b = ReBrktsfinal(2,ii);            
                %                 Sa = S11_S21(a,Lh,Xh,Bh,BC,opts);
                                if AdjustXRXLforLambda
                                    [xL_is_alphaXR] = alpha_XR(real(LEvalGrid(:,1)),imag(LEvalGrid(:,1)),abs(Bh(1)),abs(Bh(end)));        
                                    XR = ( Xh(end) - Xh(1))/(1 + xL_is_alphaXR);
                                    XL = - xL_is_alphaXR*XR;          
                                else
                                    XR = ( Xh(end) - Xh(1))/(1 + 1);
                                    XL = - 1*XR;
                                end
                                Nh = length(Bh);
                                Xh = linspace(XL, XR, Nh);
%                                 Rfactor = 2*10^(-4);
%                                 RfactorOrigninal = 2*10^(-4);
                                XhTmp = Xh;
                                BhTmp = Bh;
                                BfieldResolve = 0;
                                if SMatrixMethod == 1
                                    ErrGood = 0;
                                    XhTmp = Xh;
                                    BhTmp = Bh;
                                    BfieldResolve = 0;
                                    while ~ErrGood
                                        [Sa, SFlag, XhTmp, BhTmp,  RelPctErr, ~, ~, ~, ~,~,~] = S_Adaptive_Complete(XhTmp, BhTmp, BfieldResolve, NFracSS, a, eps_AS, Rfactor);
                                        %[S_Adapt1121, ErrFlag, Xuniform, b2nsymm,  RelPctErr, ErrStep, StepsTaken, S_Adaptive_Time, msg_S_terse, msg_S_verbose, SA_ERR]
                                        Err_ratio = abs(RelPctErr)/(100*abs(eps_AS_original));
                                         if Err_ratio > 1
                                            fprintf(['Just before bracket iteration: Re(S11) = ',s11string,' with (Max. err)/tolerance = ',num2str(abs(RelPctErr)),'/',num2str((100*abs(eps_AS_original))),' --> (Max. err)/tolerance = ',num2str(Err_ratio),' > 1 and this will be repeated with a smaller tolerance\r']);                

                                            ErrGood = 0;

                                            eps_AS_tmp = eps_AS_tmp/min((Err_ratio),5)*0.25;
                                            if eps_AS_tmp < 10*(-2)*eps_AS_original
                                                eps_AS_tmp = eps_AS_original;
                                                Rfactor = Rfactor*0.1;
                                                BfieldResolve = 1;
                                            end


                                        elseif ( Err_ratio > 0.7 ) && ( Err_ratio <= 1.0 )
                                            fprintf(['Just before bracket iteration: Re(S11) =  ',s11string,' with (Max. err)/tolerance = ',num2str(abs(RelPctErr)),'/',num2str((100*abs(eps_AS_original))),' --> (Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(a),', Avg # steps =  ', num2str(mean(TotalSteps(ii))),', Number of grid points = ',num2str(length(Bh)),'\r']);
                                            eps_AS_tmp = eps_AS_tmp*0.5;
                                            ErrGood = 1;
                                        elseif ( Err_ratio >= 0.5 ) && ( Err_ratio <= 0.7 )
                                            fprintf(['Just before bracket iteration: Re(S11) = ',s11string,' with (Max. err)/tolerance = ',num2str(abs(RelPctErr)),'/',num2str((100*abs(eps_AS_original))),' --> (Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(a),', Avg # steps =  ', num2str(mean(TotalSteps(ii))),', Number of grid points = ',num2str(length(Bh)),'\r']);
                                            eps_AS_tmp = eps_AS_tmp*0.75;
                                            ErrGood = 1;
                                        elseif Err_ratio < 0.5
                                            fprintf(['Just before bracket iteration: Re(S11) = ',s11string,' with (Max. err)/tolerance = ',num2str(abs(RelPctErr)),'/',num2str((100*abs(eps_AS_original))),' --> (Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(a),', Avg # steps =  ', num2str(mean(TotalSteps(ii))),', Number of grid points = ',num2str(length(Bh)),'\r']);
                                            eps_AS_tmp = eps_AS_tmp*min(2/Err_ratio,3);
                                            ErrGood = 1;
                                        end       
                                    end


                                else
                                    [Sa, SFlag] = S11_S21(a,Lh,Xh,Bh,BC,opts);
                                end
                                fa = sign(real(Sa(1)));
                                ga = sign(real(1i*Sa(2)));
                                SignVerify = 0;
                                k = 0;
                                RepeatK = 0;
                                Rfactor = RfactorOriginal;
%                                 RfactorOrigninal = 2*10^(-4);
                                XhTmp = Xh;
                                BhTmp = Bh;
                                BfieldResolve = 0;
                                Err_ratio_old = -1;
                                c_old = -1;
                                while k < NBSIter && ~SignVerify%actual iteration step of Bisection Method
                                    k = k+1;

                                    c = (a+b)/2; %1/2 way between points bracketing a sign change
                %                     Sc = S11_S21(c,Lh,Xh,Bh,BC,opts);
                                    [xL_is_alphaXR] = alpha_XR(real(c),imag(c),abs(Bh(1)),abs(Bh(end)));        
                                    XR = ( Xh(end) - Xh(1))/(1 + xL_is_alphaXR);
                                    XL = - xL_is_alphaXR*XR;
                                    Nh = length(Bh);
                                    Xh = linspace(XL, XR, Nh);
                                    if SMatrixMethod == 1
                                        [Sc, SFlag,  XhTmp, BhTmp,  RelPctErr, ~, ~, ~, ~,~,~] = S_Adaptive_Complete(XhTmp, BhTmp, BfieldResolve, NFracSS, c, eps_AS, Rfactor);

                                        %[S_Adapt1121, ErrFlag, Xuniform, b2nsymm,  RelPctErr, ErrStep, StepsTaken, S_Adaptive_Time, msg_S_terse, msg_S_verbose, SA_ERR] =
                                         % *******************************************************
                                         % *    Error control for Bracket Method
                                         % *******************************************************

                                        Err_ratio = abs(RelPctErr)/(100*abs(eps_AS_original));
                                        s11string = sprintf('%4.4E',real(Sc(1)));

                                        if Err_ratio > 1
                                            fprintf(['For k = , ',num2str(k),'/',num2str(NBSIter),', Re(S11) = ',s11string,' with (Max. err)/tolerance = ',num2str(abs(RelPctErr)),'/',num2str((100*abs(eps_AS_original))),' --> (Max. err)/tolerance = ',num2str(Err_ratio),' > 1 and this will be repeated with a smaller tolerance\r']);                
                                            k = k - 1; 
                                            RepeatK = 1;


                                            eps_AS_tmp = eps_AS_tmp/min((Err_ratio),5)*0.25;
                                            if (eps_AS_tmp < (10*(-2)*eps_AS_original) ) || ( abs(1 - Err_ratio/Err_ratio_old) < 0.05 )
                                                eps_AS_tmp = eps_AS_original;
                                                Rfactor = Rfactor*0.1;
                                                BfieldResolve = 1;
                                            end

                                            Err_ratio_old = Err_ratio;
                            %                 if eps_AS_original > eps_AS_tmp
                            %                     eps_AS_tmp = eps_AS_tmp/(Err_ratio)*.5;
                            %                 else
                            %                     eps_AS_tmp = eps_AS_original/(Err_ratio)*.92653589793;
                            %                 end

                                        elseif ( Err_ratio > 0.7 ) && ( Err_ratio <= 1.0 )
                                            fprintf(['For k = , ',num2str(k),'/',num2str(NBSIter),', Re(S11) = ',s11string,' with (Max. err)/tolerance = ',num2str(abs(RelPctErr)),'/',num2str((100*abs(eps_AS_original))),' --> (Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(c),', Avg # steps =  ', num2str(mean(TotalSteps(ii))),', Number of grid points = ',num2str(length(Bh)),'\r']);
                                            eps_AS_tmp = eps_AS_tmp*0.5;
                                            BfieldResolve = 0;
                                            RepeatK = 0;
                                            Err_ratio_old = -1;
                                        elseif ( Err_ratio >= 0.5 ) && ( Err_ratio <= 0.7 )
                                            fprintf(['For k = , ',num2str(k),'/',num2str(NBSIter),', Re(S11) = ',s11string,' with (Max. err)/tolerance = ',num2str(abs(RelPctErr)),'/',num2str((100*abs(eps_AS_original))),' --> (Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(c),', Avg # steps =  ', num2str(mean(TotalSteps(ii))),', Number of grid points = ',num2str(length(Bh)),'\r']);
                                            eps_AS_tmp = eps_AS_tmp*0.75;
                                            BfieldResolve = 0;
                                            RepeatK = 0;
                                            Err_ratio_old = -1;
                                        elseif Err_ratio < 0.5
                                            fprintf(['For k = , ',num2str(k),'/',num2str(NBSIter),', Re(S11) = ',s11string,' with (Max. err)/tolerance = ',num2str(abs(RelPctErr)),'/',num2str((100*abs(eps_AS_original))),' --> (Max. err)/tolerance = ',num2str(Err_ratio),', lambda = ',num2str(c),', Avg # steps =  ', num2str(mean(TotalSteps(ii))),', Number of grid points = ',num2str(length(Bh)),'\r']);
                                            eps_AS_tmp = eps_AS_tmp*min(2/Err_ratio,3);
                                            BfieldResolve = 0;
                                            RepeatK = 0;
                                            Err_ratio_old = -1;
                                        end       


                                    else
                                        [Sc, SFlag] = S11_S21(c,Lh,Xh,Bh,BC,opts);
                                    end
                                    c_old = c;
                                    fc = sign(real(Sc(1)));
                                    gc = sign(real(1i*Sc(2)));
                                    if fc*fa < 0 && ~RepeatK
                                        b = c; % if so fc must have the same sign as b, so narrow the bracket
                                        if ga*gc > 0
                                            SignVerify = 1;
                                        else
                                            SignVerify = 0;
                                        end
                                    else
                                        if ~RepeatK 
                                            a = c; % if not, fc must have the same sign as a, so narrow the bracket

                                        end
                                    end
                                end
                                ReEvFinal(ii) = real(a+b)/2;% return the 1/2 way point of the last bracket
    %                             Signs1p(ev) = IS21FineBrkt(S11BrktZeros(1))./S11BrktSlope(S11BrktZeros(1));
                                ['Real Eigenvalue calculation: at ', num2str(ii), ' of ', num2str(NBrkts)]
                            end
                        end 

                        % Evaluate Signs for final real eigenvalues
                        for ii = 1:NBrkts
                            if Signs1p(ii) == 0
                                l1 = ReEvFinal(ii);
                %                 Dl = 0.0001;
                %                 s21 = 1i*imag(S11_S21(l1,Lh,Xh,Bh,BC,opts))*[0,1]';
                %                 s11p = real(S11_S21(l1 + reDel,Lh,Xh,Bh,BC,opts))*[1,0]';
                %                 s11m = real(S11_S21(l1 - reDel,Lh,Xh,Bh,BC,opts))*[1,0]';
                                if AdjustXRXLforLambda
                                    [xL_is_alphaXR] = alpha_XR(real(LEvalGrid(:,1)),imag(LEvalGrid(:,1)),abs(Bh(1)),abs(Bh(end)));        
                                    XR = ( Xh(end) - Xh(1))/(1 + xL_is_alphaXR);
                                    XL = - xL_is_alphaXR*XR;          
                                else
                                    XR = ( Xh(end) - Xh(1))/(1 + 1);
                                    XL = - 1*XR;
                                end
                                Nh = length(Bh);
                                Xh = linspace(XL, XR, Nh);
                                if SMatrixMethod == 1
                                    [S21, SFlag, ~, ~,  ~, ~, ~, ~, ~,~,~] = S_Adaptive_Complete(Xh, Bh, 0, NFracSS, l1, eps_AS, Rfactor);
                                    s21 = 1i*imag(transpose(S21))*[0,1]';
                                    [S11p, SFlag, ~, ~,  ~, ~, ~, ~, ~,~,~] = S_Adaptive_Complete(Xh, Bh, 0, NFracSS, l1 + reDel, eps_AS, Rfactor);
                                    s11p = real(transpose(S11p))*[1,0]';
                                    [S11m, SFlag,  ~, ~,  ~, ~, ~, ~, ~,~,~] = S_Adaptive_Complete(Xh, Bh, 0, NFracSS, l1 - reDel, eps_AS, Rfactor);
                                    s11m = real(transpose(S11m))*[1,0]';
                                else
                                    [S21, SFlag] = S11_S21(l1,Lh,Xh,Bh,BC,opts);
                                    s21 = 1i*imag(S21)*[0,1]';
                                    [S11p, SFlag] = S11_S21(l1 + reDel,Lh,Xh,Bh,BC,opts);
                                    s11p = real(S11p)*[1,0]';
                                    [S11m, SFlag] = S11_S21(l1 - reDel,Lh,Xh,Bh,BC,opts);
                                    s11m = real(S11m)*[1,0]';
                                end
                                Ds11 = real(s11p - s11m)/reDel;
                                Signs1p(ii) = real(real(1i*s21)/Ds11);
                            end
                        end
                    end
                end
            

                if OldBracketMethod
                    %% Finalize the eigenvalue brackets found above
                    if NBrkts > 0 %Pair down the array to match the actual # of points
                        ReBrktsfinal = zeros(2,NBrkts);
                        for ii = 1:NBrkts
                            ReBrktsfinal(:,ii) = ReBrkts(:,ii);
                        end   
                    end

                    % ReEvFinal(1) = {'No real eigenvalues found.'};%default case valid if NBrkts = 0
                    ReEvFinal = [];%default case valid if NBrkts = 0

                    %% Background on Bisection Method used to refine real eigenvalues
                    % Apply the Bisection method to the bracketting l-values recorded in 'ReEvlist'
                    %to improve precision. 
                    % The Bisection Method:
                    % Suppose la and lb are listed sequentially in ReEvlist, then if
                    % fa = sign(real(S11(la)) and fb = sign(real(S11(lb)), it must be
                    % that fa*fb < 0
                    %Take c = (a+b)/2, then there are two cases ...
                    % 1) if fa*fc < 0, then the zero must be between a and c, so narrow 
                    %the bracket with b = c and repeat.
                    % 2) else fa and fc have the same sign and the zero must be between
                    % c and b, then update a = c
                    % Iterate a fixed number of times or until targeted precision is
                    % reached. The method will work if there is just a single zero
                    % between each pair of bracket points. So it is possible it could
                    % be blind to multiple zeroes ... and so overlook zeros.

                    %% Refine Real Eigenvalues, if any, and evaluate signs of corresponding solitons
                    NReEvs = NBrkts;
                    Signs1p = 0;
                    if NBrkts > 0 % if at least one pair of bracketting points is found
                                  % then refine eigenvalues and evaluate their signs

                        ReEvFinal = zeros(1,NBrkts); % the # real EV's = # of pairs of bracketting points  
                        Signs1p = zeros(1,NBrkts);

                        % Evaluate final real eigenvalues using the Bisection Method
                        for ii = 1:NBrkts
                            a = ReBrktsfinal(1,ii);
                            b = ReBrktsfinal(2,ii);            
            %                 Sa = S11_S21(a,Lh,Xh,Bh,BC,opts);
                            [xL_is_alphaXR] = alpha_XR(real(a),imag(a),abs(Bh(1)),abs(Bh(end)));        
                            XR = ( Xh(end) - Xh(1))/(1 + xL_is_alphaXR);
                            XL = - xL_is_alphaXR*XR;
                            Nh = length(Bh);
                            Xh = linspace(XL, XR, Nh);
                            if SMatrixMethod == 1
                                [Sa, SFlag, ~, ~,  ~, ~, ~, ~, ~,~,~] = S_Adaptive_Complete(Xh, Bh, 0, NFracSS, a, eps_AS, Rfactor);
                                % [S_Adapt1121, ErrFlag, Xuniform, b2nsymm,  RelPctErr, ErrStep, StepsTaken, S_Adaptive_Time, msg_S_terse, msg_S_verbose, SA_ERR]
                            else
                                [Sa, SFlag] = S11_S21(a,Lh,Xh,Bh,BC,opts);
                            end
                            fa = sign(real(Sa(1)));
                            for k = 1:NBSIter %actual iteration step of Bisection Method
                                c = (a+b)/2; %1/2 way between points bracketing a sign change
            %                     Sc = S11_S21(c,Lh,Xh,Bh,BC,opts);
                                [xL_is_alphaXR] = alpha_XR(real(c),imag(c),abs(Bh(1)),abs(Bh(end)));        
                                XR = ( Xh(end) - Xh(1))/(1 + xL_is_alphaXR);
                                XL = - xL_is_alphaXR*XR;
                                Nh = length(Bh);
                                Xh = linspace(XL, XR, Nh);
                                if SMatrixMethod == 1
                                    [Sc, SFlag,  ~, ~,  ~, ~, ~, ~, ~,~,~] = S_Adaptive_Complete(Xh, Bh, 0, NFracSS, c, eps_AS);
                                else
                                    [Sc, SFlag] = S11_S21(c,Lh,Xh,Bh,BC,opts);
                                end
                                fc = sign(real(Sc(1)));
                                if fc*fa < 0
                                    b = c; % if so fc must have the same sign as b, so narrow the bracket
                                else
                                    a = c; % if not, fc must have the same sign as a, so narrow the bracket
                                end
                            end
                            ReEvFinal(ii) = real(a+b)/2;% return the 1/2 way point of the last bracket
                            ['Real Eigenvalue calculation: at ', num2str(ii), ' of ', num2str(NBrkts)]
                        end 

                        % Evaluate Signs for final real eigenvalues
                        for ii = 1:NBrkts
                            l1 = ReEvFinal(ii);
            %                 Dl = 0.0001;
            %                 s21 = 1i*imag(S11_S21(l1,Lh,Xh,Bh,BC,opts))*[0,1]';
            %                 s11p = real(S11_S21(l1 + reDel,Lh,Xh,Bh,BC,opts))*[1,0]';
            %                 s11m = real(S11_S21(l1 - reDel,Lh,Xh,Bh,BC,opts))*[1,0]';
                            [xL_is_alphaXR] = alpha_XR(real(l1),imag(l1),abs(Bh(1)),abs(Bh(end)));        
                            XR = ( Xh(end) - Xh(1))/(1 + xL_is_alphaXR);
                            XL = - xL_is_alphaXR*XR;
                            Nh = length(Bh);
                            Xh = linspace(XL, XR, Nh);
                            if SMatrixMethod == 1
                                [S21, SFlag, ~, ~,  ~, ~, ~, ~, ~,~,~] = S_Adaptive_Complete(Xh, Bh, 0, NFracSS, l1, eps_AS, Rfactor);
                                s21 = 1i*imag(transpose(S21))*[0,1]';
                                [S11p, SFlag, ~, ~,  ~, ~, ~, ~, ~,~,~] = S_Adaptive_Complete(Xh, Bh, 0, NFracSS, l1 + reDel, eps_AS, Rfactor);
                                s11p = real(transpose(S11p))*[1,0]';
                                [S11m, SFlag,  ~, ~,  ~, ~, ~, ~, ~,~,~] = S_Adaptive_Complete(Xh, Bh, 0, NFracSS, l1 - reDel, eps_AS, Rfactor);
                                s11m = real(transpose(S11m))*[1,0]';
                            else
                                [S21, SFlag] = S11_S21(l1,Lh,Xh,Bh,BC,opts);
                                s21 = 1i*imag(S21)*[0,1]';
                                [S11p, SFlag] = S11_S21(l1 + reDel,Lh,Xh,Bh,BC,opts);
                                s11p = real(S11p)*[1,0]';
                                [S11m, SFlag] = S11_S21(l1 - reDel,Lh,Xh,Bh,BC,opts);
                                s11m = real(S11m)*[1,0]';
                            end
                            Ds11 = real(s11p - s11m)/reDel;
                            Signs1p(ii) = real(real(1i*s21)/Ds11);
                        end
                    end
                end
            
            
%             ReEvFinal
%             Signs1p
%             NBrkts
        end
%         
        
        
    end

    function [ReEvFinal, Signs1p, NReEvs] = RealEVs_Grock(LREval, NEval, NBSIter, Lh, Xh, Bh, BC, opts)
        %(LREval, NREval, NBSIter, Lh, Xh, Bh, opts)
        %Return list of real eigenvalues in range within LREval
        %NEval specifies the resolution
        NBrkts = 0; %counter for # of pairs of bracket points
        LEvalGrid = linspace(LREval(1),LREval(2),NEval);
        reDel = 0.1*(LREval(2) - LREval(1))/NEval;
        %*****************************
        %S11(l) evalutated on the real axis for 0 < l < Bo must be real valued
        %We will use the Bisection method which relies on knowing the
        %locations, ln, such that s11(ln)*s11(ln+1) < 0.
        %the first loop below sweeps through a fairly fine grid of points
        %and records one value of ln every time the sign is observed to change
        %See: http://www.math.niu.edu/~dattab/MATH435.2013/ROOT_FINDING.pdf
        ReBrkts = zeros(2,NEval);%store bracket points in pairs (a, b)        
        SFlag = 1;
%         Bh = Bh*exp(1i*0.5*pi);
        [S, SFlag] = S11_S21(LEvalGrid(:,1),Lh,Xh,Bh,BC, opts);
        if SFlag < 1
            ['Error in evaluation of Scattering Data: Psiminus']
            return
        end
        s11 = S(1);
        signold = sign(real(s11));
%         signold
        Sdebug = zeros(2,200);
        for ii = 2:NEval
            [S, SFlag] = S11_S21(LEvalGrid(:,ii),Lh,Xh, Bh, BC,opts);
            if SFlag < 1
                ['Error in evaluation of Scattering Data: Psiminus']
                % ReEvFinal(1) = {'Error in evaluation of Scattering Data: Psiminus.'};%default case valid if NBrkts = 0
                ReEvFinal(1) = [];%default case valid if NBrkts = 0
                Signs1p = 0;
                NReEvs = 0;
                return
            end
            s11 = S(1);
            Sdebug(:,ii) = transpose(S);
            signnew = sign(real(s11));
%             signnew
            if signnew ~= signold
                NBrkts = NBrkts + 1;
                ReBrkts(1,NBrkts) = LEvalGrid(:,ii-1);%record la of this bracket pair   
                ReBrkts(2,NBrkts) = LEvalGrid(:,ii);%record lb of this bracket pair 
                signold = signnew; 
            end
%             ['evaluating signs ', num2str(ii), ' of ', num2str(NEval)]
            ['s11 = ',num2str(real(transpose(S(1)))) ', i(s21) =  ',num2str(real(1i*transpose(S(2)))),': at ', num2str(ii), ' of ', num2str(NEval)]

        end
%          ReBrkts
%          NBrkts

        if NBrkts > 0 %Pair down the array to match the actual # of points
            ReBrktsfinal = zeros(2,NBrkts);
            for ii = 1:NBrkts
                ReBrktsfinal(:,ii) = ReBrkts(:,ii);
            end   
        end
        
        % ReEvFinal(1) = {'No real eigenvalues found.'};%default case valid if NBrkts = 0
        ReEvFinal(1) = [];%default case valid if NBrkts = 0
        
        %Apply the Bisection method to the bracketting l-values recorded in 'ReEvlist'
        %to improve precision. 
        % The Bisection Method:
        % Suppose la and lb are listed sequentially in ReEvlist, then if
        % fa = sign(real(S11(la)) and fb = sign(real(S11(lb)), it must be
        % that fa*fb < 0
        %Take c = (a+b)/2, then there are two cases ...
        % 1) if fa*fc < 0, then the zero must be between a and c, so narrow 
        %the bracket with b = c and repeat.
        % 2) else fa and fc have the same sign and the zero must be between
        % c and b, then update a = c
        % Iterate a fixed number of times or until targeted precision is
        % reached. The method will work if there is just a single zero
        % between each pair of bracket points. So it is possible it could
        % be blind to multiple zeroes ... and so overlook zeros.
        NReEvs = NBrkts;
        Signs1p = 0;
        if NBrkts > 0 %if at least one pair of bracketting points is found
            ReEvFinal = zeros(1,NBrkts); % the # real EV's = # of pairs of bracketting points  
            Signs1p = zeros(1,NBrkts);
            for ii = 1:NBrkts
                a = ReBrktsfinal(1,ii);
                b = ReBrktsfinal(2,ii);            
                Sa = S11_S21(a,Lh,Xh,Bh,BC,opts);
                fa = sign(real(Sa(1)));
                for k = 1:NBSIter %actual iteration step of Bisection Method
                    c = (a+b)/2; %1/2 way between points bracketing a sign change
                    Sc = S11_S21(c,Lh,Xh,Bh,BC,opts);
                    fc = sign(real(Sc(1)));
                    if fc*fa < 0
                        b = c; % if so fc must have the same sign as b, so narrow the bracket
                    else
                        a = c; % if not, fc must have the same sign as a, so narrow the bracket
                    end
                end
                ReEvFinal(ii) = real(a+b)/2;% return the 1/2 way point of the last bracket
                ['Real Eigenvalue calculation: at ', num2str(ii), ' of ', num2str(NBrkts)]
            end 
            
            for ii = 1:NBrkts
                l1 = ReEvFinal(ii);
%                 Dl = 0.0001;
                s21 = 1i*imag(S11_S21(l1,Lh,Xh,Bh,BC,opts))*[0,1]';
                s11p = real(S11_S21(l1 + reDel,Lh,Xh,Bh,BC,opts))*[1,0]';
                s11m = real(S11_S21(l1 - reDel,Lh,Xh,Bh,BC,opts))*[1,0]';
                Ds11 = real(s11p - s11m)/reDel;
                Signs1p(ii) = real(real(1i*s21)/Ds11);
            end
        end
%         NBrkts
%         ReEvFinal
%         Signs1p
%         
        
        
    end

    function [Psi, PsiFlag] = Psiminus(l,Lh,Xh,Bh,Bo, opts, Psiminus_Debug)
        %   Returns Psi = [Psi11, Psi21] given eigenvalue, l, length, L,
        %   B-field, B4, and grid, X4.
        %  Integration of Psi ... as step towards evaluating scattering
        %  coefficients Phi = Psi.T.R -->Psix = Dpsi.Psi where ... Dpsi =
        %  Inv(T.R).[D.T.R + i*l*z*T.sigma3.R] which simplifies to Dpsi =
        %  l/z*[iA, B*Exp(-2ilzx);Bb*B*Exp(+2ilzx),-iA]
        %      with A = bo(Re(b) - bo); B = (l*del - i*bo^2*Im(b)/(l -
        %      z))*Exp(-2ilzx);
        %                               Bb = (l*Conj(del) + i*bo^2*Im(b)/(l
        %                               - z))*Exp(+2ilzx) del = b - bo
        %Asymptotics of Psi
        PsiFlag = 1;
        Psi0 = [1, 0];range = Lh;% Lh = [EVL, EVR] & is either full range or truncated to match RABC section of profile
        
        ['From Psiminus: before ode']
        sprintf('l = %4.4f',l)
        
        [x4,y4] = ode113(@(x4,y4)derivPsi(x4,y4,Xh,Bh,Bo, l),range,Psi0,opts);
%         [x4,y4] = ode15s(@(x4,y4)derivPsi(x4,y4,Xh,Bh,Bo, l),range,Psi0,opts);
        
        ['From Psiminus: after ode']
        sprintf('l = %4.4f, y4(end,1) = %4.4f, length(y4(:,1)) = %4.0f & range = %4.4f',l, y4(end,1),length(y4(:,1)),x4(end) - x4(1))
        if (x4(end) - x4(1)) < (range(2) - range(1))
            PsiFlag = -1;
            Psi = Psi0;
            ['Error in Psiminus evaluation.']
            return
        end
        
%         Psiminus_Debug = 1;
        if Psiminus_Debug == 1
            PsiFig = figure('Color', [1, 1, 1],'units','normalized','outerposition',[(NMonitors -1 +.1) 0.1 0.8 0.8]);
                cla();
                PXmin = x4(1); PXmax = x4(end);PXrange = PXmax - PXmin;
                PZmax = 1.5*max(abs(y4(:,1)));
                plot(x4,real(y4(:,1)));
                hold('on');
                plot(x4,imag(y4(:,1)));
                hold('on');

                plot(x4,abs(y4),'x','MarkerIndices', 1:round(N/100):length(X4),'Color',[0 , 0, 1]);  % |B| is blue

                hold('on');
                text(PXmin + 0.1*PXrange,0.9*PZmax,sprintf('lambda = %4.4f',l));
                text(PXmin + 0.1*PXrange,0.7*PZmax,sprintf('Left click to proceed.'));
    %             text(EVL*0.7,0.92*zmax,sprintf('Trial ID = %s',TrialID));
    %             text(EVL*0.9,0.85*zmax,sprintf('time = %4.4f',TimePerLine*(LineNumber - 1)));
    %             text(EVL*0.9,0.75*zmax,sprintf('line # = %4.0f',LineNumber));

                axis([PXmin PXmax -PZmax PZmax]);
                drawnow(); 
                w = waitforbuttonpress;
                if w ==0
                    disp('Waiting for click.')
                else
                    disp('Moving on ... ')
                end
                
                 
            close(PsiFig);
            
                
        end
        
        Psi = y4(end,:)  ;
        
        
        
        
        
       
    end

    function [Psi, PsiFlag] = Psiminus2(l,Lh,Xh,Bh,Bo, opts, Psiminus_Debug)
        %%   Returns Psi = [Psi11, Psi21] given eigenvalue, l, length, L,
        %   B-field, B4, and grid, X4.
        %   Integration of Psi ... using unalterned scattering problem:
        %   Psi(+/-) --> W(+/-).Exp[-iLambda(+/-).sigma3*x]
        %
        %   W(+/-) = {{l+z(+/-), z(+/-) - l},{i*conj(b(+/-)),-i*conj(b(+/-))}}
        %   Lambda(+/-) = l*Sqrt(l^2-abs(b(+/-))^2)
        %
        %   d(Psi)/dx = {{-i*l^2, l*b(x)},{l*conj(b), +i*l^2}}.Psi
        %
        %   S = inv(Psi(+)).Psi(-) for x --> +Inf
        %
        %   S = Exp[iLambda(+).sigma3*x].inv(W(+)).Psi(-)(x) for x --> + Inf
        
        %% Asymptotics of Psi 
        PsiFlag = 1;
        
        Xm = Lh(1,1);
        ExpXm = [ exp(-1i*l*z(l,Bo)*Xm), 0; 0, exp(1i*l*z(l,Bo)*Xm) ];        
        Wm = [l+z(l,Bo), z(l,Bo) - l;1i*conj(Bo),-1i*conj(Bo)];
        Psi0 = Wm * ExpXm * [1, 0]';range = Lh;% Lh = [EVL, EVR] & is either full range or truncated to match RABC section of profile
        
        if Psiminus_Debug == 1
            ['From Psiminus: before ode']
            sprintf('l = %4.4f',l)
        end
%         opts4 = odeset('RelTol',1e-6,'AbsTo',1e-8);
        [x4,y4] = ode113(@(x4,y4)derivPsi2(x4,y4,Xh,Bh,Bo, l,Lh),range,Psi0,opts);
%         opts5 = odeset('RelTol',1e-8,'AbsTo',1e-10);
%         [x5,y5] = ode113(@(x5,y5)derivPsi2(x5,y5,Xh,Bh,Bo, l,Lh),range,Psi0,opts5);
%         [x4,y4] = ode15s(@(x4,y4)derivPsi(x4,y4,Xh,Bh,Bo, l),range,Psi0,opts);
        
        if Psiminus_Debug == 1
            ['From Psiminus: after ode']
            sprintf('l = %4.4f, y4(end,1) = %4.4f, length(y4(:,1)) = %4.0f & range = %4.4f',l, y4(end,1),length(y4(:,1)),x4(end) - x4(1))
        end
        if (x4(end) - x4(1)) < (range(2) - range(1))
            PsiFlag = -1;
            Psi = Psi0;
            ['Error in Psiminus evaluation.']
            return
        end
        
%         Psiminus_Debug = 1;
        if Psiminus_Debug == 1
            PsiFig = figure('Color', [1, 1, 1],'units','normalized','outerposition',[(NMonitors -1 +.1) 0.1 0.8 0.8]);
                cla();
                PXmin = x4(1); PXmax = x4(end);PXrange = PXmax - PXmin;
                PZmax = 1.5*max(abs(y4(:,1)));
                plot(x4,real(y4(:,1)));
                hold('on');
                plot(x4,imag(y4(:,1)));
                hold('on');

                plot(x4,abs(y4),'x','MarkerIndices', 1:round(N/100):length(X4),'Color',[0 , 0, 1]);  % |B| is blue

                hold('on');
                text(PXmin + 0.1*PXrange,0.9*PZmax,sprintf('lambda = %4.4f',l));
                text(PXmin + 0.1*PXrange,0.7*PZmax,sprintf('Left click to proceed.'));
    %             text(EVL*0.7,0.92*zmax,sprintf('Trial ID = %s',TrialID));
    %             text(EVL*0.9,0.85*zmax,sprintf('time = %4.4f',TimePerLine*(LineNumber - 1)));
    %             text(EVL*0.9,0.75*zmax,sprintf('line # = %4.0f',LineNumber));

                axis([PXmin PXmax -PZmax PZmax]);
                drawnow(); 
                w = waitforbuttonpress;
                if w ==0
                    disp('Waiting for click.')
                else
                    disp('Moving on ... ')
                end
                
                 
            close(PsiFig);
            
                
        end
        
        Psi = y4(end,:)  ;    
       
    end

    function [PhiML, PhiFlag] = PhiMinus(l,Xm,BB, PhiMinus_Debug)
        %%   Returns PhiMinus = Yp(m)*exp(-iLambda*x*sigma3) given eigenvalue, l, length, x = +L(-L) and B-field,
        %   B-field, B4, and grid, X4.
        %
        %   Yp(m) = {{            exp(-iLambda*x)      ,   -ib/(l+z)*exp(-iLambda*x)  }
        %           ,{ +i*conj(b)/(l+z)*exp(+iLambda*x),   exp(+iLambda*x)            }}
        %
        %   where x = -L and
        %   Lambda(+/-) = l*Sqrt(l^2-abs(b(+/-))^2)
        
        %% Asymptotics of Psi 
        PhiFlag = 1;
        LMBDA = l*z(l,BB);    
        PhiML = [       exp(-1i*LMBDA*Xm)             ,   -1i*BB/(l+z(l,BB))*exp(-1i*LMBDA*Xm)  ; ...
                   +1i*conj(BB)/(l+z(l,BB))*exp(+1i*LMBDA*Xm),       exp(+1i*LMBDA*Xm)                      ];
        
        
        if PhiMinus_Debug == 1
            ['From PhiMinus: after assignment of PhiML']
            sprintf('l = %4.4f',l)
        end
       
    end

    function [PhiPL, PhiFlag] = PhiPlus(l,Xm,BB, PhiPlus_Debug)
        %%   Returns PhiPlus = Inverse(   Yp*exp(-iLambda*x*sigma3) )given eigenvalue, l, length, x = +L and B-field,
        %   B-field, B4, and grid, X4.
        %
        %   Yp(m) = {{            exp(-iLambda*x)      ,   -ib/(l+z)*exp(-iLambda*x)  }
        %           ,{ +i*conj(b)/(l+z)*exp(+iLambda*x),   exp(+iLambda*x)            }}
        %
        %   where x = -L and
        %   Lambda(+/-) = l*Sqrt(l^2-abs(b(+/-))^2)
        
        %% Asymptotics of Psi 
        PhiFlag = 1;
        LMBDA = l*z(l,BB);    
        PhiPL = 1/2*[ ( 1 + l/z(l,BB) )*exp(+1i*LMBDA*Xm)     ,   1i*BB/(z(l,BB))*exp(-1i*LMBDA*Xm)  ; ...
                      -1i*conj(BB)/(z(l,BB))*exp(+1i*LMBDA*Xm),   ( 1 + l/z(l,BB) )*exp(-1i*LMBDA*Xm)                      ];
        
        
        if PhiPlus_Debug == 1
            ['From PhiPlus: after assignment of PhiML']
            sprintf('l = %4.4f',l)
        end
       
    end

    function [S1121, SFlag] = S11_S21(l,Lh,Xh,Bh,BC, opts)
        %   Returns [s11, s21] given eigenvalue, l, length, L, B-field, B4,
        %   and grid, X4, and direction of B-field to the left, theta_M,
        %   and to the right, theta_P.
        %  First, determine the value of Psi_minus ...        
        
        PsiFlag = 1;
        SFlag = 1;
        PsiDebug = 0;
        if PsiDebug == 1
            ['From S1121']
            [l Xh(1) Xh(end)]
        end
%         [Psi,PsiFlag] = Psiminus(l,Lh,Xh,Bh,BC(1),opts, PsiDebug);
        Xh_Symmetric = linspace(Lh(1),Lh(2),length(Bh));
        [Psi,PsiFlag] = Psiminus2(l,Lh,Xh_Symmetric,Bh,BC(1),opts, PsiDebug);
        if PsiFlag < 1
            ['Error in Psiminus evaluation.']
            SFlag = -1;
%             return
        end
        
        %Psiminus(l,Lh,Xh,Bh,Bo, opts)
        
        %Next evaluate the matrix, Mg, which connects Psi to the Scattering Matrix.
        bm = BC(1); bp = BC(2);
%         Mg = Mgamma(bm, bp, l, Lh);
        Mg = Mgamma2(bm, bp, l, Lh);
        Stemp = Mg*transpose(Psi);         
        S1121 = transpose(Stemp);
%         S1121
    end

    function Mg = Mgamma(bm, bp, l, Lh)
        g = conj(bm/bp);
        Xp = Lh(2);
        
        wmp = l + z(l,bm); wmm = l - z(l,bm);
        wpp = l + z(l,bp); wpm = l - z(l,bp);
        
        ExpLp = exp(1i*l*z(l,bp)*Xp);
        ExpLm = exp(1i*l*z(l,bm)*Xp);
        
        W11 = wmp - g*wpm; W12 = -(wmm - g*wpm);
        W21 = wmp - g*wpp; W22 = -( wmm - g*wpp);
        
        Aw = W11*ExpLp/ExpLm;   Bw = W12*ExpLp*ExpLm;       
        Cw = W21/(ExpLp*ExpLm); Dw = W22*ExpLm/ExpLp;
        
        Mg = (l/z(l,bp))*[Aw, Bw; Cw, Dw];
    end

    function Mg = Mgamma2(bm, bp, l, Lh)
        zp = z(l,bp);
        bps = conj(bp);
        WpInv = [ -1i*bps, -(zp - l);-1i*bps, zp + l]*1i/(2*bps*zp);
        
        Xp = Lh(2);
        ExpXp = [ exp(1i*l*z(l,bp)*Xp), 0; 0, exp(-1i*l*z(l,bp)*Xp) ];
        
        Mg = ExpXp*WpInv;
    end

    function dY_dx = derivPsi(x,Y,x4,b,Bo, l)
        bb = interp1(x4,b,x);
        dY_dx = Dpsi(l,bb,Bo,x)*Y;
    end

    function dY_dx = derivPsi2(x,Y,x4,b,Bo, l, Lh)
        bb = interp1(x4,b,x,'spline');
        dY_dx = Dpsi2(l,bb,Bo,x)*Y;
    end

    function out = Dpsi2(l, b, bo, x)
        out = (alpha/dispersion)*[-1i*l^2,l*b;l*conj(b),1i*l^2];
    end

    function out = Dpsi(l, b, bo, x)
        out = [1i*A(b,bo),B(l,b,bo,x);Bstar(l,b,bo,x),-1i*A(b,bo)]*l/z(l,bo);
    end

    function out = A(b, bo)
%         out = bo*(real(b) - bo);
        out = b*real(bo) - 1i*bo*imag(b) - abs(bo)^2;
    end

    function out = B(l,b, bo,x)
%         out = (-l*real(b - bo) - 1i*z(l,bo)*imag(b))*(exp(+2i*l*z(l,bo)*x));
        out = 1i/(2*conj(bo))*(-(b - bo)*conj(bo)^2 - conj(b - bo)*(l - z(l,bo))^2).*(exp(2i*l*z(l,bo)*x));
    end

    function out = Bstar(l,b, bo,x)
%         out = (-l*real(b - bo) + 1i*z(l,bo)*imag(b))*(exp(-2i*l*z(l,bo)*x));
        out = 1i/(2*conj(bo))*((b - bo)*conj(bo)^2 + conj(b - bo)*(l + z(l,bo))^2).*(exp(-2i*l*z(l,bo)*x));
    end

    function out = z(in,B0)
        out = sqrt(in^2 - abs(B0)^2);
        if imag(out) < 0
            out = -out;
        end
    end

    function mouseMove(object, eventdata)
        CPos = get(gca, 'CurrentPoint');
       
        title(gca, num2str(CPos(1,1) + 1i*CPos(1,2)));
    end
    
    function[idMax, idSL, idSR, xShockL, xShockR, Iz] = IS_ShockBoundaries2(bb,x,bL, basymptote, ntest, tol)
        
    %% ARGUMENTS: bb, x, bL, basymptote, ntest, tol
        % bb --> section of B-field with IntermediatE Shock on LHS preceeded
        % by section with constant B-field = ( bL, 0)
        %
        % x --> array of position values of same length as bb
        %
        % bL ---> as described in bb above
        %
        % basymptote --> it is assumed that the z-component of the B-field
        % goes to this constant asymptotic value to the right of the shock.
        % It is should be evaluated prior to call as basymptote = b_L*sin(theta_kbw)
        %
        %ntest --> number of indices over which the running average of Bz
        %is evaluated   
        
    %% OUTPUTS: idMax, idSL, idSR, xShockL, xShockR, Iz
        % idMax --> the location of the maximum value of the z-component of
        % the B-field. For a dispersive Intermediate Shock the z-component
        % eventually settles into an oscillatory wave similar in appearance
        % to a damped sine-wave. prior to settling to this steady-state
        % this is not a well defined measure of shock's location.
        %
        % idSL --> index of the LHS of Shock
        % idSR --> index of the RHS of Shock
        % xShockL --> postion of the LHS of Shock
        % xShockR --> position of the RHS of Shock
        
        %evaluate the index of a 'fixed' location within the shock
        %Note: this location is only fixed once the shock reaches
        %steady-state
        
   
    %% Find LHS of Shock: idSL, xShockL assumes shock has a constant B-field on its LHS
        % with B = (bL, 0)
        LHS = abs(abs(bb-bL)) <0.001; 
        %this returns an array of length(bb) with logical '1' in every 
        %location where the condition is satisfied and a '0' otherwise
        
        SATL = find(LHS ==1); %returns an array containing an order list of indices where the condition above holds
        idSL = SATL(end); % returns the rightmost index where the condition holds
        xShockL = x(idSL); % returns the position of the left edge of the shock
     
     %% Find first peak of Im(bb): idMax starting point for search of RHS
        [val, idMax] = max(abs(imag(bb)));
 
    %% Find RHS of Shock: assumes Im(b) has general form of damped sign wave towards asymptote
        %Assuming Coplanar shock and imag(Bh) --> 0 to right of shock
        Nb = length(bb);
        idSR = Nb;
        
        % Check for monotonic decrease to asymptote
        BRasymptoteSuccess = 0;
        minDeltaBz = 10;
        for idb = idMax:Nb
            DeltaBz = abs(imag(bb(idb)) - basymptote);
            if DeltaBz < minDeltaBz
                minDeltaBz = DeltaBz;
            else
%                 'Not monotonic'
                break;
            end
            if minDeltaBz < tol*abs(b_L)
                idSR = idb;
%                 xShockR = x(idSR);
                BRasymptoteSuccess = 1;
                break;
            end
        end
%         [val,idYmin] = min(Y);
        
        if ~BRasymptoteSuccess
            for idb = idMax:Nb
                Bzavg = sum(imag(bb(idb-ntest:idb)))/ntest;
                if ( ( imag(bb(idb-1)) - basymptote )*( imag(bb(idb)) - basymptote )) < 0 && (abs(Bzavg - basymptote ) < tol*abs(b_L))
                    idSR = idb;
                    break;
                end
            end
        end
        xShockR = x(idSR);
        
        dxx = (x(2) - x(1));
        Iz = sum(imag(bb(idSL:idSR)))*dxx;
    end

    function[idMax, idSL, idSR, xShockL, xShockR, Iz] = IS_ShockBoundaries(bb,x,bL, basymptoteL, basymptoteR, ntest, tol)
        
    %% ARGUMENTS: bb, x, bL, basymptote, ntest, tol
        % bb --> section of B-field with IntermediatE Shock on LHS preceeded
        % by section with constant B-field = ( bL, 0)
        %
        % x --> array of position values of same length as bb
        %
        % bL ---> as described in bb above
        %
        % basymptote --> it is assumed that the z-component of the B-field
        % goes to this constant asymptotic value to the right of the shock.
        % It is should be evaluated prior to call as basymptote = b_L*sin(theta_kbw)
        %
        %ntest --> number of indices over which the running average of Bz
        %is evaluated   
        
    %% OUTPUTS: idMax, idSL, idSR, xShockL, xShockR, Iz
        % idMax --> the location of the maximum value of the z-component of
        % the B-field. For a dispersive Intermediate Shock the z-component
        % eventually settles into an oscillatory wave similar in appearance
        % to a damped sine-wave. prior to settling to this steady-state
        % this is not a well defined measure of shock's location.
        %
        % idSL --> index of the LHS of Shock
        % idSR --> index of the RHS of Shock
        % xShockL --> postion of the LHS of Shock
        % xShockR --> position of the RHS of Shock
        
        %evaluate the index of a 'fixed' location within the shock
        %Note: this location is only fixed once the shock reaches
        %steady-state
     %% Find first peak of Im(bb): idMax starting point for search of RHS
        [val, idMax] = max(abs(imag(bb)));
    %% Find LHS of Shock: 
    %     1st) assume Im(b) has monotonic decrease to asymptote to the left
    %     2nd) see if Im(b) has oscillatory approach to asymptote
        % Check for monotonic decrease to asymptote on LHS
        
        BLasymptoteSuccess = 0;
        minDeltaBz = 10;       
        idSL = 1;
        Lexception = 0;
        for idb = 1:(idMax - (ntest - 1))
            idx = idMax-(idb-1);
%             DeltaBz = abs(imag(bb(idx)) - basymptoteL);
%             
%             if DeltaBz < minDeltaBz
%                 minDeltaBz = DeltaBz;
%                 Lexception = 0;
%             else
%                 Lexception = Lexception +1';
%             end
%             
%             if Lexception > 50
%                 msg = 'Not monotonic'
%                 minDeltaBz = 10;
%             end
            
           
            idxLeft = idx-(ntest - 1);
            idxRight = idx + (ntest - 1);
            
            DeltaBzAvgL = sum(imag(bb(idxLeft:idx)))/ntest - basymptoteL;
            DeltaBzAvgR = sum(imag(bb(idx:idxRight)))/ntest - basymptoteL;
            
            SGNBzAvgL = sign(DeltaBzAvgL);
            SGNBzAvgR = sign(DeltaBzAvgR);
            
            bbzslopeR = (imag(bb(idxRight)) - imag(bb(idx)))/(x(idxRight) - x(idx));
            bbzslopeL = (imag(bb(idx)) - imag(bb(idxLeft)))/(x(idx) - x(idxLeft));
            
            SGNCrossing = SGNBzAvgL*SGNBzAvgR; % both should have the same sign if not evaluated at a crossing
            SGNslopes = sign(bbzslopeR*bbzslopeL);% Should be positive if approaching an asymptote from either above or below.
                                                  % Will be negative if coincidentally at min or max
            if (abs(DeltaBzAvgL) < tol*abs(b_L)) &&  (abs(DeltaBzAvgR) < tol*abs(b_L))  && (SGNCrossing > 0) && (SGNslopes > 0)
                idSL = idx;
                BLasymptoteSuccess = 1;
                break;
            end
        end
        
        BLoscillatory = 0;
        if ~BLasymptoteSuccess
            idSL = 1;
            for idb = 1:(idMax - ntest + 1)
                idxLeft = idMax-(idb-1)-(ntest - 1);
                idxRight = idMax-(idb-1);
                BzAvg = sum(imag(bb(idxLeft:idxRight)))/ntest;
                signLhsidb = sign(imag(bb(idxLeft)) - basymptoteL);
                signLhsidbPlusone = sign(abs(bb(idxLeft+1)) - basymptoteL);
                if (( signLhsidb*signLhsidbPlusone ) < 0) && (abs(BzAvg - basymptoteL ) < tol*abs(bL))
                    idSL = idb;
                    BLoscillatory = 1;
                    break;
                end
            end
        elseif ~BLoscillatory && ~BLasymptoteSuccess
            idSL = floor(idMax/2);
            msg = 'No well defined LHS for Shock: ... going with midpoint'
        end
        
        xShockL = x(idSL);
        
       
 
    %% Find RHS of Shock: assumes Im(b) has general form of damped sign wave towards asymptote
        %Assuming Coplanar shock and imag(Bh) --> 0 to right of shock
        Nb = length(bb);
        idSR = floor((Nb + idMax)/2);
        Rexception = 0;
        
        % Check for monotonic decrease to asymptote on LHS
        BRasymptoteSuccess = 0;
        minDeltaBz = 10;
        for idb = (idMax + (ntest - 1)):(Nb - (ntest - 1))
%             DeltaBz = abs(imag(bb(idb)) - basymptoteR);
%             if DeltaBz < minDeltaBz
%                 minDeltaBz = DeltaBz;
%                 Rexception = 0;
%             else
%                 Rexception = Rexception +1;
%             end
%             
%             if Rexception > 50
%                 msg = 'Not monotonic'
% %                 break;
%             end
            
            idxLeft = idb - (ntest - 1);
            idxRight = idb + (ntest - 1);
             
            DeltaBzAvgL = sum(imag(bb(idxLeft:idb)))/ntest - basymptoteR;
            DeltaBzAvgR = sum(imag(bb(idb:idxRight)))/ntest - basymptoteR;
            
            SGNBzAvgL = sign(DeltaBzAvgL);
            SGNBzAvgR = sign(DeltaBzAvgR);
            
            bbzslopeR = (imag(bb(idxRight)) - imag(bb(idb)))/(x(idxRight) - x(idb));
            bbzslopeL = (imag(bb(idb)) - imag(bb(idxLeft)))/(x(idb) - x(idxLeft));
            
            SGNCrossing = SGNBzAvgL*SGNBzAvgR; % both should have the same sign if not evaluated at a crossing
            SGNslopes = sign(bbzslopeR*bbzslopeL);% Should be positive if approaching an asymptote from either above or below.
                                                  % Will be negative if coincidentally at min or max
           
            if (abs(DeltaBzAvgL) < tol*abs(b_L)) &&  (abs(DeltaBzAvgR) < tol*abs(b_L))  && (SGNCrossing > 0) && (SGNslopes > 0)
                idSR = idb;
                BRasymptoteSuccess = 1;
                break;
            end
        end

        BRoscillatory = 0;        
        if ~BRasymptoteSuccess
            for idb = idMax:(Nb - (ntest - 1) )
                idxLeft = idb;
                idxRight = idb + (ntest -1);
                Bzavg = sum( imag( bb( idxLeft:idxRight ) ) )/ntest;
                if ( ( imag(bb(idb-1)) - basymptoteR )*( imag(bb(idb)) - basymptoteR )) < 0 && (abs(Bzavg - basymptoteR ) < tol*abs(b_L))
                    idSR = idb;
                    BRoscillatory = 1;
                    break;
                end
            end
            elseif ~BRoscillatory && ~BRasymptoteSuccess
                idSR = floor((Nb + idMax)/2);
                msg = 'No well defined RHS for Shock: ... going with floor((Nb + IdMax)/2)'
        end
        xShockR = x(idSR);
        
        dxx = (x(1) - x(2)); % dxx defined as '-' to be consistent with convention in KBW paper in integrating across the shock from right to left
        Iz = sum(imag(bb(idSL:idSR)))*dxx;
    end

    function IS_Shock_Profile_Label(FigHandle, Zmin, Zmax, idOS, idSL, idSR, xSL, xSR, x, b)
        
        % Invoke Figure
        figure(FigHandle);
        
        %*********************************************
        % Plot vertical lines to Left and Right of Shock 
        % & horizontal line connecting them
        %*********************************************
        plot([xSL xSR],[0.8*Zmin, 0.8*Zmin],'Color',[0 0 0]); % Horizontal line underneath Shock
        hold('on');
        
        plot([xSL xSL],[0.8*Zmin, 0.8*Zmax],'Color',[0 0 0]); % Vertical Line on LHS of Shock
        hold('on');
        
        plot([xSR xSR],[0.8*Zmin, 0.8*Zmax],'Color',[0 0 0]); % Vertical Line on RHS of Shock
        hold('on');

        %*********************************************
        % Mark LHS of Shock
        %*********************************************
        
        XBzL = x(idSL); YBzL = imag(b(idSL));             % 'X' marker intersection of LHS with Im(b)
        plot(XBzL,YBzL,'x','MarkerFaceColor','red','MarkerSize',5);
        hold('on');

        idBzTxtL = idSL - floor(0.2*idOS); % Text label, value of Im(b) at LHS
        if idBzTxtL < 1
            idBzTxtL = 1;
        end
        XBzTxtL = x(idBzTxtL); YBzTxtL = imag(b(idSL))-0.1*Zmax;
        text(XBzTxtL,YBzTxtL,sprintf('Im(b) = %s ',num2str(imag(b(idSL)))),'Color',[1 0 0],'horizontalAlignment', 'right');
        hold 'on';

        XBmagL = x(idSL); YBmagL = abs(b(idSL));          % 'X' marker intersection of LHS with Abs(b)
        plot(XBmagL,YBmagL,'x','MarkerFaceColor','blue','MarkerSize',5);
        hold 'on';

        idBmagTxtL = idSL - floor(0.2*idOS);            % Text label, value of Abs(b) at LHS
        XBmagTxtL = x(idSL); YBzTxtL = abs(b(idSL))+ 0.1*Zmax;
        text(XBmagTxtL,YBzTxtL,sprintf('Abs(b) = %s ',num2str(abs(b(idSL)))),'Color',[0 0 1],'horizontalAlignment', 'right');
        hold 'on';

        %*********************************************
        % Mark RHS of Shock
        %*********************************************

        XBzR = x(idSR); YBzR = imag(b(idSR));
        plot(XBzR,YBzR,'x','MarkerFaceColor','red','MarkerSize',5); % 'X' marker intersection of RHS with Im(b)
        hold 'on';

        idBzTxtR = idSR + floor(0.2*idOS);              % Text label, value of Im(b) at RHS
        XBzTxtR = x(idBzTxtR); YBzTxtR = imag(b(idSR))-0.1*Zmax;
        text(XBzTxtR,YBzTxtR,sprintf('Im(b) = %s ',num2str(imag(b(idSR)))),'Color',[1 0 0]);
        hold 'on';

        XBmagR = x(idSR); YBmagR = abs(b(idSR));          % 'X' marker intersection of RHS with Abs(b)
        plot(XBmagR,YBmagR,'x','MarkerFaceColor','blue','MarkerSize',5);
        hold 'on';

        idBmagTxtR = idSR + floor(0.2*idOS);            % Text label, value of Abs(b) at RHS
        XBmagTxtR = x(idBmagTxtR); YBmagTxtR = abs(b(idSR))+ 0.1*Zmax;
        text(XBmagTxtR,YBmagTxtR,sprintf('Abs(b) = %s ',num2str(abs(b(idSR)))),'Color',[0 0 1]);
        hold 'on';
        
        idCRTxtR = idSR + floor(0.2*idOS);            % Text label, value of Abs(b) at RHS
        XCRTxtR = x(idCRTxtR); YCRTxtR = abs(b(idSR))+ 0.2*Zmax;
        text(XCRTxtR,YCRTxtR,sprintf('|Compression Ratio| = %s ',num2str( abs(b(idSR))/abs(b(idSL)) )),'Color',[0 0 1]);
        hold 'on';

    end

    function Figure_Label(FigHandle, Zmax, idOS, idL, idR, x, b)
        
            % Invoke Figure
            figure(FigHandle);
            
            %*********************************************
            % Mark Far LHS of Figure
            %*********************************************           

            XBzFL = x(idL); YBzFL = imag(b(idL));                         % 'X' marker intersection of Far LHS with Im(b)
            plot(XBzFL,YBzFL,'x','MarkerFaceColor','red','MarkerSize',5);
            hold 'on';
            
            idBzTxtFL = idL + floor(0.2*idOS);                      % Text label, value of Im(b) at Far LHS
            XBzTxtFL = x(idBzTxtFL); YBzTxtFL = imag(b(idL))-0.15*Zmax;
            text(XBzTxtFL,YBzTxtFL,sprintf('Im(b) = %s ',num2str(imag(b(idL)))),'Color',[1 0 0]);
            hold 'on';
            
            XBmagFR = x(idL); YBmagFR = abs(b(idL));                      % 'X' marker intersection of Far LHS with Abs(b)
            plot(XBmagFR,YBmagFR,'x','MarkerFaceColor','blue','MarkerSize',5);
            hold 'on';
             
            idBmagTxtFL = idL + floor(0.2*idOS);                    % Text label, value of Abs(b) at Far LHS
            XBmagTxtFL = x(idBmagTxtFL); YBmagTxtFL = abs(b(idL)) - 0.15*Zmax;
            text(XBmagTxtFL,YBmagTxtFL,sprintf('Abs(b) = %s ',num2str(abs(b(idL)))),'Color',[0 0 1]);
            hold 'on';
            
            %*********************************************
            % Mark Far RHS of Figure
            %*********************************************
            
            XBzFR = x(idR); YBzFR = imag(b(idR));                         % 'X' marker intersection of Far RHS with Im(b)
            plot(XBzFR,YBzFR,'x','MarkerFaceColor','red','MarkerSize',5);
            hold 'on';
            
            idBzTxtFR = idR - idOS;                                 % Text label, value of Im(b) at Far RHS
            XBzTxtFR = x(idBzTxtFR); YBzTxtFR = imag(b(idR))-0.15*Zmax;
            text(XBzTxtFR,YBzTxtFR,sprintf('Im(b) = %s ',num2str(imag(b(idR)))),'Color',[1 0 0],'horizontalAlignment', 'right');
            hold 'on';
            
            XBmagFR = x(idR); YBmagFR = abs(b(idR));                      % 'X' marker intersection of Far RHS with Abs(b)
            plot(XBmagFR,YBmagFR,'x','MarkerFaceColor','blue','MarkerSize',5);
            hold 'on';
             
            idBmagTxtFR = idR - idOS;                    % Text label, value of Abs(b) at Far RHS
            XBmagTxtFR= x(idBmagTxtFR); YBmagTxtFR = abs(b(idR)) - 0.15*Zmax;
            text(XBmagTxtFR,YBmagTxtFR,sprintf('Abs(b) = %s ',num2str(abs(b(idR)))),'Color',[0 0 1],'horizontalAlignment', 'right');
            hold 'on';

 
        
    end

    function [UU] = MatrixExp(lambda1, Bsel, ddx)
        LMDA1 = lambda1*z(lambda1, abs(Bsel));
        UU11 = cos(LMDA1*ddx) - 1i*lambda1*sin(LMDA1*ddx)/z(lambda1, abs(Bsel));   UU12 = Bsel*sin(LMDA1*ddx)/z(lambda1, abs(Bsel)); 
        UU21 = conj(  Bsel  )*sin(LMDA1*ddx)/z(lambda1, abs(Bsel));       UU22 = cos(LMDA1*ddx) + 1i*lambda1*sin(LMDA1*ddx)/z(lambda1, abs(Bsel));
                           
        UU = [UU11, UU12; UU21, UU22];
        
    end

    function [ppi] = Grock_Pi_N(StepSize, Xuni, Spline_N, bSpline, lmda, bb)
        %% deprecated to S_Adaptive
        %
        % Set up call to 'Grock_Pi_N'
        % Given: StepSize, dxx, Spline_N, bSpline AND lambda, bb
        % Evaluate Pi_N = U(N)*U(N-1)*U(N-2)* ... *U(2)*U(1)
        % where U(j) is Matrix_Exp(i*lambda*zeta*M)
        % with M(B,lambda) defining the 'scattering problem' ODE -->
        % d/dx(Phi) = M*Phi with Phi being a 2x2 solution matrix
        dxx = Xuni(2) - Xuni(1);
        DDxStep = StepSize*dxx;
        ppi = [1,0; 0,1];           
        npi = length(bb);

        for pidx = 1: StepSize: (npi - 1)
            bb_idx = bb(pidx);

            % If a higher, uniform resolution is required across the entire
            % B-field, it is much faster to compute this once at
            % higher memory cost PRIOR to iteration, so set b8spline = 0.
            % It may be more efficient to interpolate sections of the
            % B-field, if this is required, in sections. However,the embedded
            % interpolation below would still be inefficient as it requires
            % the field to be interpolated through a region of width,
            % 'Spline_N', at each point in a problematic section thereby
            % requiring 'Spline_N' times more computations than is
            % necessary.
            if ((pidx - Spline_N/2) > 0) && ((pidx + Spline_N/2 - 1) < npi) && bSpline
                Xpi8 = linspace(Xuni(pidx - Spline_N/2), Xuni(pidx + Spline_N/2 - 1),Spline_N*8); 
                bb8 = interp1(Xuni(pidx - Spline_N/2: (  pidx + Spline_N/2 -1)), bb(pidx - Spline_N/2: (  pidx + Spline_N/2 -1)), Xpi8); 
                [val, iiIntIdx] = min(abs(Xpi8 - Xuni(pidx)));
                bb_idx = bb8(iiIntIdx);
            end

            % Evaluate Propagation Matrix
            ppTmp = MatrixExp(lmda, bb_idx, DDxStep);
            ppi = ppTmp*ppi;

            % Output to Command Window for debugging:
            % fprintf(['index  = ',num2str(Pidx),'--> det(pi) =',num2str(det(pTmp)),'\n\n']);
        end
        
    end

    function [ppiOne, ppiSS, ErrPi11] = Grock_Pi_N_TwoStep(StepSize, Xuni, lmda, bb)
        %% deprecated to S_Adaptive
        %
    %% Uses basic 'Grock_Pi_N' to return Pi_N evaluated at StepSize = 1
    % AS WEll as incremental difference between U evaluated at
    % mod(idx,StepSize) = 0 with StepSize*dx for its distance
    %
    %  dx   dx   dx   dx   dx   dx   dx   dx   dx
    %  1,   2,   3,   4,   5,   6,   7,   8,   9,   ...
    %                 /\                  /\
    %                 ||                  ||
    % So, U(idx4 = '1',4*dx) 'should' equal U(1,dx)*U(2,dx)*U(3,dx)*U(4,dx)
    % 
    % Generally, U(idx*StepSize, StepSize*dx) = U( idx -(StepSize - 1),dx)*U(idx + 1 -(StepSize - 1),dx)* ... *U(idx,dx)
    % So, this function returns Pi_N evaluated both ways AND ALSO evaluates
    % the incremental differences between single StepSize U values and the
    % corresponding StepSize factors.
    % The purpose then is to graph the magnitude of these differences over
    % a B-field graph to see where the error accumulates the most.
    % Presumably sections with smaller radii of curvature are part of the
    % error, but areas with steep gradients may be as, or more important.
    % Insight will be helpful in crafting an error control algorithm.
    
        %% Set up call to 'Grock_Pi_N'
        % Given: StepSize, dxx, Spline_N, bSpline AND lambda, bb
        % Evaluate Pi_N = U(N)*U(N-1)*U(N-2)* ... *U(2)*U(1)
        % where U(j) is Matrix_Exp(i*lambda*zeta*M)
        % with M(B,lambda) defining the 'scattering problem' ODE -->
        % d/dx(Phi) = M*Phi with Phi being a 2x2 solution matrix
        dxx = Xuni(2) - Xuni(1);
        DDxStep = StepSize*dxx;
        ppiOne_0 = [1,0; 0,1];
        PPiSS_0 = [1,0; 0,1];
        ppiOne_Factor = [1,0; 0,1];
        npi = length(bb);
        ErrPi11 = zeros(1,npi);
        ppiOne = ppiOne_0;
        ppiSS = PPiSS_0;

        pidx = 1;
        while pidx < npi
            pidx = pidx + 1;
            bb_idx = bb(pidx);

            % Evaluate Propagation Matrix
            ppTmp = MatrixExp(lmda, bb_idx, dxx);
            ppiOne = ppTmp*ppiOne;
            ppiOne_Factor = ppTmp*ppiOne_Factor;
            
            if mod(pidx,StepSize) == 0
                ppTmp = MatrixExp(lmda, bb_idx, DDxStep);
                diff = ppTmp - ppiOne_Factor;
                
                for erridx = pidx-(StepSize-1):pidx
                    ErrPi11(pidx-(StepSize-1):pidx) = diff(1,1);
                end
                
                ppiOne_Factor = ppiOne_0;
                ppiSS = ppTmp*ppiSS;
            end

            % Output to Command Window for debugging:
            % fprintf(['index  = ',num2str(Pidx),'--> det(pi) =',num2str(det(pTmp)),'\n\n']);
        end
        
    end

    function [xL_is_alphaXR] = alpha_XR(ReLmda,ImLmda,bL,bR)
        %% This is an effort to balance the exponentials that arise in implementing 
        %  the Left-Right asymptotics of the scattering functions which are 
        %  required to evaluate the scattering matrix.
        %  Previously the range had been made symmetric with x: -L --> +L
        %  However Exp(-i*lambda*sqrt(l^2 - bo^2)*(-L)) will generally NOT
        %  balance with Exp(i*lambda*sqrt(l^2 - b1^2)*(L))
        %  This is because the terms, 
        %  sqrt(lambda^2 - bL^2) and sqrt(lambda^2 - bR^2), differ
        %  significantly for lambda --> bR
        %
        %  Because of this, if XR is a location on the RHS of the profile,
        %  then the choice for XL that leads to these exponential terms
        %  being balanced is XL = - (xL_is_alphaXR)*XR with xL_is_alphaXR
        %  defined below.
        %
        %  To relate these to the range of the profile as selected by the
        %  user, use XR - XL = XR*( 1 + xL_is_alphaXR) = range which
        %  uniquely determines XR given the range and the value of
        %  xL_is_alphaXR as determined below.
        
        %% Evaluation of xL_is_alphaXR as used in XL = xL_is_alphaXR* XR
        xL_is_alphaXR = sqrt(abs(bR)^2 + ImLmda^2 - ReLmda^2 +sqrt(4*ImLmda^2*ReLmda^2 + (-abs(bR)^2 - ImLmda^2 + ReLmda^2)^2))/sqrt(abs(bL)^2 + ImLmda^2 - ReLmda^2 + sqrt(4*ImLmda^2*ReLmda^2 + (-abs(bL)^2 - ImLmda^2 + ReLmda^2)^2));        
    end

    function [ErrMsg, pctErr, ErrStep, StepsTaken, S_N, S_N_Big, ErrFlag] = S_Adaptive(StepSize2,StepSize, ASTOL, eps_global, Xuni, lmda, bb)
    %% Returns Pi_N using an adaptive step size to enforce maximum relative error constraint 
    % ErrStep = (S_BigStep - S_half_BigStep)/S_BigStep for each step
    % The incremental difference between U evaluated at
    % mod(idx,StepSize) = 0 with StepSize*dx for its distance
    %
    %  dx   dx   dx   dx   dx   dx   dx   dx   dx
    %  1,   2,   3,   4,   5,   6,   7,   8,   9,   ...
    %                 /\                  /\
    %                 ||                  ||
    % So, U(idx4 = '1',4*dx) 'should' equal U(1,dx)*U(2,dx)*U(3,dx)*U(4,dx)
    % 
    % Generally, U(idx*StepSize, StepSize*dx) = U( idx -(StepSize - 1),dx)*U(idx + 1 -(StepSize - 1),dx)* ... *U(idx,dx)
    % So, this function evaluates factors of U between fixed sections of the profile
    % at 2^m steps of 2^(n - m) AND at 2^(m+1) steps of 2^(n - m -1) starting from m = 0
    % It then compares the relative error between these two evaluations and iterates increasing
    % the number of steps (m++) until the relative error is less than the maximum allowed.
    
    %% Set Up

        dxx = Xuni(2) - Xuni(1);
        pi_new = [1,0; 0,1];
        PN = [1,0; 0,1];
        PNBig = [1,0; 0,1];
        ND = length(bb) - 2;% drops the last two points of B-field resulting in ND = 2^n
        STEP = 2*StepSize(1); % attempt to make a big step of this size (2^n with n < log2(Nd))
        STEP2 = 1 + StepSize2(1); % attempt to make a big step of this size (2^n with n < log2(Nd))
        StepsTaken = zeros(1,ND/STEP);
        ErrStep = zeros(1,ND/STEP);
        ErrMsg = ['Relative error exceed at X = '];
        ErrMsggood = ['Relative error not exceeded at any interation.'];
        ErrFlag = 0;
        
   % Evaluate asymptotics for Phi(-) and Phi(+)
        % Used to evaluate Scattering Matrix, S, as Phi(-) = Phi(+)*S
        PhiMinus_Debug = 0;
        PhiPlus_Debug = 0;
        [PhiML, PhiFlag] =    PhiMinus(lmda ,Xuni(1)   , bb(1)  , PhiMinus_Debug);
        [PhiPL_Inv, PhiFlag] = PhiPlus(lmda ,Xuni(end) , bb(end), PhiPlus_Debug);      
        
    %% Begin Error-Controlled Evaluation
        
        for StepIdx = 1:(ND/STEP) % assuming that ND is 2^N and STEP is 2^n with n < N
            UBigStep = MatrixExp(lmda, bb(StepIdx*( 2^STEP2 ) ), ( 2^STEP2 )*dxx);
            diff = 10;
            errorMet = 1;
            bbg2 = bb(1 + STEP*( StepIdx - 1) :STEP*StepIdx); % picking out the section of the B-field corresponding to the current 'STEP'  
            bbg2Iterp = bb(1 + STEP*( StepIdx - 1) :STEP*StepIdx);
            
            for st = 1:length(StepSize2)
                BigStep = STEP;
                LittleStep = 2^StepSize2(st);
                Xg2 = Xuni(1 + STEP*( StepIdx - 1) :STEP*StepIdx); % a section of the grid spanning the current 'BigStep', but at original resolution

                if StepSize2(st) < 0 % this is reached if the resolution of the original B-field was not sufficient to fall below the maximum error tolerance
                    BigStep = 2^(-StepSize2(st))*STEP; 
                   % adjusting # steps of 'BigStep' to match higher resolution
                   % for example if StepSize(st) = -2, the step size is
                   % 2^(-2) = 1/4 of the original grid spacing requiring
                   % the # grid points spanning the same region of x to increase 4-fold
                   % which is implemented by BigStep = 2^(-StepSize(st))*STEP = 2^(-(-2))*STEP = 4*STEP
                    LittleStep = 1;
                    % Of course, as the resolution is only increased the
                    % minimum amount, the first increase will be to double
                    % and in this case, and in all other higher resolution
                    % cases, the finest increment will be used which is a
                    % stepsize of 1.
                    Xg1 = Xuni(1 + STEP*( StepIdx - 1) :STEP*StepIdx); % a section of the original grid using the original defined 'STEP'
                    Xg2 = linspace(Xg1(1),Xg1(end),BigStep);% a grid spanning current 'Step', with resolution increased to 'BigStep'
                    bbg2 = interp1(Xg1,bbg2Iterp,Xg2); % resolving the B-field of current section to the new, finer resolution.
                end
                Unew = pi_new;
                StepsTaken(StepIdx) = 2^StepSize2(st);
                for factor = 1:(2^( STEP2 - StepSize2(st) ) ) % assuming (STEP/StepSize(st)) = 2^m with m = 1, 2, 3, 4, ... MAX
                    idxbbg2 = BigStep - (factor - 1)*( LittleStep );
                    if (BigStep > STEP) && ( idxbbg2 == 1 ) && (StepSize2(st) < -1) && 0
                        FactorMsg = sprintf([char(955),' = ',num2str(lmda),', StepIdx = ',num2str(StepIdx),'/',num2str(ND/STEP),', st = ',num2str(st),', BigStep = ',num2str(BigStep),', idxbbg2 = ',num2str(idxbbg2)]);
%                         FactorMsg
                    end
                    Unew = MatrixExp(lmda, bbg2(  BigStep - (factor - 1)*( LittleStep )  ), ( 2^StepSize2(st) )*dxx)*Unew;
                end
                [PhiMLtmp, PhiFlag]  = PhiMinus(lmda ,Xg2(1)       , bbg2(1)      , PhiMinus_Debug);
                [PhiPL_Invtmp, PhiFlag] = PhiPlus( lmda ,Xg2(BigStep) , bbg2(BigStep), PhiPlus_Debug );  
                S_big = PhiPL_Invtmp*(UBigStep)*PhiMLtmp;
                S_new = PhiPL_Invtmp*(Unew)*PhiMLtmp;
                diff = (S_big - S_new)./S_big;
                [Drow,Dcol] = find( abs(diff) == max(max(abs(diff))) );
%                 if Drow == 1 || Dcol == 1
%                     ['Whaaat?']
%                 end
%                 maxErrElement = diff(Drow,Dcol);
                diffR = real( diff(Drow(1),Dcol(1))     );
                diffI = imag(  diff(Drow(1),Dcol(1))  );
                
%                 if (~islogical(abs(diffR) > ASTOL)) || (~islogical(abs(diffI) > ASTOL))
%                     ['Whaat?']
%                 end
                
                if ( abs(diffR) > ASTOL ) || ( abs(diffI) > ASTOL )
                    if st < length(StepSize2)
                        UBigStep = Unew;
                    end
                    errorMet = -1;                    
                elseif ( abs(diffR) < ASTOL ) && ( abs(diffI) < ASTOL ) && ( sqrt(diffR^2 + diffI^2) < ASTOL )
                    errorMet = +1;
                    break; % break out of 'st' loop and begin next BIGSTEP                 
                end 
                if st == length(StepSize2)
                    ErrFlag = 1; % This is only reached after the smallest step size is evaluated AND the error tolerance has failed.
%                     FactorMsg
                end
            end
            if errorMet == -1
                ErrMsg = [ErrMsg,' ', num2str(Xuni(StepIdx*STEP)),', '];
            end
            ErrStep(StepIdx) = diffR + 1i*diffI;           
            PNBig = UBigStep*PNBig;
            PN = Unew*PN;
%             UBigStep = Unew;
        end
        
        PN=MatrixExp(lmda, bb(length(bb)), dxx)*MatrixExp(lmda, bb(length(bb) - 1), dxx)*PN;
        PNBig=MatrixExp(lmda, bb(length(bb)), dxx)*MatrixExp(lmda, bb(length(bb) - 1), dxx)*PNBig;
        S_N = PhiPL_Inv*PN*PhiML;
        S_N_Big = PhiPL_Inv*PNBig*PhiML;
        pctErrTmp = (S_N_Big - S_N)./S_N_Big*100;
        pctErr = max([pctErrTmp(1,1), pctErrTmp(1,2), pctErrTmp(2,1), pctErrTmp(2,2) ]);
        if ErrFlag == 0
            ErrMsg = ErrMsggood;
        elseif ErrFlag == 1
            if ( abs(pctErr) < (eps_global*100) )
                ErrFlag = 0;
                ErrMsg = ['Relative error exceeded for some iterations, but global error satisfied.']; 
            end
        end
        
    end

    function [S_Adapt1121, ErrFlag, Xuniform, b2nsymm,  RelPctErr, ErrStep, StepsTaken, S_Adaptive_Time, msg_S_terse, msg_S_verbose, SA_ERR] = S_Adaptive_Complete(Xas, Bas, FieldPrep, NFractionSTEPS, lambda_as, eps_as, Rfactor)
 %           [Sval(ii,jj), SFlag, Xh, Bh,  RelPctErr, ErrStep, StepsTaken, SA_Time, msg_S_terse, msg_S_verbose,SA_ERR] =                       S_Adaptive_Complete(Xh, Bh, BF_Prep, NFracSS, l, eps_AS);   
        % Xas and Bas       --> the symmetrized range and B-field selected previously by user
        % FieldPrep         --> 1 if Xas and Bas need to be modified
        %                   --> 0 if Xas and Bas DO NOT need to be modified
        % NFractionSTEPS    --> the number of fractional steps to be attempted
        % lambda_as         --> the eigenvalue used in evaluation. Note, this must be the actual value, not the fraction of the background field
        % eps_as            --> the fractional overall error. 
        %
        % For example, eps_as =.5/100 will produce S_Adapt with S11 having a relative error no more than 0.5 percent
        
        %% Input Values
        Xh = Xas;
        b2nsymm = Bas;        
        NFractionalSteps = NFractionSTEPS; % maximum number of times to double resolution of a 'step' in effort to meet error tolerance requirement for that step.       
        lambda = lambda_as;  % actual, not as fraction of background field
        eps_Divisor = 1;
        
        %% Check that B-field is adequately resolved and prepared
       
        % Symmetric grid
        Lh_symmetric = [-1, 1]*(  Xh(end) - Xh(1)  )/2;           % switch to symmetric boundary conditions
        Nd = length(b2nsymm);
        Xuniform = Xas;
        if FieldPrep
            DoThis = 1;
            Nd_old = Nd;
            Xsym = linspace(  Lh_symmetric(1), Lh_symmetric(2), Nd  );  % 'original' symmetric grid
            Xuniform = linspace(  Lh_symmetric(1), Lh_symmetric(2), Nd  );

            % Evaluate the minimum radius of curvature and set dX = epsilon*(Min Radius)
            NRadTst = floor(Nd/1);
            RR = ones(1,NRadTst)*100;
            XRadTst = linspace(Xuniform(1),Xuniform(end),NRadTst);
            b2nsymm_RadTst = interp1(Xuniform,b2nsymm,XRadTst,'spline');
            [val,imMx]=max(imag(b2nsymm_RadTst));
            RidxMin = imMx - floor(0.05*NRadTst);
            RidxMax = imMx + floor(0.05*NRadTst);
            for idxR = RidxMin:RidxMax
                xR = transpose(XRadTst(idxR - 4:idxR + 4));
                yR = transpose(imag(b2nsymm_RadTst(idxR - 4:idxR + 4)));
                AR = [-2*xR -2*yR ones(length(xR),1)];
                xR = AR\-(xR.^2+yR.^2);
                xo=xR(1);
                yo=xR(2);
                RR(idxR) = sqrt(  xo.^2 + yo.^2  - xR(3));
            end       
            [Rmin, RminIdx] = min(RR);
            dRx = Rfactor*abs(RR(RminIdx));
            dXh = Xh(2) - Xh(1);
            debugRadiusTest = 0;
            if debugRadiusTest
                RRgraph = RR;
                for irdx = 1:NRadTst
                    if abs(RR(irdx)) > 20
                        RRgraph(irdx) = 20;
                    end
                end
                RRgraph = RRgraph/20;
                figure;
                plot(XRadTst,real(b2nsymm_RadTst),'color', [0 1 0]);hold 'on';
                plot(XRadTst,imag(b2nsymm_RadTst),'color', [1 0 0]);hold 'on';
                plot(XRadTst,abs(b2nsymm_RadTst),'color', [0 0 1]);hold 'on';
                plot([XRadTst(imMx) XRadTst(imMx)],[-2 2]);
                plot(XRadTst,RRgraph,'color', [0 0 0]);
            end

            % Find Nd, # grid points, so that dX is small enough to resolve
            % sharpest feature in B-field profile
            PowerNd = floor(log2(length(Bh)) + 1); % default to original resolution, except adjust to smallest 2^n + 2 => Nh
            if (dRx < dXh) && dRx > 0
                factorR = dXh/dRx;
                PowerNd = floor(log2(factorR*length(Bh)) + 1);
            end
            ResolutionBoost = 1;
            Nd = 2^(PowerNd + ResolutionBoost)+2;   % so that Nd - 1 = 2^n + 1

            % Create a uniform grid with Nd points with XL and XR chosen to
            % balance exponential effects arising from boundary conditions
            % whose dominant terms go as Exp[I*lambda*sqrt(lambda^2 -bo^2)*x]. 
            % As the bo is the abs. val. of the field and |b(-L)| = b_left
            % and |b(+L)| = b_right with b_left/b_right > 1 across a shock,
            % if xL = -xR, then the effects of applying bc's will be
            % exponentially asymmetric especially for lambda --> b_right.
            % The relation between XL and XR through the function,
            % alpha_XR, makes these exponential effects balanced and,
            % hopefully, allows the numerical estimation of eigenvalues to
            % be less problematic.
            
            if ( Nd_old > 2^16 ) && (Nd > 2^20 )
                Qmsg = ['Nd = ',num2str(Nd_old),' & will increase to ', num2str(Nd),'!'];
                % Construct a questdlg with three options
                choice = questdlg('Proceed?', ...
                Qmsg, ...
                'Proceed','Do Not Do This','Do Not Do This');
                % Handle response
                switch choice
                    case 'Proceed'
                        disp([choice ' Proceeding.'])
                        DoThis = 1;
                    case 'Do Not Do This'
                        disp([choice ' Keeping old resolution!'])
                        DoThis = 0;

                end
                
            end
            if DoThis == 1
                [xL_is_alphaXR] = alpha_XR(real(lambda),imag(lambda),abs(Bh(1)),abs(Bh(end)));        
                XR = ( Xh(end) - Xh(1))/(1 + xL_is_alphaXR);
                XL = - xL_is_alphaXR*XR;
                Nh = length(Bh);
                Xsym = linspace(XL, XR, Nh);
                Xuniform = linspace(XL, XR, Nd);
    %             Xuniform = linspace(  Lh_symmetric(1), Lh_symmetric(2), Nd  );
                b2nsymm = interp1(Xsym,Bh,Xuniform,'spline'); % User's Selection interpolated to fit symmetric grid with (2^n + 2) data points
            end
        end

        %% Background on Propagator Solution Method
            % Recursively evaluate Quaternion state
            % First evaluate initial state:
            % With:
            % U1 = Exp(M1*Dx)
            % M1 = ( -i*lambda^2     , lambda*b(x1)  )
            %      (lambda*b_star(x1), +i*lambda^2   )
            % then: psi(x1 + DX) = U1*psi(x1)
            % and iterating ...
            % psi(xn + Dx) = Pn*Psi(x1)
            % with Pn =  U(n-1)*U(n-2)* ... *U(1)
            %
            % Quaternion Representation:
            % U(k) = Exp(M(k)*Dx) = Exp(i*Lk*n(k).q) = Cos(Lk) + Sin(Lk)*n(k).q
            % where Lk = lambda*zeta(k) with zeta(k) = sqrt(lambda^2 - abs(bk)^2 )
            % n(k) = (lambda*Re(b(k)),-lambda*Im(b(k)),-i*lambda^2)/(i*lambda*zeta(k)) 
            % is a unit vector defining the axis of rotation
            %
            % with q = (I, J, K) with I, J and K being quaternions which could be defined as i times
            % the respective Pauli spin matrices or, formally, by their
            % definition: I^2 = J^2 = K^2 = I*J*K = -1, so that IJ = K, JK = I,
            % KI = J and JI = -IJ etc, that is I, J and K anticommute
            %    

        %% Prepare Boundary Conditions and Adaptive Step Parameters       

            % Evaluate asymptotics for Phi(-) and Phi(+)
            % Used to evaluate Scattering Matrix, S, as Phi(-) = Phi(+)*S
            PhiMinus_Debug = 0;
            PhiPlus_Debug = 0;
            [PhiML, PhiFlag] =    PhiMinus(lambda ,Xuniform(1)   , b2nsymm(1)  , PhiMinus_Debug);
            [PhiPL_Inv, PhiFlag] = PhiPlus(lambda ,Xuniform(end) , b2nsymm(end), PhiPlus_Debug);

            NumberSteps = 2^floor(log2(length(b2nsymm))/2 -1);
            dimOfStepSize = log2(floor((Nd+1)/NumberSteps))-1;% 'StepSize' holds the biggest TWO STEP to be taken
            StepSize = zeros(1,dimOfStepSize);
            for iss = 1:( 1 + dimOfStepSize)
                StepSize(iss) = 2^(dimOfStepSize - ( iss - 1) );
            end

            for iss = 1:( 1 + NFractionalSteps + dimOfStepSize)
                StepSize2(iss) = (dimOfStepSize - ( iss - 1) );
            end

            step = 2*StepSize(1);
            step2 = 1 + StepSize2(1);

            mu = real(lambda);nu = imag(lambda);
            
        %% Set error tolerance per step
            if ( mu >= 0.2 && nu >= 0.2) || (mu >= 1 && nu >= 0.1) %&& sqrt(mu^2 + nu^2) > 0.1
                AdaptiveStepTolerance = sqrt( 1 - ((1 - eps_as^2))^(2/NumberSteps)  )/eps_Divisor; 
                %% A value of eps_Divisor = 50 helps speed up evaluation for
                % the range of complex lambda used in the condition above.
                % For eps = 0.5/100 and NumberSteps = 512 a value of
                % eps_Divisor = 455 makes this 'efficient' tolerance
                % equivalent to the stricter tolerance after the else. 
                % It should be noted that the appropriateness of this
                % cut-off, or the value of eps_Divisor, may not work for
                % profiles other than that used in development.
            else
                AdaptiveStepTolerance = ( (1 + eps_as))^(1/NumberSteps) - 1;   

                %% Notes on Evaluation of AdaptiveStepTolerance
                % Initial observations show that complex eigenvalues lead to much faster convergence. 
                % The reason for this is not comopletely clear. There may well be exceptions.
                % This seems to provide an upper bound on overall error.
                % Assumes the relation between error of a single step, eps_n, is
                % related to the error of multiple factors to be ...
                % (1 + eps_n)^NSteps = 1 + eps
                % Given a desired maximum error, eps, AdaptiveStepTolerance as evaluated after the
                % 'else' gives the corresponding maximum error per step.
                %
                % NOTE: if the signs of eps_n are 1/2 positive and 1/2 negative, then the resulting overall error
                % will go as ( 1 - eps_n^2)^(NSteps/2) giving eps_n = (+/-) sqrt( (1 + eps)^(2/NSteps) - 1 )
                %
                % Example: if eps = 0.5/100 and NSteps = 512
                % if all eps_n > 0 then eps_n = 9.743e-06
                % if 1/2 eps_n > 0 and 1/2 are < 0, then eps_n = (+/-) 0.004
                % ==> Presumably then the more rapid convergence for complex eigenvalues is related to
                % the more random alternation of signs of step-wise errors. It is not clear why this would be so. 
                % Given this, it is not known how to build this into a more efficient evaluation
                % of eps_n
            end
            
        %% Call S_Adaptive
            StartTime = cputime;   
            [ErrMessage, RelPctErr, ErrStep, StepsTaken, S_Adapt, S_Adapt_Big, SA_ERR] = S_Adaptive(StepSize2,StepSize, AdaptiveStepTolerance, eps_as, Xuniform, lambda, b2nsymm);
            S_Adaptive_Time = cputime-StartTime;
            S_Adapt1121 = S_Adapt*[1,0]';
            ErrFlag = prod(isnan(S_Adapt1121)<1) > 0; % if S_Adapt1121 all numeric, then isnan(S_Adapt1121) returns 0's
                                                   % if so,isnan(S_Adapt1121)<1, returns 1's
                                                   % then prod(isnan(S_Adapt1121)<1)returns 1 if 
                                                   % -> all entries are numeric and 0 if even a single entry is NaN,
                                                   % so, prod(isnan(S_Adapt1121)<1) > 0 
                                                   % returns logical 1 if every entry is numeric 
                                                   % and logical 0 if there is an NaN entry
                                                   % (to be consistent with twisted usage with RealEV's, all is good
            
        %% Prepare Diagnostic Messages
            maxErr = max(ErrStep);
%             clc;
            msg_S_adaptive1 = sprintf('*****************************************************************************************\n\n');
            msg_S_adaptive1b = sprintf(['For ',char(955), ' = ',num2str(lambda),':\n\n']);
            msg_S_adaptive2 = sprintf(['Final relative percent error in S  = ',num2str( RelPctErr  ),' time = ',num2str(S_Adaptive_Time),newline,newline]);
            msg_S_adaptive3 = sprintf(['Min. (step size, #steps) = ( ',num2str(  min(StepsTaken)  ),', ',num2str(sum(StepsTaken==min(StepsTaken))),')  ']);
            msg_S_adaptive4 = sprintf(['Max. (step size, #steps) = ( ',num2str(  max(StepsTaken)  ),', ',num2str(sum(StepsTaken==max(StepsTaken))),')\n\n']); 
            msg_S_adaptive5 = sprintf(['Avg. step size = ',num2str(  mean(StepsTaken)  ),', prod(1 + ErrStep)-1 = ',num2str(prod(1 + ErrStep)-1),'\n\n']);  
            msg_S_adaptive6 = sprintf(['Percent Error from iterative errors = 100*( prod(1 + ErrStep) -1 ) = ',num2str(  100*( prod(1 + ErrStep) -1 )  ),'\n\n']);   
            msg_S_adaptive7 = sprintf(['Mean error in ErrStep = ',num2str( 10^(floor(log10(abs(real(mean( ErrStep )))))) ),'(',num2str(10^(-floor(log10(abs(real(mean( ErrStep ))))))*mean(ErrStep)),...
                '), Standard deviation of error = ',num2str( 10^(floor(log10(abs(std(real( ErrStep )))))) ),'{',num2str( 10^(-floor(log10(abs(std(real( ErrStep ))))))*std(real(ErrStep)) ),' + i(',num2str( 10^(-floor(log10(abs(std(real( ErrStep ))))))*std(imag(ErrStep)) ),')},\n\n']);
            msg_S_adaptive8 = sprintf(['# iterative errors with real part (above/below) zero = (',num2str(  sum(real(ErrStep) > 0 )),', ',num2str(  sum(real(ErrStep) < 0 )),') ... ',num2str(abs((sum(real(ErrStep) < 0) - length(ErrStep)/2))/(sqrt(length(ErrStep))/2)),'*sigma \n\n']);
            msg_S_adaptive9 = sprintf(['# iterative errors with imaginary part (above/below) zero = (',num2str(  sum(imag(ErrStep) > 0 )),', ',num2str(  sum(imag(ErrStep) < 0 )),') ... ',num2str(abs((sum(imag(ErrStep) < 0) - length(ErrStep)/2))/(sqrt(length(ErrStep))/2)),'*sigma \n\n']);
            msg_S_adaptive10 = sprintf([ErrMessage,'\n\n']);

            msg_S_adaptive11 = sprintf('*****************************************************************************************\n\n');
            msg_S_verbose = sprintf([msg_S_adaptive1,msg_S_adaptive1b,msg_S_adaptive2,msg_S_adaptive3,msg_S_adaptive4,msg_S_adaptive5,msg_S_adaptive6,msg_S_adaptive7,msg_S_adaptive8,msg_S_adaptive9,msg_S_adaptive10,msg_S_adaptive11]);

            msg_S_terse = sprintf([msg_S_adaptive1b,msg_S_adaptive10]);
        
    end
%*************************************************************
%       END FUNCTION DEFINITIONS

%*************************************************************
%% Set UP
%   timing
        t = cputime;   
    %close all figures that may be open
    close all;
    thisFile = mfilename();
    ['From: ',thisFile,' Begin Set Up ']
    
    BStatus = 1;
        
    %Folder Paths Relative to Archive path
        TrialFolderRel = ['\Trial ',TrialNumber];%e.g. '\Trial 11'

    %Absolute Folder Paths
        TrialFolderAbs = [ArchiveBase,TrialFolderRel];%e.g. 'F:\ ...\Trial 11\'
        
        FigRel = ['\Figures'];
        JpgRel = ['\JPEG'];
        
        FigAbs = [TrialFolderAbs,FigRel];%e.g. 'F:\ ...\Trial 11\'
        JpgAbs = [TrialFolderAbs,JpgRel];%e.g. 'F:\ ...\Trial 11\'
        
        ProfileFigBase = [FigAbs,'\Profile_T',TrialNumber,'_R',RunNumber,'_L'];%e.g. 'F:\ ...\Trial 11\'
        ProfileJpgBase = [JpgAbs,'\Profile_T',TrialNumber,'_R',RunNumber,'_L'];%e.g. 'F:\ ...\Trial 11\'
        
        HodoFigBase = [FigAbs,'\Hodo_T',TrialNumber,'_R',RunNumber,'_L'];%e.g. 'F:\ ...\Trial 11\'
        HodoJpgBase = [JpgAbs,'\Hodo_T',TrialNumber,'_R',RunNumber,'_L'];%e.g. 'F:\ ...\Trial 11\'
        
        AviFolderRel = [TrialFolderRel,'\avi ',TrialNumber];%e.g. '\Trial 11\avi 11' folder
        AviFolderAbs = [ArchiveBase,AviFolderRel];%e.g. 'F:\ ...\Trial 11\avi 11'
        
%         if strfind(ArchiveBase, 'Kayla') > 0
%             BFieldFile = [TrialFolderAbs,'\Bfields\Bfield_',TrialNumber,'_',RunNumber];
%         else
%             BFieldFile = [TrialFolderAbs,'\Bfields\Bfield_',TrialNumber,'_',RunNumber,'.txt'];%e.g. 'F:\ ...\Trial 11\'      
%         end
        
        
        BFieldFile = [TrialFolderAbs,'\Bfields\Bfield_',TrialNumber,'_',RunNumber,'.txt'];%e.g. 'F:\ ...\Trial 11\'   
        BfieldMatfile = [TrialFolderAbs,'\Bfields\Bfield_',TrialNumber,'.mat'];%e.g. 'F:\ ...\Trial 11\'   
        
        %Functions to invoke:
        BoolConKeys = BoolCon.keys;
        BoolConValues = BoolCon.values;
        
        boolDevelopment = BoolCon('boolDevelopment');
        boolSpeedCheck = BoolCon('boolSpeedCheck');
        
        boolPeakOscAnalysis = 0;
        POAIdx = find(contains(BoolConKeys,'boolPeakOscAnalysis'),1);
        if POAIdx > 0
            boolPeakOscAnalysis = BoolCon(BoolConKeys{POAIdx});
        end
                        
        boolIz = BoolCon('boolIz');
        boolBgraph = BoolCon('boolBgraph');
        boolHodogram = BoolCon('boolHodogram');
        boolVideo = BoolCon('boolVideo');
        if boolVideo == 2
            boolVideo = 1;
            boolmp4 = 1;
        else
            boolmp4 = 0;
        end
        boolHodogramVideo = BoolCon('boolHodogramVideo');
        boolREV = BoolCon('boolREV');
        EVZoom = BoolCon('EVZoom');
        boolEnergy = BoolCon('boolEnergy');
        boolFE = BoolCon('boolFE');
        boolConservedQuantities = BoolCon('boolConservedQuantities');
        boolRoot = BoolCon('boolRoot');
        boolCEV = BoolCon('boolCEV');
        boolCEVgraph = BoolCon('boolCEVgraph');
        boolRAD = BoolCon('boolRAD');
        boolRADgraph = BoolCon('boolRADgraph');
        NMonitors = BoolCon('NMonitors');
        evaluateKBWShockBoundaries = 0; %default value. Reset when parameters are read from excel file ... if ... KBW profile
        CompressionRatio = 0.8586;
        BZasymptoteR = 0;
        BZasymptoteL = 0;
        Ntest = 10;
        tol = 0.002;
        
        LRWin = EvalCon('LRWin');
        LREval = EvalCon('LREval');
        NREval = EvalCon('NREval');
        NBSIter = EvalCon('NBSIter');
        zinit = EvalCon('zinit');
        LComplexRWin = EvalCon('LComplexRWin');
        LComplexEval = EvalCon('LComplexEval');
        NCEval = EvalCon('NCEval');
%% User Declared Parameters
        

        
%         boolREV = 1; % set to 1 to evaluate real eigenvalues using
%         EVZoom = 1; %  = '+1' if user is to chose a smaller section of B-field for Real EV evaluation
%         boolRoot = 1; %Use Newton's Method to estimate complex eigenvalue
%         boolCEV = 0; % evaluate complex eigenvalues
%         boolCEVgraph = 0; % plot sign matrix for s11 in complex plane
%         boolRAD = 0; % evaluate continuous scattering data
%         boolRADgraph = 0; % plot continous scattering data
        %Default options passed onto ODE45
%         opts = odeset('RelTol',1e-6,'AbsTo',1e-8);
        opts = odeset('RelTol',1e-6,'AbsTo',1e-8);
        
        theta_M = 0;%default values for B-field direction upstream, theta_M, and downstream, theta_P
        theta_P = 0;% changed to archived values for Helical profile
            
%         %Evaluate real eigenvalues  
%         LRWin = [0.01, 0.99]; %future development: allow graphical display of region w/o evaluating entirely
%         LREval = [.7, 0.99]; %range of real values to study
%         NREval = 100;             % number of points in range to evaluate
%         NBSIter = 10;    % '6' appears satisfactory and should return lr with accuracy of ...
%                         %     (lmax - lmin)/NREval/2^6 = (1 - 0)/50/64 ~ +/- 0.0003
        
%         %Evaluate complex eigenvalues
%         LComplexRWin = [0.3, 0.99; -.9, 0.9]; % future development ... as with LRWin
%         LComplexEval = [0.3, 0.99; -.02, 0.7]; %range of real values to study ... [Real Min, Real Max; Im Min, Im Max]
%         NCEval = 20;  % Max # of intervals to evaluate along either Re(lambda) or Im(lambda) directions
        
        %Parameters used in evaluation of continuous scattering data
        LReBranchEval = [1.01, 2]; %range along lambda_real ... end point chosen so Rho --> 0 there
        NReBranch = 500;             % # of data points to evaluate ... sufficient to resolve function
        RBeps = 0.0001;             % offset above real axis
        LImBranchEval = [0.8, 0.01];  %range along imaginary axis ... the 1st value should be high enough that Rho --> 0
        NImBranch = 500;             % # of points to evaluate along imaginary axis ... sufficient to resolve function
        IBeps = 0.0001;             % offset of contour to right of imaginary axis
        
        %Declarations
        CEVList = cell(1,3);
        ReEVs = zeros(1,1);
        RealEV_List = [];

%% Read B-field Line Deprecated  
% ['From: ',thisFile,' Begin Read B-field Line ...  ']
%     [B4, BFlines] = FindBfieldLines(BfieldMatfile, BFieldFile, LineNumber, N);

%     if exist(BfieldMatfile,'file') && 0
%         ['From: ',thisFile,' Begin Read B-field Line. Trying Matfile Archive  ']
%         BFMat = matfile(BfieldMatfile);
%         BfieldNames = fieldnames(BFMat);
%         idxRun = strfind(BfieldNames,'Run');
%         if str2num(RunNumber) <= length(BFMat.PStr(:,1))
%             runN = ['Run',RunNumber];
%         end
%         if sum(~cellfun('isempty',strfind(BfieldNames,runN))) > 0
%             if LineNumber <= length(BFMat.(runN)(:,1))              
%                 B4 = BFMat.(runN)(LineNumber,:);
%             else
%                 ['Line # not in Matfile']
%                 length(BFMat.(runN)(:,1))
%                 return;
%             end
%         else
%             ['Run # not in file/nr']
%             ['Field Names:']
%             fieldnames(BFMat)
%             return;
%         end
%     else 
%         %Read all lines from Bfield file for given TrialNumber and RunNumber
%         ['From: ',thisFile,' Begin Read B-field Line. Trying text file archive']
%             BFh = fopen(BFieldFile);
%             if BFh < 1
%                 [BFieldFile, ' does not exist.']
%                 return
%             end
%             BFline = fgetl(BFh);
%             BFlines = cell(0,1);
%             NumOfLines = 0;
%             while ischar(BFline)
% %                 BFlines{end+1,1} = BFline;
%                 BFline = fgetl(BFh);
%                 if strfind(BFline,'LineStart')>0
%                     NumOfLines = NumOfLines + 1;
%                     BFlines{end+1,1} = fgetl(BFh);
% %                     clc;
%                     ['From: ',thisFile,' verifying B-field Line # ',num2str(NumOfLines)]
%                 end
%             end
%             fclose(BFh);
% 
%         %Find Bfield profile for Line n  
%             LineStartN = ['LineStart',num2str(LineNumber)];
%             numlines = length(BFlines); %returns number of lines in BFh
%             IndexBF = strfind(BFlines, LineStartN); %returns line # of flag before Bfield Line n
%             lineN = find(not(cellfun('isempty',IndexBF))) + 1; %gives line # of nth Bfield Line
%             if LineNumber < (NumOfLines + 1)
% %                 clc;
%                 ['From: ',thisFile,' loading B-field Line # ',num2str(LineNumber)]
%                 temp = BFlines(LineNumber);
%                 C = textscan(temp{1},'%f');
%                 B4 = C{1}.'; %B is an N/4 complex valued, double precision array
%             else
%                 ['Line # ', num2str(LineNumber),' is greater than the ', num2str(NumOfLines), ' lines in the file.' ]
%                 return;
%             end
%     end
%     if length(B4)> (N/4)
%         B4temp = B4(1:round(length(B4)/(N/4)):end);
%         B4 = B4temp;
%     end
%             B4 = B4*exp(1i*0.2*pi)    ;
%% Read Parameters from Excel File
% clc;
['From: ',thisFile,' Begin reading parameters ...  ']
    %Read Excel File with Trial and Run Parameters
        % Extract PStr: the structure with members such as Physics, Numerics etc.
        % with each member containing {keys, values} of parameters describing a given trial & run  
        BFMatWorks = 1; % Assuming BfieldMatfile exists and PSTr is a variable
        if exist(BfieldMatfile,'file')
            ['From: ',thisFile,' Begin reading parameters from Matfile Archive  ']
            BFMat = matfile(BfieldMatfile);
            varList = who(BFMat);
            
            if ismember('PStr',varList) 
                PStr = BFMat.PStr(str2num(RunNumber),1);
            else
                BFMatWorks = 0;
            end
        end
        
        if ~BFMatWorks
            ['From: ',thisFile,' Begin reading parameters from Excel Archive  ']
            Efile = [TrialFolderAbs,'\*.xlsx'];%wild card for files in Trial folder with '.xls' extension
            RelXLfile = dir(Efile);% returns, presumably, the single file with '.xlsx' extension
            XLfile = [TrialFolderAbs,'\',RelXLfile.name];%absolute address of excel file with trial parameters

            %Read sheet 'RunNumber' of the Excel file 'XLfile' and create
            %'PStr' ... a  Parameter structure defined as:
            % Pstr = struct('Dcon',Diagnostics, 'Fcon', Folder,'GCon',
            % Graphics, 'NCon',Numerics,'PCon',Physics,'TCon',Trial)
            % To find the value of 'L' stored in the Numerics container use:
            % L = Pstr.NCon('L') ... for example

            [PStr] = DNLS_Excel_Read_v3(XLfile, RunNumber);
        end

        %********************************************************************
        % Extract Parameter Containers by category 
        Physics     = PStr.PCon;
        Prfl        = PStr.PrflCon;
        Numerics    = PStr.NCon;
        Graphics    = PStr.GCon;
        Trial       = PStr.TCon;
        Diagnostics = PStr.DCon;
        Profile = Prfl('Profile');
                 
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
        TrialNumber = PStr.TCon('TrialNumber');% might not be in Trial container? ... but it is in the excel sheet???
        RunNumber = PStr.TCon('RunNumber');% might not be in Trial container? ... but it is in the excel sheet???
        TrialName = PStr.TCon('TrialName');
        TimeStamp = PStr.TCon('TimeStamp');   
        if strcmp(Profile, 'double_1psoliton')
            speed1 = Prfl('speed1');
            speed2 = Prfl('speed2');
            speed = speed1;
            if speed2 > speed1
                speed = speed2;
            end
        else
                    speed = Prfl('speed');
        end


        LineMax = PStr.DCon('timesteps');
        
        %Load Physics parameters based on Profile type
        EVL = -L; EVR = +L; % default to full profile
        if EVZoom && boolREV && (  LineNumber > LNstart  )
            % if evaluating lines iteratively, for the 2nd and subsequent lines
            % use LHS and RHS limits found on the first line evaluated
            EVL = XLeft; EVR = XRight; 
        end

        % update LHS and RHS for next lines
        EVLeft_return = EVL; EVRight_return = EVR;
        EVLeft_return = EVL; EVRight_return = EVR;

%% Read B-field Line  
    ['From: ',thisFile,' Begin Read B-field Line ...  ']
    [B4, BFlines] = FindBfieldLines(BfieldMatfile, BFieldFile, LineNumber, N);
    B4end = B4;
    if EVZoom && boolREV && (  LNend > LNstart  ) && (  LineNumber ==  LNstart )
        [B4, BFlines] = FindBfieldLines(BfieldMatfile, BFieldFile, LNend, N); 
        [B4end, BFlines] = FindBfieldLines(BfieldMatfile, BFieldFile, LineNumber, N);         
    end
    
    X4 = linspace(-L,L,length(B4));
    dx4 = X4(2) - X4(1);
        
    if strcmp(Profile,'KBW') || strcmp(Profile, 'Generated13IS') || contains(Profile, 'HybridIS') || sum(strcmp(Prfl.keys,'CR13'))
        if sum(strcmp(Prfl.keys,'theta_kbw')) 
            theta_kbw = Prfl('theta_kbw');
            bkbw = Prfl('b_kbw');
            xkbw = Prfl('Xkbw');
        elseif sum(strcmp(Prfl.keys,'theta_IS')) 
            theta_kbw = Prfl('theta_IS');
            bkbw = abs(Prfl('CR13'));
            xkbw = Prfl('xIS');
        end
        
        b_L = Prfl('b_L');
        
        evaluateKBWShockBoundaries = 1;
        Prfl        = PStr.PrflCon;
        if sum(strcmp(Prfl.keys,'Compression_Ratio')) < 1
            CompressionRatio = 0.31415926;
        else
            CompressionRatio = Prfl('Compression_Ratio');     
        end
        
        if (sum(strcmp(Prfl.keys,'CR13')) < 1) && (sum(strcmp(Prfl.keys,'CR23')) < 1)
            CR13 = -0.31415926;
            CR23 = 1-CR13;
        else
            CR13 = Prfl('CR13'); 
            CR23 = Prfl('CR23'); 
        end 
            
        TrialData = split(TrialID,'_');
        TrialTxt = sprintf('TrialID %s_{%s}',TrialData{1},TrialData{2});
                  
        HBflag = 1;
        idL1 = 1;idR1 = length(B4);
           
            
        HBfield = figure('Color', [1, 1, 1],'units','normalized','outerposition',[(NMonitors -1 +.1) 0.1 0.8 0.8]);
        cla();
        xIdxOffsetTxt = floor(0.05*(idR1 - idL1));
        if evaluateKBWShockBoundaries
            %Following assumes Shock is in LHS of profile
            [idRelMax, idxShockL, idxShockR, XShockL, XShockR, Izcalc] = IS_ShockBoundaries(B4(1:floor(length(B4)/2)), X4(1:floor(length(B4)/2)), b_L, BZasymptoteL, BZasymptoteR, Ntest, tol);
            
            IS_Shock_Profile_Label(HBfield, zmin, zmax, xIdxOffsetTxt, idxShockL, idxShockR, XShockL, XShockR, X4, B4)
            
        end
        
        Figure_Label(HBfield, zmax, xIdxOffsetTxt, idL1, idR1, X4, B4)
        
        plot(X4(idL1:idR1),real(B4(idL1:idR1)),'Color',[0 , 1, 0]);
        hold('on');
        plot(X4(idL1:idR1),imag(B4(idL1:idR1)),'Color',[1 , 0, 0]);
        hold('on');
        plot(X4(idL1:idR1),abs(B4(idL1:idR1)),'Color',[0 , 0, 1]);
        hold('on');
        
        xL1 = X4(idL1); xR1 = X4(idR1);
        X0 = xL1 + 0.1*(xR1 - xL1)/2; X1 = xL1 + 0.3*(xR1 - xL1)/2;
        Y1 = 0.94*zmax; Y2 = 0.92*zmax; Y3 = 0.85*zmax; Y4 = 0.75*zmax;

        text(X0,Y1,sprintf('Profile = %s',Profile));
        text(X1,Y2,sprintf('Trial ID = %s',TrialTxt));
        text(X0,Y3,sprintf('time = %4.4f',TimePerLine*(LineNumber - 1)));
        text(X0,Y4,sprintf('line # = %4.0f',LineNumber));


        axis([X4(idL1) X4(idR1) zmin zmax]);
        drawnow(); 
            
            
            
            
            choice = questdlg('Pick Compression Ratio:', ...
            'Compression ratio:', ...
            'Use |CR13| < = 0.5?','Use |CR23| = > 0.5?','No.','No.');
            % Handle response
            switch choice
            case 'Use |CR13| < = 0.5?'
                CompressionRatio = CR13;               
                disp([choice 'Using CR13'])
            case 'Use |CR23| = > 0.5?'
                CompressionRatio = CR23;
                disp([choice 'Using CR23'])
            case 'No.'
                disp('Moving on.')
            end
            close all;
        
        
        BZasymptoteR = abs(CompressionRatio)*b_L*sin(theta_kbw);
        BZasymptoteL = 0;
    end
        
    if strcmp(Profile,'Helical') % EVL and EVR are intended to select just the 'helical' section of profile
            EVL = Prfl('EVL');
            EVR = Prfl('EVR');
            helix_angle = Prfl('helix_angle');
            theta_M = 0;
            theta_P = helix_angle;
        end

%% Select Section of Profile
['From: ',thisFile,' Select Section of Profile to Examine  ']
        if (EVZoom || boolIz) && (  LNstart == LNend  )
%             IzZoom1 = sum(imag(B4(1:1070)))*dx4; % idx = 1070 was found in one trial to be the right edge of the IS
                                                % NOTE: for non-coplanar
                                                % profiles, bz != 0 to the
                                                % right of the IS, so Iz is
                                                % NOT well defined ...
                                                % without a fixed upper
                                                % limit
            TrialData = split(TrialID,'_');
            TrialTxt = sprintf('TrialID %s_{%s}',TrialData{1},TrialData{2});
                      
            HBflag = 1;
            idL1 = 1;idR1 = length(B4);
            while HBflag == 1                
                
            HBfield = figure('Color', [1, 1, 1],'units','normalized','outerposition',[(NMonitors -1 +.1) 0.1 0.8 0.8]);
            cla();
            xIdxOffsetTxt = floor(0.05*(idR1 - idL1));
            if evaluateKBWShockBoundaries
                %Following assumes Shock is in LHS of profile
                [idRelMax, idxShockL, idxShockR, XShockL, XShockR, Izcalc] = IS_ShockBoundaries(B4(1:floor(length(B4)/2)), X4(1:floor(length(B4)/2)), b_L, BZasymptoteL, BZasymptoteR, Ntest, tol);
                
                IS_Shock_Profile_Label(HBfield, zmin, zmax, xIdxOffsetTxt, idxShockL, idxShockR, XShockL, XShockR, X4, B4)
                
            end
            
            Figure_Label(HBfield, zmax, xIdxOffsetTxt, idL1, idR1, X4, B4)
            
            plot(X4(idL1:idR1),real(B4(idL1:idR1)),'Color',[0 , 1, 0]);
            hold('on');
            plot(X4(idL1:idR1),imag(B4(idL1:idR1)),'Color',[1 , 0, 0]);
            hold('on');
            plot(X4(idL1:idR1),abs(B4(idL1:idR1)),'Color',[0 , 0, 1]);
            hold('on');
            
            xL1 = X4(idL1); xR1 = X4(idR1);
            X0 = xL1 + 0.1*(xR1 - xL1)/2; X1 = xL1 + 0.3*(xR1 - xL1)/2;
            Y1 = 0.94*zmax; Y2 = 0.92*zmax; Y3 = 0.85*zmax; Y4 = 0.75*zmax;

            text(X0,Y1,sprintf('Profile = %s',Profile));
            text(X1,Y2,sprintf('Trial ID = %s',TrialTxt));
            text(X0,Y3,sprintf('time = %4.4f',TimePerLine*(LineNumber - 1)));
            text(X0,Y4,sprintf('line # = %4.0f',LineNumber));
            txalert = text(L*0,0.9*zmax,sprintf('Left click to select region boundaries'));
            txalert.FontSize = 14;
            txalert.FontWeight = 'bold';

            axis([X4(idL1) X4(idR1) zmin zmax]);
            drawnow(); 
                
                
                [x,y] = ginput(2);
                x = sort(x);
                xL1 = x(1);%xL2 = x(2);
                xR1 = x(2);%xR2 = x(4);
                figure(HBfield);
                plot([xL1 xL1],[-1.5 1.5],'Color',[1 0 0]);
                hold('on');
%                 plot([xL2 xL2],[-1.5 1.5],'Color',[1 0 0]);
%                 hold('on');
                plot([xR1 xR1],[-1.5 1.5],'Color',[0 0 1]);
                hold('on');
%                 plot([xR2 xR2],[-1.5 1.5],'Color',[0 0 1]);
%                 hold('on');

                axis([X4(idL1) X4(idR1) zmin zmax]);
                drawnow();

                % Construct a questdlg with three options
                choice = questdlg('Selection correct?', ...
                'Select Left-Right boundary of shock:', ...
                'Yes','Try Again','Zoom In','Get Me Outa Here');
                % Handle response
                switch choice
                    case 'Try Again'
                        disp([choice ' Try, try again.'])
                        HBflag = 1;
                    case 'Yes'
                        disp([choice ' Woot!'])
                        HBflag = 2;
                    case 'Zoom In'
                        disp([choice ' Woot!'])
                        tmpL1 = abs(X4 - xL1);
                        tmpR1 = abs(X4 - xR1);

                        [val, idL1] = min(tmpL1);
                        [val, idR1] = min(tmpR1);

                        EVL = X4(idL1);EVR = X4(idR1);
                        HBflag = 1;
                    case 'Get Me Outa Here'
                        disp('Prepare for Transport.')
                        HBflag = -1;
                end
                        close(HBfield);
            end
    

            tmpL1 = abs(X4 - xL1);
            tmpR1 = abs(X4 - xR1);

            [val, idL1] = min(tmpL1);
            [val, idR1] = min(tmpR1);

            EVL = X4(idL1);EVR = X4(idR1);
            
        end
        if EVZoom && boolREV && (  LNend > LNstart  ) && (  LineNumber ==  LNstart )
            TrialData = split(TrialID,'_');
            TrialTxt = sprintf('TrialID %s_{%s}',TrialData{1},TrialData{2});
                      
            HBflag = 1;
            idL1 = 1;idR1 = length(B4);
            while HBflag == 1                
                
            HBfield = figure('Color', [1, 1, 1],'units','normalized','outerposition',[(NMonitors -1 +.1) 0.1 0.8 0.8]);
            cla();
            xIdxOffsetTxt = floor(0.05*(idR1 - idL1));
            
            Figure_Label(HBfield, zmax, xIdxOffsetTxt, idL1, idR1, X4, B4)
            
            plot(X4(idL1:idR1),real(B4(idL1:idR1)),'Color',[0 , 1, 0]);
            hold('on');
            plot(X4(idL1:idR1),imag(B4(idL1:idR1)),'Color',[1 , 0, 0]);
            hold('on');
            plot(X4(idL1:idR1), abs(B4(idL1:idR1)),'Color',[0 , 0, 1]);
            hold('on');

            plot(X4(idL1:idR1),real(B4end(idL1:idR1)),'Color',[0 , 1, 0],'LineStyle', '--');
            hold('on');
            plot(X4(idL1:idR1),imag(B4end(idL1:idR1)),'Color',[1 , 0, 0],'LineStyle', '--');
            hold('on');
            plot(X4(idL1:idR1), abs(B4end(idL1:idR1)),'Color',[0 , 0, 1],'LineStyle', '--');
            hold('on');
            
            xL1 = X4(idL1); xR1 = X4(idR1);
            X0 = xL1 + 0.1*(xR1 - xL1)/2; X1 = xL1 + 0.3*(xR1 - xL1)/2;
            Y1 = 0.94*zmax; Y2 = 0.92*zmax; Y3 = 0.85*zmax; Y4 = 0.75*zmax;

            text(X0,Y1,sprintf('Profile = %s',Profile));
            text(X1,Y2,sprintf('Trial ID = %s',TrialTxt));
            text(X0,Y3,sprintf('time = %4.4f',TimePerLine*(LineNumber - 1)));
            text(X0,Y4,sprintf('line # = %4.0f',LineNumber));
            txalert = text(L*0,0.9*zmax,sprintf('Left click to select region boundaries'));
            txalert.FontSize = 14;
            txalert.FontWeight = 'bold';

            axis([X4(idL1) X4(idR1) zmin zmax]);
            drawnow(); 
                
                
                [x,y] = ginput(2);
                x = sort(x);
                xL1 = x(1);%xL2 = x(2);
                xR1 = x(2);%xR2 = x(4);
                figure(HBfield);
                plot([xL1 xL1],[-1.5 1.5],'Color',[1 0 0]);
                hold('on');
%                 plot([xL2 xL2],[-1.5 1.5],'Color',[1 0 0]);
%                 hold('on');
                plot([xR1 xR1],[-1.5 1.5],'Color',[0 0 1]);
                hold('on');
%                 plot([xR2 xR2],[-1.5 1.5],'Color',[0 0 1]);
%                 hold('on');

                axis([X4(idL1) X4(idR1) zmin zmax]);
                drawnow();

                % Construct a questdlg with three options
                choice = questdlg('Selection correct?', ...
                'Select Left-Right boundary of shock:', ...
                'Yes','Try Again','Zoom In','Get Me Outa Here');
                % Handle response
                switch choice
                    case 'Try Again'
                        disp([choice ' Try, try again.'])
                        HBflag = 1;
                    case 'Yes'
                        disp([choice ' Woot!'])
                        HBflag = 2;
                    case 'Zoom In'
                        disp([choice ' Woot!'])
                        tmpL1 = abs(X4 - xL1);
                        tmpR1 = abs(X4 - xR1);

                        [val, idL1] = min(tmpL1);
                        [val, idR1] = min(tmpR1);

                        EVL = X4(idL1);EVR = X4(idR1);
                        HBflag = 1;
                    case 'Get Me Outa Here'
                        disp('Prepare for Transport.')
                        HBflag = -1;
                end
                        close(HBfield);
            end
    

            tmpL1 = abs(X4 - xL1);
            tmpR1 = abs(X4 - xR1);

            [val, idL1] = min(tmpL1);
            [val, idR1] = min(tmpR1);

            EVL = X4(idL1);EVR = X4(idR1);
            EVLeft_return = EVL; EVRight_return = EVR;
            
        end
                
        B0 = PStr.PCon('bo');
        nu = PStr.PCon('nu');
        if sum(strcmp(PStr.PCon.keys,'alpha')) < 1
            alpha = 1;
        else
            alpha = PStr.PCon('alpha');     
        end
        if sum(strcmp(PStr.PCon.keys,'dispersion')) < 1
            dispersion = 1;
        else
            dispersion = PStr.PCon('dispersion');     
        end
%         alpha = PStr.PCon('alpha');
%         dispersion = PStr.PCon('dispersion');
        
        [val, idL] = min(abs(EVL - X4));
        [val, idR] = min(abs(X4 - EVR));
        
%         if idL > 1 && EVL < X4(idL)
%             idL = idL - 1;
%         end
        
%         if idR < length(X4) && EVR > X4(idR)
%             idR = idR + 1;
%         end
%         theta_M = atan(imag(B4(idL))/real(B4(idL)));
%         theta_P = atan(imag(B4(idR))/real(B4(idR)));
        Bh = B4(idL:idR);
        Nh = 1 + idR - idL;
        Xh = X4(idL:idR);
        Lh = [EVL, EVR];
        

        
        if 0
            % Find N = 2^n s.t. N > length(Bh)
            PowerNd = floor(log2(length(Bh)) + 1);
            Nd = 2^PowerNd+2;   % so that Nd - 1 = 2^n + 1
            PowBoost = 4;
            N4d = 2^(PowerNd+PowBoost)+2;   % so that Nd - 1 = 2^n + 1

            % Create uniform, symmetric grid with Nd points
            Lh_symmetric = [-1, 1]*(  Lh(2) - Lh(1)  )/2;
            Xsym = linspace(  Lh_symmetric(1), Lh_symmetric(2), Nh  );  % 'original' symmetric grid
            Xuniform = linspace(  Lh_symmetric(1), Lh_symmetric(2), Nd  );

            b2nsymm = interp1(Xsym,Bh,Xuniform); % User's Selection interpolated to fit symmetric grid with (2^n + 2) data points

            X4uniform = linspace(  Lh_symmetric(1), Lh_symmetric(2), N4d  );

            b2T4nsymm = interp1(Xsym,Bh,X4uniform,'spline'); % User's Selection interpolated to fit symmetric grid with (2^n + 2) data points

            Bh = b2T4nsymm;
            Xuniform = X4uniform;
            Nd = N4d;
        end


        

        BC = [Bh(1), Bh(end)];
        bomin = min(abs([Bh(1), Bh(end)]));
%         LRWin = [0.01, 0.99*bomin]; %future development: allow graphical display of region w/o evaluating entirely
%         LREval = [0.01, 0.99*bomin]; %range of real values to study
%% Evaluate Iz in selected region of Profile 
['From: ',thisFile,' Select Section of Profile to Evaluate Iz  ']
        if boolIz
%             IzZoom1 = sum(imag(B4(1:1070)))*dx4; % idx = 1070 was found in one trial to be the right edge of the IS
                                                % NOTE: for non-coplanar
                                                % profiles, bz != 0 to the
                                                % right of the IS, so Iz is
                                                % NOT well defined ...
                                                % without a fixed upper
                                                % limit
            TrialData = split(TrialID,'_');
            TrialTxt = sprintf('TrialID %s_{%s}',TrialData{1},TrialData{2});
            idIzL = idL1;idIzR = idR1;
            xIdxMid = floor( ( idIzL + idIzR )/2);
            xIdxMid0 = floor( ( idL1 + idR1 )/2);
            xMid0 = X4(xIdxMid0);
                      
            HBflag = 1;
%             idL1 = 1;idR1 = length(B4);
            while HBflag == 1      
                
                HBfield = figure('Color', [1, 1, 1],'units','normalized','outerposition',[(NMonitors -1 +.1) 0.1 0.8 0.8]);
                cla();
                xIdxOffsetTxt = floor(0.05*(idR1 - idL1));

                if evaluateKBWShockBoundaries
                    %Following assumes Shock is in LHS of profile
                    [idRelMax, idxShockL, idxShockR, XShockL, XShockR, Izcalc] = IS_ShockBoundaries(B4(1:floor(length(B4)/2)), X4(1:floor(length(B4)/2)), b_L, BZasymptoteL, BZasymptoteR, Ntest, tol);

                    IS_Shock_Profile_Label(HBfield, zmin, zmax, xIdxOffsetTxt, idxShockL, idxShockR, XShockL, XShockR, X4, B4)

                end

                Figure_Label(HBfield, zmax, xIdxOffsetTxt, idL1, idR1, X4, B4)

                plot(X4(idL1:idR1),real(B4(idL1:idR1)),'Color',[0 , 1, 0]);
                hold('on');
                plot(X4(idL1:idR1),imag(B4(idL1:idR1)),'Color',[1 , 0, 0]);
                hold('on');
                plot(X4(idL1:idR1),abs(B4(idL1:idR1)),'Color',[0 , 0, 1]);
                plot([X4(idIzL) X4(idIzL)],[-1.5 1.5],'Color',[1 0 0]);
                hold('on');
                plot([X4(idIzR) X4(idIzR)],[-1.5 1.5],'Color',[0 0 1]);
                hold('on');
                

                xL1 = X4(idL1); xR1 = X4(idR1);
                X0 = xL1 + 0.1*(xR1 - xL1)/2; X1 = xL1 + 0.4*(xR1 - xL1)/2;
                xIdxMid = floor( ( idIzL + idIzR )/2);
                xIzMid = X4(xIdxMid);
                Y1 = 0.94*zmax; Y2 = 0.92*zmax; Y3 = 0.85*zmax; Y4 = 0.75*zmax;Y5 = -0.8*zmax; %Y5 = 0.5*abs(B4(xIdxMid));
                IzSelectedRegion = sum(imag(B4(idIzL:idIzR)))*(-dx4); % to be consistent with KBW convention

                text(X0,Y1,sprintf('Profile = %s',Profile));
                text(X1,Y2,sprintf('Trial ID = %s',TrialTxt));
                text(X0,Y3,sprintf('time = %4.4f',TimePerLine*(LineNumber - 1)));
                text(X0,Y4,sprintf('line # = %4.0f',LineNumber));
                text(xIzMid,Y5,sprintf('Iz = %4.6f',IzSelectedRegion));
                plot([X4(idIzL) X4(idIzR)],[1.1*Y5, 1.1*Y5],'Color',[0 0 0]); % Horizontal line underneath Shock
                hold('on');
                txalert = text(xMid0,0.9*zmax,sprintf('Left click to select region to evaluate Iz'));
                txalert.FontSize = 14;
                txalert.FontWeight = 'bold';

                axis([X4(idL1) X4(idR1) zmin zmax]);
                drawnow();


                [x,y] = ginput(2);
                x = sort(x);
                xL1 = x(1);%xL2 = x(2);
                xR1 = x(2);%xR2 = x(4);
                figure(HBfield);
                plot([xL1 xL1],[-1.5 1.5],'Color',[1 0 0]);
                hold('on');
                %                 plot([xL2 xL2],[-1.5 1.5],'Color',[1 0 0]);
                %                 hold('on');
                plot([xR1 xR1],[-1.5 1.5],'Color',[0 0 1]);
                hold('on');
                %                 plot([xR2 xR2],[-1.5 1.5],'Color',[0 0 1]);
                %                 hold('on');

                axis([X4(idL1) X4(idR1) zmin zmax]);
                drawnow();

                % Construct a questdlg with three options
                choice = questdlg('Selection correct?', ...
                    'Select Iz Region:', ...
                    'Yes','Try Again','Zoom In','Get Me Outa Here');
                % Handle response
                switch choice
                    case 'Try Again'
                        disp([choice ' Try, try again.'])
                        HBflag = 1;
                    case 'Yes'
                        disp([choice ' Woot!'])
                        tmpL1 = abs(X4 - xL1);
                        tmpR1 = abs(X4 - xR1);

                        [val, idIzL] = min(tmpL1);
                        [val, idIzR] = min(tmpR1);
                        HBflag = 1;
                    case 'Zoom In'
                        disp([choice ' Woot!'])
                        tmpL1 = abs(X4 - xL1);
                        tmpR1 = abs(X4 - xR1);

                        [val, idL1] = min(tmpL1);
                        [val, idR1] = min(tmpR1);

                        EVL = X4(idL1);EVR = X4(idR1);
                        HBflag = 1;
                    case 'Get Me Outa Here'
                        disp('Prepare for Transport.')
                        HBflag = -1;
                end
                close(HBfield);
            end
    

            tmpL1 = abs(X4 - xL1);
            tmpR1 = abs(X4 - xR1);

            [val, idL1] = min(tmpL1);
            [val, idR1] = min(tmpR1);

            EVL = X4(idL1);EVR = X4(idR1);
            
        end
        B0 = PStr.PCon('bo');
        nu = PStr.PCon('nu');
        if sum(strcmp(PStr.PCon.keys,'alpha')) < 1
            alpha = 1;
        else
            alpha = PStr.PCon('alpha');     
        end
        if sum(strcmp(PStr.PCon.keys,'dispersion')) < 1
            dispersion = 1;
        else
            dispersion = PStr.PCon('dispersion');     
        end

        gamma_Landau = 0;
        if sum(strcmp(PStr.PCon.keys,'gamma_Landau')) > 0 
           gamma_Landau = PStr.PCon('gamma_Landau');     
        end
%         alpha = PStr.PCon('alpha');
%         dispersion = PStr.PCon('dispersion');
        
        [val, idL] = min(abs(EVL - X4));
        [val, idR] = min(abs(X4 - EVR));
        
%         if idL > 1 && EVL < X4(idL)
%             idL = idL - 1;
%         end
        
%         if idR < length(X4) && EVR > X4(idR)
%             idR = idR + 1;
%         end
%         theta_M = atan(imag(B4(idL))/real(B4(idL)));
%         theta_P = atan(imag(B4(idR))/real(B4(idR)));
        Bh = B4(idL:idR);
        Nh = 1 + idR - idL;
        Xh = X4(idL:idR);
        Lh = [EVL, EVR];
        

      
        

        BC = [Bh(1), Bh(end)];
        bomin = min(abs([Bh(1), Bh(end)]));
%         LRWin = [0.01, 0.99*bomin]; %future development: allow graphical display of region w/o evaluating entirely
%         LREval = [0.01, 0.99*bomin]; %range of real values to study
%% Create B-field Graph Figure/Jpg
        if boolBgraph == 1 && boolVideo ~= 1 && (~boolDevelopment)
%             BFieldGraphGroot = groot;
%             NMonitors = length(BFieldGraphGroot.MonitorPositions(:,1));
            
            ProfileFigFile = [ProfileFigBase,num2str(LineNumber),'.fig'];%e.g. 'F:\ ...\Trial 11\'
            ProfileJpgFile = [ProfileJpgBase,num2str(LineNumber),'.jpg'];%e.g. 'F:\ ...\Trial 11\'
            
            IzZoom = sum(imag(Bh))*(-dx4); % to be consistent with KBW convention
            TrialData = split(TrialID,'_');
            TrialTxt = sprintf('TrialID %s_{%s}',TrialData{1},TrialData{2});
            XDelta = 0.17*(EVR-EVL);YDelta = 0.04*(zmax - zmin);
            XLeft=EVL + 0.02*(EVR-EVL) ; YTop = 1.05*zmax + 0.04*(zmax - zmin);
            XT0 = XLeft; XT1=XLeft + 1*XDelta; XT2=XLeft + 2*XDelta; YT0 = YTop;YT1 = YTop - 1*YDelta;YT2 = YTop - 2*YDelta;YT3 = YTop - 3*YDelta;YT4 = YTop - 4*YDelta;YT5 = YTop - 5*YDelta;
            
            HBfield = figure('Color', [1, 1, 1],'units','normalized','outerposition',[(NMonitors -1 +.02) 0.1 0.6 0.6]);
            cla();

%             evaluateKBWShockBoundaries = 1;
            if evaluateKBWShockBoundaries
%                 BZ_asymptote = 0;
%                 Ntest = 10;
%                 tol = 0.01;
                [idRelMax, idxShockL, idxShockR, XShockL, XShockR, Izcalc] = IS_ShockBoundaries(Bh, Xh, b_L, BZasymptoteL, BZasymptoteR, Ntest, tol);
                plot([XShockL XShockL],[0.8*zmin, 0.8*zmax],'Color',[0 0 0]);
                hold('on');

                plot([XShockR XShockR],[0.8*zmin, 0.8*zmax],'Color',[0 0 0]);
                hold('on');
                text(XT2,YT0,sprintf('Iz between lines = %4.4f',Izcalc));
            end


            plot(Xh,real(Bh),'Color',[1 0 0]);
            hold('on');

            plot(Xh,imag(Bh),'Color',[0 1 0]);
            hold('on');

            plot(Xh,abs(Bh),'Color',[0 0 1]);
            hold('on');
            
 
    
            text(XT0,YT0,sprintf('Profile = %s',Profile));
            text(XT0,YT1,sprintf('Iz = %4.4f',IzZoom));
%             text(XT0,YT2,sprintf('\\Delta\\theta  = %4.4f rad',theta_kbw));
%             text(XT0,YT3,sprintf('b_{kbw} = %4.4f',bkbw));
            text(XT0,YT4,sprintf('line # = %4.0f',LineNumber));
            text(XT0,YT5,sprintf('time = %4.4f',TimePerLine*(LineNumber - 1)));
            
            
            text(XT1,YT0,sprintf('Trial ID = %s',TrialTxt));
            text(XT1,YT1,sprintf('dispersion = %4.4f',dispersion));
            text(XT1,YT2,sprintf('dissipation = %4.4f',nu));
            text(XT1,YT3,sprintf('alpha = %4.4f',alpha));
            
           


            axis([EVL EVR zmin zmax]);
            set(gca, 'box', 'off');
            drawnow(); 
            
            % GUI: Save Profile as JPG/Fig?
%             ProfileFile = [FigFileBase,num2str(LineNumber),
            choice = questdlg('Save Graph?', ...
            'Save Graph:', ...
            'Save as *.jpeg?','Save as *.jpg and *.fig?','No.','No.');
            % Handle response
            switch choice
            case 'Save as *.jpeg?'
                utilities(JpgAbs, JpgRel); %Creates JPEG folder if needed
                saveas(gcf,ProfileJpgFile);               
                disp([choice 'Saved as *.jpeg.'])
            case 'Save as *.jpg and *.fig?'
                utilities(JpgAbs, JpgRel); %Creates JPEG folder if needed
                utilities(FigAbs, FigRel); %Creates JPEG folder if needed
                saveas(gcf,ProfileJpgFile);
                saveas(gcf,ProfileFigFile);
                disp([choice 'Saved as *.jpg and *.fig'])
            case 'No.'
                disp('Moving on.')
            end
        end
%% Testing detection of constant phase: Downstream terminus of Intermediate Shock        
        
    if 0
            
            ProfileFigFile = [ProfileFigBase,num2str(LineNumber),'.fig'];%e.g. 'F:\ ...\Trial 11\'
            ProfileJpgFile = [ProfileJpgBase,num2str(LineNumber),'.jpg'];%e.g. 'F:\ ...\Trial 11\'
            
            IzZoom = sum(imag(Bh))*(-dx4); % to be consistent with KBW paper
            TrialData = split(TrialID,'_');
            TrialTxt = sprintf('TrialID %s_{%s}',TrialData{1},TrialData{2});
            XDelta = 0.17*(EVR-EVL);YDelta = 0.04*(zmax - zmin);
            XLeft=EVL + 0.02*(EVR-EVL) ; YTop = 1.05*zmax + 0.04*(zmax - zmin);
            XT0 = XLeft; XT1=XLeft + 1*XDelta;YT0 = YTop;YT1 = YTop - 1*YDelta;YT2 = YTop - 2*YDelta;YT3 = YTop - 3*YDelta;YT4 = YTop - 4*YDelta;YT5 = YTop - 5*YDelta;
            
         HBfield = figure('Color', [1, 1, 1],'units','normalized','outerposition',[(NMonitors -1 +.02) 0.1 0.6 0.6]);
            cla();
            plot(Xh,real(Bh),'Color',[1 0 0]);
            hold('on');
            plot(Xh,imag(Bh),'Color',[0 1 0]);
            hold('on');
            plot(Xh,abs(Bh),'Color',[0 0 1]);
            hold('on');
            plot(Xh,angle(Bh),'Color',[0 0 1]);
            hold('on');
    
            text(XT0,YT0,sprintf('Profile = %s',Profile));
            text(XT0,YT1,sprintf('Iz = %4.4f',IzZoom));
            text(XT0,YT2,sprintf('\\Delta\\theta  = %4.4f rad',theta_kbw));
            text(XT0,YT3,sprintf('b_{kbw} = %4.4f',bkbw));
            text(XT0,YT4,sprintf('line # = %4.0f',LineNumber));
            text(XT0,YT5,sprintf('time = %4.4f',TimePerLine*(LineNumber - 1)));
            
            
            text(XT1,YT0,sprintf('Trial ID = %s',TrialTxt));
            text(XT1,YT1,sprintf('dispersion = %4.4f',dispersion));
            text(XT1,YT2,sprintf('dissipation = %4.4f',nu));
            text(XT1,YT3,sprintf('alpha = %4.4f',alpha));


            axis([EVL EVR -4 4]);
            set(gca, 'box', 'off');
            drawnow(); 
            
            % GUI: Save Profile as JPG/Fig?
%             ProfileFile = [FigFileBase,num2str(LineNumber),
            choice = questdlg('Save Graph?', ...
            'Save Graph:', ...
            'Save as *.jpeg?','Save as *.jpg and *.fig?','No.','No.');
            % Handle response
            switch choice
            case 'Save as *.jpeg?'
                utilities(JpgAbs, JpgRel); %Creates JPEG folder if needed
                saveas(gcf,ProfileJpgFile);               
                disp([choice 'Saved as *.jpeg.'])
            case 'Save as *.jpg and *.fig?'
                utilities(JpgAbs, JpgRel); %Creates JPEG folder if needed
                utilities(FigAbs, FigRel); %Creates JPEG folder if needed
                saveas(gcf,ProfileJpgFile);
                saveas(gcf,ProfileFigFile);
                disp([choice 'Saved as *.jpg and *.fig'])
            case 'No.'
                disp('Moving on.')
            end
        end
%% Create Hodogram Figure/Jpg        
        if boolHodogram == 1 && boolVideo ~= 1 && (~boolDevelopment)
            
            HodoFigFile = [HodoFigBase,num2str(LineNumber),'.fig'];%e.g. 'F:\ ...\Trial 11\'
            HodoJpgFile = [HodoJpgBase,num2str(LineNumber),'.jpg'];%e.g. 'F:\ ...\Trial 11\'
            
            Hodofield = figure('Color', [1, 1, 1],'units','normalized','outerposition',[(NMonitors -1 +.65) 0.1 0.3 0.6]);
            
%             evaluateKBWShockBoundaries = 1;
            if evaluateKBWShockBoundaries
%                 BZ_asymptote = 0;
%                 Ntest = 10;
%                 tol = 0.01;
                [idRelMax, idxShockL, idxShockR, XShockL, XShockR, Izcalc] = IS_ShockBoundaries(Bh, Xh, b_L, BZasymptoteL, BZasymptoteR, Ntest, tol);
                
                XCL = real(Bh(idxShockL)); YCL = imag(Bh(idxShockL));
                plot(XCL,YCL,'o','MarkerFaceColor','green','MarkerSize',5);
                hold 'on';
                
                XCR = real(Bh(idxShockR)); YCR = imag(Bh(idxShockR));
                plot(XCR,YCR,'o','MarkerFaceColor','green','MarkerSize',5);
                hold 'on';
                           
            end
            
            XDelta = 0.17*(EVR-EVL);YDelta = 0.04*(zmax - zmin);
            XLeft=EVL + 0.02*(EVR-EVL) ; YTop = 1.05*zmax + 0.04*(zmax - zmin);
            XT0 = XLeft; XT1=XLeft + 1*XDelta; XT2=XLeft + 2*XDelta; YT0 = YTop;YT1 = YTop - 1*YDelta;YT2 = YTop - 2*YDelta;YT3 = YTop - 3*YDelta;YT4 = YTop - 4*YDelta;YT5 = YTop - 5*YDelta;
            
            
            IzZoom = sum(imag(Bh))*(-dx4); % to be consistent with KBW paper
            Bmax = 1.5*max(abs([real(Bh), imag(Bh)]));
            BH_L = Bh(1); BH_R = Bh(end);
%             cla();
            TrialData = split(TrialID,'_');
            TrialTxt = sprintf('TrialID %s_{%s}',TrialData{1},TrialData{2});
            txtL = ['bL = ',num2str(BH_L)];txtR = ['bR = ',num2str(BH_R)];
            txtCR = ['CR = ',num2str(abs(Prfl('CR13')))];
            CR = Prfl('CR13'); 
            XCR = real(BH_L)*CR; YCR = imag(BH_L);
            XCRtxt = XCR + Bmax*0.1; YCRtxt = YCR -Bmax*0.2;
            XL = real(BH_L)-Bmax*0.2;YL = imag(BH_L)-Bmax*0.05;
            XR = real(BH_R)-Bmax*0.1;YR = imag(BH_R)-Bmax*0.05;
            plot(real(Bh),imag(Bh));
            hold on;
            
            plot(XCR,YCR,'r+');
            p1 = [XCRtxt YCRtxt];                         % First Point
            p2 = [XCR YCR];                         % Second Point
            dp = p2-p1;                         % Difference
            quiver(p1(1),p1(2),dp(1),dp(2),0);
            hold on;
            YDelta = 0.05;
          
            
                % There will likely be errors thrown for profiles other than 'KBW'
                % Text display will need to be formatted and selected based
                % on profile type!!!!
                text(Bmax*(0.0), Bmax*(1.0 + YDelta),TrialTxt);
                text(Bmax*(0.0), Bmax*(0.9 + YDelta),sprintf('dispersion = %4.4f',dispersion));
                text(Bmax*(0.0), Bmax*(0.8 + YDelta),sprintf('dissipation = %4.4f',nu));
                text(Bmax*(0.0), Bmax*(0.7 + YDelta),sprintf('alpha = %4.4f',alpha));
                text(Bmax*(-0.9), Bmax*(1.0 + YDelta),sprintf('Profile:%s',Profile));
                text(Bmax*(-0.9), Bmax*(0.9 + YDelta),sprintf('Iz = %4.4f',IzZoom));
                text(Bmax*(-0.9), Bmax*(0.8 + YDelta),sprintf('\\Delta\\theta  = %4.4f rad',theta_kbw));
                text(Bmax*(-0.9), Bmax*(0.7 + YDelta),sprintf('b_{kbw} = %4.4f',bkbw));
                text(Bmax*(-0.9), Bmax*(0.6 + YDelta),sprintf('line # = %4.0f',LineNumber));
                text(Bmax*(-0.9), Bmax*(0.5 + YDelta),sprintf('time = %4.4f',TimePerLine*(LineNumber - 1)));
                text(XL, YL,txtL);
                text(XR, YR,txtR);
                text(XCRtxt, YCRtxt,txtCR);
            
%             text(Bmax*(-0.9), Bmax*(0.9),sprintf('Iz = %4.4f',IzZoom));
%             text(XL, YL,txtL);
%             text(XR, YR,txtR);
            axis([-1 +1 -1 +1]*Bmax);
            drawnow();
            
            % GUI: Save Hodo as JPG/Fig?
%             ProfileFile = [FigFileBase,num2str(LineNumber),
            choice = questdlg('Save Hodogram?', ...
            'Save Hodogram:', ...
            'Save as *.jpeg?','Save as *.jpg and *.fig?','No.','No.');
            % Handle response
            switch choice
            case 'Save as *.jpeg?'
                utilities(JpgAbs, JpgRel); %Creates JPEG folder if needed
                saveas(gcf,HodoJpgFile);               
                disp([choice 'Saved as *.jpeg.'])
            case 'Save as *.jpg and *.fig?'
                utilities(JpgAbs, JpgRel); %Creates JPEG folder if needed
                utilities(FigAbs, FigRel); %Creates JPEG folder if needed
                saveas(gcf,HodoJpgFile);
                saveas(gcf,HodoFigFile);
                disp([choice 'Saved as *.jpg and *.fig'])
            case 'No.'
                disp('Moving on.')
            end
        end  
%% Animations
    %% B-Field Graph Video
        if (boolVideo || boolPeakOscAnalysis)  && (~boolDevelopment)
            

            if EVZoom
                X3min = EVL;X3max = EVR;% <-- use if passed through EVZoom
            else
                X3min = -L;X3max = +L;
            end
            
            tmp3L = abs(X4 - (X3min));
            tmp3R = abs(X4 - (X3max));
            [val, idx3L] = min(tmp3L);
            [val, idx3R] = min(tmp3R);
            
            idx3L = idx3L + mod(idx3L,2); % there's probably an easier way to do this, but make idx3L an even number
            idx3R = idx3R + (1 - mod(idx3R,2) ); % results in idx3R being odd
            N3 = idx3R - idx3L; %with idx3L being even and idx3R being odd, this results in N3 being odd
            N3 = N3 + mod(N3,2); % This makes N3 even AND makes it match the # of data points idx3L:idx3R inclusive
            X3 = linspace(X3min,X3max,N3);
            
            TrialData = split(TrialID,'_');
            TrialTxt = sprintf('TrialID %s_{%s}',TrialData{1},TrialData{2});
            XDelta = 0.17*(EVR-EVL);YDelta = 0.04*(zmax - zmin);
            XLeft=EVL + 0.02*(EVR-EVL) ; YTop = 1.05*zmax + 0.04*(zmax - zmin);
            XT0 = XLeft; XT1=XLeft + 1*XDelta; XT2=XLeft + 2*XDelta;
            YT0 = YTop;YT1 = YTop - 1*YDelta;YT2 = YTop - 2*YDelta;YT3 = YTop - 3*YDelta;YT4 = YTop - 4*YDelta;YT5 = YTop - 5*YDelta;YT6 = YTop - 6*YDelta;

%             Vold = Prfl('speed');
            if boolSpeedCheck == 1
                xFSmiddle = zeros(length(BFlines));
                tFSmiddle = zeros(length(BFlines));
                vFSmiddle = zeros(length(BFlines));
            end
            
            if (boolVideo == 1)
                utilities(AviFolderAbs, AviFolderRel); %Creates both Trial & avi folders if needed
                FrameRate = 10;
                
                VidFileName = ['\BFvideo_',TrialID,'.avi'];
                BFVideoFile = [AviFolderAbs,VidFileName];
                writerObjBF = VideoWriter(BFVideoFile);
                writerObjBF.FrameRate = FrameRate;
                open(writerObjBF);
               
                if boolmp4 == 1
                    VidFileNamemp4 = ['\BFvideo_',TrialID,'.mp4'];
                    BFVideoFilemp4 = [AviFolderAbs,VidFileNamemp4];
                    writerObjBFmp4 = VideoWriter(BFVideoFilemp4,'MPEG-4');
                    writerObjBFmp4.FrameRate = FrameRate;
                    open(writerObjBFmp4);
                end
           
                BFAnim =figure('Color', [1, 1, 1],'units','normalized','outerposition',[(NMonitors -1 +.1) 0.1 0.8 0.8]);
                set(BFAnim,'doublebuffer','on');

                for kk = 1:length(BFlines)
                    temp = BFlines(kk*1);
                    C = textscan(temp{1},'%f');
                    Bk = C{1}.'; %B is an N/4 complex valued, double precision array
                    
                    Bk3 = Bk(idx3L:(idx3R));
                    IzZoom = sum(imag(Bk3))*(X3(1) - X3(2)); % to be consistent with KBW paper
                                        
                    figure(BFAnim);
                    clf();
                    set(gca,'nextplot','replacechildren');
                    set(gcf,'Renderer','zbuffer');
                    
                    
                    
                    %                 [val, idL] = min(abs(xkbw-0.1*L - X3));
                    %                 [val, idR] = min(abs(xkbw+0.1*L - X3));
                    %                 [val, idRelMax] = max(abs(imag(Bk3(idL:idR))));
                    %                 idMax = idRelMax + idL;
                    %                 NinterpolateSteps = 40;
                    %                 NXsteps = 20;
                    %                 Xint = linspace(X3(idMax - NXsteps),X3(idMax + NXsteps),2*NinterpolateSteps);
                    %                 Binterp = interp1(X3(idMax-NXsteps:idMax+NXsteps),Bk3(idMax-NXsteps:idMax+NXsteps),Xint);
                    %                 [val,idIntMax] = max(abs(imag(Binterp)));
                    %                 Xnew = Xint(idIntMax);
                    %                 vtmp = Vold - (Xold - Xnew)/(timestepsPerFrame*dt);
                    %                 v = Vold*.75 + 0.25*vtmp;
                    %                 Xold = Xnew;
                    %                 Vold = v;
                    %                 q = -1/2*(1+sqrt(1 + 4*(v-Prfl('b_L'))));
                    %                 evaluateKBWShockBoundaries = 1;
                    if evaluateKBWShockBoundaries
                        BZ_asymptote = 0;
                        Ntest = 10;
                        tol = 0.01;
                        [idRelMax, idxShockL, idxShockR, XShockL, XShockR, Izcalc] = IS_ShockBoundaries(Bk3, X3, b_L, BZasymptoteL, BZasymptoteR, Ntest, tol);
                        plot([XShockL XShockL],[0.8*zmin, 0.8*zmax],'Color',[0 0 0]);
                        
                        hold('on');
                        
                        plot([XShockR XShockR],[0.8*zmin, 0.8*zmax],'Color',[0 0 0]);
                        hold('on');
                        text(XT2,YT0,sprintf('Iz between lines = %4.4f',Izcalc));
                    end
                    
                    if boolSpeedCheck == 1
                        BmidShock = 1/2*(real(Bk3(1)) + real(Bk3(end)));
                        [~,idxFSmiddle]=min(abs(real(Bk3)-BmidShock));
                        xFSmiddle(kk) = X3(idxFSmiddle);
                        tFSmiddle(kk) = TimePerLine*(kk - 1);
                        
                        plot([xFSmiddle(kk) xFSmiddle(kk)],[0.8*zmin, 0.8*zmax],'Color',[0 0 1]);
                        hold('on');
                        xShockLabel = xFSmiddle(kk) + 0.02*(EVR-EVL);
                        text(xShockLabel  - 0.16*(EVR-EVL), 0.8*zmin + 3*YDelta, sprintf('Mid point of Fast Shock -->'));
                        text(xShockLabel  + 0.05*(EVR-EVL), 0.25*zmax + 1*YDelta, sprintf('Bleft  = %4.4f', abs( Bk3(1)   )  ));
                        text(xShockLabel  + 0.05*(EVR-EVL), 0.25*zmax + 0*YDelta, sprintf('Bright = %4.4f', abs( Bk3(end) )  ));
                        hold('on');
                        vFSmiddle(kk) = speed;
                        
                        if kk > 1
                            deltaSpeed =( xFSmiddle(kk) - xFSmiddle(kk-1) )/TimePerLine;
                            vFSmiddle(kk) = vFSmiddle(kk) + deltaSpeed;
                        end
                        
                        Kstart = 18;
                        if kk > Kstart
                            vFSmax    =    max(vFSmiddle(Kstart:kk));
                            vFSmean   =   mean(vFSmiddle(Kstart:kk));
                            vFSmedian = median(vFSmiddle(Kstart:kk));
                            vFSmin    =    min(vFSmiddle(Kstart:kk));
                            vFSerror  =   std( vFSmiddle(Kstart:kk) ) / sqrt( length( vFSmiddle(Kstart:kk) ));
                            
                            text(xShockLabel, 0.8*zmin + 6*YDelta, sprintf('Vcalc           = %4.4f', speed));
                            text(xShockLabel, 0.8*zmin + 5*YDelta, sprintf('Vactual - Vcalc = %4.4f', deltaSpeed));
                            text(xShockLabel, 0.8*zmin + 4*YDelta, sprintf('Vactual max     = %4.4f', vFSmax ));
                            text(xShockLabel, 0.8*zmin + 3*YDelta, sprintf('Vactual mean    = %4.4f', vFSmean ));
                            text(xShockLabel, 0.8*zmin + 2*YDelta, sprintf('Vactual median  = %4.4f', vFSmedian ));
                            text(xShockLabel, 0.8*zmin + 1*YDelta, sprintf('Vactual min     = %4.4f', vFSmin ));
                            text(xShockLabel, 0.8*zmin + 0*YDelta, sprintf('Vactual StdErr) = %4.4f', vFSerror ));
                            
                            hold('on');
                        end
                        
                    end
                    
                    plot(X3,real(Bk3),'Color',[0 , 1, 0]);
                    hold('on');
                    plot(X3,imag(Bk3),'Color',[1 , 0, 0]);
                    hold('on');
                    plot(X3,abs(Bk3),'Color',[0 , 0, 1]);
                    hold('on');
                    
                    text(XT0,YT1,sprintf('Trial ID = %s',TrialTxt));
                    if strcmp(Profile, 'KBW')
                        text(XT0,YT2,sprintf('Iz = %4.4f',IzZoom));
                        text(XT0,YT3,sprintf('\\Delta\\theta  = %4.4f rad',theta_kbw));
                        text(XT0,YT4,sprintf('b_{kbw} = %4.4f',bkbw));
                        %                     text(XT0,YT5,sprintf('line # = %4.0f',kk));
                    elseif strcmp(Profile, '1psoliton')
                        text(XT0,YT2,sprintf('b_0 = %4.4f',Prfl('b1p')));
                        text(XT0,YT3,sprintf('\\lambda = (%4.4f)%4.4f',Prfl('lambda'),Prfl('b1p')));
                        text(XT0,YT4,sprintf('sign = %4.0f',Prfl('sign')));
                        %                     text(XT0,YT5,sprintf('line # = %4.0f',kk));
                    elseif strcmp(Profile, '2psoliton')
                        text(XT0,YT2,sprintf('b_0 = %4.4f',Prfl('b2p')));
                        text(XT0,YT3,sprintf('\\lambdar = (%4.4f)%4.4f',Prfl('lr'),Prfl('b2p')));
                        text(XT0,YT4,sprintf('\\lambdai = (%4.4f)%4.4f',Prfl('li'),Prfl('b2p')));
                        %                     text(XT0,YT5,sprintf('line # = %4.0f',kk));
                    end
                    
                    
                    text(XT1,YT1,sprintf('Profile = %s',Profile));
                    text(XT1,YT2,sprintf('dispersion = %4.4f',dispersion));
                    text(XT1,YT3,sprintf('dissipation = %4.4f',nu));
                    text(XT1,YT4,sprintf('alpha = %4.4f',alpha));
                    if gamma_Landau ~= 0
                        text(XT1,YT5,sprintf('gamma_Landau = %4.4f',gamma_Landau));
                    end
                    text(XT2,YT2,sprintf('time = %4.4f/%4.1f',TimePerLine*(kk - 1),TimePerLine*length(BFlines)));
                    text(XT2,YT3,sprintf('line # = %4.0f/%4.0f',kk,length(BFlines)));
                    
                    
                    if strcmp(Profile, 'KBW')
                        text(XT2,YT1,sprintf('Compression Ratio = %4.4f',CompressionRatio));
                    end
                    
                    
                    axis([X3min X3max zmin zmax]);
                    drawnow();
                    
                    frame = getframe(gcf);
                    writeVideo(writerObjBF,frame);
                    if boolmp4 == 1
                        writeVideo(writerObjBFmp4,frame);
                    end
                    
                    
                    
                end
                close(writerObjBF);
                if boolmp4 == 1
                    close(writerObjBFmp4);
                end
            end
            
            if (boolPeakOscAnalysis == 1)                
                kkStart = LineNumber;
                POALines = length(BFlines) - kkStart;
                PLb = zeros(1, POALines);
                PCb = zeros(1, POALines);
                PRb = zeros(1, POALines);
                
                PLx = zeros(1, POALines);
                PCx = zeros(1, POALines);
                PRx = zeros(1, POALines);
                
                ISPeakOscillationFigName = 'Oscillations of central peak magnitudes';
                ISPeakOscillationFig     = figure('Name',ISPeakOscillationFigName,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',[(NMonitors -1 +.1) 0.55 0.8 0.40]);

                ISFFTOscillationFigName = 'FFT of |b|(t) of central peaks';
                ISFFTFig                 = figure('Name',ISFFTOscillationFigName ,'NumberTitle','off','Color', [1, 1, 1],'units','normalized','outerposition',[(NMonitors -1 +.1) 0.05 0.8 0.40]);

                for kk = 1:length(BFlines)
                    temp = BFlines(kk*1);
                    C = textscan(temp{1},'%f');
                    Bk = C{1}.'; %B is an N/4 complex valued, double precision array
                    
                    Bk3 = Bk(idx3L:(idx3R));
                    IzZoom = sum(imag(Bk3))*(X3(1) - X3(2)); % to be consistent with KBW paper
                    
                    if kk > kkStart
                        % Find location of:
                        % PC = central peak of |b| of oscillatory shock
                        % PL = central peak of |b| of oscillatory shock
                        % PC = central peak of |b| of oscillatory shock
                        
                        POAidx = kk - kkStart;
                        
                        K3x = [0:N3/2-1,-N3/2:-1] * pi/L; %N3 defined above to be even
                        Bk3abs = abs(Bk3);
%                         Bk3AbsPrimed = real( ifft( 1i*K3x.*fft( Bk3abs ) ) );
%                         Bk3AbsDoublePrimed = real( ifft( -K3x.*K3x.*fft( Bk3abs ) ) );
%                         Bk3AbsPrimedShift = circshift(Bk3AbsPrimed,1);
%                         Bk3PrimeSignChanges = Bk3AbsPrimed.*Bk3AbsPrimedShift;
%                         Bk3MaxMinIdxs = find(Bk3PrimeSignChanges < 0);
                        
                        [Bk3MaxMinIdxs,Bk3AbsPrimed] = FindMaxMinIndicesOriginal(Bk3);
                       
                        
                        [val, PCidx] = max( abs( Bk3 ) );
                        
                        
                        [val, MaxIdx] = min(abs(Bk3MaxMinIdxs - PCidx));
                        PLidx = Bk3MaxMinIdxs(MaxIdx - 2);
                        PRidx = Bk3MaxMinIdxs(MaxIdx + 2);
                        
                        PLb(POAidx) = Bk3abs(PLidx);
                        PCb(POAidx) = Bk3abs(PCidx);
                        PRb(POAidx) = Bk3abs(PRidx);
                        
                        PLx(POAidx) = X3(PLidx);
                        PCx(POAidx) = X3(PCidx);
                        PRx(POAidx) = X3(PRidx);
                        

                        
                    end
                    
                    
                    
                end
                
                % cla;plot(1.27);hold on; plot(1.082);hold on;plot(0.5771);hold on;plot(PLb);hold on;plot(PCb);hold on; plot(PRb);hold off;
                
                
                PLbs = smoothdata(PLb);
                PCbs = smoothdata(PCb);
                PRbs = smoothdata(PRb);
                T3 = (kkStart:1:length(BFlines)-1 )*TimePerLine;
                
                Ymax = 1.1*max(PCb);
                Ymin = 0.7*min(PRb);

                figure(ISPeakOscillationFig);
                PLbAvg = mean(PLb);
                PCbAvg = mean(PCb);
                PRbAvg = mean(PRb);
                PLinit =PLbAvg*ones(1,POALines);
                PCinit = PCbAvg*ones(1,POALines);
                PLinit = PLbAvg*ones(1,POALines);
                PRinit = PRbAvg*ones(1,POALines);
                
                [PLbMaxMinIdxs,PLbPP] = FindMaxMinIndicesOriginal(PLbs);
                [PCbMaxMinIdxs,PCbPP] = FindMaxMinIndicesOriginal(PCbs);
                [PRbMaxMinIdxs,PRbPP] = FindMaxMinIndicesOriginal(PRbs);
                
                lenL = length(PLbMaxMinIdxs);
                lenC = length(PCbMaxMinIdxs);
                lenR = length(PRbMaxMinIdxs);
                NPeaksMax = min([lenL, lenC, lenR]) - 1;
                NPeaksMax = (NPeaksMax + mod(NPeaksMax,2) )/2;
                
                [PLbMinIdxs,PLbMaxIdxs] = FindMaxMinIndices(PLb, NPeaksMax);
                [PCbMinIdxs,PCbMaxIdxs] = FindMaxMinIndices(PCb, NPeaksMax);
                [PRbMinIdxs,PRbMaxIdxs] = FindMaxMinIndices(PRb, NPeaksMax);
                
               % Linear Regression Fit to Y = Mx + B
                [XLMaxtmp, bLMaxtmp, PLbMaxCalc] = MaxMinLinearRegression(T3, PLbMaxIdxs, PLb);
                [XCMaxtmp, bCMaxtmp, PCbMaxCalc] = MaxMinLinearRegression(T3, PCbMaxIdxs, PCb);
                [XRMaxtmp, bRMaxtmp, PRbMaxCalc] = MaxMinLinearRegression(T3, PRbMaxIdxs, PRb);
                
                [XLMintmp, bLMintmp, PLbMinCalc] = MaxMinLinearRegression(T3, PLbMinIdxs, PLb);
                [XCMintmp, bCMintmp, PCbMinCalc] = MaxMinLinearRegression(T3, PCbMinIdxs, PCb);
                [XRMintmp, bRMintmp, PRbMinCalc] = MaxMinLinearRegression(T3, PRbMinIdxs, PRb);
                
                % Linear Regression Fit to B = Bo + A*exp(kt) <--> log(B - Bavg) = log(A) + kx
                YLexpMax = log(PLb - PLbAvg);
                [XLMaxtmp, bLMaxtmp, PLbMaxExpCalc] = MaxMinLinearRegression(T3, PLbMaxIdxs, YLexpMax);
                ALMax = exp(bLMaxtmp(1));
                BLMaxExpCalc = PLbAvg + ALMax*exp(bLMaxtmp(2)*XLMaxtmp);
                DY = (Ymax - Ymin)*0.02;
                idxYmid = round(NPeaksMax/2);
                Xtxt = XLMaxtmp(idxYmid);
                
                YLtoptxt  = PLbMaxCalc( idxYmid ) + DY;
                YLdowntxt = PLbMinCalc( idxYmid ) - DY;               
                Ltxtup   = [num2str(bLMaxtmp(1)),' + ',num2str(bLMaxtmp(2)),'x'];
                Ltxtdown = [num2str(bLMintmp(1)),' + ',num2str(bLMintmp(2)),'x'];
                
                YCtoptxt  = PCbMaxCalc( idxYmid ) + DY;
                YCdowntxt = PCbMinCalc( idxYmid ) - DY;               
                Ctxtup   = [num2str(bCMaxtmp(1)),' + ',num2str(bCMaxtmp(2)),'x'];
                Ctxtdown = [num2str(bCMintmp(1)),' + ',num2str(bCMintmp(2)),'x'];
                
                YRtoptxt  = PRbMaxCalc( idxYmid ) + DY;
                YRdowntxt = PRbMinCalc( idxYmid ) - DY;               
                Rtxtup   = [num2str(bRMaxtmp(1)),' + ',num2str(bRMaxtmp(2)),'x'];
                Rtxtdown = [num2str(bRMintmp(1)),' + ',num2str(bRMintmp(2)),'x'];
                              
                cla;
                text(Xtxt, YLtoptxt , Ltxtup); hold on;
                text(Xtxt, YLdowntxt, Ltxtdown); hold on;
                
                text(Xtxt, YCtoptxt , Ctxtup); hold on;
                text(Xtxt, YCdowntxt, Ctxtdown); hold on;
                
                text(Xtxt, YRtoptxt , Rtxtup); hold on;
                text(Xtxt, YRdowntxt, Rtxtdown); hold on;
                
                plot( T3,PLinit,'r--');hold on; 
                plot( T3,PCinit,'g--');hold on;
                plot( T3,PRinit,'b--');hold on;
                legend('Initial |b| Left','Initial |b| Center','Initial |b| Right');
                
                plot( T3,PLb,'r-','DisplayName','Left of center');hold on;
%                 plot( T3,PLbs,'r:','DisplayName','Left of center');hold on;
                plot( XLMaxtmp,PLb(PLbMaxIdxs),'r+','DisplayName','Left of center');hold on;
                plot( XLMintmp,PLb(PLbMinIdxs),'ro','DisplayName','Left of center');hold on;
                plot( XLMaxtmp,PLbMaxCalc,'r:','DisplayName','Left of center EXP Fit');hold on;
                plot( XLMaxtmp,BLMaxExpCalc,'rx','DisplayName','Left of center LR');hold on;
                plot( XLMintmp,PLbMinCalc,'r:','DisplayName','Left of center LR');hold on;
                
                plot( T3,PCb,'g-','DisplayName','Central Peak');hold on; 
%                 plot( T3,PCbs,'g:','DisplayName','Central Peak');hold on; 
                plot( XCMaxtmp,PCb(PCbMaxIdxs),'g+','DisplayName','Center');hold on;
                plot( XCMintmp,PCb(PCbMinIdxs),'go','DisplayName','Center');hold on;
                plot( XCMaxtmp,PCbMaxCalc,'g:','DisplayName','Center LR','LineWidth',1.5);hold on;
                plot( XCMintmp,PCbMinCalc,'g:','DisplayName','Center LR','LineWidth',1.5);hold on;
                
                plot( T3,PRb,'b-','DisplayName','Right of center');hold on;
%                 plot( T3,PRbs,'b:','DisplayName','Right of center');hold on;
                plot( XRMaxtmp,PRb(PRbMaxIdxs),'b+','DisplayName','Right of center');hold on;
                plot( XRMintmp,PRb(PRbMinIdxs),'bo','DisplayName','Right of center');hold on;
                plot( XRMaxtmp,PRbMaxCalc,'b:','DisplayName','Right of center LR');hold on;
                plot( XRMintmp,PRbMinCalc,'b:','DisplayName','Right of center LR');hold off;
                
                xlim([T3(1) T3(end)]);
                ylim([Ymin Ymax]);
              
%                 PLbMaxValues = -1*ones(1,NPeaksMax);


                %Prep FFTs for graphing
                Npoa = POALines;
                Npoa = Npoa - mod(Npoa,2); % force Npoa to be even and less than # of data points
                if Npoa < POALines
                    PLb=PLb(1:end-1); % adjust # data points to match Npoa by excluding the last point
                    PCb=PCb(1:end-1);
                    PRb=PRb(1:end-1);
                end
                FLb = fft(PLb);
                FCb = fft(PCb);
                FRb = fft(PRb);
                fs = 1/TimePerLine;
                f = (0:length(FLb)-1)*fs/length(FLb);
%                 fshift = (-Npoa/2:Npoa/2-1)*(fs/Npoa);
                fshift = [0:Npoa/2-1,-Npoa/2:-1] *fs/Npoa;
                
                [val, FCbPeakidx] = max( abs( FCb(2:end) ) );
                FCbPeakidx = FCbPeakidx + 1;
                FFTMax = max(abs([FLb(2:end),FCb(2:end),FRb(2:end)]));
                
                    
                WidthFshift = max(fshift) - min(fshift);
                Xtxt = fshift(FCbPeakidx) + 0.05*WidthFshift;
                Y1txt = FFTMax; Y2txt = 0.95*FFTMax;
                FFT1txt = ['Peak frequency = ',num2str(fshift(FCbPeakidx))];
                FFT2txt = ['Period of oscillation = ',num2str(1/fshift(FCbPeakidx))];
                
                figure(ISFFTFig);
                plot([fshift(FCbPeakidx) fshift(FCbPeakidx)],[0 1.1*FFTMax],'k--','LineWidth',1); hold on;
                text(Xtxt,Y1txt,FFT1txt); hold on;
                text(Xtxt,Y2txt,FFT2txt); hold on;
                plot(fshift(2:end),abs(FLb(2:end)),'r'); hold on;
                plot(fshift(2:end),abs(FCb(2:end)),'g'); hold on;
                plot(fshift(2:end),abs(FRb(2:end)),'b'); hold off;
                legend('FFT |b|(t) Left Peak','FFT |b|(t) Center Peak','FFT |b|(t) Right Peak');
                
                
                
                "here's hoping the arrays are filled"
            end
        end
    %% Hodogram Video
        if boolHodogramVideo && (~boolDevelopment)
                                   
            utilities(AviFolderAbs, AviFolderRel); %Creates both Trial & avi folders if needed
            FrameRate = 10;
            VidFileName = ['\BFHodvideo_',TrialID,'.avi'];
            BFHodoVideoFile = [AviFolderAbs,VidFileName];
            BFHodoAnim =figure('Color', [1, 1, 1],'units','normalized','outerposition',[(NMonitors -1 +.65) 0.1 0.3 0.6]);
            set(BFHodoAnim,'doublebuffer','on');
            writerObjBF = VideoWriter(BFHodoVideoFile);
            writerObjBF.FrameRate = FrameRate;
            open(writerObjBF);

            if EVZoom
                X3min = EVL;X3max = EVR;% <-- use if passed through EVZoom
            else
                X3min = -L;X3max = +L;
            end
            tmp3L = abs(X4 - (X3min));
            tmp3R = abs(X4 - (X3max));
            [val, idx3L] = min(tmp3L);
            [val, idx3R] = min(tmp3R);
            N3 = idx3R - idx3L;
            X3 = linspace(X3min,X3max,N3);
            
            temp = BFlines(1);
            C = textscan(temp{1},'%f');
            Bk = C{1}.'; %B is an N/4 complex valued, double precision array            
            Bh = Bk(idx3L:(idx3R-1));
            Bmax = 1.5*max(abs([real(Bh), imag(Bh)]));
            BH_L = Bh(1); BH_R = Bh(end);
            % txtL = ['bL = ',num2str(abs(BH_L))];txtR = ['bR = ',num2str(abs(BH_R))];
            XL = real(BH_L)-Bmax*0.2;YL = imag(BH_L)-Bmax*0.05;
            XR = real(BH_R)-Bmax*0.1;YR = imag(BH_R)-Bmax*0.1;
            TrialData = split(TrialID,'_');
            TrialTxt = sprintf('TrialID %s_{%s}',TrialData{1},TrialData{2});
            for kk = 1:length(BFlines)
                temp = BFlines(kk*1);
                C = textscan(temp{1},'%f');
                Bk = C{1}.'; %B is an N/4 complex valued, double precision array

                Bh = Bk(idx3L:(idx3R-1));

                figure(BFHodoAnim);
                clf();
                set(gca,'nextplot','replacechildren');
                set(gcf,'Renderer','zbuffer');
                
                IzZoom = sum(imag(Bh))*(-dx4); %to be consistent with the KBW paper

                cla();
                
%                 evaluateKBWShockBoundaries = 1;
                if evaluateKBWShockBoundaries
%                     BZ_asymptote = 0;
%                     Ntest = 10;
%                     tol = 0.01;
                    [idRelMax, idxShockL, idxShockR, XShockL, XShockR, Izcalc] = IS_ShockBoundaries(Bh, Xh, b_L, BZasymptoteL, BZasymptoteR, Ntest, tol);

                    XCL = real(Bh(idxShockL)); YCL = imag(Bh(idxShockL));
                    plot(XCL,YCL,'o','MarkerFaceColor','green','MarkerSize',5);
                    hold 'on';

                    XCR = real(Bh(idxShockR)); YCR = imag(Bh(idxShockR));
                    plot(XCR,YCR,'o','MarkerFaceColor','green','MarkerSize',5);
                    hold 'on';

                end
                
                txtL = ['b_{L} = ',num2str(abs(BH_L))];txtR = ['b_{R} = ',num2str(abs(BH_R))];

                plot(real(Bh),imag(Bh));
                
%                 text(Bmax*(0.0), Bmax*(1.0),sprintf('Trial ID = %s',TrialID));
                text(Bmax*(0.0), Bmax*(1.0),TrialTxt);
                text(Bmax*(0.0), Bmax*(0.9),sprintf('dispersion = %4.4f',dispersion));
                text(Bmax*(0.0), Bmax*(0.8),sprintf('dissipation = %4.4f',nu));
                text(Bmax*(0.0), Bmax*(0.7),sprintf('alpha = %4.4f',alpha));
                text(Bmax*(-0.9), Bmax*(1.0),sprintf('Profile:%s',Profile));                
                text(Bmax*(-0.9), Bmax*(0.9),sprintf('line # = %4.0f',kk));
                text(Bmax*(-0.9), Bmax*(0.8),sprintf('time = %4.4f',TimePerLine*(kk - 1)));
                text(Bmax*(-0.9), Bmax*(0.7),sprintf('Iz = %4.4f',IzZoom));
                if sum(strcmp(PStr.PCon.keys,'theta_kbw')) > 0                  
                    text(Bmax*(-0.9), Bmax*(0.6),sprintf('\\Delta\\theta  = %4.4f rad',theta_kbw));
                end
                if sum(strcmp(PStr.PCon.keys,'CR13')) > 0                  
                    text(Bmax*(-0.9), Bmax*(0.5),sprintf('b_{kbw} = %4.4f',abs(Prfl('CR13'))));
                end
                text(XL, YL,txtL);
                text(XR, YR,txtR);
                axis([-1 +1 -1 +1]*Bmax);
               
               
                drawnow();              
                
                frame = getframe(gcf);
                writeVideo(writerObjBF,frame);
                
            end
        end
%% Evaluate Energy, FELandau etc.
Energy = 0;
if boolEnergy && (~boolDevelopment)
    Energy = Energy_Calc(Lh,Xh,Bh,BC);
end

FELandau = 0;
if boolFE && (~boolDevelopment)
    FELandau = FELandau_Calc(Lh,Xh,Bh,BC);
end
    %% Development    
        
    if boolDevelopment
        ['About to enter development.']
%         idL = 4020;idR = 4702;
%         Bh = B4(idL:idR);
%         Xh = X4(idL:idR);
        Nh = length(Bh);
%         Lsymm = [-1/2,1/2]*( X4(idR) - X4(idL) );
%         Xsym = linspace(Lsymm(1),Lsymm(2),Nh);
        lmda = 0.28;
        PhiMinus_Debug = 0;
        PhiPlus_Debug = 0;
        
        [xL_is_alphaXR] = alpha_XR(real(lmda),imag(lmda),abs(Bh(1)),abs(Bh(end)));        
        XR = ( Xh(end) - Xh(1))/(1 + xL_is_alphaXR);
        XL = - xL_is_alphaXR*XR;
        
        Xsym = linspace(XL,XR,Nh);
        
        [PhiML_condition, PhiFlag] =    PhiMinus(lmda ,Xsym(1)   , Bh(1)  , PhiMinus_Debug);
        [PhiPL_Inv_condition, PhiFlag] = PhiPlus(lmda ,Xsym(end) , Bh(end), PhiPlus_Debug); 
        
        [PhiML, PhiFlag] =    PhiMinus(lmda ,Xh(1)   , Bh(1)  , PhiMinus_Debug);
        [PhiPL, PhiFlag] = PhiPlus(lmda ,Xh(end) , Bh(end), PhiPlus_Debug); 
        
%         [S11_21, ErrFlag, Xuniform, b2nsymm, RelPctErr, ErrStep, StepsTaken, msg_S_terse, msg_S_verbose] = S_Adaptive_Complete(Xh , Bh , 1        , 3             , 0.2576   , 0.005);
% [S_Adapt1121, ErrFlag, Xuniform, b2nsymm,  RelPctErr, ErrStep, StepsTaken, msg_S_terse, msg_S_verbose] = S_Adaptive_Complete(Xas, Bas, FieldPrep, NFractionSTEPS, lambda_as, eps_as)

        %% Graph 'Resolved' b-field
        ErrStep = 1;
        StepsTaken = 1;
        StepSize(1) = 1;
        zmin = -zmax;
        step = 1;
        [ErrMaxAbs, Iabs] = max(abs(ErrStep)); %#ok<*SHVAI>
        BfieldFig = figure('Color', [1, 1, 1],'units','normalized','outerposition',[(NMonitors -1 +.02) 0.1 0.6 0.6],'visible','off');
        clf;
        
        plot(Xuniform,real(b2nsymm),'Color',[1 0 0]);
        hold('on');

        plot(Xuniform,imag(b2nsymm),'Color',[0 1 0]);
        hold('on');

        plot(Xuniform,abs(b2nsymm),'Color',[0 0 1]);
        hold('on');
        figure(BfieldFig);

        zmax = 1.5*max(abs(b2nsymm));
        zmin = -zmax;
        

        PosErrIndices = (ErrStep>=0);
        NegErrLine = ErrStep;
        PosErrLine = ErrStep;       
        PosErrLine(~PosErrIndices) = NaN; % Set the values you don't want to get drawn to nan
        NegErrLine(PosErrIndices) = NaN;

        plot([Xuniform(Iabs*step) Xuniform(Iabs*step)],[0.8*zmin, 0.8*zmax],'Color',[0 0 1]);
        hold('on');

        plot(Xuniform(step:step:end),(PosErrLine/ErrMaxAbs),'Color',[0.2 1 0.2],'LineStyle','--');
        hold('on');
        
        plot(Xuniform(step:step:end),(NegErrLine/ErrMaxAbs),'Color',[1 0.2 0.2],'LineStyle','--');
        hold('on');
        
        plot(Xuniform(step:step:end),(StepsTaken/StepSize(1)),'.','MarkerSize',3,'MarkerEdgeColor',[0 0 0]);
        hold('on');
        
        axis([Xuniform(1) Xuniform(end) StepsTaken zmin zmax]);
        set(gca, 'box', 'off');
        drawnow();
        
        step

    end
%% Real Eigenvalues    
NReEvs = 0; NReEvs2 = 0;
ReEvs = 0; ReEvs2 = 0;
Signs1p = 0; Signs1p2 = 0;
SFlag = 1;
    if boolREV && (~boolDevelopment)
        %Evaluate real eigenvalues
        %If the # of eigenvalues in the range is comparable to NEval, some 
        %eigenvalues may well be overlooked. 
        %By construction the method used here returns just the real part.
        %The tacit assumption is that the imaginary part is 0 or at least
        %negligible.
        ['About to evaluate real eigenvalues.']
        
        NFracSS = 0;       % # of fractional steps to use in S_Adaptive_Complete
        eps_AS = (1 + 0.2/100)^1-1;  % Relative error enforced in scattering matrix output of S_Adaptive_Complete
        SMatrixMethod = 2; % 1 to use S_Adaptive_Complete
                           % 2 to use S11_S21


        [ReEvs    , Signs1p, NReEvs, Xrange, Bconditions, Xh, Bh, SFlag, s11, is21] = RealEVs(LREval, NREval, NBSIter, Xh, Bh,  NFracSS, eps_AS, SMatrixMethod, opts);
        [EnergyReEV] = Energy_Calc(Lh,Xh,Bh,abs(Bh(1)));
        [FELandauReEV] = FELandau_Calc(Lh,Xh,Bh,abs(Bh(1)));
        RealEV_List = [EnergyReEV, FELandauReEV, ReEvs];
        ;
      % [ReEvFinal, Signs1p, NReEvs, Xrange, Bconditions, SFlag] = RealEVs(LREval, NEval, NBSIter, Xh, Bh, NFracSS, eps_AS, SMatrixMethod, opts)
%         ['Real Eigenvalues: ']
%         ReEvs
%         Signs1p
% Xrange

    end
%% Complex Eigenvalues  
    boolHalley = 0;
    if boolRoot && (~boolDevelopment)
        boolNewton      = 0;
        boolHalley      = 1;
        boolHalley24    = 0;
%         zev = 0.837 + 1i*0.3;% the actual root
        zev =zinit;% the actual root
        z0 = zinit;% the initial guess
%         ['About to execute Halley 2.4:']
%         CEVList(1,1) = 0;
%         z0 = 0.4 + 1i*0.4;% initial guess
%         ['z',num2str(0),' = ', num2str(real(z0)), ' + i(',num2str(imag(z0)),')']
% 
%         del = (zev - z0)/10;
% 
%         for ni = 1:5
% 
%             foo = S11_S21(z0,Lh,Xh,Bh,B0,theta_M,theta_P,opts)*[1,0]';
%             fro = S11_S21(z0 + del,Lh,Xh,Bh,B0,theta_M,theta_P,opts)*[1,0]';
%             fpo = (fro - foo)/(del);
%             
%             y1 = z0 - foo/fpo;
%             
%             fo1 = S11_S21(y1,Lh,Xh,Bh,B0,theta_M,theta_P,opts)*[1,0]';
%             fr1 = S11_S21(y1 + del,Lh,Xh,Bh,B0,theta_M,theta_P,opts)*[1,0]';
%             fp1 = (fr1 - fo1)/(del);
%             
%             den = 2*foo*fp1^2 - fpo^2*fo1 + fpo*fp1*fo1;
%             num = 2*foo*fo1*fp1;
% 
%             z1 = y1 - num/den;
%             ['z',num2str(ni),' = ', num2str(real(z1)), ' + i(',num2str(imag(z1)),') with Dz = ', num2str(real(z1-z0)), ' + i(',num2str(imag(zev-z0)),') and percent error = ', num2str(abs((z1 - zev)/zev)*100)]
%              z0 = z1;             
%         end
%         CEVList(1,1) = z0;
%         ['Exit Halley 2.4 ... Complex eigenvalues found are:']
%         CEVList
        
        if boolNewton
            NewtTime = cputime;
            ['About to execute Newton:']
            CEVList = cell(1,3);
            Nnewt = 40; % # of iterations to use
            [CEVList] = Newton(z0,zev,Lh,Xh,Bh,BC,opts, Nnewt);
            ['Exit Newton ... cpu time = ',num2str(cputime - NewtTime)]
            CEVList
        end
        
        if boolHalley
            HalleyTime = cputime;
            ['About to execute Halley:']
            CEVList = cell(1,3);
            NHalley = 25; % # of iterations to use in Halley Method
            [CEVList, HalleyVerbose] = Halley(z0,zev,Lh,Xh,Bh,B0,theta_M,theta_P,opts,NHalley, LComplexEval);
            ['z',num2str(0),' = ', num2str(real(z0)), ' + i(',num2str(imag(z0)),')']         
            ['Exit Halley ... cpu time = ',num2str(cputime - HalleyTime)]
            ['hello']
            %CEVList
        end
        
         if boolHalley24
             Halley24Time = cputime;
            ['About to execute Halley24:']
            CEVList = cell(1,3);
            NHalley24 = 25; % # of iterations to use in Halley Method
            [CEVList] = Halley24(z0,zev,Lh,Xh,Bh,B0,theta_M,theta_P,opts,NHalley24);
            ['z',num2str(0),' = ', num2str(real(z0)), ' + i(',num2str(imag(z0)),')']         
            ['Exit Halley24 ... cpu time = ',num2str(cputime - Halley24Time)]
            %CEVList
        end
        
        
    end
    if (boolCEV || boolCEVgraph) && (~boolDevelopment)
        %Evaluate Complex Eigenvalues
        % 1) Evaluate Scattering Coefficient S11 on grid in complex plane
        % 2) Analyze signs of s11 on grid to flag regions containing zeros
        
            %Evaluate Scattering Coefficients, s11 and s21 at grid points
            ['About to evaluate Scattering Data Grid']
            BC = [Bh(1), Bh(end)];
            
            CEflag = 1;     % 1 to evaluate Complex Eigenvalues
            NFracSS = 0;       % # of fractional steps to use in S_Adaptive_Complete
            eps_AS = (1 + 0.5/100)^4-1;  % Relative error enforced in scattering matrix output of S_Adaptive_Complete
            SMatrixMethod = 1; % 1 to use S_Adaptive_Complete
                               % 2 to use S11_S21
            Rfactor = 2*10^(-4);
            
            CEcount = 0;
            FirstPass = 1;
            WayBack = 0;
            LComplexEvalNow = LComplexEval;
            FieldPrepped = 0;
            while CEflag > 0 % evaluate signs of Re(s11(l)) and Im(s11(l)) within LComplexEval [ Re(l)min,Re(l)max; Im(l)min,Iml)max]
                CEZoomflag = 1; % 1 to allow user to successively zoom on a region of the complex plane
                CEflag = 0;
%                 Sval = ScatteringDataGrid(NCEval, LComplexEvalNow,Lh, Bh, Xh,BC,opts);
                if FirstPass < 3
                    [ Sval, SAtime, S11Err, SAErrFlag, TotalSteps, FieldPrepped, Lh, Bh, Xh] = ScatteringDataGrid(NCEval, FieldPrepped, LComplexEvalNow,Lh, Bh, Xh,BC, SMatrixMethod, NFracSS, eps_AS, Rfactor, opts) ;
                    %[ Sval, SAtime, S11Err, SAErrFlag, TotalSteps, FieldPrepped, Lh, Bh, Xh] = ScatteringDataGrid(NCEval, FieldPrepped, LComplexEval,Lh, Bh, Xh, BC, SMatrixMethod, NFracSS, eps_AS, opts)
                    
%                     WayBack = 0;
                end
                [CevState,printgrid, Xr, Yr, Xg, Yg, Xb, Yb, Xw, Yw] = CevStateMatrix(LComplexEvalNow,NCEval,Sval);
                close all;

                S11ErrFig = figure('units','normalized','outerposition',[(NMonitors -1 +.00) 0.05 0.3 0.5]);%,'title','Error as fraction of Tolerance'); % create figure to display color map               
                surf(linspace(LComplexEvalNow(1),LComplexEvalNow(3),NCEval),linspace(LComplexEvalNow(2),LComplexEvalNow(4),NCEval),S11Err/(100*eps_AS));
                tS11ErrFig = title('Error/Tolerance');
                
                TotalStepsFig = figure('units','normalized','outerposition',[(NMonitors -1 +.7) 0.05 0.3 0.5]);%,'title','Total Number of Steps for Evaluation'); % create figure to display color map               
                surf(linspace(LComplexEvalNow(1),LComplexEvalNow(3),NCEval),linspace(LComplexEvalNow(2),LComplexEvalNow(4),NCEval),TotalSteps);
                tTotalStepsFig = title('Total # Steps for Evaluation');
                
                CEcount = CEcount + 1;
                ['Completed Scattering Data Grid evaluation']
                %Search through values of s11 and determine the 'state' of each value
                %on the grid:
                %   if Re(s11) > 0 and Im(s11) > 0 --> 'PP' and color = red and factor = 2
                %   if Re(s11) > 0 and Im(s11) < 0 --> 'PN' and color = green and factor = 3
                %   if Re(s11) < 0 and Im(s11) > 0 --> 'NP' and color = blue and factor = 5
                %   if Re(s11) < 0 and Im(s11) < 0 --> 'NN' and color = white and factor = 7
                % Output is matrix with 'state' of each eigenvalue
                %Note that the X,Y arrays are required for execution of printgrid
%                 if ~WayBack
%                     [CevState,printgrid, Xr, Yr, Xg, Yg, Xb, Yb, Xw, Yw] = CevStateMatrix(LComplexEval,NCEval,Sval);
%                     WayBack = 0;
%                 end
                if CEcount == 1
                    LComplexEval_original = LComplexEval;Sval_original = Sval;
                end

                if boolCEVgraph % display color map of Re(s11(lambda) vs. Im(s11(lambda)) within LComplexEval [ Re(l)min,Re(l)max; Im(l)min,Iml)max]
                    HCEV = figure('units','normalized','outerposition',[(NMonitors -1 +.3) 0.05 0.4 0.8]); % create figure to display color map
                    cla();
                    eval(sprintf(printgrid));
                    drawnow();
                    
                    CursorCheck = 1;
                    
                    while CursorCheck
                        figure(HCEV);
                        set(gcf, 'WindowButtonMotionFcn',@mouseMove);
                        CCtitle = get(gca, 'title');
                        Cvalue = str2num(CCtitle.String);
                        if real(Cvalue) > LComplexEvalNow(1,2)
                            CursorCheck = 0;
                        end
                        pause(0.1);
                    end
                    
    ['FirstPass = ',num2str(FirstPass)]
                    while CEZoomflag == 1
                        [x,y] = ginput(2);
                        x = sort(x);y = sort(y);
%                         LComplexEval_old = LComplexEvalNow;

                        figure(HCEV);
                        cla;
                        eval(sprintf(printgrid));
                        hold('on');
                        plot([x(1) x(1)],[y(1) y(2)],'Color',[1 0 0]);
                        hold('on');
                        plot([x(2) x(2)],[y(1) y(2)],'Color',[0 1 0]);
                        hold('on');
                        plot([x(1) x(2)],[y(1) y(1)],'Color',[0 0 1]);
                        hold('on');
                        plot([x(1) x(2)],[y(2) y(2)],'Color',[0 0 0]);
                        hold('on');

                        drawnow();
                         ['FirstPass = ',num2str(FirstPass)]
                        if FirstPass == 1
                            % Construct a questdlg with three options
                            choice = questdlg('Selection correct?', ...
                            'Select region in Complex EV plane to zoom in on:', ...
                            'Yes','Find Root','Try Again','Get Me Outa Here');
                            % Handle response
                            switch choice
    %                         case 'Zoom Out'
    %                             disp([choice ' Try, try again.'])
    %                             printgrid = printgrid_old; LComplexEval = LComplexEval_old;
    %                             CEflag = 1;
    %                             CEZoomflag = 0;
                            case 'Find Root'
                                boolRoot = 1;
                                zmid = 1/2*(LComplexEvalNow(1,1) + LComplexEvalNow(1,2)) + 1i/2*(LComplexEvalNow(2,1) + LComplexEvalNow(2,2));
                                Nnewt = 10; % # of iterations to use
                                [CEVList] = Newton(zmid,zmid,Lh,Xh,Bh,BC,opts, Nnewt);
                                CEVList
                                CEZoomflag = 3; % CEZoomflag ~= 1, kick out of gui loop
                                CEflag = 1; % cycle back to gui with no change to LComplexEvalNow
                                FirstPass = 3; % on new gui cycle, do not update Sval or printgrid or LComplexEvalNow
                                 ['FirstPass = ',num2str(FirstPass)]
                            case 'Yes'
                                disp([choice ' Woot!'])
                                LComplexEvalPrevious = LComplexEvalNow;
                                SvalPrevious = Sval;
                                LComplexEvalNow = [x(1), x(2); y(1), y(2)]; %range of real values to study ... [Real Min, Real Max; Im Min, Im Max]
%                                 LComplexEval_old = LComplexEval;
                                CEZoomflag = 2; % kick out of gui loop
                                CEflag = 1; % cycle back to gui
                                FirstPass = 2; % up date Sval
                                 ['FirstPass = ',num2str(FirstPass)]
                            case 'Try Again'
                                disp('Prepare for Transport.')
                                TryAgain = 1;
                                CEflag = 1; % CEflag = 1 --> cycle trhough display again, no change to LComplexEvalNow
                                FirstPass = 3; % cycle back to gui with no change to Sval or printgrid
                                CEZoomflag = 2; % kick out of gui loop
                            end
                        else
                                                        % Construct a questdlg with three options
                            choice = questdlg('Selection correct?', ...
                            'Select region in Complex EV plane to zoom in on:', ...
                            'Yes','Find Root','Zoom Out','Get Me Outa Here');
                            % Handle response
                            switch choice
                            case 'Zoom Out'
                                disp([choice ' Try, try again.'])                                                         
                                if LComplexEvalNow == LComplexEvalPrevious
                                    LComplexEvalNow = LComplexEval_original;Sval = Sval_original;
                                else
                                    LComplexEvalNow = LComplexEvalPrevious; Sval = SvalPrevious;
                                end
                                CEflag = 1; % repeat gui loop and ...
                                FirstPass = 3; % revert back to original Sval, printgrid and LComplexEval
%                                 WayBack = 1;
                                 ['FirstPass = ',num2str(FirstPass)]
                                
                                CEZoomflag = 2; % kick out of gui loop
                            case 'Find Root'
                                boolRoot = 1;
                                zmid = 1/2*(LComplexEvalNow(1,1) + LComplexEvalNow(1,2)) + 1i/2*(LComplexEvalNow(2,1) + LComplexEvalNow(2,2));
                                Nnewt = 10; % # of iterations to use
                                FirstPass = 3;
                                 ['FirstPass = ',num2str(FirstPass)]
                                [CEVList] = Newton(zmid,zmid,Lh,Xh,Bh,BC,opts, Nnewt);
                                CEVList
                                CEZoomflag = 3; % kick out of gui loop and end Complex evaluation
                            case 'Yes'
                                disp([choice ' Woot!'])
                                LComplexEvalPrevious = LComplexEvalNow;
                                SvalPrevious = Sval;
                                LComplexEvalNow = [x(1), x(2); y(1), y(2)]; %range of real values to study ... [Real Min, Real Max; Im Min, Im Max]
                                CEZoomflag = 2;
                                FirstPass = 2; % kick out and reevaluate printgrid, Sval with updated LComplexEvalNow
                                 ['FirstPass = ',num2str(FirstPass)]
                                CEflag = 1; % start cycle over with updated values
                            case 'Try Again'
                                disp('Prepare for Transport.')
                                TryAgain = 1;
                                FirstPass = 4; % no up date, keep last values of printgrid, LComplexEvalNow and Sval
                                 ['FirstPass = ',num2str(FirstPass)]
                                 CEZoomflag = 2;
                                CEflag = 1;
                            end
                            
                        end

%                         % Construct a questdlg with three options
%                         choice = questdlg('Selection correct?', ...
%                         'Select region in Complex EV plane to zoom in on:', ...
%                         'Yes','Find Root','Zoom Out','Get Me Outa Here');
%                         % Handle response
%                         switch choice
%                         case 'Zoom Out'
%                             disp([choice ' Try, try again.'])
%                             printgrid = printgrid_old; LComplexEval = LComplexEval_old;
%                             CEflag = 1;
%                             CEZoomflag = 1;
%                         case 'Find Root'
%                             zmid = 1/2*(LComplexEval(1,1) + LComplexEval(1,2)) + 1i/2*(LComplexEval(2,1) + LComplexEval(2,2));
%                             Nnewt = 10; % # of iterations to use
%                             [CEVList] = Newton(zmid,zmid,Lh,Xh,Bh,BC,opts, Nnewt);
%                             CEVList
%                             CEZoomflag = 3;
%                         case 'Yes'
%                             disp([choice ' Woot!'])
%                             LComplexEval = [x(1), x(2); y(1), y(2)]; %range of real values to study ... [Real Min, Real Max; Im Min, Im Max]
%                             CEZoomflag = 2;
%                             CEflag = 1;
% %                         case 'Try Again'
% %                             disp('Prepare for Transport.')
% %                             TryAgain = 4;
% %                             CEflag = 0;
%                         end
                    end
                end
            end
         
            
            %Search through the 'state matrix' and record indices of all those points
            % that lie at the center of a 3x3 grid which contains all of the four
            % states; remove multiple entries for the same eigenvalue; 
            %Store resulting complex eigenvalues in the 1xN,complex valued array
%             CEVList = ComplexEVFinder(LComplexEval,NCEval,CevState);
%             CEVList
            
            
    end
%% Radiation   
    if boolRAD && (~boolDevelopment)
        %Evaluate continuous scattering data, Rho = s21/s11 along branch
        %cuts in 1st quadrant
        % RhoR --> values above real branch
        % RhoI --> values to right of imaginary branch
        bomin = min(abs([Bh(1), Bh(end)]));
        LReBranchEval = [1.01*bomin, 1.2];
        [LReBvals, RhoR, AbsRhoR, LImBvals, RhoI, AbsRhoI] = Rho(NReBranch, NImBranch, LReBranchEval, LImBranchEval, RBeps,IBeps, Lh, Bh, Xh,B0,theta_M,theta_P,opts); 
        
        if boolRADgraph
            HRAD = figure('units','normalized','outerposition',[(NMonitors -1 +.75) 0.0 0.25 0.5]);
            subplot(2,1,1);
            plot(LReBvals, real(RhoR),'--go',LReBvals, imag(RhoR),':r+',LReBvals, AbsRhoR,'-b');
            title('rho: real branch');hold on;

            subplot(2,1,2);
            plot(LImBvals, real(RhoI),'--go',LImBvals, imag(RhoI),':r+',LImBvals, AbsRhoI,'-b');
            title('rho: imaginary branch');drawnow;
        end
    end
%% Update BoolCon
% In case boolRoot set = 1 through user choice via graphing of complex
% plane
    boolKeys    = {'title'       , 'boolBgraph', 'boolHodogram', 'boolVideo',...
        'EVZoom', 'boolEnergy', 'boolREV', 'REV_ErrFlag','boolRoot', 'zinit', 'boolCEV', ...
        'boolCEVgraph', 'boolRAD', 'boolRADgraph', 'boolHalley'};
    boolValues  = {'User Choices',  boolBgraph ,  boolHodogram , ...
        boolVideo, EVZoom, boolEnergy, boolREV,SFlag,boolRoot, zinit, boolCEV, ...
        boolCEVgraph, boolRAD, boolRADgraph, boolHalley};
    BoolCon     = containers.Map(boolKeys, boolValues);
%% Output   
    %clc;
    debugDNLS_Eigenvalues = 0;
    if debugDNLS_Eigenvalues==1 
        ['About to write eigenvalue results to file:'] 
    end
    if (boolREV || boolCEV || boolRoot || boolEnergy) && (~boolDevelopment)
        msgHalleyFinal = 'boolHalley not executed.';
        if boolRoot && boolHalley
            msgHalleyFinal = HalleyVerbose;
        end
        DNLS_EigenValues(PStr, BoolCon, EvalCon, ArchiveBase, TrialFolderRel, ...
            TrialNumber, RunNumber, LineNumber, ...
            Energy, Signs1p, NReEvs, ReEvs, LREval, CEVList, debugDNLS_Eigenvalues, ...
            msgHalleyFinal, Xh, Bh, dispersion, alpha);
        % if boolREV && ()
        %     Real_EV_Excel_Archive(PStr,ReEvs,t);
        % end
    end
    if debugDNLS_Eigenvalues==1 
        ['Finished call to write eigenvalues to file.']
    end
end