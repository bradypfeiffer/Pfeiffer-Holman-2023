function [PSTR, growth, decay, timestep] = DNLS_v4(profile, PStr)
%% New in v4:
% A symmetric Fourier split step method is implemented following
% Ralf Deiterding and Stephen W. Poole,
% 'Robust split-step Fourier methods for simulating propagation of
% ultra-short pulses in single- and two-mode optical communication fibers,'
% (https://www.dropbox.com/s/fb3lf4203a8ualt/Deiterding_Poole_Robust_split-step_Fourier_methods.pdf?)
% Specifically their Section 3.2 gives sufficient detail for implementation
% of the linear step in Fourier space.
% Section 3.3 provides a detailed discussion of implementation of the
% nonlinear step with the benefit of a stability bound relating dt, dx and
% wave amplitude and, likely, gradients.
%% ToDo Archive    

% eventually use these, or similar, parameters for diagnostics and control    
        growth = 1.2345; 
        decay = 0.98765;
    %   timing
        t = cputime;   
        clear b;
    
    %This 1D version handles aliasing by setting the initial value of N with intent of N/4
    %providing sufficient resolution, then setting high-k elements of FFT
    %to zero each time an FFT is computed. However this is not necessary
    %within the leapfrog method in FFT space as there is no possibility of
    %introducing a non-zero FFT element.
    
    % ToDo:  
    %-1) --> Create dual 'Writer' & 'Reader' Codes with efficient I/O method 
    %       DONE        'Writer' = Find b(x,t) given b(x, t = 0) and archive b(x,t)
    %       IN PROCESS  'Reader' = Read in archived b(x,t) and evaluate eigenvalues etc.
    %
    %0)     --> Create a test for "success" or instability and flag to
    %       abort trial before fatal error occurs.
    %       -->  Post the target number of iterations ... or time steps?
    %       ... to take for a given trial
    %       --> in post diagnostics create test to adjust 'Courant' towards
    %       stability.
    %
    %1) DONE --> Create a class 'PStr' to pass trial parameters back
    %and forth between 'Call_DNLS_vx' and 'DNLS_vx'.
       
    %As time allows:
    %1)     --> Matlab handles all the 'N' k-element calculations in one command,
    %but only N/4 of these are necessary. It is worth checking if
    %implementing this 'time stepping' line in a for-loop over the
    %non-trivial k-elements might speed up the execution.
    
    %2)     --> Consider helpful annotations such as: total time, time step
    %number, integration parameters (N, L, Courant), initial profile
    %parameters (1p soliton, bright/dark, eigenvalue), other?
    
    %3)     --> Are there more efficient methods than Leap Frog? Surely
    %there are, but this is not mission critical at the moment.
    %tic
    
    %4)     --> there is a significant time lag on the first 'drawnow()'
    %call (measured by 'tic' & 'toc')

    %********************************************************************
%% Function Definitions
    function[idMax, idSL, idSR, xShockL, xShockR, Iz,OscLHS] = IS_ShockBoundaries_old(bb,x,bL, basymptote, ntest, tol)
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
            % ntest --> number of indices over which the running average of Bz
            % is evaluated
            
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
            
        %% Find first peak of Im(bb): idMax starting point for search for both LHS & RHS
            % Seems to be useful if the KBW profile is near to being coplanar
            [val, idMax] = max(abs(imag(bb)));
        %% Find LHS of Shock: idSL, xShockL 
        
            % Estimate number of oscillations on LHS of Shock          
            OscLHS = 0;
            DeltaLHS = abs(bb(1)) - (bL);
            signLhsOld = sign(DeltaLHS);
            for ilhs = 2:idMax
                DeltaLHS = abs(bb(ilhs)) - bL;
                signLhsNew = sign(DeltaLHS);
                if ((signLhsNew*signLhsOld) < 0) && (abs(DeltaLHS) > 5*tol) % not counting crossings possibly due to numerical error
                    signLhsOld = signLhsNew;
                    OscLHS = OscLHS + 1;
                end
            end
            
            if OscLHS < 2       
                %assumes shock has a constant B-field on its LHS with B = (bL, 0)
                %and does not have oscillations approaching Shock from LHS
                
                npasses = 0;
                SATL=0;
                while length(SATL) < floor(idMax/10)
                    npasses = npasses + 1;
                    LHS = abs(real(bb(1:idMax)-bL)) < npasses*tol; %Returns
                    %this returns an array of length idMax with logical '1' in every 
                    %location where real(b) > bL '0' otherwise
    %                 LHSNEG = abs(real(bb(1:idMax)-bL)) < tol; %returns '1' where b is below bL                
                    SATL = find(LHS ==1); %returns an array containing an order list of indices where the condition above holds
                end
                idSL = SATL(end); % returns the rightmost index where the condition holds
                xShockL = x(idSL); % returns the position of the left edge of the shock
            elseif OscLHS > 1
                Nb = length(bb);
                idSL = 1;
                for idb = 1:(idMax - ntest + 1)
                    idxLeft = idMax-(idb-1)-(ntest - 1);
                    idxRight = idMax-(idb-1);
                    BmagAvg = sum(abs(bb(idxLeft:idxRight)))/ntest;
                    signLhsidb = sign(abs(bb(idb)) - bL);
                    signLhsidbPlusone = sign(abs(bb(idb+1)) - bL);
                    if ( signLhsidb*signLhsidbPlusone) < 0 && (abs(BmagAvg - bL ) < tol*abs(bL))
                        idSL = idb;
                        break;
                    end
                end
                xShockL = x(idSL);
            end            
        %% Find RHS of Shock: assumes Im(b) has general form of damped sign wave towards asymptote
            %Assuming Coplanar shock and imag(Bh) --> 0 to right of shock
            Nb = length(bb);
            idSR = Nb;
            for idb = idMax:Nb
                Bzavg = sum(imag(bb(idb-ntest:idb)))/ntest;
                if ( ( imag(bb(idb-1)) - basymptote )*( imag(bb(idb)) - basymptote )) < 0 && (abs(Bzavg - basymptote ) < tol*abs(b_L))
                    idSR = idb;
                    break;
                end
            end
            xShockR = x(idSR);

            dxx = (x(2) - x(1));
        %% Evaluate Iz
            Iz = sum(imag(bb(idSL:idSR)))*dxx;
    end
    function[idMax, idSL, idSR, xShockL, xShockR, Iz] = IS_ShockBoundaries_v2(bb,x,bL, basymptoteL, basymptoteR, ntest, tol)
        
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
        for idb = 1:(idMax - 1)
            idx = idMax-(idb-1);
            DeltaBz = abs(imag(bb(idx)) - basymptoteL);
            
            if DeltaBz < minDeltaBz
                minDeltaBz = DeltaBz;
                Lexception = 0;
            else
                Lexception = Lexception +1';
            end
            
            if Lexception > 50
                msg = 'Not monotonic'
                break;
            end
            
            if minDeltaBz < tol*abs(b_L) && (idx > ntest)
                BzAvg = sum(abs(imag(bb(idx - (ntest - 1)):idx) - basymptoteL))/ntest;
                if abs(BzAvg) < tol*abs(b_L)
                    idSL = idx;
                    BLasymptoteSuccess = 1;
                else
                    msg = 'Not monotonic'
                end
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
        for idb = idMax:Nb
            DeltaBz = abs(imag(bb(idb)) - basymptoteR);
            if DeltaBz < minDeltaBz
                minDeltaBz = DeltaBz;
                Rexception = 0;
            else
                Rexception = Rexception +1;
            end
            
            if Rexception > 5
                msg = 'Not monotonic'
                break;
            end
            
            if minDeltaBz < tol*abs(b_L)
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
        for idb = (ntest - 0):(idMax - (ntest - 1))
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

%% Extract Parameters from PStr
    % 'Numerics', 'Physics', 'Graphics', 'Trial', 'Driver' and 'Folder' Containers
    Numerics = PStr.NCon;
    Prfl = PStr.PrflCon;
    Physics  = PStr.PCon;
    Graphics = PStr.GCon;
    Trial    = PStr.TCon;
    Folder   = PStr.FCon;
    Driver   = PStr.DrvrCon;
    
        Profile = Prfl('Profile');
        PrflKeys = Prfl.keys;
        PrflValues = Prfl.values;
    
        N = Numerics('N');
        L = Numerics('L');
        dt = Numerics('dt');
        dx = Numerics('dx');
        if sum(strcmp(Numerics.keys, 'Lva0')) > 0
            a0 = Numerics('Lva0');
        else
            a0 = Numerics('Lva0');
        end
        if sum(strcmp(Numerics.keys, 'L2nd')) > 0
            c0 = Numerics('L2nd');
        else
            c0 = Numerics('L2nd');
        end
        B0 = Physics('bo');
        nu = Physics('nu');
        gamma_Landau = Physics('gamma_Landau');
        dispersion = Physics('dispersion');
        alpha = Physics('alpha');
        
       
        TimestepAttempted = Numerics('TimestepsAttempted');
        timestepsPerFrame = Numerics('timestepsPerFrame');
        TotalTime         = Numerics('TotalTime');
        
        Numericskeys   = Numerics.keys;
        Numericsvalues = Numerics.values;
        NLinesSavedIdx = find(contains(Numericskeys,'NumberLinesSaved'),1);
        NumberLinesSaved = round( TimestepAttempted/timestepsPerFrame +.5);
        if NLinesSavedIdx > 0
            NumberLinesSaved = Numerics(Numericskeys{NLinesSavedIdx});
        end

        
        boolShockFit = 0; %default value. Reset when parameters are read from excel file ... if ... KBW profile
        CompressionRatio = 0.8586;
        BZasymptoteR = 0;
        BZasymptoteL = 0;
        Ntest = 10;
        tol = 0.01;
        
        IsDriven = Driver('IsDriven');
        if IsDriven == 1
            DriverAmplitude = Driver('Amplitude');
            K_driver = Driver('K');
            X_driver = Driver('X');
            W_Driver = Driver('W');
            
            DriverRange = @(x, X, K, W) 1/2*( tanh( (x + X)/W) - tanh( (x - X)/W) );
            
            Circular_Driver = @(x, t, v, A, X, K, W) A*DriverRange(x, X, K, W).*exp( 1i*K*(x - v*t));
        end
        if strcmp(Profile, '1psoliton')        
            lambda = Prfl('lambda');
            onePsign = Prfl('sign');
        end

        if strcmp(Profile, 'double_1psoliton')
            lambda1 = Prfl('lambda1');
            lambda2 = Prfl('lambda2');
            onePsign1 = Prfl('sign1');
            onePsign2 = Prfl('sign2');
        end

        if strcmp(Profile, '2psoliton')
            lr = Prfl('lr');
            li = Prfl('li');
        end
        if strcmp(Profile, 'Helical')
            %{'title', 'Profile','W_AS', 'W_0', 'L_AS', 'K_Helix', 'helix_angle','xL', 'xR','EVL', 'EVR'};
            W_AS = Prfl('W_AS');
            W_0 = Prfl('W_0');
            L_AS = Prfl('L_AS');
            K_Helix = Prfl('K_Helix');
            helix_angle = Prfl('helix_angle');
            theta_p = Prfl('theta_p');
            xL = Prfl('xL');
            xR = Prfl('xR');
            EVL = Prfl('EVL');
            EVR = Prfl('EVR');
        end
        if strcmp(Profile, 'FastShock')
        end
        if ~isempty(strfind(Profile,'Hybrid')) && ~isempty(strfind(Profile,'KBW'))
            idxLtemp = Prfl('IdxISLeft');
            idxRtemp = Prfl('IdxISRight');
            idxISDelta = idxRtemp - idxLtemp;
            idxL = idxLtemp - round(idxISDelta*0.3);
            idxR = idxRtemp + round(idxISDelta*0.1);
        end
        if contains(Profile,'KBW')
            idxL = 1;
            idxR = round(N/2);
            Prflkeys = Prfl.keys;
            Prflvalues = Prfl.values;
            ShockFitIdx = find(contains(Prflkeys,'ShockFit'),1);

            if ShockFitIdx > 0
                boolShockFit = Prfl(Prflkeys{12});
            else
                 boolShockFit = 0;
            end
%             boolShockFit = Prfl('boolShockFit');
            theta_kbw = Prfl('theta_kbw');
            bkbw = Prfl('b_kbw');
            b_L = Prfl('b_L');
            xkbw = Prfl('Xkbw');
            CompressionRatio = abs(Prfl('CR13'));
            BZasymptoteR = CompressionRatio*b_L*sin(theta_kbw);
            BZasymptoteL = 0;
            xkbw = Prfl('Xkbw');
            b_kbw = Prfl('b_kbw');
            b_L = Prfl('b_L');
        end
        if strcmp(Profile,'Generated13IS')

            idxL = 1;
            idxR = round(N/2);
            boolShockFit = Prfl('boolShockFit');
            theta_kbw = Prfl('theta_IS');
            b_kbw = abs(Prfl('CR13'));
            b_L = Prfl('b_L');
            xIS = Prfl('xIS');
            CompressionRatio = Prfl('CR13');
            BZasymptoteR = abs(CompressionRatio)*b_L*sin(theta_kbw);
            BZasymptoteL = 0;
        end
        TimestepStart = 0;
        LineNumberStart = 1;
        PrflMatches = strfind(PrflKeys,'ProfileNew');
        PrflContinuationBool = any( vertcat( PrflMatches{:} ) );
        if PrflContinuationBool
            TimestepStart = Numerics('TimestepStart');
            LineNumberStart = Numerics('LineNumberStart');
            
            
        end
        
        
        if strcmp(Profile, 'double_1psoliton')
            v1 = Prfl('speed1');
            v2 = Prfl('speed2');
            v = v1;
            if v2 > v1
                v = v2;
            end
        else
            v = Prfl('speed');
        end
        
        exact = Graphics('DisplayExactSolution');
        zmax  = Graphics('zmax');
        zmin  = Graphics('zmin');
        SystemGraphics = groot;
        NMonitors = length(SystemGraphics.MonitorPositions(:,1)); % returns the number of connected monitors.
        
        TrialID = Trial('TrialID');
        TrialName = Trial('TrialName');
        TimeStamp = Trial('TimeStamp');
        RunNumber = Trial('RunNumber');
        
        BFfile = Folder('Bfieldfile');
        CreateMatfile =  Folder('CreateMatfile');
        if CreateMatfile == 1            
            BFfileMat = Folder('BfieldMatfile');
            BFMat = matfile(BFfileMat,'Writable',true);
        end
        % Check if BFfileMat exists, delete if so
%         if exist(BFfileMat,'file')
%             delete(BFfileMat);
%         end
%         % Create BFfileMat
        
%         BFMat.PStr(str2num(RunNumber),1) = PStr;
%         save BFMat;
%% Initialize FFT    
    X = linspace(-L,L,N);
    Kx = [0:N/2-1,-N/2:-1] * pi/L;
    b = profile(Prfl, Physics, X, 0);
    
%     bp = ifft(Kx.*fft(b));
%     [val, idL] = min(abs(xkbw-0.1*L - X));
%     [val, idR] = min(abs(xkbw+0.1*L - X));
%     [val, idMax] = max(abs(bp(idL:idR)));
    
    if length(b) > N
        btemp = b(1:round(length(b)/N):end);
        b = btemp;
    end
    
    
    
    nanidx = isnan(b);      %Soliton Profiles may return NaN values for large x --> nanidx is cell array with [0] if numeric and [1] if NaN
    NZnan = find(nanidx > 0);
    if (~isempty(strfind(Profile, '1psoliton'))  || ~isempty(strfind(Profile, '2psoliton'))) && ~isempty(NZnan)  
        % if profile is a soliton and some values are NaN
        NZnan = find(nanidx > 0);
        if ~isempty(strfind(Profile, '1psoliton')) 
            %peak a max (min) if bright (dark) soliton
            if onePsign > 0
                [bMax, bcIdx] = max(abs(b)); %bcIdx is index of peak value of bright soliton
            else
                [bMin, bcIdx] = min(abs(b)); %bcIdx is index of min value of dark soliton
            end
        else
            [bMax, bcIdx] = max(abs(b)); %bcIdx is index of peak value of 2p-soliton
        end

        PrflResolution = (length(find(abs(b) > 1.05*B0)) + length(find(abs(b) < 0.95*B0))); % est. of # of samples across soliton
        if ((bcIdx + PrflResolution) < N) && (NZnan(1) > (bcIdx + PrflResolution)) % Check that indices are valid and consistent
            BGoodbc = find(abs(abs(b(bcIdx + PrflResolution:NZnan(1))) - B0)< 0.001); 
            % find locations to right of soliton peak that do not match asymptotics
            b(BGoodbc(end) + bcIdx:NZnan(1)) = B0;
        end
    end
    
    
    b(nanidx) = B0;         %Set b = bo which is the asymptotic value for solitons
 
    
    if IsDriven == 1
        drv = Circular_Driver(X, 0, v, DriverAmplitude, X_driver, K_driver, W_Driver);
        Fdrv = fft(drv);
    end
    fb = fft(b);
    fb(N/8+1:end-N/8) = 0;  %Zero out the middle 3/4 of the FFT components 

    %Create fbNew using split-step method given fb
    b = ifft(fb);
    FHilbert = -1i*sign(Kx).*fft(  abs(b).^2  );
    NLHilbert = ifft(  FHilbert  );
    nonlinearStuff = Kx.*fft((alpha*abs(b).^2 - v + gamma_Landau*NLHilbert ).*b);%agrees with earlier version for alpha = 0 and v = B0^2
    nonlinearStuff(N/8+1:end-N/8) = 0;
    
    piece = 1i * dt * 0.5 * Kx .*Kx * (dispersion + 1i*nu);
    if IsDriven == 1
        fbNew = ( fb - 1i * dt * nonlinearStuff) .* ((1+piece)./(1-piece))+dt*Fdrv./(1 - piece);
    else
        fbNew = ( fb - 1i * dt * nonlinearStuff) .* ((1+piece)./(1-piece));        
    end
%     Follows from the split-step method.

    fbOld = fb;
    fb = fbNew;
    
    
    X4 = linspace(-L,L,N/4);
    elapsed = num2str(cputime - t); t = cputime;
    stars = '*************************************\r\n';
    us = '_____________________________________\r\n';
    nl = '\r\n';
%% Timestepping   
    close all;
    debugging = 1;
    figureCoordinates = [(NMonitors - 1 + 0.1) 0.1 .8 .8];
    if NMonitors == 1 && debugging == 1
        figureCoordinates = [(NMonitors - 1 + 0.51) 0.1 .48 .8];
    end
    figureHandle = figure('Color', [1, 1, 1],'units','normalized','outerposition',figureCoordinates);
    % Be sure timestep, TimestepAttempted, linenum and TotalTime are
    % updated prior to timestepping
    timestep =  TimestepStart;
    linenum = LineNumberStart - 1;
    idxMid = round(length(b)/2*0.7);

     
        if  contains(Profile,'KBW')  ||  boolShockFit          

            KBW_tol = 0.001;
            ntest = 10;%(bb,x,bL, basymptote, ntest, tol)
            %[idMax, idSL, idSR, xShockL, xShockR, Iz] = IS_ShockBoundaries(bb,x,bL, bZasymptoteL, bZasymptoteR, ntest, tol)
            [idMax, IdxShockL, IdxShockR, XShockL, xShockR, IZ_KBW_init] = IS_ShockBoundaries(b(1:floor(length(b)/2)), X(1:floor(length(b)/2)), b_L, BZasymptoteL, BZasymptoteR, ntest, KBW_tol);
%             [val, idL] = min(abs(xkbw-0.1*L - X));
%             [val, idR] = min(abs(xkbw+0.1*L - X));
%             [val, idMax] = max(abs(imag(b(idL:idR))));
            Xold = X(idMax);
            Vold = v;
            cstheta = cos(theta_kbw);
            sgntheta = sign(cstheta);
            q = b_kbw/b_L*sgntheta;
%             Varray(1:10) = v;
%             Vweight(1:10) = 1;
%             vIdx = 1;
%             for i = 1:10
%                 Vweight(i)=sqrt(1.0/(1.0*i));
%             end
        end
    NoWorries = 1;
    TrialData = split(TrialID,'_');
    TrialTxt = sprintf('TrialID %s_{%s}',TrialData{1},TrialData{2});
%     timeSoFar = cputime;
    timeSoFar = rem(now,1);
%     timeSoFar = now;
    secondsPerday = 24*3600;
    FirstTimeNumber = -1;
    NFrames = round( (TimestepAttempted - timestep)/timestepsPerFrame - 0.5 );
    elapsedtimeTrac = zeros(1,timestepsPerFrame);
    timePerLineTrac = zeros(1, NFrames);
%     timingCorrection = 1.56;% 4.013 for N = 76,428 and timestepsPerFrame = 12,1312
                             % 1.56 for N = 2^17 and timestepsPerFrame =  12721
    while ishandle(figureHandle) && timestep <= (TimestepAttempted ) && NoWorries
        if timestep == (TimestepAttempted )
             "hello"
         end
        FrameIdx = mod(timestep,timestepsPerFrame) + 1; % FrameIdx = 1 if timestep = n*timestepsPerFrame for n = 0, 1, 2, ...
                                                        % ... so it counts from 1 to timestepsPerFrame
%         elapsedTime = (cputime - timeSoFar)/timingCorrection; timeSoFar = cputime;
        % now returns the date in the form xxxxx.yyyyyy where xxxxx is a
        % counting number giving the date and yyyyyy is the decimal
        % fraction of a day into the current
        if rem(now,1) < timeSoFar
            elapsedTimeSeconds = ( timeSoFar - rem(now,1) )*secondsPerday; 
        else            
            elapsedTimeSeconds = ( rem(now,1) - timeSoFar )*secondsPerday;
        end
        timeSoFar = rem(now,1);
%         timeSoFar = now;
        elapsedtimeTrac(FrameIdx) = elapsedTimeSeconds;
        if FrameIdx > 1
            temptimeremaining = sum(elapsedtimeTrac(2:FrameIdx))/(FrameIdx-1)*(timestepsPerFrame-FrameIdx);
        else
            temptimeremaining = 0;
        end
        if mod(timestep,floor(timestepsPerFrame/100)) == 1
            clc;
            if temptimeremaining > 60
                timemsg = fprintf('Est. time left to complete this line = %4.2f minutes.\n',temptimeremaining/60);
            else 
                timemsg = fprintf('Est. time left to complete this line = %4.2f seconds.\n',temptimeremaining);
            end
%             ['Est. time left to complete this line = ',num2str(temptimeremaining/60)]
        end

        if FrameIdx == 1
            FirstTimeNumber = FirstTimeNumber + 1;
            linenum = linenum + 1;
%             Bline = num2str(linenum);
            figure(figureHandle);
            b = ifft(fbOld);
            
           
            
            cla();
%             if contains(Profile,'KBW')
%                 bp = ifft(Kx.*fft(b));
%                 [val, idL] = min(abs(xkbw-0.1*L - X));
%                 [val, idR] = min(abs(xkbw+0.1*L - X));
%                 [val, idMax] = max(abs(bp(idL:idR)));
%                 plot([idMax idMax],[zmin, zmax],'Color',[1 0 0]);
%                 hold('on');
%             end
            plot(X,real(b), 'Color', [1, 0, 0]); % By is green
            hold('on');
            plot(X,imag(b), 'Color', [0,1,0]);  % Bz is red
            hold('on');
            plot(X,abs(b) , 'Color',[0 , 0, 1]);  % |B| is blue
            hold('on');
            if strcmp(Prfl('Profile'),'Helical')
%                 figureHandle;
                xL = Prfl('xL'); xR = Prfl('xR');
                EVL = Prfl('EVL'); EVR = Prfl('EVR');
                plot([xL xL],[-1.5 1.5],'Color',[1 0 0]);
                hold('on');
                plot([xR xR],[-1.5 1.5],'Color',[1 0 0]);
                hold('on');
                plot([EVL EVL],[-1.5 1.5],'Color',[0 1 0]);
                hold('on');
                plot([EVR EVR],[-1.5 1.5],'Color',[0 1 0]);
                hold('on');
            end
           
            if exact ==1
                plot(X,real(profile(Prfl,X,timestep*dt)));
                hold('on');
                plot(X,imag(profile(Prfl,X,timestep*dt)));
                hold('on');
                plot(X,abs(profile(Prfl,X,timestep*dt)));
                hold('on');
            end
            
            if  contains(Profile,'KBW') || boolShockFit
                
                %% Identify leading edge of profile
          
                [idMax, idxL, idxR, XShockL, XShockR, IZ_KBW] = IS_ShockBoundaries(b(1:floor(length(b)/2)), X(1:floor(length(b)/2)), b_L, BZasymptoteL, BZasymptoteR, ntest, KBW_tol);
                
                
%                 Izmid = sum(imag(b(idxL:idxMid)))*dx;

%                 [val, idL] = min(abs(xkbw-0.1*L - X));
%                 [val, idR] = min(abs(xkbw+0.1*L - X));
%                 [val, idRelMax] = max(abs(imag(b(idL:idR))));
% %                 idMax = idRelMax + idL;
%                 NinterpolateSteps = 40;
%                 NXsteps = 20;
%                 Xint = linspace(X(idMax - NXsteps),X(idMax + NXsteps),2*NinterpolateSteps);
%                 bb = interp1(X(idMax-NXsteps:idMax+NXsteps),b(idMax-NXsteps:idMax+NXsteps),Xint);
%                 [val,idIntMax] = max(abs(imag(bb)));
                Xnew = X(idMax);
                vtmp = Vold - (Xold - Xnew)/(timestepsPerFrame*dt);
                v = Vold*.75 + 0.25*vtmp;
                Xold = Xnew;
                Vold = v;
%                 v_scaled_by_alpha = v*(alpha);
                lambda1 = v/(alpha*Prfl('CR13')^2) - 1;
                r1 = sqrt(lambda1 + 0.25);
                q13 = 1/(-1/2 - r1);
                q23 = -1 - q13;
                
                  %Find LHS of Shock: assuming 'downstream' extends to N/2
%                 LHS = abs(abs(b-b_L)) <0.001;
%                 SATL = find(LHS(1:floor(N/2)) ==1);
%                 IdxShockL = SATL(end);
%                 XShockL = X(IdxShockL);

%               Find RHS of Shock: 
%               As the Intermediate Shock must carry the entirety of the
%               field rotation, then the downstream foot of the
%               Intermediate Shock is identified as the location where the
%               magnetic field direction approaches the full rotation of
%               the initial profile. The algorithm begins the search of the
%               magnetic field immediately to the right after the maximum
%               of Bz, which should be located within the shock.

%                 RHS = abs(angle(b(idMax:idMax+floor(0.2*N))) - theta_kbw)< 0.002*abs(theta_kbw);
%                 SATR = find(RHS == 1);
%                 IdxShockR = SATR(1)+idMax;
%                 XShockR = X(IdxShockR);

%                 IZ_KBW = sum(imag(b(IdxShockL:IdxShockR)))*dx;
                
                plot([XShockL XShockL],[0.8*zmin, 0.5*zmax],'Color',[0 0 0]);
                hold('on');
                plot([XShockR XShockR],[0.8*zmin, 0.5*zmax],'Color',[0 0 0]);
                hold('on');
                plot([XShockL XShockR],[0.8*zmin, 0.8*zmin],'Color',[0 0 0]);
                hold('on');
                
%                 plot([Xnew Xnew],[zmin, zmax],'Color',[1 0 0]);
%                 hold('on');
%                 plot([X(idxMid) X(idxMid)],[-1.5 1.5],'Color',[0 1 0]);
%                 hold('on');
%                 plot([X(idxL) X(idxL)],[-1.5 1.5],'Color',[0 1 0]);
%                 hold('on');
%                 if timestep == 0
% %                     Iz = sum(imag(b(idxL:idxR)))*dx;
% %                     KBW_tol = 0.01;
% %                     ntest = 10;%(bb,x,bL, basymptote, ntest, tol)
%                     %[idMax, idSL, idSR, xShockL, xShockR, Iz] = IS_ShockBoundaries(bb,x,bL, basymptote, ntest, tol)
%                     [idMax, idSL, idSR, XShockL, XShockR, Iz] = IS_ShockBoundaries(b(1:floor(length(b)/2)), X(1:floor(length(b)/2)), b_L, BZasymptoteL, BZasymptoteR, ntest, KBW_tol);
%                 end
%                 plot([X(idxL) X(idxL)],[-1.5 1.5],'Color',[0 0 1]);
%                 hold('on');
%                 plot([X(idxR) X(idxR)],[-1.5 1.5],'Color',[0 0 1]);
%                 hold('on');
                
                pIZmid = text(-L*0.9,0.65*zmax,sprintf('Iz initial # = %4.4f',IZ_KBW_init));
                pIZmid.Color = [0 1 0];
                pIZIS = text(-L*0.9,0.55*zmax,sprintf('Iz # = %4.4f',IZ_KBW));
                pIZIS.Color = [0 0 1];
                pspeed = text(-L*0.9,0.45*zmax,sprintf('IS shock speed = %4.4f',v));
                pspeed.Color = [0 0 0];
                pthetakbw = text(-L*0.9,0.35*zmax,sprintf('\\theta_{kbw} = %4.4f rad',theta_kbw));
                pthetakbw.Color = [0 0 0];
                pCR = text(-L*0.9,0.25*zmax,sprintf('1-3 IS compression ratio = %4.4f',q13));
                pCR.Color = [0 0 0];
                pLambda = text(-L*0.9,0.15*zmax,sprintf('2-3 IS compression ratio = %4.4f',q23));
                pLambda.Color = [0 0 0];
                text(X(idMax)-.1*L,0.85*zmin,sprintf('shock front: X = %4.2f',Xnew));
                
            end
            
            text(-L*0.9,0.95*zmax,sprintf('profile = %s',Profile));
            text(-L*0.9,0.85*zmax,sprintf('time = %4.3f/%4.3f',timestep*dt,TotalTime));
            text(-L*0.9,0.75*zmax,sprintf('line # = %4.0f/%4.0f',linenum, NumberLinesSaved));
            %(-L*0.9,0.65*zmax) reserved for profile specific parameters
            %(-L*0.9,0.55*zmax) reserved for profile specific parameters
            
            text(-L*0.4,0.95*zmax,sprintf('Trial ID = %s',TrialTxt));  
            text(-L*0.4,0.85*zmax,sprintf('alpha = %4.4f',alpha));
            text(-L*0.4,0.75*zmax,sprintf('dissipation = %4.4f',nu));
            text(-L*0.4,0.65*zmax,sprintf('gamma\\_Landau = %4.4f',gamma_Landau));
            text(-L*0.4,0.55*zmax,sprintf('dispersion = %4.4f',dispersion));
            
            
            text(L*0.1,0.95*zmax,sprintf('N # = %4.0f',N));            
            text(L*0.1,0.85*zmax,sprintf('dx # = %4.4f',dx));
            
%             if FirstTimeNumber == 1
%                 elapsedTime = cputime - timeSoFar; timeSoFar = cputime;
%                 timePerLineTrac(FirstTimeNumber) = elapsedTime;
%                 TimePerLine = elapsedTime;
%             end
            if FirstTimeNumber > 0
%                 elapsedTime = cputime - timeSoFar; timeSoFar = cputime;
%                 FTNidx = mod(FirstTimeNumber - 1,10) + 1;
%                 elapsedtimeTrac(FTNidx) = elapsedTime;
%                 if FirstTimeNumber > 10
%                     TimePerLine = mean(elapsedtimeTrac);
%                 else
%                     TimePerLine = elapsedTime;
%                 end
                timesPerstep = nonzeros(elapsedtimeTrac);
%                 TimePerLine = sum(timesPerstep);
%                 TimePLErr = std(timesPerstep)*sqrt(length(timesPerstep));
                TimePerLine = sum(timesPerstep);
                timePerLineTrac(FirstTimeNumber) = TimePerLine;
                timePerLineTracTmp=nonzeros(timePerLineTrac);
                
                if FirstTimeNumber > 5
                    TimePerLine = mean(timePerLineTracTmp);
                    TimePLErr = std(timePerLineTracTmp)/length(timePerLineTracTmp);
                else
%                     TimePerLine = sum(timesPerstep);
                    TimePLErr = std(timesPerstep)*sqrt(length(timesPerstep));
                end
                
                text(L*0.4,0.95*zmax,sprintf('Time per Line = (%4.3f +/- %4.3f) s',TimePerLine,TimePLErr));
                linesRemaining = NumberLinesSaved - linenum;
                text(L*0.4,0.85*zmax,sprintf('Lines Remaining = %4.0f', linesRemaining ));
                timeRemaining = linesRemaining*TimePerLine;
                hour = 3600;
                day = 24*hour;
                minnet = 60;
                if timeRemaining < minnet
                    text(L*0.4,0.75*zmax,sprintf('Estimated time remaining = (%4.1f +/- %4.2f) s', linesRemaining*TimePerLine, linesRemaining*TimePLErr ));
                elseif timeRemaining < hour
                    text(L*0.4,0.75*zmax,sprintf('Estimated time remaining = %4.0f s = (%4.2f +/- %4.2f)min ', linesRemaining*TimePerLine, linesRemaining*TimePerLine/60,linesRemaining*TimePLErr/60));
                elseif timeRemaining < day
                    text(L*0.4,0.75*zmax,sprintf('Estimated time remaining = %4.1f min = (%4.2f +/- %4.2f)hr',linesRemaining*TimePerLine/60, linesRemaining*TimePerLine/60/60, linesRemaining*TimePLErr/60/60));
                elseif timeRemaining > day
                    text(L*0.4,0.75*zmax,sprintf('Estimated time remaining = %4.1f hr = (%4.3f +/- %4.3f)days', linesRemaining*TimePerLine/60/60, linesRemaining*TimePerLine/60/60/24, linesRemaining*TimePLErr/60/60/24));
                    text(L*0.4,0.65*zmax,sprintf('Really ought to reconsider this ... '));
                    text(L*0.4,0.55*zmax,sprintf('or pursue a hobby. Maybe start learning to play the piano?'));
                end
            end

            if IsDriven
                plot([-X_driver -X_driver],[-0.8*zmax 0.5*zmax],'Color',[0 0 0]);
                plot([X_driver X_driver],[-0.8*zmax 0.5*zmax],'Color',[0 0 0]);
                p1 = [-X_driver -0.78*zmax]; p2 = [X_driver -0.78*zmax]; dp = p2-p1;
                quiver(p1(1),p1(2),dp(1),dp(2),0,'color', [0 0 0],'MaxHeadSize',1/(2*X_driver));
                quiver(p2(1),p2(2),-dp(1),-dp(2),0,'color', [0 0 0],'MaxHeadSize',1/(2*X_driver));
                text((p1(1) + p2(1))/2,p1(2) - zmax*0.05,sprintf('Driver Range'),'HorizontalAlignment','center');
                
                text(L*0.1,0.75*zmax,sprintf('Driver: Amp.  = %4.4f',DriverAmplitude),'HorizontalAlignment','left');
                text(L*0.1,0.68*zmax,sprintf('            K     = %4.4f',K_driver),'HorizontalAlignment','left');
                text(L*0.1,0.61*zmax,sprintf('            range = +/- %4.4f',X_driver),'HorizontalAlignment','left');
                text(L*0.1,0.54*zmax,sprintf('            Width = %4.4f',W_Driver),'HorizontalAlignment','left');
            else
                text(L*0.1,0.75*zmax,sprintf('No Driver'));
            end

            
            
            
            if strcmp(Profile, '1psoliton')
%                 text(-L*0.8,0.8*zmax,sprintf('profile = %s',Profile));
                text(-L*0.9,0.65*zmax,sprintf('lambda = %4.4f and sign = %4.0f',lambda,onePsign));
            end
            if strcmp(Profile, '2psoliton')
%                 text(-L*0.8,0.8*zmax,sprintf('profile = %s',Profile));
                text(-L*0.9,0.65*zmax,sprintf('lambda = %4.4f + i(%4.4f)',lr,li));
            end
            
           % text(-L*0.8,0.7,sprintf('time = %4.4f',timestep*dt));
%             text(-L*0.8,0.7,sprintf('step # = %4.4f',timestep));
            axis([-L L zmin zmax]);
            drawnow();
            
            
            Stemp1 = [nl,'LineStart',num2str(linenum),nl];
            Bfile = fopen(BFfile,'a');
            fprintf(Bfile,Stemp1);
            fclose(Bfile);
%             b4 = complex(real(b(1:4:length(b))), imag(b(1:4:length(b))));
            b4 = b(1:4:length(b)) + 1i*10^(-25);
            dlmwrite(BFfile,b4,'precision','%.6f','-append','delimiter',' ');
            if CreateMatfile == 1
                BFMat.(['Run',RunNumber])(linenum,1:length(b4)) = b4;
                clear BFMat;
            end
%             save BFMat;             
        end
%          if timestep == (TimestepAttempted )
%              "hello"
%          end

        timestep = timestep + 1    ;

            b = ifft(fb);
            FHilbert = -1i*sign(Kx).*fft(  abs(b).^2  );
            NLHilbert = ifft(  FHilbert  );
            nonlinearStuff = Kx.*fft(  (alpha*abs(b).^2 + gamma_Landau*NLHilbert ).*b  );
            
            nonlinearStuffLF = - fft((alpha*abs(b).^2  + gamma_Landau*NLHilbert ).*b);
            linearstuffLF =(a0*v + Kx.*(dispersion +1i*nu).*c0).*fb;
%             nonlinearStuff = Kx.*fft((abs(b).^2 - v/alpha).*b);

            nonlinearStuff(N/8+1:end-N/8) = 0;
            nonlinearStuffLF(N/8+1:end-N/8) = 0;
            linearstuffLF(N/8+1:end-N/8) = 0;
            
            piece = 1i * dt * Kx .* ((dispersion + 1i*nu)*Kx); % i*Dt*k^2 for v = 0 and nu = 0
            pieceAB = 1i*Kx.*(v*(1 - a0) + Kx.*(dispersion + 1i*nu*(1 - c0))  ).*dt;
            
            if IsDriven == 1
                if gamma_Landau ~= 0
                    fbNew = fbOld .* ((1+piece)./(1-piece)) - nonlinearStuff .* ((2i*dt)./(1-piece)) + dt*Fdrv./(1 - piece);
                else
                    fbNew = fbOld .* ((1+piece)./(1-piece)) - nonlinearStuff .* ((2i*dt)./(1-piece)) + dt*Fdrv./(1 - piece);
                end
                
            else
                if gamma_Landau ~= 0
%                     fbNew = fbOld .* ((1+piece)./(1-piece)) - nonlinearStuff .* ((2i*dt)./(1-piece)) ;
                    fbNew = fbOld.* ((1+pieceAB)./(1-pieceAB)) + 2i*Kx.*dt.*(  linearstuffLF + nonlinearStuffLF  )./(1-pieceAB);
                else
                    fbNew = fbOld .* ((1+piece)./(1-piece)) - nonlinearStuff .* ((2i*dt)./(1-piece)) ;
                end
            end
            
            
            fbOld = fb;
            fb = fbNew;
            
            nanidx = isnan(b);
            if any(nanidx(:) ==1)
                thisIsNotGood = 1;
                growth = ['B-field evaluated to NAN at step = ',num2str(timestep)];
                growth
                NoWorries = 0;
            end
%             timestep
%             ishandle(figureHandle) && timestep <= (TimestepAttempted - 1) && NoWorries
%             if timestep > 105
%                 msg = 'let''s see what''s happenging.'
%             end
    end
    
    if  contains(Profile,'KBW') && boolShockFit
        Prfltitle = 'Profile Parameters';
        Iz = IZ_KBW_init;
        Wkbw=Prfl('Wk');WR=Prfl('Wr');theta_L=Prfl('theta_L');theta_kbw=Prfl('theta_kbw');theta_R=Prfl('theta_R');b_L=Prfl('b_L');b_kbw=Prfl('b_kbw');b_R=Prfl('b_R');xR=Prfl('Xr');
        PrflKeys   = {'title',    'Profile', 'speed'                , 'CR13', 'CR23', 'Wk'  , 'Wr' , 'theta_L', 'theta_kbw', 'Iz', 'theta_R' , 'b_L', 'b_kbw', 'b_R','Xkbw','Xr', 'boolKbwShockFit'};
        PrflValues = { Prfltitle,  Profile ,  v     ,  q13     q23  ,  Wkbw ,  WR  ,  theta_L ,  theta_kbw ,  Iz ,  theta_R  ,  b_L ,  b_kbw ,  b_R , xkbw , xR ,  boolShockFit };
        Prfl       = containers.Map(PrflKeys,PrflValues);
    end
    
    if  contains(Profile,'Generated13IS') && boolShockFit
        Prfltitle = 'Profile Parameters';
        Iz = IZ_KBW_init;
        CR13=Prfl('CR13');WR=Prfl('Wr');theta_L=Prfl('theta_L');theta_IS=Prfl('theta_IS');theta_R=Prfl('theta_R');b_L=Prfl('b_L');b_R=Prfl('b_R');xR=Prfl('Xr');
        PrflKeys   = {'title',    'Profile', 'speed', 'CR13' , 'CR23'       , 'Wr' , 'theta_L', 'theta_IS', 'xIS', 'Iz'         , 'theta_R' , 'b_L', 'b_R','Xr', 'boolShockFit'};
        PrflValues = { Prfltitle,  Profile ,  v     ,  CR13  , (-1 - CR13)  ,  WR  ,  theta_L ,  theta_IS ,  xIS ,  IZ_KBW_init ,  theta_R  ,  b_L ,  b_R , xR ,  boolShockFit };
        Prfl       = containers.Map(PrflKeys,PrflValues);
    end
    
    if isKey(Prfl,'ProfileNew')
    end
    
    PSTR = struct('TCon',Trial,'NCon',Numerics,'PCon',Physics,...
                      'PrflCon', Prfl, 'GCon',Graphics,'FCon',Folder, 'DrvrCon',Driver); 
    
%      if ishghandle(figureHandle)
%        close(figureHandle);
%      end
end