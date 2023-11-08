function DNLS_Bfield_Footer(file, Diagnostics, Folder, Graphics, Numerics, Physics, Trial, k, time)

    elapsed = num2str(cputime - time); t = cputime;
 
    stars = '*************************************\r\n';
    us = '_____________________________________\r\n';
    nl = '\r\n';

%     TrialLine = [stars,'Trial Identifiers:\r\n',us,...
%         'Time Stamp =      ',num2str(Trial('TimeStamp')),nl,'Trial Name =      ',num2str(Trial('TrialName')),nl,'Trial Objective = ',num2str(Trial('TrialObjective')),nl,nl];
% 
% 
%     PhysicsLine = [stars,'Physics Parameters:\r\n',us,...
%         'bo =     ',num2str(Physics('bo')),nl,'lambda = ',num2str(Physics('lambda')),nl,'nu =     ',num2str(Physics('nu')),nl,'sign =   ',num2str(Physics('sign')),nl,nl];
% 
%     NumericsLine = [stars,'Numeric Parameters:\r\n',us,...
%         'Courant = ',num2str(Numerics('Courant')),nl,'L = ',num2str(Numerics('L')),nl,'N = ',num2str(Numerics('N')),nl,'dt = ',num2str(Numerics('dt')),nl,'dx = ',num2str(Numerics('dx')),nl,nl];

    DiagnosticsLine = [nl,stars,'Diagnostic Parameters:\r\n',us,...
        'CPU Time = ',num2str(Diagnostics('CpuTimeUsed')),nl,'Time Steps Achieved = ',num2str(Diagnostics('timesteps')),nl,'Time Steps Attempted = ',num2str(Diagnostics('TimestepAttempted')),nl,'growth = ',num2str(Diagnostics('growth')),nl,'decay = ',num2str(Diagnostics('decay')),nl,nl];

    Bfile = fopen(file,'a');
%     fprintf(Bfile,TrialLine);
%     fprintf(Bfile,PhysicsLine);
%     fprintf(Bfile,NumericsLine);
    fprintf(Bfile,DiagnosticsLine);
    fprintf(Bfile,[nl,'eof']);
    fclose(Bfile);
    elapsed = num2str(cputime - t); t = cputime;
end