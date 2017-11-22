function success = xgenlattice(lambda,lambdau,Elim,Klim,Nwig,Nsec,Fodo,Taper,Delay,Root,Line)
% XGENLATTICE  Generate a lattice file for Genesis  
% success = xgenlattice(lambda,Elim,Klim,Nwig,Nsec,Fodo,Taper,Root,Line)
%
% lambda  : resonant wavelength in m
% lambdau : undulator period in m
% Elim    : [Emin, Emax] - Array for energy range in GeV
% Klim    : [Kmin, Kmax] - Array for energy range in GeV
% Nwig    : Number of undulator periods
% Nsec    : Number of section (undulator modules)
% Fodo    : [k1,Lq,Ld,Lfodo] = Array for focusing strength, quad length,
%           k1    - Focusing strength in 1/m
%           Lq    - Quadrupole length in m
%           Ld    - Drift between module and quadrupole m
%           Lfodo - Fodo cell length (2 modules) in m
% Taper   : [z,a,b,c] = Array for taper
%           z - taper start in m
%           a - general linear taper with K=K0*(1+a*z)
%           b,c - K=K0*(1+b*(z-z0)^c) for z>z0
% Root    : root file name, adding the extention ".lat" for file name
% Line    : name of final line to be generated
%
% returns the resonant energy or -1 for case of no solution


success = -1;


Emin = Elim(1);
Emax = Elim(2);
Kmin = Klim(1);
Kmax = Klim(2);


% Step one - calculating resonance condition. Higher energy is always
% preferred

gam=Emax/0.511e-3;   % Lorentz factor

Ktmp=(lambda*2*gam*gam/lambdau-1)*2;  % this is K^2 and should be a positive number


if Ktmp < 0
    K = Kmin;
    gam=sqrt(0.5*lambdau/lambda*(1+0.5*K*K));
else
    K=sqrt(Ktmp);
    if K>Kmax
        K=Kmax;
        gam=sqrt(0.5*lambdau/lambda*(1+K*K*0.5));
    end
end

Etmp=gam*0.511e-3;
if ((Etmp < Emin) || (Etmp > Emax))
       fprintf('No resonance condition for lambda = %e with given limits in energy and K\n',lambda);
       return;
end

fprintf('Wavelength (nm): %f\n', lambda*1e9);
fprintf('K-Parameter    : %f\n', K);
fprintf('Energy (GeV)   : %f\n', gam*0.511*1e-3);


% step two - generating the file


Lund=Nwig*lambdau;
qgrad=Fodo(1);
Lq=Fodo(2);
Ld1=Fodo(3);
Ldrift=0.5*Fodo(4)-Lund;

if (Lq+Ld1) > Ldrift
      fprintf('Quadrupole and Drift not fitting into FODO lattice\n');
      return;
end
Ld2 = Ldrift-Lq-Ld1;

fid = fopen(sprintf('%s.lat',Root),'w');

fprintf(fid,'D1: DRIFT = { l = %f };\n',Ld1);
fprintf(fid,'D2: DRIFT = { l = %f };\n',Ld2);
fprintf(fid,'QF: QUADRUPOLE = { l = %f, k1= %f };\n',Lq,qgrad);
fprintf(fid,'QD: QUADRUPOLE = { l = %f, k1= %f };\n',Lq,-qgrad);


nchic=length(Delay);

Lb=0.25*0.5*Ld1;
Lc=0.5*0.5*Ld1;

for i=1:nchic
    fprintf(fid,'C%2.2d: CHICANE = { l = %f,lb = %f, ld = %f, delay = %e };\n',i,Ld1,Lb,Lc,Delay(i));
end

 ztap=Taper(1);
 atap=Taper(2);
 btap=Taper(3);
 ctap=Taper(4);

 z=0;
 for i = 1:Nsec
     if z>ztap
        Kloc=K*(1+atap*z+btap*(z-ztap)^ctap);     
     else
        Kloc=K*(1+atap*z);     
     end
     z=z+0.5*Fodo(4);
     fprintf(fid,'U%2.2d: UNDULATOR = { lambdau=%f,nwig=%d,aw=%f};\n',i,lambdau,Nwig,Kloc/sqrt(2));
end
 
 fprintf(fid,'\n%s: LINE={',Line);
 
 for i=1:Nsec
    fprintf(fid,'U%2.2d,',i);
    if i<nchic
        fprintf(fid,'C%2.2d,',i);
    else
        fprintf(fid,'D1,');
    end
    if (mod(i,2)==1)
        fprintf(fid,'QF');
    else
        fprintf(fid,'QD');
    end
    if (Ld2>0)
        fprintf(fid,',D2');
    end
    if (i<Nsec)
        fprintf(fid,',\n\t');
    else
        fprintf(fid,'};\n');
    end
 end

fclose(fid);

success=gam;
return











end