function xgeninput(lambda, gamma, lambdau, sample, Beam , Field, Root, Line, zmatch, seed)
% XGENINPUT  Generate an inputfile for Genesis  
% success = xgenlattice(lambda,gamma,lambdau,sample,Beam,Field,Root,Line,zmatch,seed)
%
% lambda  : resonant wavelength in m
% gamma   : resonant energy
% lambdau : undulator period in m = basic integration step
% sample  : sample rate and factor for integration step
% Beam    : [current, emit, espread, blen]
%            current  - current in A
%            emit     - emittance in m
%            espread  - RMS spread of Lorent factor
%            blen     - bunch length. If zero or negative, steady-state
%            bpad     - padding of zero current before the bunch
% Field   : [power, dgrid, ngrid, w0]
%            power     - radiation power in W
%            dgrid     - grid size from center to edge
%            ngrid     - number of grid points per dimension
%            w0        - initial mode size  in m
% Root    : root file name. Input file has extension '.in'
% Line    : lattice file = Root+'.lat'
% zmatch  : Fodo cell length for matching
% seed    : seed base for shot noise


% nothing to return



fid=fopen(strcat(Root,'.in'),'w');
fprintf(fid,'&setup\n');
fprintf(fid,'rootname=%s\n',Root);
fprintf(fid,'lattice=%s.lat\n',Root);
fprintf(fid,'beamline=%s\n',Line);
fprintf(fid,'lambda0=%e\n',lambda);
fprintf(fid,'gamma0=%f\n',gamma);
fprintf(fid,'delz=%f\n',lambdau*sample);
fprintf(fid,'shotnoise=1\n');
fprintf(fid,'seed=%d\n',seed);
fprintf(fid,'&end\n\n');
        
fprintf(fid,'&lattice\nzmatch=%f\n&end\n\n',zmatch);

tlen=Beam(4);
tpad=Beam(5);
cur =Beam(1);

if tlen > 0
    fprintf(fid,'&profile_step\nlabel=current\nc0=%f\ns_start=0\ns_end=%e\n&end\n\n',cur,tlen);
    fprintf(fid,'&time\nslen=%e\nsample=%f\n&end\n\n',tlen+tpad,sample);
end    

power=Field(1);
dg = Field(2);
ng = floor(Field(3));
w0= Field(4);
fprintf(fid,'&field\npower=%e\ndgrid=%e\nngrid=%d\nwaist_size=%e\n&end\n\n',power,dg,ng,w0);
    


dgam=Beam(3);
emit=Beam(2);
fprintf(fid,'&beam\n');
if (tlen>0)
    fprintf(fid,'current=@current\n');
else
    fprintf(fid,'current=%f\n',cur);
end

fprintf(fid,'delgam=%f\nex=%e\ney=%e\n&end\n\n',dgam,emit,emit);
fprintf(fid,'&track\n&end\n');
fclose(fid);

 
   
end














