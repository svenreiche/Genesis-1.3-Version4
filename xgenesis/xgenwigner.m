function [ output_args ] = xgenwigner( ref )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global xgenstat

[dat, lab]=xgenreaddataset('intensity-nearfield');

dat1=dat{1};
dims=size(dat1);
ns=dims(1);
nz=dims(2);

if (ns<2)
    fprintf('No time-dependent data set available\n');
    return;
end


[dat, lab]=xgenreaddataset('phase-nearfield');
dat2=dat{1};


z=xgenstat.zplot;
       
dz=abs(z-ref);
[zmin,idx]=min(dz);
fprintf('Output for closest data point at z = %f m\n',z(idx));


signal=sqrt(dat1(:,idx)).*exp(1i*dat2(:,idx));


n0=ns;
s=((1:ns)-1)*xgenstat.ds;
f0=1/xgenstat.sref;
df=0.5/xgenstat.ds;
E0=1240e-9*f0;
dE=1240e-9*df;
f=((1:ns)-1)*dE/ns-0.5*dE+E0;




sig0=(1:(3*n0))*0;
for i=1:n0
    sig0(i+n0)=signal(i);
end
    
wig=zeros(n0,n0);

up=round((0:(n0-1))/2);
down=round((1:n0)/2);
for it=1:n0
    t=it+n0;   % index
    sig1=sig0(t+up); 
    sig2=sig0(t-down);
    sig=sig1.*conj(sig2);
    spec=abs(fftshift(fft(sig)));
    wig(:,it)=spec;

end

imagesc(s,f,wig);

%surface(f,s,wig,'EdgeColor','none')
%ax=gca;
%ax.YDir='normal';
