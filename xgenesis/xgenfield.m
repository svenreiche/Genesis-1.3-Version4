function xgenfield(file)

ns=h5read(file,'/slicecount');
dg=h5read(file,'/gridsize');
lam=h5read(file,'/wavelength');
ds=h5read(file,'/slicespacing');

ngrid=151;
power=(1:ns)*0.1;
inten=(1:(ngrid*ngrid))*0;

i=2350;
field1=sprintf('/slice%6.6d/field-real',i);
field2=sprintf('/slice%6.6d/field-imag',i);
fre=h5read(file,field1)*1.0;
fim=h5read(file,field2)*1.0;




for i=1:ns
    field1=sprintf('/slice%6.6d/field-real',i);
    field2=sprintf('/slice%6.6d/field-imag',i);
    fre=h5read(file,field1)*1.0;
    fim=h5read(file,field2)*1.0;
    sum1=sum(fre.^2);
    sum2=sum(fim.^2);
%    intloc=fre.^2+fim.^2;
%    inten=inten+intloc;
    power(i)=sum1+sum2;
end

figure(1)
plot(power)