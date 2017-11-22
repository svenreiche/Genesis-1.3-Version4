function xgeninit(file)

global xgeninfo xgenfile xgenstat

xgenfile=file;
xgeninfo=h5info(file);

z=xgenreaddataset('/Lattice/z');
dz=xgenreaddataset('/Lattice/dz');
nz=xgenreaddataset('/Global/nzout');
sref=xgenreaddataset('/Global/lambdaref');
ds=xgenreaddataset('/Global/sample');

xgenstat.z=z{1};
xgenstat.dz=dz{1};
xgenstat.nz=nz{1};
xgenstat.sref=sref{1};
xgenstat.ds=ds{1}*sref{1};
