function xgeninit(file)

global xgeninfo xgenfile xgenstat

xgenfile=file;
xgeninfo=h5info(file);

z=xgenreaddataset('/Lattice/z'); % read z and zplot
sref=xgenreaddataset('/Global/lambdaref');
ds=xgenreaddataset('/Global/sample');

xgenstat.z=z{1};
xgenstat.zplot=z{2};
xgenstat.sref=sref{1};
xgenstat.ds=ds{1}*sref{1};
