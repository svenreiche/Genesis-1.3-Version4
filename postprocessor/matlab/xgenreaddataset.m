function [dat label]=xgenreaddataset(dset)


global xgeninfo xgenfile
dat={};
label={};

for i=1:length(xgeninfo.Groups)
    for j=1:length(xgeninfo.Groups(i).Datasets)
        field=sprintf('%s/%s',xgeninfo.Groups(i).Name,xgeninfo.Groups(i).Datasets(j).Name);
        if (length(regexpi(field,dset))>0)
            dd=h5read(xgenfile,field);
            if ~isempty(strfind(xgeninfo.Groups(i).Name,'Lattice'))
                dd=dd';
            end
            ndd=size(dd);
            fprintf('Dataset found: %s [%dx%d]\n',field,ndd(1),ndd(2))
            dd(isnan(dd))=0;
            dat{end+1}=dd;
            label{end+1}=field;
        end
    end
end

end

