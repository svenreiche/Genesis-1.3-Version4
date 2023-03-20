function ret = xgenplot(field,mode,ref)
% XGENPLOT generate a plot with the information specified
% xgenplot(field,mode,ref)
%
% returns in an cell array the data which has been plotted
%
% field  - regular expression of a dataset name in the hdf output file.
%          It will plot all datasets which matches the string
% mode   - various methods to plot the data
%           'normal' - along undulator, requires s-position as 'ref'
%           'profile'- along beam frame, requires z-position as 'ref'
%           'mean'   - along undulator, using the mean value at each step
%           'max'    - along undulator, using the max value at each step
%           'rms'    - along undulator, RMS value only useful for power and spectrum
%           'weighted' - along the undulator, usign a weighted mean calculation
%                        using power profile for field properties or
%                        current for beam properties
%           '2d'     - Image along beam frame in x and along undulator in y
%           '2dnorm' - same as '2d' but normalized at each step along
%                      undulator
%
    if nargin<3
        ref=0;
    end    
    if nargin<2
        mode='normal';
    end
    
    if isempty(strfind(field,'spectrum'))
        [dat, lab]=xgenreaddataset(field);
    else
        field1=strrep(field,'spectrum','intensity');
        [dat, lab]=xgenreaddataset(field1);
    end
    
    if isempty(dat)
        fprintf('XGENESIS - ERROR: Invalid Field - no dataset found\n');
        return
    end
    
    if ~isempty(strfind(field,'spectrum'))
        for i=1:length(dat)
            field2=strrep(lab{i},'intensity','phase');
            datphase=xgenreaddataset(field2);
            sig=sqrt(dat{i}).*exp(1i*datphase{1});
            spec=abs(fftshift(fft(sig),1));
            dat{i}=spec;
            lab{i}=strrep(lab{i},'intensity','spectrum');
        end
    end
    
    if ~isempty(strfind(mode,'rms'))
        for i=1:length(dat)
            dat{i}=xgenrms(dat{i},lab{i});
            lab{i}='rms length';
        end
    end
    
    if ~isempty(strfind(mode,'weight'))
        for i=1:length(dat)
            dat{i}=xgenweight(dat{i},lab{i});
            lab{i}=sprintf('%s - weighted',lab{i});
        end
    end
    
    ret = {};
    for i=1:length(dat)
        [x,y]=xgensingleplot(dat{i},lab{i},mode,ref);
        ret{end+1}={x,y};
        
        hold on
    end
 
    

    legend(lab)
    hold off
end

    
function [x]=xgenweight(dat,lab)

global xgenstat

    dims=size(dat);
    ns=dims(1);
    nz=dims(2);

    s=1:ns;
    x=1:nz;
    
    % getting the weighting distribution
    if ~isempty(strfind(lab,'Beam'))
        wdat=xgenreaddataset('/Beam/current');
        w=wdat{1};
        for i=1:nz
            prof=dat(:,i);
            x(i)=sum(prof'.*w')/sum(w');           
        end
    else
        idx=strfind(lab,'/');
        llab=sprintf('%spower',lab(1:idx(2)));
        wdat=xgenreaddataset(llab);
        w=wdat{1};
        for i=1:nz
            prof=dat(:,i);
            ww=w(:,i);
            x(i)=sum(prof'.*ww')/sum(ww');           
        end

    end
    

    
    return
end

function [x]=xgenrms(dat,lab)

global xgenstat

    dims=size(dat);
    ns=dims(1);
    nz=dims(2);

    s=1:ns;
    x=1:nz;
    if isempty(strfind(lab,'spectrum'))
        for i=1:nz
            prof=dat(:,i);
            xmean=sum(prof'.*s)/sum(prof);
            xsig=sum(prof'.*s.*s)/sum(prof);
            x(i)=sqrt(xsig-xmean*xmean);

        end 
        x=x*xgenstat.ds/299792458.0;
    else
        
        for i=1:nz
            prof=dat(:,i);
            c=max(prof)*0.01;
            idx=prof>c;
            xmean=sum(prof(idx)'.*s(idx))/sum(prof(idx));
            xsig=sum(prof(idx)'.*s(idx).*s(idx))/sum(prof(idx));
            x(i)=sqrt(xsig-xmean*xmean);
        end
        df=1/xgenstat.ds/ns;
        dE=1240e-9*df;
        x=x*dE;
    end        
         
    
    return
end
    

function [x,y]=xgensingleplot(dat,label,mode,ref)

    global xgenstat

    if ~isempty(strfind(label,'Lattice'))
        x=xgenstat.z;
        y=dat';
        plot(x,y);
        xlabel('z (m)');
        return;
    end

    isspec=0;
    if (~isempty(strfind(label,'spectrum')))
        isspec=1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % process data
    
    if strcmp(mode,'mean')
        dat=mean(dat);
    end

    if strcmp(mode,'max')
        dat=max(dat);
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prepare the plot
    
    x=[];
    y=[];

    dims=size(dat);
    ns=dims(1);
    nz=dims(2);
    z=xgenstat.zplot;    
    s=((1:ns)-1)*xgenstat.ds;
    zlab='z (m)';
    slab='s (m)';
    
    if isspec~=0
        f0=1/xgenstat.sref;
        df=1/xgenstat.ds;
        E0=1240e-9*f0;
        dE=1240e-9*df;
        s=(1:ns)*dE/ns-0.5*dE+E0;
        slab='E_{ph} (eV)';
    end
    

    if (nz<2) && (ns<2)
        fprintf('ERROR - Dataset is a single value\n');
        return;
    end
        
    if (ns<2)
      x=z;
      y=dat';
      plot(x,y)
      xlabel(zlab);
      return;
    end
    
    if (nz<2)
        x=s;
        y=dat;
        plot(x,y);
        xlabel(slab);
        return
    end

    if (strcmp(mode,'2d')||strcmp(mode,'2dnorm')) 
        y=dat;
        if (strcmp(mode,'2dnorm'))
          for i=1:nz
            y(:,i)=y(:,i)/mean(dat(:,i));      
          end
        end  
        x=s;
        %imagesc(s,z,y');
        %meshz(s,z,y');
        surface(s,z,y','EdgeColor','none')
        ax=gca;
        ax.YDir='normal';
        return
    end
    
    if (strcmp(mode,'profile'))
        dz=abs(z-ref);
        [zmin,idx]=min(dz);
        x=s;
        y=dat(:,idx);
        fprintf('Output for closest data point at z = %f m\n',z(idx));
        plot(x,y)
        xlabel(slab);
        return;
    end
    
    if (strcmp(mode,'normal'))
        ds=abs(s-ref);
        [zmin,idx]=min(ds);
        x=z;
        y=dat(idx,:)';
        fprintf('Output for closest data point at s = %e micron\n',s(idx));
        plot(x,y)
        xlabel(zlab);
       return;
    end
    
    fprintf('ERROR - Unknown plotting mode: %s\n',mode);
    return;
end
