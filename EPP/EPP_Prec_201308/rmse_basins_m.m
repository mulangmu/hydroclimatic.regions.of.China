%**************************************************************
%The cor for provinces
%**************************************************************
clear 
%figure;
%provinces={'HAINAN','HEILONGJIANG','NEIMENGGU','XINJIANG','JILIN','LIAONING','HEBEI','GANSU','BEIJING','SHANXI',...
%           'TIANJIN','QINGHAI','SHANNXI','NINGXIA','SHANDONG','XIZHANG','HENAN','JIANGSHU','ANHUI','SHICHUANG','HUBEI',...
%           'SHANGHAI2','SHANGHAI','ZHEJIANG','JIANGXI','HUNAN','YUNNAN','GUIZHOU','FUJIAN','TAIWAN','GUANGXI','GUANGDONG'}
basins={'1','2','3','4','5','6','7','8','9','10',...
           '11','12','13','14','15','16','17'} 

files=dir('F:\Thesis\CFS version 2\201308\文章修改\EPP_Prec_201308\basins17_0\001\*.txt');
%print(files.name)
for n = 1:length(files)
    
    figure(n)
    set(gcf,'unit','normalized','position',[0.1,0.2,0.35,0.35]);
    f = ['F:\Thesis\CFS version 2\201308\文章修改\EPP_Prec_201308\basins17_0\001\',files(n).name];
    fid = fopen(f, 'r');  
    k=0;
    s1 = fscanf(fid,'%s',1); 
    while ~feof(fid)
        k=k+1;
        xy =  fscanf(fid,'%d',2); 
        ncol=xy(1);
        ncow=xy(2);
        x1=zeros(ncol,ncow);
        for i=1:ncol   
            for j=1:ncow  
                x1(i,j) = fscanf(fid,'%f',1); 
            end
        end
        if (strcmp(s1,'cor'))
            x11=x1;
            Title=basins(n);
            s1=Title{1}
            contourf(1:ncow,1:ncol, x1,'DisplayName',s1,'LineStyle','none');
            caxis([0,1]);  %if (kk==11)  colorbar; end
            title(s1,'FontSize',26,'FontWeight','bold');
            xlabel('Lead time','FontSize',20,'FontWeight','bold');
            ylabel('Date','FontSize',20,'FontWeight','bold');
            set(gca,'xtick',[1 2 3 4 5 6 7 ],'fontsize',20,'FontWeight','bold');
            set(gca,'ytick',[100 200 300],'fontsize',20,'FontWeight','bold')
        end
        s1 = fscanf(fid,'%s',1);
    end
    fclose(fid);
    ncow=12;ncol=365;
end

