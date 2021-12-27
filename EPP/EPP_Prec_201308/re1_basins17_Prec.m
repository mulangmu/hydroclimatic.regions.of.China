%**************************************************************
%The cor ccor rhoopt trccor
%**************************************************************
clear ;
basins={'1','2','3','4','5','6','7','8','9','10',...
           '11','12','13','14','15','16','17'} 

figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.84,0.62]);
for kk=1:17 %1:6
  %if (kk~=12) 
    if (kk<10) 
     f=['/mnt/LVM/langyang/EA/EPP_Prec_201308/basins17_0/001/00',num2str(kk),'.txt'];
    else
      f=['/mnt/LVM/langyang/EA/EPP_Prec_201308/basins17_0/001/0',num2str(kk),'.txt'];   
    end
    fid = fopen(f, 'r');  
    k=0;
    s1 = fscanf(fid,'%s',1); 
    while ~feof(fid)
      k=k+1;
      %s1 = fscanf(fid,'%s',1); 
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
        i=mod(kk,6);
        j=fix((kk-0.5)/6);
        if (i==0) i=6; end
        subplot(3,6,kk);
        if (i==1) set(gca,'position',[0.06,0.75-0.33*j,0.10,0.20]);   end;   
        if (i==2) set(gca,'position',[0.22,0.75-0.33*j,0.10,0.20]);   end; 
        if (i==3) set(gca,'position',[0.38,0.75-0.33*j,0.10,0.20]);   end;
        if (i==4) set(gca,'position',[0.54,0.75-0.33*j,0.10,0.20]);   end;
        if (i==5) set(gca,'position',[0.70,0.75-0.33*j,0.10,0.20]);   end;
        if (i==6) set(gca,'position',[0.86,0.75-0.33*j,0.10,0.20]);   end;
        %s1= strtrim(Prov(kk,1:12)) ;
        Title=basins(kk);
        s1=Title{1}
        contourf(1:ncow,1:ncol, x1,'DisplayName',s1,'LineStyle','none');
        caxis([0,1]);  %if (kk==11)  colorbar; end
        title(s1,'FontSize',24,'FontWeight','bold');
        xlabel('Lead time','FontSize',15,'FontWeight','bold');
        ylabel('Date','FontSize',15,'FontWeight','bold');
        set(gca,'xtick',[1 2 3 4 5 6 7 ],'fontsize',13,'FontWeight','bold');
        set(gca,'ytick',[100 200 300],'fontsize',13,'FontWeight','bold')
      end      
      s1 = fscanf(fid,'%s',1); 
    end  
    fclose(fid); 
  %end
end
h=colorbar('FontSize',15);
 %h=colorbar;
% get(h, 'Position')
set(h,'Position',[0.90    0.75-0.33*2    0.02    0.2])
caxis([0 1])

