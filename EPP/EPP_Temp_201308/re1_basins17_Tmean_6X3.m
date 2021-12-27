%**************************************************************
%The cor ccor rhoopt trccor
%**************************************************************
clear ;
basins={'1','2','3','4','5','6','7','8','9','10',...
           '11','12','13','14','15','16','17'} 

figure;
%set(gcf,'unit','normalized','position',[0.1,0.2,0.84,0.62]);
set(gcf,'unit','normalized','position',[0.1,0.1,0.65,1.5]);

%subplot(6,3,17);

for kk=1:17 %1:6
  %if (kk~=12) 
    if (kk<10) 
     f=['/mnt/LVM/langyang/EA/EPP_Temp_201308/Tmean/EPP_Tmean/basins17_0/001/00',num2str(kk),'.txt'];
    else
      f=['/mnt/LVM/langyang/EA/EPP_Temp_201308/Tmean/EPP_Tmean/basins17_0/001/0',num2str(kk),'.txt'];    
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
        i=mod(kk,3);
        j=fix((kk-0.5)/3);
        if (i==0) i=3; end
        subplot(6,3,kk);
        if (i==1) set(gca,'position',[0.12,0.86-0.16*j,0.20,0.1]);   end;   
        if (i==2) set(gca,'position',[0.40,0.86-0.16*j,0.20,0.1]);   end; 
        if (i==3) set(gca,'position',[0.68,0.86-0.16*j,0.20,0.1]);   end;
       
        %if (j==1) set(gca,'position',[0.04,0.75-0.33*i,0.12,0.2]);   end;   
        %if (j==2) set(gca,'position',[0.4,0.75-0.33*i,0.12,0.2]);   end; 
        %if (j==3) set(gca,'position',[0.8,0.75-0.33*i,0.12,0.2]);   end;
        %if (j==4) set(gca,'position',[1.2,0.75-0.33*i,0.12,0.2]);   end;
        %if (j==5) set(gca,'position',[1.6,0.75-0.33*i,0.12,0.2]);   end;
        %if (j==6) set(gca,'position',[2.0,0.75-0.33*i,0.12,0.2]);   end;
        
        %PPP=get(H1,'pos'); 
        %PPP(1)=PPP(1)-0.00;
        %PPP(2)=PPP(2)-0.00;
        %PPP(3)=PPP(3)-0.04;   
        %PPP(4)=PPP(4)-0.02;    
        %set(H1,'pos',PPP)
        
        %subplot(6,3,kk,'replace');
        %subplot('position',[PPP(1)+0.02,PPP(2)+0.02,PPP(3)-0.08,PPP(4)-0.02],'fontsize',8);
        %s1= strtrim(Prov(kk,1:12)) ;
        Title=basins(kk);
        s1=Title{1}
        contourf(1:ncow,1:ncol, x1,'DisplayName',s1,'LineStyle','none');
        caxis([0,1]);  %if (kk==11)  colorbar; end
        title(s1,'FontSize',27,'FontWeight','bold');
        xlabel('Lead time','FontSize',20,'FontWeight','bold');
        ylabel('Date','FontSize',20,'FontWeight','bold');
        set(gca,'xtick',[1 2 3 4 5 6 7 ],'fontsize',17,'FontWeight','bold');
        set(gca,'ytick',[100 200 300],'fontsize',17,'FontWeight','bold')
      end
     
      s1 = fscanf(fid,'%s',1); 
    end  
    fclose(fid); 
  %end
end
 %k=caxis;
 %caxis(k)
 %colorbar
 %colorbar('ytick',[0,0.25,0.5,0.75,1.0],'yticklabel',{'0.0','0.25','0.5','0.75','1.0'});
 h=colorbar('FontSize',20);
 %h=colorbar;
% get(h, 'Position')
set(h,'Position',[0.77    0.91-0.17*5    0.03    0.1])
caxis([0 1])
