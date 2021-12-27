;2012-12-19 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


begin

;###################################################################################
;path="/mnt/LVM/data/CFS/CFSv2/prate/"
;fpath=systemfunc("ls "+path+"*")

dir  = "/mnt/LVM/data/CFS/CFStmp2m"
cmd  = "find " + dir + "/* -type f -print"
filepath = systemfunc(cmd)
sqsort(filepath) ;Sorts a singly dimensioned arrays of strings.
       
;asciiwrite("/mnt/LVM/langyang/CFSv2/Prate/ex/filepath.txt",filepath )
;filepath有8176个grb2文件

Tavmax=new((/17,2044,280/),float,1e+20) ;/basins,dates,forecast days/
Tavmin=new((/17,2044,280/),float,1e+20) ;/basins,dates,forecast days/
Tempav=new((/17,2044,280/),float,1e+20) ;/basins,dates,forecast days/
Tmax=new((/17,2044,280/),float,1e+20) ;/basins,dates,forecast days/
Tmin=new((/17,2044,280/),float,1e+20) ;/basins,dates,forecast days/
b = addfile("/mnt/LVM/langyang/CFSv2/Basins17/basins17.nc","r") ; read in basins mask basemap file
lsm = b->basins(0:95,0:192)
d = addfile("/mnt/LVM/langyang/CFSv2/Basins17/DEM.nc","r") ; read in basins mask basemap file
DEM = d->DEM(0:95,0:192)
do i =0,dimsizes(filepath)-2, 4
   q=addfile(filepath(i),"r")
   print(filepath(i))
   Temp=q->TMP_P0_L103_GGA0(0:1119,0:95,0:192)
   TempDEM=Temp
   TempDEM=0
   do leadT=0,1119,1
      TempDEM(leadT,:,:)=Temp(leadT,:,:)+DEM*0.006
   end do
   do p_l =0 , 16
     Temp_mask = mask(TempDEM,lsm.eq.(p_l+1),True)        ; mask out provinces number points from Prate			      
     copy_VarMeta(TempDEM,Temp_mask)
     ;;; prepare
     Temp_mask@units="K"
     Temp_mask@long_name="Temperature"
     Temp_AvgLon=dim_avg_n_Wrap(Temp_mask,2)
     Temp_t=dim_avg_n_Wrap(Temp_AvgLon,1)
     Temp_av_max=Temp_t
     Temp_av_min=Temp_t
     Temp_av=Temp_t
     Temp_max = dim_max_n(Temp_mask, (/1,2/)) 
     Temp_min = dim_min_n(Temp_mask, (/1,2/)) 
     
     do t = 0, 1119,4
        Temp_av_max(t)= max((/Temp_av_max(t),Temp_av_max(t+1),Temp_av_max(t+2),Temp_av_max(t+3)/))
        Temp_av_min(t)= min((/Temp_av_min(t),Temp_av_min(t+1),Temp_av_min(t+2),Temp_av_min(t+3)/))
        Temp_av(t)= (Temp_av(t)+ Temp_av(t+1)+ Temp_av(t+2)+ Temp_av(t+3))/4
        Temp_max(t)= max((/Temp_max(t),Temp_max(t+1),Temp_max(t+2),Temp_max(t+3)/))
        Temp_min(t)= min((/Temp_min(t),Temp_min(t+1),Temp_min(t+2),Temp_min(t+3)/))
     end do
     do j = 0, 279
        Tavmax(p_l,i/4,j) =  Temp_av_max(j*4)
        Tavmin(p_l,i/4,j) =  Temp_av_min(j*4)
        Tempav(p_l,i/4,j) =  Temp_av(j*4)
        Tmax(p_l,i/4,j) =  Temp_max(j*4)
        Tmin(p_l,i/4,j) =  Temp_min(j*4)
     end do
   end do
 
end do
;print(x(:,277:279))
;如何把矩阵保存出来为2_00.txt文件
  N    = 2044    ; number of rows
  M    = 280     ; number of columns
  fWidth = 10    ; specify the format width
  fDec   = 0     ; specify the number to the right of decimal point
  fmtx   = M + "f" + fWidth +".2"              ; fmtx="50f10.3"
  
  do p_n = 0 , 16
    opt = True
    opt@fout = "/mnt/LVM/langyang/CFSv2/Basins17/Temp_Basins/Tavmax/"+(p_n+1) + ".txt"  ;output filename
    opt@title  = "1982  1   1   2009  12  31  280  TempavMax"
    write_matrix (Tavmax(p_n,:,:), fmtx, opt)
    opt=False
  end do

  do p_n = 0 , 16
    opt = True
    opt@fout = "/mnt/LVM/langyang/CFSv2/Basins17/Temp_Basins/Tavmin/"+(p_n+1) + ".txt"  ;output filename
    opt@title  = "1982  1   1   2009  12  31  280  TempavMin"
    write_matrix (Tavmin(p_n,:,:), fmtx, opt)
    opt=False
  end do

  do p_n = 0 , 16
    opt = True
    opt@fout = "/mnt/LVM/langyang/CFSv2/Basins17/Temp_Basins/Tempav/"+(p_n+1) + ".txt"  ;output filename
    opt@title  = "1982  1   1   2009  12  31  280  Tempav"
    write_matrix (Tempav(p_n,:,:), fmtx, opt)
    opt=False
  end do

  do p_n = 0 , 16
    opt = True
    opt@fout = "/mnt/LVM/langyang/CFSv2/Basins17/Temp_Basins/Tmax/"+(p_n+1) + ".txt"  ;output filename
    opt@title  = "1982  1   1   2009  12  31  280  TempMax"
    write_matrix (Tmax(p_n,:,:), fmtx, opt)
    opt=False
  end do

  do p_n = 0 , 16
    opt = True
    opt@fout = "/mnt/LVM/langyang/CFSv2/Basins17/Temp_Basins/Tmin/"+(p_n+1) + ".txt"  ;output filename
    opt@title  = "1982  1   1   2009  12  31  280  TempMin"
    write_matrix (Tmin(p_n,:,:), fmtx, opt)
    opt=False
  end do

end
