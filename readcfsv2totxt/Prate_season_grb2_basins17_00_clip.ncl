load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


begin

;###################################################################################
;path="/mnt/LVM/data/CFS/CFSv2/prate/"
;fpath=systemfunc("ls "+path+"*")

dir  = "/home/data/CFSv2/prate"
cmd  = "find " + dir + "/* -type f -print"
filepath = systemfunc(cmd)
sqsort(filepath) ;Sorts a singly dimensioned arrays of strings.
       
;asciiwrite("/mnt/LVM/langyang/CFSv2/Prate/ex/filepath.txt",filepath )
;filepath有8176个grb2文件
 
x=new((/17,2044,280/),float,1e+20) ;/basins,dates,forecast days/
b = addfile("/home/langyang/CFSv2/Basins17/Prate/basins17.nc","r") ; read in basins mask basemap file
lsm = b->basins(0:95,0:192)

do i =0,dimsizes(filepath)-2, 4
   q=addfile(filepath(i),"r")
   print(filepath(i))
   Prate=q->PRATE_P0_L1_GGA0(0:1119,0:95,0:192)

   do p_l =0 ,16
     Prate_mask = mask(Prate,lsm.eq.(p_l+1),True)        ; mask out basins number points from Prate			
     ;Prate_mask = mask(Prate,(lsm.ne.2),False)
     ;print(Prate_mask(1,:,:))        
     copy_VarMeta(Prate,Prate_mask)
     ;;; prepare
     Prate_mm=Prate_mask*3600*6*100
     copy_VarMeta(Prate_mask,Prate_mm)
     Prate_mm@units="0.01mm"
     Prate_mm@long_name="Precipitation"
     P_AvgLon=dim_avg_n_Wrap(Prate_mm,2)
     P_t=dim_avg_n_Wrap(P_AvgLon,1)
     ;print(P_t(0:10))
     do t = 0, 1119,4
       P_t(t)=P_t(t)+P_t(t+1)+P_t(t+2)+P_t(t+3)
     end do
     do j = 0, 279
       x(p_l,i/4,j) = P_t(j*4)
     end do
   end do
 
end do
;print(x(:,277:279))
;如何把矩阵保存出来为2_00.txt文件
  N    = 2044    ; number of rows
  M    = 280     ; number of columns
  fWidth = 10    ; specify the format width
  fDec   = 0     ; specify the number to the right of decimal point
  fmtx   = M + "f" + fWidth +".0"              ; fmtx="50f10.3"
  do p_n = 0 , 16
    opt = True
    opt@fout = "/home/langyang/CFSv2/Basins17/Prate/Prate_Basins17_results/"+(p_n+1) + ".txt"  ;output filename
    opt@title  = "1982  1   1   2009  12  31  280 "
    write_matrix (x(p_n,:,:), fmtx, opt)
    opt=False
  end do

end
