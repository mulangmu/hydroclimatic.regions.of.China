!This fortran program is main program of Ensemble Pre-processor. 

 
!Observed data of USA
!undef -999.9
!xdef 464 linear -124.9375 0.125
!ydef 224 linear   25.0625 0.125
!prcp 0 99 daily total precipitation (mm/day)
!tmax 0 99 daily max temperature (C)
!tmin 0 99 daily min temperature (C)
!wind 0 99 daily average wind speed (m/s)

program PreProcessP
use PrecStru 

implicit none
 
character(c100) :: tvalue,inputfile, outputfile,Weightfile,BToAfile,ControlFile
character(c1000) :: tempc
type(GFSdata) ::  GFSdata1
type(GFSdata) ::  GFSdata2
type(Obsdata) ::  Obsdata1
type(ASCFile) ::  MyASC,obw
type(CRelEffP) ::  CRE
type(CanonicalE) ::  CanE,CanE1,CanE2,CanE3,CanE4
type(stat24) ::  stat24a
type(stat06) ::  stat06a
type(Ssort) ::  sort1
type(Sgammcf) ::  gammcf1
type(Sthreshold) ::  threshold1
type(Scgauss) ::  cgauss1
type(Scaldat) ::  caldat1
type(Sbeta) :: beta1  
type(Svtrans) :: vtrans1  
type(Sbivar) :: bivar1  
type(Sfcstparm) :: fp1  
type(Sfcst2ensem) :: fe1  
type(Scrps) ::  crps1
type(SreadData) :: rData1 
type(SgenerateCaE) ::  gCanE1 ,gCanE2 
type(SreadEppPara) ::  rEPara1
type(SEPP) ::  epp,epp1,epp2
type(CFSData) ::  CFSData1,CFSData2,CFSData3,CFSData4,CFSData5

integer  i,j,k,ii,jj,l,ll,flagC,flagOs ,status,ifile

integer :: iyear, imonth,iday,ipd,cyear
integer :: irec,nday,BeginYear,EndYear
integer    months(12)
real   w(464,224,4),sump  
real(r8),dimension(464,224,14) :: daily
real(r8)  zavg,zstd,xx,xval,gasdev, GSET
integer*4 irno,ISET
integer julday,jday1,cday
real(r8) tempf,GAMMLN, gammp
real(r8) weibb,weibcv,gausspi,betacf,betai
real(r8) betapi2,vtrans,weibp,binormp2,frho,findrho
 
real(r8) x(365),xts(365),erfcc
real(r8)  x1,y1,yy,wx1,wy1,wx2,wy2
real x2,crps,cumprob
real,pointer ::  rdata(:,:,:)
 
real(r8),pointer ::   try(:,:)          
character(8) :: Tdate


  !read control file
  call getarg (1,ControlFile)
  epp.CoFile=ControlFile
  epp1.CoFile=ControlFile
  CanE.CoFile=ControlFile
  CanE1.CoFile=ControlFile
  CanE2.CoFile=ControlFile
  CanE3.CoFile=ControlFile
 
!***************************************************************************************************** 
!0. input parameter
!***************************************************************************************************** 
  open(unit=10,file=ControlFile)
    do 
      read(10,*,end=10)tempc,tvalue
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('flagC')) read(tvalue,*) flagC 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('flagOs')) read(tvalue,*) flagOs 

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('BeginYear')) read(tvalue,*) GFSdata1.BeginYear
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('EndYear')) read(tvalue,*) GFSdata1.EndYear

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('OldGFSdata1.ncols')) read(tvalue,*) GFSdata1.ncols
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('OldGFSdata1.nrows')) read(tvalue,*) GFSdata1.nrows
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('OldGFSdata1.leadT')) read(tvalue,*) GFSdata1.leadT
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('OldGFSdata1.cellsize')) read(tvalue,*) GFSdata1.cellsize
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('OldGFSdata1.xllcorner')) read(tvalue,*) GFSdata1.xllcorner
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('OldGFSdata1.yllcorner')) read(tvalue,*) GFSdata1.yllcorner
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('NewGFSdata2.ncols')) read(tvalue,*) GFSdata2.ncols
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('NewGFSdata2.nrows')) read(tvalue,*) GFSdata2.nrows
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('NewGFSdata2.leadT')) read(tvalue,*) GFSdata2.leadT
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('NewGFSdata2.cellsize')) read(tvalue,*) GFSdata2.cellsize
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('NewGFSdata2.xllcorner')) read(tvalue,*) GFSdata2.xllcorner
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('NewGFSdata2.yllcorner')) read(tvalue,*) GFSdata2.yllcorner

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Weightfile')) read(tvalue,*) Weightfile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('GFSdata1.BName')) read(tvalue,*) GFSdata1.BName
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('GFSdata2.BName')) read(tvalue,*) GFSdata2.BName
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('GFSdata1.MName')) read(tvalue,*) GFSdata1.MName

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('BeginYear')) read(tvalue,*) Obsdata1.BeginYear
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('EndYear')) read(tvalue,*) Obsdata1.EndYear
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Ngrids')) read(tvalue,*) Obsdata1.SumN
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('originalFile')) read(tvalue,*) Obsdata1.BName
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('6hourFile')) read(tvalue,*) Obsdata1.BName06
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('obw')) read(tvalue,*) Obsdata1.obw
 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Obsdata1.Tim')) read(tvalue,*) Obsdata1.Tim
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Obsdata1.CuVa')) read(tvalue,*) Obsdata1.CuVa
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Obsdata1.ncols')) read(tvalue,*) Obsdata1.ncols
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Obsdata1.nrows')) read(tvalue,*) Obsdata1.nrows
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Obsdata1.NVara')) read(tvalue,*) Obsdata1.NVara
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('Obsdata1.undef')) read(tvalue,*) Obsdata1.undef

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('xllcorner')) read(tvalue,*) MyASC.xllcorner
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('yllcorner')) read(tvalue,*) MyASC.yllcorner
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('cellsize')) read(tvalue,*) MyASC.cellsize

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('BToAfile')) read(tvalue,*) BToAfile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.CFile')) read(tvalue,*) CFSData1.CFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.DFile')) read(tvalue,*) CFSData1.DFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.OutFile')) read(tvalue,*) CFSData1.OutFile 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.EName')) read(tvalue,*) CFSData1.EName 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.leadT')) read(tvalue,*) CFSData1.leadT
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.BeginYear')) read(tvalue,*) CFSData1.BeginYear
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.EndYear')) read(tvalue,*) CFSData1.EndYear
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.ifile')) read(tvalue,*) CFSData1.ifile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.MName')) read(tvalue,*) CFSData1.MName
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.ihour'))   read(tvalue,*) CFSData1.ihour

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData1.WFile')) read(tvalue,*) CFSData1.WFile

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.DFile')) read(tvalue,*) CFSData2.DFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.OutFile')) read(tvalue,*) CFSData2.OutFile 
 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.leadT')) read(tvalue,*) CFSData2.leadT
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.BeginYear')) read(tvalue,*) CFSData2.BeginYear
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.EndYear')) read(tvalue,*) CFSData2.EndYear
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.ifile')) read(tvalue,*) CFSData2.ifile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.MName')) read(tvalue,*) CFSData2.MName
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.EName')) read(tvalue,*) CFSData2.EName 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.nparint')) read(tvalue,*) CFSData2.nparint 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.x1')) read(tvalue,*) CFSData2.x1
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.x2')) read(tvalue,*) CFSData2.x2
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.y1')) read(tvalue,*) CFSData2.y1
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData2.y2')) read(tvalue,*) CFSData2.y2
 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE1.flagPo')) read(tvalue,*) CanE1.flagPo
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE1.leadT')) read(tvalue,*) CanE1.leadT
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE1.ObFile')) read(tvalue,*) CanE1.ObFile

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE1.SiFile')) read(tvalue,*) CanE1.SiFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE1.ReFile')) read(tvalue,*) CanE1.ReFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE1.DoFile')) read(tvalue,*) CanE1.DoFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE1.minwin')) read(tvalue,*) CanE1.minwin
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE1.maxwin')) read(tvalue,*) CanE1.maxwin
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE1.BeginYear ')) read(tvalue,*) CanE1.BeginYear 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE1.EndYear')) read(tvalue,*) CanE1.EndYear 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE1.Events')) then
        read(tvalue,*) CanE1.Events 
        allocate(CanE1.Estart(CanE1.Events))
        allocate(CanE1.Estop(CanE1.Events))
      end if 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE1.Estart'))read(10,*)(CanE1.Estart(i),i=1,CanE1.Events)  
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE1.Estop'))read(10,*)(CanE1.Estop(i),i=1,CanE1.Events)  
   


      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.flagPo')) read(tvalue,*) CanE2.flagPo
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.leadT')) read(tvalue,*) CanE2.leadT
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.nparint')) read(tvalue,*) CanE2.nparint
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.x1')) read(tvalue,*) CanE2.x1
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.x2')) read(tvalue,*) CanE2.x2
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.y1')) read(tvalue,*) CanE2.y1
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.y2')) read(tvalue,*) CanE2.y2
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.ObFile')) read(tvalue,*) CanE2.ObFile

      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.SiFile')) read(tvalue,*) CanE2.SiFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.ReFile')) read(tvalue,*) CanE2.ReFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.DoFile')) read(tvalue,*) CanE2.DoFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.minwin')) read(tvalue,*) CanE2.minwin
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.maxwin')) read(tvalue,*) CanE2.maxwin
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.BeginYear ')) read(tvalue,*) CanE2.BeginYear 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.EndYear')) read(tvalue,*) CanE2.EndYear 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.Events')) then
        read(tvalue,*) CanE2.Events 
        allocate(CanE2.Estart(CanE2.Events))
        allocate(CanE2.Estop(CanE2.Events))
      end if 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.Estart'))read(10,*)(CanE2.Estart(i),i=1,CanE2.Events)  
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CanE2.Estop'))read(10,*)(CanE2.Estop(i),i=1,CanE2.Events)  


      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.ObFile')) read(tvalue,*) CFSData3.ObFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.DFile')) read(tvalue,*) CFSData3.DFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.OutFile')) read(tvalue,*) CFSData3.OutFile 
 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.leadT')) read(tvalue,*) CFSData3.leadT
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.BeginYear')) read(tvalue,*) CFSData3.BeginYear
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.EndYear')) read(tvalue,*) CFSData3.EndYear
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.ifile')) read(tvalue,*) CFSData3.ifile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.MName')) read(tvalue,*) CFSData3.MName
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.EName')) read(tvalue,*) CFSData3.EName 
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.x1')) read(tvalue,*) CFSData3.x1
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.x2')) read(tvalue,*) CFSData3.x2
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.y1')) read(tvalue,*) CFSData3.y1
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData3.y2')) read(tvalue,*) CFSData3.y2


      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData5.CFile')) read(tvalue,*) CFSData5.CFile
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData5.MName')) read(tvalue,*) CFSData5.MName
      if (tempc(1:1).ne.'$'.and.trim(tempc).eq.trim('CFSData5.WFile')) read(tvalue,*) CFSData5.WFile

    
  end do 
10 close(10) 
 


  select case (flagC)
! **************************************************
! 00. debug subroutines and functions
! **************************************************
  case(0) 
    rEPara1.PaFile="E:\EPP4\DATA\Para\EP\464.162.txt"
    rEPara1.Events=14
    call ReadEpara(rEPara1)

    write(*,*) 'GM 0.4',GAMMLN(0.4d0)
  write(*,*) 'GM 1',GAMMLN(1.0d0)
  write(*,*) 'GM 10',GAMMLN(10d0)
  write(*,*) 'GM 100',GAMMLN(100.d0)
  allocate(sort1.indx(10))
  allocate(sort1.a(10))
    sort1.n=10
    sort1.a(1:10)=(/9.,10.,1.,3.,7.,11.,2.,4.,8.,5./)
    call indexx(sort1)
    write(*,*)sort1.a
    write(*,*)sort1.indx
    call sort(sort1)
    write(*,*)sort1.a
    gammcf1.a=0.2
    gammcf1.x=01.   
    call gcf(gammcf1) 
    write(*,*)'gcf:', gammcf1.gammcf,gammcf1.gln
    call gser(gammcf1)
    write(*,*)'gser:', gammcf1.gammcf,gammcf1.gln
    write(*,*)'gammp:',  gammp(gammcf1.a,gammcf1.x +100), gammp(gammcf1.a,gammcf1.x)
    write(*,*)'weibcv(2.)',weibcv(2.0),'weibb(weibcv(2.0))',weibb(weibcv(2.0))
    threshold1.obs=>sort1.a
    threshold1.nobs=10
    threshold1.pop_fraction=0.97
    call threshold (threshold1)
    write(*,*) 'pthresh=', threshold1.pthresh  
    write(*,*) 'gausspi=0.95 0.9', gausspi(0.95) , gausspi(0.9) , gausspi(0.05) , gausspi(0.1) 
    cgauss1.u=0.8
    cgauss1.u0=1.645
    cgauss1.nv=10
    cgauss1.rho=0.6
    cgauss1.vmin=0
    cgauss1.vmax=1
    allocate(cgauss1.vval(cgauss1.nv))
    allocate(cgauss1.cpdfv(cgauss1.nv))
    allocate(cgauss1.ccdfv(cgauss1.nv))
    call cgauss(cgauss1)
    write(*,*)'cgauss standard normal deviate v', cgauss1.vval
    write(*,*)'cgauss probability density', cgauss1.cpdfv
    write(*,*)'cgauss cumulative distribution', cgauss1.ccdfv
    write(*,*)'julday(2,23,1996)',julday(2,28,1997),julday(3,1,1997)
    caldat1.julian=julday(2,23,1996)
    call caldat(caldat1)
  
  write(*,*)'y-m-d',caldat1.iyyy,caldat1.mm,caldat1.id 
  beta1.a=0.1
  beta1.b=0.9
  beta1.f=0.5

  write(*,*)'betacf(a,b,x)', betacf(beta1)
  write(*,*)'betai(a,b,x)', betai(beta1)
  beta1.a=0.5
  beta1.b=0.5
  beta1.f=0.6

  write(*,*)'betapi2(a,b,x)',betapi2(beta1)
    write(*,*)'weibp:',weibp(0.5,0.5,0.6)
  vtrans1.avg=0.5
  vtrans1.cv=0.9
  vtrans1.pz=0.3
  vtrans1.obs=0.6
  vtrans1.ivopt=0
  write(*,*)'vtrans0:', vtrans(vtrans1)
  vtrans1.ivopt=1
  write(*,*)'vtrans1:', vtrans(vtrans1)
  vtrans1.ivopt=2
  write(*,*)'vtrans2:', vtrans(vtrans1)
  vtrans1.ivopt=03
  write(*,*)'vtrans3:', vtrans(vtrans1)
  vtrans1.ivopt=4
  write(*,*)'vtrans4:', vtrans(vtrans1)
  vtrans1.ivopt=5
  write(*,*)'vtrans5:', vtrans(vtrans1)
  vtrans1.ivopt=6
  write(*,*)'vtrans6:', vtrans(vtrans1)
  vtrans1.ivopt=7
  write(*,*)'vtrans7:', vtrans(vtrans1)
    !bivar1=Sbivar(ax,bx,cx,tol,rho22,p11,p12,p21,p22,umin,u0,umax,vmin,v0,vmax,rho)
    bivar1=Sbivar(0.02,0.5,0.98,0.001,0.2,0.3,0.4,0.5,0.6,0,0.3,0.9,0.1,0.2,0.8,0.7)
  write(*,*)'frho(bivar1):', frho(bivar1) 
  write(*,*)'binormp2(bivar1):', binormp2(bivar1)
  write(*,*)'findrho(bivar1):',  findrho(bivar1)
  call estrho3 (bivar1)
    write(*,*)'estrho3 (bivar1):', bivar1.rho
    fp1.nobs=25
    fp1.pthresh1=1.9
    fp1.pthresh2=1.9
    fp1.iobstran=1
    fp1.ifcsttran=1
    fp1.cor_weight=0.5
    fp1.iopt_rho=1
    !1 cor        ! correlation coefficient between obs and fcst (includes zeros)
    !2 rhoopt     ! correlation coefficient from transformed forecasts and obs
    !3 w*cor + (1. - w)*rhoopt
    !4 ccor       ! correlation coefficient between wet obs and fcsts (>0)
    !5 trccor     ! correlation coefficient between transformed wet obs and fcsts (>0)
    allocate(fp1.obsx(fp1.nobs))
    allocate(fp1.fcstx(fp1.nobs)) 
  fp1.obsx(1:25)=(/0.8,0.9,1.8,2.1,2.4,2.7,3.,3.3,3.6,3.9,4.2,4.5,4.8,5.1,5.4,5.7,6.,6.3,6.6,6.9,7.2,7.5,7.8,8.1,8.4/)
  fp1.fcstx(1:25)=(/0.,2.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.,10.5,11.,11.5,12./)

    call epp_fcstparms(fp1) 

    write(*,*)'cavgobsx',fp1.cavgobsx
    write(*,*)'cavgfcstx',fp1.cavgfcstx
    write(*,*)'avgobs',fp1.avgobs
    write(*,*)'popobs',fp1.popobs
    write(*,*)'cavgobs',fp1.cavgobs
    write(*,*)'ccvobs',fp1.ccvobs
    write(*,*)'avgfcst',fp1.avgfcst
    write(*,*)'popfcst',fp1.popfcst
    write(*,*)'cavgfcst',fp1.cavgfcst
    write(*,*)'ccvfcst',fp1.ccvfcst
    write(*,*)'rho_est',fp1.rho_est
    write(*,*)'rmsfcst',fp1.rmsfcst
    write(*,*)'effnsfcst',fp1.effnsfcst
    write(*,*)'ratio',fp1.ratio
    write(*,*)'eqts_fcst',fp1.eqts_fcst
    write(*,*)'cor',fp1.cor
    write(*,*)'ccor',fp1.ccor
    write(*,*)'trccor',fp1.trccor
    write(*,*)'rhoopt',fp1.rhoopt
    write(*,*)'nvalcal',fp1.nvalcal
    write(*,*)'avgobscal',fp1.avgobscal
    write(*,*)'stdobscal',fp1.stdobscal
    write(*,*)'avgfcstcal',fp1.avgfcstcal
    write(*,*)'stdfcstcal',fp1.stdfcstcal

    x=0
  do i=1,365,3
    x(i)=1.0*i
  end do
    call fill (3,x,xts)
    !write(*,*)xts
    x2=0.1 
    x1=x2
    write(*,*)'erfcc(0.3)',erfcc(x1)
     



    fe1.fcst=10.d0
    fe1.ifcst=1
    fe1.pthresh_fcst=0.5d0
    fe1.fpop=0.6d0
    fe1.fcavg=5.d0
    fe1.fccv=1.1d0

    fe1.iobs=1
    fe1.pthresh_obs=0.5d0
    fe1.obspop=0.6d0
    fe1.obscavg=6.d0
    fe1.obsccv=1.2d0
    fe1.rho=0.8d0
    fe1.npp=5
    allocate(fe1.pp(fe1.npp+1))
    call fcst2ensem(fe1)

    write(*,*)'genavg ',  fe1.genavg 
    write(*,*)'genpop ',  fe1.genpop 
    write(*,*)'gencavg ', fe1.gencavg 
    write(*,*)'genccv ',  fe1.genccv 

    write(*,*)'vmax ',     fe1.vmax 
    write(*,*)'exprobmin ',fe1.exprobmin 
    write(*,*)'pmin ',     fe1.pmin 
    write(*,*)'pmax ',     fe1.pmax 
    write(*,*)'fe1.pp ',   fe1.pp


    crps1.n=10
    crps1.x0=03.0
  allocate(crps1.xx(crps1.n))
  allocate(try(10,10))
    crps1.xx(1:10)=(/0.,2.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5/)
  crps1.n=10
    write(*,*)'crps(crps1)', crps(crps1)
    try(1:10,1)=(/0.,2.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5/)
    crps1.xx=try(:,1)
    write(*,*)'crps(crps1)', crps(crps1)
    write(*,*)'cumprob(crps1)', cumprob(crps1)



!***************************************************************************************************** 
!1. inverse distance weighting method for interpolation. 
     !GFS precipitation data 2.5 by 2.5 was interpolated to 0.125 by 0.125 
!***************************************************************************************************** 
  case(1) 
    GFSdata1.ndays=julday(1,1,GFSdata1.EndYear)-julday(1,1,GFSdata1.BeginYear)+1
    i=GFSdata1.EndYear-GFSdata1.BeginYear+1
    allocate(GFSdata1.Gdatat(GFSdata1.ncols,GFSdata1.nrows,i,12,31,GFSdata1.leadT))   ! GFS [Lon,lat,lt] 
    allocate(GFSdata2.Gdata(i,12,31,GFSdata2.leadT))   ! GFS [Lon,lat,lt]   
    !2.1calculate the Inverse distance weighting
    GFSdata1.flagOs=flagOs
    GFSdata2.flagOs=flagOs

    MyASC.ncols=GFSdata2.ncols 
    MyASC.nrows=GFSdata2.nrows 
    MyASC.NODATA_value=-9999
    MyASC.xllcorner=GFSdata2.xllcorner
    MyASC.yllcorner=GFSdata2.yllcorner
    MyASC.cellsize=GFSdata2.cellsize

    allocate(MyASC.MValue(MyASC.ncols,MyASC.nrows))
    MyASC.flagAF=1
    MyASC.MValue=0

    if (flagOs .eq. 2)   call slashchange(GFSdata1.MName)
    call  ReadASC(obw,GFSdata1.MName)
    ! USA 0.125: X1=X0+i*0.125, Y1=Y0+j*0.125    
    ! GFS 2.5:   X2=X0+i*2.5,   Y2=Y0+j*2.5
    do i=1,GFSdata2.ncols
      do j=1,GFSdata2.nrows
     
        ii=int((GFSdata2.xllcorner-GFSdata1.xllcorner+(i-1)*GFSdata2.cellsize)/GFSdata1.cellsize)+1
        jj=int((GFSdata2.yllcorner-GFSdata1.yllcorner+(j-1)*GFSdata2.cellsize)/GFSdata1.cellsize)+1
        ii=max(1,ii)
        jj=max(1,jj)  
        ii=min(ii,GFSdata1.ncols-1)
        jj=min(jj,GFSdata1.nrows-1)  
        sump=0      
        w(i,j,1)= (GFSdata2.xllcorner-GFSdata1.xllcorner+(i-1)*GFSdata2.cellsize-(ii-1)*GFSdata1.cellsize)*&
          (GFSdata2.xllcorner-GFSdata1.xllcorner+(i-1)*GFSdata2.cellsize-(ii-1)*GFSdata1.cellsize)+&
          (GFSdata2.yllcorner-GFSdata1.yllcorner+(j-1)*GFSdata2.cellsize-(jj-1)*GFSdata1.cellsize)*&
          (GFSdata2.yllcorner-GFSdata1.yllcorner+(j-1)*GFSdata2.cellsize-(jj-1)*GFSdata1.cellsize)
        w(i,j,1)=1/w(i,j,1)
        sump=sump+w(i,j,1)
        w(i,j,2)= (GFSdata2.xllcorner-GFSdata1.xllcorner+(i-1)*GFSdata2.cellsize-(ii)*GFSdata1.cellsize)*&
          (GFSdata2.xllcorner-GFSdata1.xllcorner+(i-1)*GFSdata2.cellsize-(ii)*GFSdata1.cellsize)+&
          (GFSdata2.yllcorner-GFSdata1.yllcorner+(j-1)*GFSdata2.cellsize-(jj-1)*GFSdata1.cellsize)*&
          (GFSdata2.yllcorner-GFSdata1.yllcorner+(j-1)*GFSdata2.cellsize-(jj-1)*GFSdata1.cellsize)
        w(i,j,2)=1/w(i,j,2)
        sump=sump+w(i,j,2)
        w(i,j,3)= (GFSdata2.xllcorner-GFSdata1.xllcorner+(i-1)*GFSdata2.cellsize-(ii)*GFSdata1.cellsize)*&
          (GFSdata2.xllcorner-GFSdata1.xllcorner+(i-1)*GFSdata2.cellsize-(ii)*GFSdata1.cellsize)+&
          (GFSdata2.yllcorner-GFSdata1.yllcorner+(j-1)*GFSdata2.cellsize-(jj)*GFSdata1.cellsize)*&
          (GFSdata2.yllcorner-GFSdata1.yllcorner+(j-1)*GFSdata2.cellsize-(jj)*GFSdata1.cellsize)
        w(i,j,3)=1/w(i,j,3)
        sump=sump+w(i,j,3)
        w(i,j,4)= (GFSdata2.xllcorner-GFSdata1.xllcorner+(i-1)*GFSdata2.cellsize-(ii-1)*GFSdata1.cellsize)*&
          (GFSdata2.xllcorner-GFSdata1.xllcorner+(i-1)*GFSdata2.cellsize-(ii-1)*GFSdata1.cellsize)+&
          (GFSdata2.yllcorner-GFSdata1.yllcorner+(j-1)*GFSdata2.cellsize-(jj)*GFSdata1.cellsize)*&
          (GFSdata2.yllcorner-GFSdata1.yllcorner+(j-1)*GFSdata2.cellsize-(jj)*GFSdata1.cellsize)
        w(i,j,4)=1/w(i,j,4)
        sump=sump+w(i,j,4)
        do k=1,4
          w(i,j,k)=w(i,j,k)/sump
        end do  
      end do
    end do
    if (flagOs .eq. 2)   call slashchange(Weightfile)
    open(10,file=Weightfile,form='unformatted',access='direct',recl=464*224*4)   !open (10,file=outputfile,form='unformatted')
      write(10,rec=1)w 
    close(10)

    BeginYear=GFSdata1.BeginYear
    EndYear=GFSdata1.EndYear
  
    !2.2. read GFS precipitation data         
    call ReadGFS(GFSdata1)
    if (flagOs .eq. 2)   call slashchange(GFSdata2.BName)

    !2.3. IDW interpolate GFS to 1/8 degree
    do i=GFSdata2.ncols,1,-1
    if (flagOs .eq. 2) then
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'mkdir ',trim(GFSdata2.BName),i             
      else
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'md ',trim(GFSdata2.BName),i   
      end if
      call system(outputfile)
      do j=1,GFSdata2.nrows
        if (abs(obw.MValue(i,j)-obw.NODATA_value)>0.01) then
        print*,'x,y=         ',i,j 
        write(outputfile,'(2(a,i3.3),a)')trim(GFSdata2.BName),i,'\',j,'.txt'        
        if (flagOs .eq. 2)   call slashchange(outputfile)
        open(20,file=trim(outputfile))!,status='unknown',form='unformatted',access='direct',recl=464*224*14
          !write(tempc,'(5a,400(2a,i4))')'Year',char(9),'Month',char(9),'Day', (char(9),'F',l,l=1,GFSdata2.leadT) 
          !call delblank(tempc)
          !write (20,'(a)')trim(tempc)
          write (20,'(i4,5a,i4,4a)') &
              BeginYear,char(9),'1',char(9),'1',char(9),EndYear,char(9),'12',char(9),'1'
            
      ! 'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12','F13','F14'
          do iyear=BeginYear,EndYear  
            cyear=iyear-BeginYear+1   
            if((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0)then
              months = (/31,29,31,30,31,30,31,31,30,31,30,31/)
            else
              months = (/31,28,31,30,31,30,31,31,30,31,30,31/)
            endif
            do imonth=1,12
              do iday=1,months(imonth)
                ii=int((GFSdata2.xllcorner-GFSdata1.xllcorner+(i-1)*GFSdata2.cellsize)/GFSdata1.cellsize)+1
                jj=int((GFSdata2.yllcorner-GFSdata1.yllcorner+(j-1)*GFSdata2.cellsize)/GFSdata1.cellsize)+1
                ii=max(1,ii)
                jj=max(1,jj)  
                ii=min(ii,GFSdata1.ncols-1)
                jj=min(jj,GFSdata1.nrows-1)  
         
                do l=1,GFSdata2.leadT
                  GFSdata2.Gdata(cyear,imonth,iday,l)=GFSdata1.Gdatat(ii,jj,cyear,imonth,iday,l)*w(i,j,1)+&
                  GFSdata1.Gdatat(ii+1,jj,cyear,imonth,iday,l)*w(i,j,2)+&
                  GFSdata1.Gdatat(ii+1,jj+1,cyear,imonth,iday,l)*w(i,j,3)+&
                  GFSdata1.Gdatat(ii,jj+1,cyear,imonth,iday,l)*w(i,j,4)
                end do
                MyASC.MValue(i,j)=MyASC.MValue(i,j)+GFSdata2.Gdata(cyear,imonth,iday,1)
                !2.4. output 0.125 degree GFS data
                !write(tempc,'(i4,a1,i2.2,a1,i2.2,365(a,i7))')&
                !  iyear,char(9),imonth,char(9),iday,(char(9),int(GFSdata2.Gdata(cyear,imonth,iday,l)*100),l=1,GFSdata2.leadT)
                write(tempc,'(365(a,i7))')(char(9),int(GFSdata2.Gdata(cyear,imonth,iday,l)*100),l=1,GFSdata2.leadT)
                 

                call delblank(tempc)
                write (20,'(a)')trim(tempc)

                !write(20,'(i4,a1,i2.2,a1,i2.2,14(TR1,f6.2))')&
                !  iyear,'-',imonth,'-',iday,(GFSdata2.Gdata(cyear,imonth,iday,l),l=1,GFSdata2.leadT)
              end do
            end do
          end do
        close(20)
          MyASC.MValue(i,j)=MyASC.MValue(i,j)/(GFSdata1.EndYear-GFSdata1.BeginYear+1)
        else
          MyASC.MValue(i,j)=MyASC.NODATA_value
        end if

      end do
    end do
    write(outputfile,'(a,a4,i4,i4,a4)')trim(GFSdata2.BName),'Mean',&
    Obsdata1.BeginYear,Obsdata1.EndYear,'.asc'
    if (flagOs .eq. 2)   call slashchange(outputfile)
    call WriteASC(MyASC,outputfile)


    deallocate(MyASC.MValue)   
    deallocate(obw.MValue)   

    deallocate(GFSdata1.Gdatat)
    deallocate(GFSdata2.Gdata)
!***************************************************************************************************** 
!2.read binary format daily map data and write 6 hourly data x.y file (sequential binary) .     
!The output precipitation units in the x.y are in 'mm' 
!***************************************************************************************************** 
  case(2)
    !2.1 read the observed data 
    zavg = .0001
    zstd = .00001
    irno = 129874633
    ISET= 0
    k=Obsdata1.EndYear-Obsdata1.BeginYear+1
    if (flagOs .eq. 2)   call slashchange(Obsdata1.obw)
    call  ReadASC(obw,Obsdata1.obw)
    allocate(Obsdata1.Obdata(Obsdata1.ncols,Obsdata1.nrows/Obsdata1.Tim,k,12,31))  !in order to save memory of computer
    allocate(Obsdata1.Obdata06(k,12,31,4))  
    !Obsdata1.CuVa=1
    Obsdata1.flagOs=flagOs

    MyASC.ncols=Obsdata1.ncols 
    MyASC.nrows=Obsdata1.nrows 
    MyASC.NODATA_value=-9999
    allocate(MyASC.MValue(MyASC.ncols,MyASC.nrows))
    MyASC.flagAF=1 
    MyASC.MValue=0
    if (flagOs .eq. 2) call   slashchange(Obsdata1.BName06)
    !=========================================================== 
  !precipitation 
    !=========================================================== 
    if (Obsdata1.CuVa.eq.1) then
    do k= 1,Obsdata1.Tim
    Obsdata1.k=k 
    call ReadObdata(Obsdata1)
    do i= 1, Obsdata1.ncols 
      if (flagOs .eq. 2) then
        write(outputfile,'(2a,i3.3)')'mkdir ',trim(Obsdata1.BName06),i         
      else
        write(outputfile,'(2a,i3.3)')'md ',trim(Obsdata1.BName06),i   
      end if
      call system(outputfile)
      do j= 1, Obsdata1.nrows/Obsdata1.Tim
        if (abs(Obsdata1.Obdata(i,j,1,1,1)-Obsdata1.undef)>0.01) then
          print *,'x=',i,'  y=',j+(k-1)*Obsdata1.nrows/Obsdata1.Tim
          do iyear=Obsdata1.BeginYear,Obsdata1.EndYear  
            cyear=iyear-Obsdata1.BeginYear+1
            if((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0)then
              months = (/31,29,31,30,31,30,31,31,30,31,30,31/)
            else
              months = (/31,28,31,30,31,30,31,31,30,31,30,31/)
            endif
            do imonth=1,12 
              do iday=1,months(imonth)
                xx=0
                if (obw.MValue(i,j+(k-1)*Obsdata1.nrows/Obsdata1.Tim)>0.001) &
                  Obsdata1.Obdata(i,j,cyear,imonth,iday)=&
                  Obsdata1.Obdata(i,j,cyear,imonth,iday)/obw.MValue(i,j+(k-1)*Obsdata1.nrows/Obsdata1.Tim)
                do ipd=1,4
                  xval = zavg + zstd*gasdev(irno,GSET,ISET)!get a random data
                  Obsdata1.Obdata06(cyear,imonth,iday,ipd)=Obsdata1.Obdata(i,j,cyear,imonth,iday)/4.0+xval 
                  xx=xx+xval
                end do  !ipd
                Obsdata1.Obdata(i,j,cyear,imonth,iday)=Obsdata1.Obdata(i,j,cyear,imonth,iday)+xx
                MyASC.MValue(i,j+(k-1)*Obsdata1.nrows/Obsdata1.Tim)=&
                  MyASC.MValue(i,j+(k-1)*Obsdata1.nrows/Obsdata1.Tim)+Obsdata1.Obdata(i,j,cyear,imonth,iday)
              end do !iday
            end do !imonth     
          end do !iyear
          MyASC.MValue(i,j+(k-1)*Obsdata1.nrows/Obsdata1.Tim)=MyASC.MValue(i,j+(k-1)*Obsdata1.nrows/Obsdata1.Tim)/&
                    (Obsdata1.EndYear-Obsdata1.BeginYear+1)
          !2.2 Output the 6 hours observed data (a point to a file)
          Obsdata1.x=i 
          Obsdata1.y=j+(k-1)*Obsdata1.nrows/Obsdata1.Tim
          call WriteObdata06(Obsdata1)

          !2.3 compute map06 statistics
          Obsdata1.x=i 
          Obsdata1.y=j 
          stat24a.pthresh=0.25
          stat24a.BeginYear=Obsdata1.BeginYear
          stat24a.EndYear=Obsdata1.EndYear
          stat24a.Obdata=>Obsdata1.Obdata(i,j,:,:,:)
          call mapstats24 (stat24a) 
          Obsdata1.avg(5,:)=stat24a.avg
          Obsdata1.pop(5,:)=stat24a.pop
          Obsdata1.cavg(5,:)=stat24a.cavg
          Obsdata1.ccv(5,:)=stat24a.ccv
          Obsdata1.npos(5,:)=stat24a.npos
          Obsdata1.nobs(5,:)=stat24a.nobs

          stat06a.pthresh=0.25
          stat06a.BeginYear=Obsdata1.BeginYear
          stat06a.EndYear=Obsdata1.EndYear
          stat06a.Obdata06=>Obsdata1.Obdata06
          call mapstats06 (stat06a) 
          Obsdata1.avg(1:4,:)=stat06a.avg
          Obsdata1.pop(1:4,:)=stat06a.pop
          Obsdata1.cavg(1:4,:)=stat06a.cavg
          Obsdata1.ccv(1:4,:)=stat06a.ccv
          Obsdata1.npos(1:4,:)=stat06a.npos
          Obsdata1.nobs(1:4,:)=stat06a.nobs

          Obsdata1.x=i 
          Obsdata1.y=j+(k-1)*Obsdata1.nrows/Obsdata1.Tim

          call Writemapstats(Obsdata1)
        else
          MyASC.MValue(i,j+(k-1)*Obsdata1.nrows/Obsdata1.Tim)=MyASC.NODATA_value
        endif  !-999.9
      end do !j
    end do !i
    end do !k
   
    write(outputfile,'(2a,i4,i4,a4)')trim(Obsdata1.BName06),'Mean',&
      Obsdata1.BeginYear,Obsdata1.EndYear,'.asc'
    if (flagOs .eq. 2)   call slashchange(outputfile)
    call WriteASC(MyASC,outputfile)
    end if !precipitation

    !=========================================================== 
  !Temperature  
    !=========================================================== 
    if (Obsdata1.CuVa.eq.2.or.Obsdata1.CuVa.eq.3 ) then
    do k=1,Obsdata1.Tim
    Obsdata1.k=k 
 
    call ReadObdata(Obsdata1)
    do i=1, Obsdata1.ncols 
      if (flagOs .eq. 2) then
        write(outputfile,'(2a,i3.3)')'mkdir ',trim(Obsdata1.BName06),i         
      else
        write(outputfile,'(2a,i3.3)')'md ',trim(Obsdata1.BName06),i   
      end if
      call system(outputfile)
      do j=1, Obsdata1.nrows/Obsdata1.Tim
        if (abs(Obsdata1.Obdata(i,j,1,1,1)-Obsdata1.undef)>0.01) then
          print *,'x=',i,'  y=',j+(k-1)*Obsdata1.nrows/Obsdata1.Tim
          do iyear=Obsdata1.BeginYear,Obsdata1.EndYear  
            cyear=iyear-Obsdata1.BeginYear+1
            if((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0)then
              months = (/31,29,31,30,31,30,31,31,30,31,30,31/)
            else
              months = (/31,28,31,30,31,30,31,31,30,31,30,31/)
            endif
            do imonth=1,12 
              do iday=1,months(imonth)

                MyASC.MValue(i,j+(k-1)*Obsdata1.nrows/Obsdata1.Tim)=MyASC.MValue(i,j+(k-1)*Obsdata1.nrows/Obsdata1.Tim)+&
                              Obsdata1.Obdata(i,j,cyear,imonth,iday)
                Obsdata1.Obdata(i,j,cyear,imonth,iday)=Obsdata1.Obdata(i,j,cyear,imonth,iday)+273.16

              end do !iday
            end do !imonth     
          end do !iyear
          MyASC.MValue(i,j+(k-1)*Obsdata1.nrows/Obsdata1.Tim)=MyASC.MValue(i,j+(k-1)*Obsdata1.nrows/Obsdata1.Tim)/&
                    (Obsdata1.EndYear-Obsdata1.BeginYear+1)/365
          !2.2 Output the 6 hours observed data (a point to a file)
          Obsdata1.x=i 
          Obsdata1.y=j+(k-1)*Obsdata1.nrows/Obsdata1.Tim 
          call WriteObdataT(Obsdata1)

          !2.3 compute map06 statistics
          Obsdata1.x=i 
          Obsdata1.y=j 
          stat24a.pthresh=0.25
          stat24a.BeginYear=Obsdata1.BeginYear
          stat24a.EndYear=Obsdata1.EndYear
          stat24a.Obdata=>Obsdata1.Obdata(i,j,:,:,:)
          call mapstats24 (stat24a) 
          Obsdata1.avg(5,:)=stat24a.avg
          Obsdata1.pop(5,:)=stat24a.pop
          Obsdata1.cavg(5,:)=stat24a.cavg
          Obsdata1.ccv(5,:)=stat24a.ccv
          Obsdata1.npos(5,:)=stat24a.npos
          Obsdata1.nobs(5,:)=stat24a.nobs
 
          Obsdata1.x=i 
          Obsdata1.y=j+(k-1)*Obsdata1.nrows/Obsdata1.Tim 

          call WritemapstatsT(Obsdata1)
        else
          MyASC.MValue(i,j+(k-1)*Obsdata1.nrows/Obsdata1.Tim)=MyASC.NODATA_value
        endif  !-999.9
      end do !j
    end do !i
    end do !k
   
    write(outputfile,'(2a,i4,i4,a4)')trim(Obsdata1.BName06),'Mean',&
      Obsdata1.BeginYear,Obsdata1.EndYear,'.asc'
    if (flagOs .eq. 2)   call slashchange(outputfile)
    call WriteASC(MyASC,outputfile)
    end if !Temperature max


    deallocate(MyASC.MValue)   
    deallocate(Obsdata1.Obdata06) 
    deallocate(Obsdata1.Obdata)
!***************************************************************************************************** 
!3.  read binary file ,write ascII file 
!                  
!***************************************************************************************************** 
  case(3)
    !call  BinToAdata(BToAfile)
    call BinToAdata06(BToAfile)

!***************************************************************************************************** 
!4. Calculate the Nash-Sutcliffe efficiency  and correlation coefficient of Canonical Events
!                  
!***************************************************************************************************** 
  case(4)
 
    if (flagOs .eq. 2)   call slashchange(CanE1.ObFile)
    if (flagOs .eq. 2)   call slashchange(CanE1.SiFile)
    if (flagOs .eq. 2)   call slashchange(CanE1.ReFile)
    if (flagOs .eq. 2)   call slashchange(CanE1.DoFile)
    CanE1.flagOs=flagOs
    call  ReadASC(obw,CanE1.DoFile)
    CanE1.win=CanE1.minwin
    i=max(CanE1.leadT,CanE1.Events) 
    allocate(CanE1.r(365,i))
    allocate(CanE1.DC(365,i))
    allocate(CanE1.Bias(365,i))
    allocate(CanE1.RMSE(365,i))

    CRE.TimeSum=(CanE1.EndYear-CanE1.BeginYear+1)*(CanE1.win+1)
    allocate(CRE.yc(CRE.TimeSum))
    allocate(CRE.yy(CRE.TimeSum))
    jday1=JULDAY(1,1,CanE1.BeginYear)
    CanE1.ndays=JULDAY(12,31,CanE1.EndYear)-jday1+1
    allocate(CanE1.ob(CanE1.ndays,CanE1.leadT))
    allocate(CanE1.si(CanE1.ndays,CanE1.leadT))
    allocate(CanE1.Eob(CanE1.ndays,CanE1.Events))
    allocate(CanE1.Esi(CanE1.ndays,CanE1.Events))
    rData1.flagOs=flagOs
    rData1.ObFile=CanE1.ObFile
    rData1.SiFile=CanE1.SiFile

    do ii=1, obw.ncols  

    if (flagOs .eq. 2) then
        write(outputfile,'(2a,i3.3)')'mkdir ',trim(CanE1.ReFile),ii         
      else
        write(outputfile,'(2a,i3.3)')'md ',trim(CanE1.ReFile),ii   
      end if
    call system(outputfile)

      do jj=1, obw.nrows 
      if (abs(obw.MValue(ii,jj)-obw.NODATA_value)>0.01) then
        print *,'x=',ii,'  y=',jj 
        CanE1.x=ii
        CanE1.y=jj
        !4.1 read observed data and forecasting data
        !call readCaE(CanE1) 
        !rData1.SiFile=CanE1.SiFile
        rData1.x =CanE1.x 
        rData1.y=CanE1.y
        rData1.BeginYear=CanE1.BeginYear
        rData1.EndYear=CanE1.EndYear
        rData1.leadT=CanE1.leadT
        rData1.ndays=CanE1.ndays
        rData1.ob=>CanE1.ob
        rData1.si=>CanE1.si
      rData1.nparint=1
        call readDataob(rData1)
        call readDatafcst(rData1)

        !4.2  statistical information between observed and forecasting
        if (CanE1.flagPo.eq.1) then
        do j=1,365
          do i=1,CanE1.leadT
            l=0
            do iyear=CanE1.BeginYear,CanE1.EndYear
              cday=JULDAY(1,1,iyear)+j-jday1
              do k=max(1,cday-CanE1.win/2),min(CanE1.ndays,cday+CanE1.win/2)
                l=l+1
                CRE.yy(l)=CanE1.ob(k,i)
                CRE.yc(l)=CanE1.si(k,i)
              end do
            end do
            CRE.TimeSum=l
            call CRelEff(CRE)
            CanE1.r(j,i)=CRE.r
            CanE1.DC(j,i)=max(0.0,CRE.DC)
            CanE1.Bias(j,i)=CRE.Bias
            CanE1.RMSE(j,i)=CRE.RMSE
          end do 
        end do      
        call WritestatOF(CanE1)
        end if
        !4.3 Generate Canonical Events 


        gCanE1.leadT=CanE1.leadT
        gCanE1.ndays=CanE1.ndays
        gCanE1.Events=CanE1.Events
        gCanE1.ob=>CanE1.ob
        gCanE1.Estart=>CanE1.Estart
        gCanE1.Estop=>CanE1.Estop
        gCanE1.Eob=>CanE1.Eob
        call generateCaE(gCanE1)

        gCanE1.ob=>CanE1.si
        gCanE1.Eob=>CanE1.Esi                 
        call generateCaE(gCanE1)

        !call generateCaE(CanE1)
        !4.4  statistical information between observed and forecasting events
        do j=1,365
          do i=1,CanE1.Events
            l=0
            do iyear=CanE1.BeginYear,CanE1.EndYear
              cday=JULDAY(1,1,iyear)+j-jday1
              do k=max(1,cday-CanE1.win/2),min(CanE1.ndays,cday+CanE1.win/2)
                l=l+1
                CRE.yy(l)=CanE1.Eob(k,i)
                CRE.yc(l)=CanE1.Esi(k,i)
              end do
            end do
            CRE.TimeSum=l
            call CRelEff(CRE)
            CanE1.r(j,i)=CRE.r
            CanE1.DC(j,i)=max(0.0,CRE.DC)
            CanE1.Bias(j,i)=CRE.Bias
            CanE1.RMSE(j,i)=CRE.RMSE
          end do 
        end do      
        call WritestatOFE(CanE1)

        end if
      end do
    end do


    deallocate(CanE1.Eob)
    deallocate(CanE1.Esi)
    deallocate(CanE1.ob)
    deallocate(CanE1.si)
    deallocate(CRE.yc)
    deallocate(CRE.yy)
    deallocate(CanE1.r)
    deallocate(CanE1.DC)
    deallocate(CanE1.Bias)
    deallocate(CanE1.RMSE)




!***************************************************************************************************** 
!5.from historical observed data and forecast to get parameter of EPP model
!                  GFS Data
!***************************************************************************************************** 
  case(5)
  CanE.CoFile=ControlFile
  call ReadCEpara(CanE)
    if (flagOs .eq. 2) then
      call slashchange(CanE.ObFile)
      call slashchange(CanE.SiFile)
      call slashchange(CanE.ReFile)
      call slashchange(CanE.DoFile)
    end if
    CanE.flagOs=flagOs
    call  ReadASC(obw,CanE.DoFile)  
    jday1=JULDAY(1,1,CanE.BeginYear)
    CanE.ndays=JULDAY(12,31,CanE.EndYear)-jday1+1
    call NewspaceCE(CanE)
     
    rData1.flagOs=flagOs
    rData1.ObFile=CanE.ObFile
    rData1.SiFile=CanE.SiFile
    rData1.BeginYear=CanE.BeginYear
    rData1.EndYear=CanE.EndYear
    rData1.leadT=CanE.leadT
    rData1.ndays=CanE.ndays
    rData1.ob=>CanE.ob
    rData1.si=>CanE.si
    rData1.nparint=1

    do ii=141,141 !3,3 !464,464 ! obw.ncols,1,-1    
    if (flagOs .eq. 2) then
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'mkdir ',trim(CanE.ReFile),ii             
      else
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'md ',trim(CanE.ReFile),ii   
      end if
    call system(outputfile)

      do jj=55,55 !186,186 !162,162 ! obw.nrows,1,-1
        if (abs(obw.MValue(ii,jj)-obw.NODATA_value)>0.01) then
          write(*,'(a2,i4,tr1,a2,i4)')'x=',ii,'y=',jj 
          CanE.x=ii
          CanE.y=jj
          !5.1 read observed data and forecasting data
          rData1.x=CanE.x 
          rData1.y=CanE.y
          call readDataob(rData1)
          if (flagOs .eq. 2)   then
            write(rData1.SiFile,'(a,2(i3.3,a1),i3.3,a4)')trim(CanE.SiFile),CanE.x,'/',CanE.x,'.',CanE.y,'.txt'          
          else
            write(rData1.SiFile,'(a,2(i3.3,a1),i3.3,a4)')trim(CanE.SiFile),CanE.x,'\',CanE.x,'.',CanE.y,'.txt'  
          end if

          call readDatafcst(rData1)
          !call readCaE(CanE) 
          !5.2 from historical observed data and forecast to get parameter of EPP model
          call EppPara(CanE)
          !5.3 Write  paramater information of EPP
          call WriteCEpara(CanE)

        end if
      end do
    end do

    call deletspaceCE(CanE)

 

 

!***************************************************************************************************** 
!6.hindcast---GFS
!                  
!***************************************************************************************************** 
  case(6)
    epp.CoFile=ControlFile
    call ReadHpara(epp)
    if (flagOs .eq. 2) then
      call slashchange(epp.ObFile)
      call slashchange(epp.SiFile)
      call slashchange(epp.PaFile)
      call slashchange(epp.DoFile)
      call slashchange(epp.ReFile)
    end if
    !write(tempc,'(tr1,3a)')"  sdfa ada",char(9)," bg a mon th er  "
  !call delblank(tempc)
    epp.flagOs=flagOs
    call  ReadASC(obw,epp.DoFile)  
    jday1=JULDAY(1,1,epp.BeginYear)
    epp.iyear=epp.BeginYear   
    epp.imonth=1  
    epp.iday=1   

    epp.ndays=JULDAY(12,31,epp.EndYear)-jday1+1
    call NewspaceEP(epp)
     
    rEPara1.Events=epp.Events
    rEPara1.pthresh_obs=>epp.pthresh_obs
    rEPara1.pthresh_fcst=>epp.pthresh_fcst
    rEPara1.popobs=>epp.popobs
    rEPara1.cavgobs=>epp.cavgobs
    rEPara1.ccvobs=>epp.ccvobs
    rEPara1.popfcst=>epp.popfcst
    rEPara1.cavgfcst=>epp.cavgfcst
    rEPara1.ccvfcst=>epp.ccvfcst
    rEPara1.rho_est=>epp.rho_est
    
    rData1.flagOs=flagOs 
    rData1.EpFile=epp.EpFile
    rData1.ObFile=epp.ObFile
    rData1.SiFile=epp.SiFile
    rData1.BeginYear=epp.BeginYear
    rData1.EndYear=epp.EndYear
    rData1.BeginYearEns=epp.BeginYearEns
    rData1.EndYearEns=epp.EndYearEns
    rData1.leadT=epp.leadT
    rData1.ndays=epp.ndays
    rData1.ep=>epp.ep
    rData1.ob=>epp.ob
    rData1.si=>epp.si
    rData1.nparint=1
    do ii= 141,141  !464,464 !3,3  ! obw.ncols,1,-1    
      if (flagOs .eq. 2) then
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'mkdir ',trim(epp.ReFile),ii             
      else
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'md ',trim(epp.ReFile),ii   
      end if
      call system(outputfile)


      do jj=55,55 !162,162 !186,186  ! obw.nrows,1,-1  
      !do ii=464,464 ! obw.ncols,1,-1    
      !  do jj=162,162 ! obw.nrows,1,-1
        if (abs(obw.MValue(ii,jj)-obw.NODATA_value)>0.01) then
          write(*,'(a2,i4,tr1,a2,i4)')'x=',ii,'y=',jj 
          epp.x=ii
          epp.y=jj

          !6.1 read observed data and forecasting data
          rData1.x=epp.x 
          rData1.y=epp.y
          if (flagOs .eq. 2)   then
            write(rData1.SiFile,'(a,2(i3.3,a1),i3.3,a4)')trim(epp.SiFile),epp.x,'/',epp.x,'.',epp.y,'.txt'          
          else
            write(rData1.SiFile,'(a,2(i3.3,a1),i3.3,a4)')trim(epp.SiFile),epp.x,'\',epp.x,'.',epp.y,'.txt'  
          end if

          call readDataob(rData1)
          call readDataEp(rData1)
          call readDatafcst(rData1)

          !6.1 read parameter of EPP
          if (flagOs .eq. 2)   then
            write(rEPara1.PaFile,'(a,i3.3,a1,i3.3,a4)')trim(epp.PaFile),rData1.x,'/',rData1.y,'.txt'        
          else
            write(rEPara1.PaFile,'(a,i3.3,a1,i3.3,a4)')trim(epp.PaFile),rData1.x,'\',rData1.y,'.txt' 
          end if
       
          call ReadEpara(rEPara1)
          epp.nparint1 =1
          call Hindcast(epp)
          !5.2 from historical observed data and forecast to get parameter of EPP model
          !call EppPara(epp)
          !5.3 Write  paramater information of EPP
          !call WriteCEpara(epp)
          call Writehincast(epp)

        end if
      end do
    end do

    call deletspaceEP(epp)
!***************************************************************************************************** 
!7.CFS data manage 
!                  
! inverse distance weighting method for interpolation. 
! CFS precipitation data 0.938 (gaussian weights) was interpolated to 0.125 by 0.125 
!***************************************************************************************************** 
  case(7)
 
   if (flagOs .eq. 2)   call slashchange(CFSData1.CFile) 
  
   call ReadCFScoord(CFSData1)
   allocate(CFSdata1.MValue00(CFSdata1.leadT*4,CFSData1.ncols,CFSData1.nrows))    
   allocate(CFSdata1.MValue06(CFSdata1.leadT*4,CFSData1.ncols,CFSData1.nrows))    
   allocate(CFSdata1.MValue12(CFSdata1.leadT*4,CFSData1.ncols,CFSData1.nrows))    
   allocate(CFSdata1.MValue18(CFSdata1.leadT*4,CFSData1.ncols,CFSData1.nrows))    
   allocate(CFSdata1.MValue1(CFSdata1.leadT,CFSData1.ncols,CFSData1.nrows))    
   allocate(CFSdata1.PValue(CFSdata1.leadT))    
   allocate(CFSdata1.PIValue(CFSdata1.leadT))    
   
  
    CFSData1.ndays=julday(12,31,CFSData1.EndYear)-julday(1,1,CFSData1.BeginYear)+1
    i=GFSdata1.EndYear-GFSdata1.BeginYear+1
    !7.1calculate the Inverse distance weighting
    CFSData1.flagOs=flagOs
 
    call  calweight(CFSdata1)
    obw=CFSData1.obw
    allocate(CFSData1.MValueN(CFSData1.leadT*4,obw.ncols,obw.nrows))    

    BeginYear=CFSdata1.BeginYear
    EndYear=CFSdata1.EndYear

    !7.2. IDW interpolate CFS to 1/8 degree
    do iyear=BeginYear,EndYear  
      !cyear=iyear-BeginYear+1   
      CFSdata1.iyear=iYear
      jday1=JULDAY(1,1,CFSdata1.iYear)
      do cday=1,366  !1,365
        CFSdata1.cday=cday        
        caldat1.julian=jday1+CFSdata1.cday-1
        call caldat(caldat1)  
        imonth=caldat1.mm
        iday=caldat1.id
        !CFSdata1.ihour=0
        !7.3. read CFS precipitation data   
        !call ReadCFS4(CFSdata1)
        call ReadCFS(CFSdata1)
        if (CFSdata1.alive ) then
          CFSdata1.iyear=iyear
          CFSdata1.imonth=imonth
          CFSdata1.iday=iday
          call CFSinterpolate(CFSData1)
          if (trim(CFSdata1.EName).eq.trim('tmp2m')) then
            CFSdata1.EName= trim('tmin2m')
            do l=1,CFSdata1.leadT
              CFSData1.MValue(l,:,:)=CFSData1.MValue1(l,:,:)
            end do     
            call CFSinterpolate(CFSData1)       
            CFSdata1.EName= trim('tmp2m')
          end if
        end if !file exist
      end do
    end do
    !do i=1,obw.ncols
    !  do j=1,obw.nrows   
    !    if (abs(obw.MValue(i,j)-obw.NODATA_value)>0.01) then
    !      MyASC.MValue(i,j)=MyASC.MValue(i,j)/(CFSdata1.EndYear-CFSdata1.BeginYear+1)
    !    else
    !      MyASC.MValue(i,j)=MyASC.NODATA_value
    !    end if
    !  end do
    !end do
              
    !write(outputfile,'(a,a4,i4,i4,a4)')trim(CFSData1.OutFile),'Mean',&
    !CFSData1.BeginYear,CFSData1.EndYear,'.asc'
    !if (flagOs .eq. 2)   call slashchange(outputfile)
    !call WriteASC(MyASC,outputfile)


     
    deallocate(obw.MValue)   

    deallocate(CFSdata1.PValue)
    deallocate(CFSdata1.PIValue)
    deallocate(CFSdata1.MValue00)
    deallocate(CFSdata1.MValue06)
    deallocate(CFSdata1.MValue12)
    deallocate(CFSdata1.MValue18)
    deallocate(CFSdata1.MValue1)

    deallocate(CFSdata1.MValueN)

  


!***************************************************************************************************** 
!8.CFS data manage  CFS data manage from date form to point form
!                  
! 
!***************************************************************************************************** 
  case(8)
 

    if (CFSData2.nparint.eq.5 ) then
      CFSdata2.monthday(1,1:8)=(/7,1,6,11,16,21,26,31/)
      CFSdata2.monthday(2,1:8)=(/5,5,10,15,20,25,0,0/)
      CFSdata2.monthday(3,1:8)=(/6,2,7,12,17,22,27,0/)
      CFSdata2.monthday(4,1:8)=(/6,1,6,11,16,21,26,0/)
      CFSdata2.monthday(5,1:8)=(/7,1,6,11,16,21,26,31/)
      CFSdata2.monthday(6,1:8)=(/6,5,10,15,20,25,30,0/)
      CFSdata2.monthday(7,1:8)=(/6,5,10,15,20,25,30,0/)
      CFSdata2.monthday(8,1:8)=(/6,4,9,14,19,24,29,0/)
      CFSdata2.monthday(9,1:8)=(/6,3,8,13,18,23,28,0/)
      CFSdata2.monthday(10,1:8)=(/6,3,8,13,18,23,28,0/)
      CFSdata2.monthday(11,1:8)=(/6,2,7,12,17,22,27,0/)
      CFSdata2.monthday(12,1:8)=(/6,2,7,12,17,22,27,0/)
    end if
    if (CFSData2.nparint.eq.1 ) then
      CFSdata2.monthday(1,:)=(/31,1:31/)  
      CFSdata2.monthday(2,:)=(/28,1:31/)  
      CFSdata2.monthday(3,:)=(/31,1:31/)  
      CFSdata2.monthday(4,:)=(/30,1:31/)  
      CFSdata2.monthday(5,:)=(/31,1:31/)  
      CFSdata2.monthday(6,:)=(/30,1:31/)  
      CFSdata2.monthday(7,:)=(/31,1:31/)  
      CFSdata2.monthday(8,:)=(/31,1:31/)  
      CFSdata2.monthday(9,:)=(/30,1:31/)  
      CFSdata2.monthday(10,:)=(/31,1:31/)  
      CFSdata2.monthday(11,:)=(/30,1:31/)  
      CFSdata2.monthday(12,:)=(/31,1:31/)  
	end if
    allocate(CFSdata2.PValue(CFSdata2.leadT))    
    CFSData2.flagOs=flagOs
    if (flagOs .eq. 2)   call slashchange(CFSData2.MName)
    if (flagOs .eq. 2)   call slashchange(CFSData2.OutFile)
    call  ReadASC(obw,CFSData2.MName)   
    CFSData2.obw=obw
    BeginYear=CFSData2.BeginYear
    EndYear=CFSData2.EndYear
    if (CFSData2.nparint.eq.5)  CFSdata2.ndays=73*(EndYear-BeginYear+1)
    if (CFSData2.nparint.eq.1)  then 
      CFSdata2.ndays=julday(12,31,EndYear)-julday(1,1,BeginYear)+1
      !366*(EndYear-BeginYear+1)
    end if
    allocate(CFSdata2.PLValue(CFSdata2.ndays,obw.nrows,CFSdata2.leadT))    
    MyASC.ncols=obw.ncols 
    MyASC.nrows=obw.nrows 
    MyASC.NODATA_value=-9999
    MyASC.xllcorner=obw.xllcorner
    MyASC.yllcorner=obw.yllcorner
    MyASC.cellsize=obw.cellsize

    allocate(MyASC.MValue(MyASC.ncols,MyASC.nrows))
    MyASC.flagAF=1
    MyASC.MValue=0

    CFSdata2.ncols=obw.ncols
    CFSdata2.nrows=obw.nrows
    !allocate(CFSdata2.MValue(CFSdata2.leadT,CFSData2.ncols,CFSData2.nrows))    

    !do i=1,CFSdata2.ncols !,1,-1
    !  if (flagOs .eq. 2) then
    !    write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'mkdir ',trim(CFSdata2.OutFile),i             
    !  else
    !    write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'md ',trim(CFSdata2.OutFile),i   
    !  end if
    !  call system(outputfile)
    !end do
    CFSdata2.Crow=1 
    CFSdata2.Crow1=1 
    CFSdata2.nstep=0
    do i=CFSData2.x1,CFSdata2.x2   !CFSdata2.ncols,CFSdata2.ncols*0.75,-1
      if (flagOs .eq. 2) then
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'mkdir ',trim(CFSdata2.OutFile),i             
      else
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'md ',trim(CFSdata2.OutFile),i   
      end if
      call system(outputfile)
      CFSdata2.ix=i
      CFSdata2.cday=0
      CFSdata2.Crow1=CFSdata2.Crow  
      do iyear=BeginYear,EndYear
        Print *,'x=',i,'    Year=',iyear  
        CFSdata2.iYear=iyear 
        do imonth=1,12
          CFSdata2.imonth=imonth
          k=CFSdata2.monthday(imonth,1)
          if(((mod(iyear,4)==0.and.mod(iyear,100)/=0).or.mod(iyear,400)==0) &  !leap year
           .and.imonth.eq.2.and.CFSData2.nparint.eq.1) k=k+1
          do iday=1,k
            CFSdata2.cday=CFSdata2.cday+1
            !read CFS data               
            CFSdata2.iday=CFSdata2.monthday(imonth,iday+1)
            !call ReadDCFS(CFSdata2) 
            CFSdata2.Crow=CFSdata2.Crow1           
            call ReadDCFSb(CFSdata2)   
			!print *,'j,k',j,CFSdata2.Crow                              
          end do
        end do
      end do
       
      do j=1,CFSdata2.nrows
        if (abs(obw.MValue(i,j)-obw.NODATA_value)>0.01) then
        print*,'x,y=         ',i,j 
        write(outputfile,'(2(a,i3.3),a)')trim(CFSdata2.OutFile),i,'\',j,'.txt'        
        if (flagOs .eq. 2)   call slashchange(outputfile)
        CFSdata2.iy=j
        !if (CFSdata2.PLValue(1,1,1).gt.-1)
        open(20,file=trim(outputfile)) 
          write (20,'(i4,5a,i4,5a,i3)') BeginYear,char(9),'1',char(9),'1',&
            char(9),EndYear,char(9),'12',char(9),'31',char(9),CFSdata2.leadT
          CFSdata2.cday=0  
          do iyear=BeginYear,EndYear
            !Print *,'xy=',i,j,'    Year=',iyear               
            do imonth=1,12  
              k=CFSdata2.monthday(imonth,1)
              if(((mod(iyear,4)==0.and.mod(iyear,100)/=0).or.mod(iyear,400)==0) &  !leap year
                .and.imonth.eq.2.and.CFSData2.nparint.eq.1)k=k+1			             
              do iday=1,k !CFSdata2.monthday(imonth,1)  
                CFSdata2.cday=CFSdata2.cday+1       
                MyASC.MValue(i,j)=MyASC.MValue(i,j)+CFSdata2.PLValue(CFSdata2.cday,j,1)
				if (trim(CFSData2.EName).eq.trim('prate')) then
                  write(tempc,'(365(a,i7))')(char(9),CFSdata2.PLValue(CFSdata2.cday,j,l),l=1,CFSdata2.leadT) 
                else
                  write(tempc,'(365(a,i7))')(char(9),CFSdata2.PLValue(CFSdata2.cday,j,l)-27316,l=1,CFSdata2.leadT) 
                endif
				 				                  
                call delblank(tempc)
                write (20,'(a)')trim(tempc)
              end do
            end do
          end do
        close(20)
        !end if
        MyASC.MValue(i,j)=MyASC.MValue(i,j)/(CFSdata2.EndYear-CFSdata2.BeginYear+1)
        else
          MyASC.MValue(i,j)=MyASC.NODATA_value
        end if
      end do
    end do

    MyASC.MValue=MyASC.MValue/(CFSdata2.EndYear-CFSdata2.BeginYear+1)
    write(outputfile,'(a,a4,i4,i4,a4)')trim(CFSData1.OutFile),'Mean',&
      CFSData1.BeginYear,CFSData1.EndYear,'.asc'
    if (flagOs .eq. 2)   call slashchange(outputfile)
    call WriteASC(MyASC,outputfile)

    deallocate(MyASC.MValue)   
    deallocate(obw.MValue)   

    deallocate(CFSdata2.PValue)
    deallocate(CFSdata2.PLValue)
    !deallocate(CFSdata2.MValue)

!***************************************************************************************************** 
!9. Calculate the Nash-Sutcliffe efficiency  and correlation coefficient of Canonical Events-CFS
!                  
!***************************************************************************************************** 
  case(9)
    if (flagOs .eq. 2)   call slashchange(CanE2.ObFile)
    if (flagOs .eq. 2)   call slashchange(CanE2.SiFile)
    if (flagOs .eq. 2)   call slashchange(CanE2.ReFile)
    if (flagOs .eq. 2)   call slashchange(CanE2.DoFile)
    CanE2.flagOs=flagOs
    call  ReadASC(obw,CanE2.DoFile)
    CanE2.win=CanE2.minwin  !/CanE2.nparint
    i=max(CanE2.leadT,CanE2.Events) 
    allocate(CanE2.r(365,i))
    allocate(CanE2.DC(365,i))
    allocate(CanE2.Bias(365,i))
    allocate(CanE2.RMSE(365,i))
    CRE.TimeSum=(CanE2.EndYear-CanE2.BeginYear+1)*(CanE2.win+1)
    allocate(CRE.yc(CRE.TimeSum))
    allocate(CRE.yy(CRE.TimeSum))
    jday1=JULDAY(1,1,CanE2.BeginYear)
    CanE2.ndays=JULDAY(12,31,CanE2.EndYear)-jday1+1
    CanE2.ndays1= (CanE2.EndYear- CanE2.BeginYear+1)*73
    allocate(CanE2.ob(CanE2.ndays+CanE2.leadT,CanE2.leadT))
    allocate(CanE2.si(CanE2.ndays,CanE2.leadT))
    allocate(CanE2.Eob(CanE2.ndays,CanE2.Events))
    allocate(CanE2.Esi(CanE2.ndays,CanE2.Events))
    rData1.flagOs=flagOs
    rData1.ObFile=CanE2.ObFile
    rData1.SiFile=CanE2.SiFile
    rData1.nparint=5
    do ii=CanE2.x1,CanE2.x2 !1, 50 !obw.ncols  
      if (flagOs .eq. 2) then
        write(outputfile,'(2a,i3.3)')'mkdir ',trim(CanE2.ReFile),ii         
      else
        write(outputfile,'(2a,i3.3)')'md ',trim(CanE2.ReFile),ii   
      end if
      call system(outputfile)
      rData1.BeginYear=CanE2.BeginYear
      rData1.EndYear=CanE2.EndYear
      rData1.leadT=CanE2.leadT
      rData1.ndays=CanE2.ndays
      rData1.ndays1=CanE2.ndays1
      rData1.ob=>CanE2.ob
      rData1.si=>CanE2.si
      do jj=CanE2.y1,CanE2.y2 !1, obw.nrows 
      if (abs(obw.MValue(ii,jj)-obw.NODATA_value)>0.01) then
        print *,'x=',ii,'  y=',jj 
        CanE2.x=ii
        CanE2.y=jj
        !9.1 read observed data and forecasting data
        !call readCaE(CanE2) 
        !rData1.SiFile=CanE2.SiFile
        rData1.x =CanE2.x 
        rData1.y=CanE2.y
        call readDataob(rData1)
        call readDatafcst(rData1)
        !9.2  statistical information between observed and forecasting
        if (CanE2.flagPo.eq.1) then
        do j=1,365,CanE2.nparint
          do i=1,CanE2.leadT
            l=0
            ll=0
            do iyear=CanE2.BeginYear,CanE2.EndYear
              cday=JULDAY(1,1,iyear)+j-jday1    
        if(((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0).and.j.gt.59) &  !leap year
                 cday=cday+1               
              do k=max(1,cday-CanE2.win/2),min(CanE2.ndays,cday+CanE2.win/2),CanE2.nparint
                l=l+1
                CRE.yy(l)=CanE2.ob(k,i)
                CRE.yc(l)=CanE2.si(k,i)
              end do
              !cday=(iyear-CanE2.BeginYear)*(365/CanE2.nparint)+j/CanE2.nparint+1              
              !do k=max(1,cday-CanE2.win/2/CanE2.nparint),min(CanE2.ndays1,cday+CanE2.win/2/CanE2.nparint) 
                !CRE.yy(l)=CanE2.ob(k,i)
              !  ll=ll+1 
              !  CRE.yc(ll)=CanE2.si(k,i)
              !end do
            end do
            CRE.TimeSum=l
            call CRelEff(CRE)
            CanE2.r(j,i)=CRE.r
            CanE2.DC(j,i)=max(0.0,CRE.DC)
            CanE2.Bias(j,i)=CRE.Bias
            CanE2.RMSE(j,i)=CRE.RMSE
          end do 
        end do  
      do i=1,CanE2.leadT
          call fill(CanE2.nparint,CanE2.r(:,i),CanE2.r(:,i))   
          call fill(CanE2.nparint,CanE2.DC(:,i),CanE2.DC(:,i))   
          call fill(CanE2.nparint,CanE2.Bias(:,i),CanE2.Bias(:,i))   
          call fill(CanE2.nparint,CanE2.RMSE(:,i),CanE2.RMSE(:,i))   
      end do     
        call WritestatOF(CanE2)
        end if
        !9.3 Generate Canonical Events 
        gCanE2.leadT=CanE2.leadT
        gCanE2.ndays=CanE2.ndays
        gCanE2.Events=CanE2.Events
        gCanE2.Estart=>CanE2.Estart
        gCanE2.Estop=>CanE2.Estop
        gCanE2.ob=>CanE2.ob
        gCanE2.Eob=>CanE2.Eob
        call generateCaE(gCanE2)
        !gCanE2.ndays=CanE2.ndays1
        gCanE2.ob=>CanE2.si
        gCanE2.Eob=>CanE2.Esi                 
        call generateCaE(gCanE2)
        !9.4  statistical information between observed and forecasting events
        do j=1,365,CanE2.nparint
          do i=1,CanE2.Events
            l=0
            ll=0
            do iyear=CanE2.BeginYear,CanE2.EndYear
              cday=JULDAY(1,1,iyear)+j-jday1

        if(((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0).and.j.gt.59) &  !leap year
                 cday=cday+1               
              do k=max(1,cday-CanE2.win/2),min(CanE2.ndays,cday+CanE2.win/2),CanE2.nparint
                l=l+1
                CRE.yy(l)=CanE2.Eob(k,i)
                CRE.yc(l)=CanE2.si(k,i)
              end do
              !cday=(iyear-CanE2.BeginYear)*(365/CanE2.nparint)+j/CanE2.nparint+1              
              !do k=max(1,cday-CanE2.win/2/CanE2.nparint),min(CanE2.ndays1,cday+CanE2.win/2/CanE2.nparint) 
              !  !CRE.yy(l)=CanE2.ob(k,i)
              !  ll=ll+1 
              !  CRE.yc(ll)=CanE2.Esi(k,i)
              !end do 
            end do
            CRE.TimeSum=l
            call CRelEff(CRE)
            CanE2.r(j,i)=CRE.r
            CanE2.DC(j,i)=max(0.0,CRE.DC)
            CanE2.Bias(j,i)=CRE.Bias
            CanE2.RMSE(j,i)=CRE.RMSE
          end do 
        end do  
      do i=1,CanE2.leadT
          call fill(CanE2.nparint,CanE2.r(:,i),CanE2.r(:,i))   
          call fill(CanE2.nparint,CanE2.DC(:,i),CanE2.DC(:,i))   
          call fill(CanE2.nparint,CanE2.Bias(:,i),CanE2.Bias(:,i))   
          call fill(CanE2.nparint,CanE2.RMSE(:,i),CanE2.RMSE(:,i))   
      end do     
        
        call WritestatOFE(CanE2)

        end if
      end do
    end do


    deallocate(CanE2.Eob)
    deallocate(CanE2.Esi)
    deallocate(CanE2.ob)
    deallocate(CanE2.si)
    deallocate(CRE.yc)
    deallocate(CRE.yy)
    deallocate(CanE2.r)
    deallocate(CanE2.DC)
    deallocate(CanE2.Bias)
    deallocate(CanE2.RMSE)

 

!***************************************************************************************************** 
!10.from historical observed data and forecast to get parameter of EPP model  for CFS+GFS
!                  
!***************************************************************************************************** 
  case(10)
    CanE3.CoFile=ControlFile
    call ReadCEparaCG(CanE3)
    if (flagOs .eq. 2) then
      call slashchange(CanE3.ObFile)
      call slashchange(CanE3.SiFile)
      call slashchange(CanE3.SiFile1)
      call slashchange(CanE3.ReFile)
      call slashchange(CanE3.DoFile)
    end if
    CanE3.flagOs=flagOs
    call  ReadASC(obw,CanE3.DoFile)  
    jday1=JULDAY(1,1,CanE3.BeginYear)
    CanE3.ndays=JULDAY(12,31,CanE3.EndYear)-jday1+1
    call NewspaceCE(CanE3)
     
    rData1.flagOs=flagOs
    rData1.ObFile=CanE3.ObFile
    rData1.SiFile=CanE3.SiFile
    rData1.BeginYear=CanE3.BeginYear
    rData1.EndYear=CanE3.EndYear
    rData1.leadT=CanE3.leadT
    rData1.ndays=CanE3.ndays
    rData1.ob=>CanE3.ob
    rData1.si=>CanE3.si
    rData1.nparint=1
    rData1.flagT=CanE3.flagT
    do ii=CanE3.x1,CanE3.x2  !3,3 !464,464 ! obw.ncols,1,-1    
      if (flagOs.eq.2) then
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'mkdir ',trim(CanE3.ReFile),ii             
      else
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'md ',trim(CanE3.ReFile),ii   
      end if
      call system(outputfile)

      do jj=CanE3.y1,CanE3.y2 !186,186 !162,162 ! obw.nrows,1,-1
        if (abs(obw.MValue(ii,jj)-obw.NODATA_value)>0.01) then
          write(*,'(a2,i4,tr1,a2,i4)')'x=',ii,'y=',jj 
          CanE3.x=ii
          CanE3.y=jj
          !5.1 read observed data and forecasting data
          rData1.x=CanE3.x 
          rData1.y=CanE3.y   
          rData1.leadT=CanE3.leadT
          if (CanE3.flagT.eq.0) then
            call readDataob(rData1)
          else
            call readDataobT(rData1)        
          end if
          !read CFS
          rData1.si=-99.0
          rData1.nparint=CanE3.nparint
          rData1.SiFile=CanE3.SiFile1  
          rData1.leadT=CanE3.leadT
          call readDatafcst(rData1)
          !if (CanE3.flagT.ne.0)  rData1.si=rData1.si-200.0
          if (CanE3.flagCG.eq.1) then
            !read GFS
            rData1.nparint=0
            rData1.SiFile=CanE3.SiFile 
            rData1.leadT=CanE3.leadT1
            call readDatafcst(rData1)
          end if 

          !5.2 from historical observed data and forecast to get parameter of EPP model
          k=CanE3.nparint
          CanE3.nparint=5
          call EppPara(CanE3)
          CanE3.nparint=k
          !5.3 Write  paramater information of EPP
          call WriteCEparaN(CanE3)

        end if
      end do
    end do

    call deletspaceCE(CanE3)

 !***************************************************************************************************** 
!11.hindcast for CFS
!                  
!***************************************************************************************************** 
  case(11)
    epp1.CoFile=ControlFile
    call ReadCHpara(epp1)
    if (flagOs .eq. 2) then
      call slashchange(epp1.EpFile)
      call slashchange(epp1.ObFile)
      call slashchange(epp1.SiFile)
      call slashchange(epp1.SiFile1)
      call slashchange(epp1.PaFile)
      call slashchange(epp1.DoFile)
      call slashchange(epp1.ReFile)
    end if
    !write(tempc,'(tr1,3a)')"  sdfa ada",char(9)," bg a mon th er  "
    !call delblank(tempc)
    epp1.flagOs=flagOs
    call  ReadASC(obw,epp1.DoFile)  
    jday1=JULDAY(1,1,epp1.BeginYear)
	epp1.iyear=epp1.BeginYear   
    epp1.imonth=1  
    epp1.iday=1   

    epp1.ndays=JULDAY(12,31,epp1.EndYear)-jday1+1
    epp1.ndays1=JULDAY(12,31,epp1.EndYearEns)-JULDAY(1,1,epp1.BeginYearEns)+1
    call NewspaceEP(epp1)
     
    rEPara1.Events=epp1.Events
    rEPara1.pthresh_obs=>epp1.pthresh_obs
    rEPara1.pthresh_fcst=>epp1.pthresh_fcst
    rEPara1.popobs=>epp1.popobs
    rEPara1.cavgobs=>epp1.cavgobs
    rEPara1.ccvobs=>epp1.ccvobs
    rEPara1.popfcst=>epp1.popfcst
    rEPara1.cavgfcst=>epp1.cavgfcst
    rEPara1.ccvfcst=>epp1.ccvfcst
    rEPara1.rho_est=>epp1.rho_est

    rData1.flagT=epp1.flagT
    rData1.flagOs=flagOs
    rData1.EpFile=epp1.EpFile
    rData1.ObFile=epp1.ObFile
    rData1.SiFile=epp1.SiFile
    rData1.BeginYear=epp1.BeginYear
    rData1.EndYear=epp1.EndYear
    rData1.BeginYearEns=epp1.BeginYearEns
    rData1.EndYearEns=epp1.EndYearEns
    rData1.leadT=epp1.leadT
    rData1.ndays=epp1.ndays
    rData1.ep=>epp1.ep
    rData1.ob=>epp1.ob
    rData1.si=>epp1.si
    rData1.nparint=1
    do ii= epp1.x1,epp1.x2 ! 141,141  !464,464 !3,3  ! obw.ncols,1,-1    
      if (flagOs .eq. 2) then
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'mkdir ',trim(epp1.ReFile),ii             
      else
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'md ',trim(epp1.ReFile),ii   
      end if
      call system(outputfile)

      do jj= epp1.y1,epp1.y2 !55,55 !162,162 !186,186  ! obw.nrows,1,-1  
      !do ii=464,464 ! obw.ncols,1,-1    
      !  do jj=162,162 ! obw.nrows,1,-1
        if (abs(obw.MValue(ii,jj)-obw.NODATA_value)>0.01) then
          write(*,'(a2,i4,tr1,a2,i4)')'x=',ii,'y=',jj 
          epp1.x=ii
          epp1.y=jj

          !6.1 read observed data and forecasting data
          rData1.x=epp1.x 
          rData1.y=epp1.y
          rData1.nparint=epp1.nparint
          rData1.leadT=epp1.leadT
          if (epp1.flagT.eq.0) then
            call readDataob(rData1)
            call readDataEp(rData1)
          else
            call readDataobT(rData1)        
            call readDataEpT(rData1)
          end if

          !call readDataob(rData1)

          !read CFS
          rData1.si=-99.0
          rData1.nparint=epp1.nparint
          rData1.SiFile=epp1.SiFile1  
          rData1.leadT=epp1.leadT
          call readDatafcst(rData1)
          !read GFS
          if (epp1.flagCG.eq.1) then
            !read GFS
            rData1.nparint=0
            rData1.SiFile=epp1.SiFile 
            rData1.leadT=epp1.leadT1
            call readDatafcst(rData1)
          end if 
 
   

          !6.1 read parameter of epp1
          if (flagOs .eq. 2)   then
            write(rEPara1.PaFile,'(a,i3.3,a1,i3.3,a4)')trim(epp1.PaFile),rData1.x,'/',rData1.y,'.txt'        
          else
            write(rEPara1.PaFile,'(a,i3.3,a1,i3.3,a4)')trim(epp1.PaFile),rData1.x,'\',rData1.y,'.txt' 
          end if
       
          call ReadEpara(rEPara1)

          call Hindcast(epp1)


          !call Writehincast(epp1)
          call Writehincastb(epp1)

        end if
      end do
    end do

    call deletspaceEP(epp1)



!***************************************************************************************************** 
!12.CFS members correlation
!                  
! 
!***************************************************************************************************** 
  case(12)
 

    !CFSdata3.monthday(:,1)=(/7,5,6,6,7,6,6,6,6,6,6,6/)
    CFSdata3.monthday(1,1:8)=(/7,1,6,11,16,21,26,31/)
    CFSdata3.monthday(2,1:8)=(/5,5,10,15,20,25,0,0/)
    CFSdata3.monthday(3,1:8)=(/6,2,7,12,17,22,27,0/)
    CFSdata3.monthday(4,1:8)=(/6,1,6,11,16,21,26,0/)
    CFSdata3.monthday(5,1:8)=(/7,1,6,11,16,21,26,31/)
    CFSdata3.monthday(6,1:8)=(/6,5,10,15,20,25,30,0/)
    CFSdata3.monthday(7,1:8)=(/6,5,10,15,20,25,30,0/)
    CFSdata3.monthday(8,1:8)=(/6,4,9,14,19,24,29,0/)
    CFSdata3.monthday(9,1:8)=(/6,3,8,13,18,23,28,0/)
    CFSdata3.monthday(10,1:8)=(/6,3,8,13,18,23,28,0/)
    CFSdata3.monthday(11,1:8)=(/6,2,7,12,17,22,27,0/)
    CFSdata3.monthday(12,1:8)=(/6,2,7,12,17,22,27,0/)
    allocate(CFSdata3.PValue(CFSdata3.leadT))    
    CFSdata3.flagOs=flagOs
    if (flagOs .eq. 2)   call slashchange(CFSdata3.MName)
    if (flagOs .eq. 2)   call slashchange(CFSdata3.OutFile)
    call  ReadASC(obw,CFSdata3.MName)   
    CFSdata3.obw=obw
    BeginYear=CFSdata3.BeginYear
    EndYear=CFSdata3.EndYear
    CFSdata3.ndays=73*(EndYear-BeginYear+1)
    allocate(CFSdata3.PLValue(CFSdata3.ndays,5,CFSdata3.leadT))    
    allocate(CFSdata3.cor(5,5,CFSdata3.leadT+2))    
    rData1.ndays=(EndYear-BeginYear+1)*366
    allocate(rData1.ob(rData1.ndays+CFSdata3.leadT,CFSdata3.leadT))    
    jday1=JULDAY(1,1,CFSdata3.BeginYear)

    CFSdata3.ncols=obw.ncols
    CFSdata3.nrows=obw.nrows
    !allocate(CFSdata3.MValue(CFSdata3.leadT,CFSdata3.ncols,CFSdata3.nrows))    


    CFSdata3.Crow=1 
    CFSdata3.Crow1=1 
    CFSdata3.nstep=0
    rData1.flagOs=flagOs
    rData1.ObFile=CFSData3.ObFile
    rData1.nparint=5

    rData1.BeginYear=CFSdata3.BeginYear
    rData1.EndYear=CFSdata3.EndYear
    rData1.leadT=CFSdata3.leadT
  
    CRE.TimeSum=CFSdata3.ndays*CFSdata3.leadT
    allocate(CRE.yc(CRE.TimeSum))
    allocate(CRE.yy(CRE.TimeSum))



    do i=CFSdata3.x1,CFSdata3.x2   !CFSdata3.ncols,CFSdata3.ncols*0.75,-1
      if (flagOs .eq. 2) then
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'mkdir ',trim(CFSdata3.OutFile),i             
      else
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'md ',trim(CFSdata3.OutFile),i   
      end if
      call system(outputfile)

      do j=CFSdata3.y1,CFSdata3.y2
        if (abs(obw.MValue(i,j)-obw.NODATA_value)>0.01) then
          print*,'x,y=         ',i,j  
          ! read 4 members of CFS
          CFSdata3.ix=i
          CFSdata3.iy=j
          CFSdata3.cday=0
          CFSdata3.Crow1=CFSdata3.Crow  
          do iyear=BeginYear,EndYear
            Print *,'x=',i,'    Year=',iyear  
            CFSdata3.iYear=iyear 
            do imonth=1,12
              CFSdata3.imonth=imonth
              do iday=1,CFSdata3.monthday(imonth,1)
                CFSdata3.cday=CFSdata3.cday+1
                !read CFS data               
                CFSdata3.iday=CFSdata3.monthday(imonth,iday+1)          
                CFSdata3.Crow=CFSdata3.Crow1           
                call ReadDCFSb1(CFSdata3) 	                           
              end do
            end do
          end do 
          ! read observed data      
          rData1.x =CFSdata3.ix
          rData1.y=CFSdata3.iy
          call readDataob(rData1)
          ! allocate(CFSdata3.PLValue(CFSdata3.ndays,5,CFSdata3.leadT))    
          l=0
          do iyear=CFSdata3.BeginYear,CFSdata3.EndYear
            do k=1,365,5  
              cday=JULDAY(1,1,iyear)+k-jday1
              if(((mod(iyear,4)==0.and.mod(iyear,100)/=0).or. mod(iyear,400)==0).and.k.gt.59) &  !leap year
                 cday=cday+1  
              l=l+1              			   
              CFSdata3.PLValue(l,5,:)=int(rData1.ob(cday,:)*100)
            end do
          end do

          ! compute correlation  allocate(CFSdata3.cor(5,5,CFSdata3.leadT+2))   
          ! CFSdata3.PLValue(CFSdata3.ndays,5,CFSdata3.leadT))   	
          CFSdata3.cor(:,:,CFSdata3.leadT+1)=0		  	     
          do k= 1,CFSdata3.leadT
            do ii=1,5
			  do l=1,CFSdata3.ndays
                CRE.yy(l)=CFSdata3.PLValue(l,ii,k)
              end do
              do jj=ii,5
			    do l=1,CFSdata3.ndays
                  CRE.yc(l)=CFSdata3.PLValue(l,jj,k)
                end do
                CRE.TimeSum=CFSdata3.ndays
                call CRelEff(CRE)
                CFSdata3.cor(ii,jj,k)=CRE.r
                CFSdata3.cor(ii,jj,CFSdata3.leadT+1)=CFSdata3.cor(ii,jj,CFSdata3.leadT+1)+abs(CRE.r)
              end do
            end do 
          end do
          CFSdata3.cor(:,:,CFSdata3.leadT+1)=CFSdata3.cor(:,:,CFSdata3.leadT+1)/CFSdata3.leadT
          do ii=1,5
            do l=1,CFSdata3.ndays*CFSdata3.leadT
              ll=mod(l,CFSdata3.leadT)
			  k=int(l/CFSdata3.leadT)+1
              if (ll.eq.0) then 
                ll=CFSdata3.leadT
                k=k-1
              end if
              CRE.yy(l)=CFSdata3.PLValue(k,ii,ll)
            end do
            do jj=ii,5
			  do l=1,CFSdata3.ndays*CFSdata3.leadT
                ll=mod(l,CFSdata3.leadT)
			    k=int(l/CFSdata3.leadT)+1
                if (ll.eq.0) then 
                  ll=CFSdata3.leadT
                  k=k-1
                end if
                CRE.yc(l)=CFSdata3.PLValue(k,jj,ll)                 
              end do
              CRE.TimeSum=CFSdata3.ndays*CFSdata3.leadT
              call CRelEff(CRE)
              CFSdata3.cor(ii,jj,CFSdata3.leadT+2)=CRE.r        
            end do
          end do 

          do k= 1,CFSdata3.leadT+2
            do ii=1,5
              do jj=ii+1,5
                CFSdata3.cor(jj,ii,k)=CFSdata3.cor(ii,jj,k)
              end do
            end do 
          end do

          !output the correlation            
          write(outputfile,'(2(a,i3.3),a)')trim(CFSdata3.OutFile),i,'\',j,'.txt'        
          if (flagOs .eq. 2)   call slashchange(outputfile)
      
          open(20,file=trim(outputfile))  
            do k=1,CFSdata3.leadT+2
              write (tempc,'(a,i3,a,i3,a,i3)') trim('cor'),k,char(9),5,char(9),5 
              call delblank(tempc)
              write (20,'(a)')trim(tempc)
	   
              do ii=1,5   
                write (tempc,'(365(a,f8.2))') (char(9),CFSdata3.cor(ii,jj,k),jj=1,5)   
                call delblank(tempc)
                write (20,'(a)')trim(tempc)             
              end do 
            end do  
          close(20)       
        else
          
        end if
      end do
    end do



 
    deallocate(obw.MValue)   
    deallocate(rData1.ob)   
    deallocate(CRE.yc)
    deallocate(CRE.yy)

    deallocate(CFSdata3.PValue)
    deallocate(CFSdata3.PLValue)
    deallocate(CFSdata3.cor)
    !deallocate(CFSdata2.MValue)




!***************************************************************************************************** 
!13. operational 
!
!
!***************************************************************************************************** 

  case(13)
    epp2.CoFile=ControlFile
    call ReadCOPpara(epp2,CFSData4)
    if (flagOs .eq. 2)   call slashchange(CFSData4.CFile) 
  
    call ReadCFScoord(CFSData4)
    allocate(CFSData4.MValue00(CFSData4.leadT*4,CFSData4.ncols,CFSData4.nrows))    
    allocate(CFSData4.MValue06(CFSData4.leadT*4,CFSData4.ncols,CFSData4.nrows))    
    allocate(CFSData4.MValue12(CFSData4.leadT*4,CFSData4.ncols,CFSData4.nrows))    
    allocate(CFSData4.MValue18(CFSData4.leadT*4,CFSData4.ncols,CFSData4.nrows))    
    allocate(CFSData4.PValue(CFSData4.leadT))    
    allocate(CFSData4.PIValue(CFSData4.leadT))    
    
    !13.1calculate the Inverse distance weighting
    CFSData4.flagOs=flagOs
    call  calweight(CFSdata4)
    obw=CFSData4.obw
    allocate(CFSData4.MValueN(CFSData4.leadT,obw.ncols,obw.nrows))    

    if (trim(CFSdata4.EName).eq.'tmp2m') then
      allocate(CFSData4.MValue1(CFSData4.leadT,CFSData4.ncols,CFSData4.nrows))    
      allocate(CFSData4.MValueN1(CFSData4.leadT,obw.ncols,obw.nrows))    
    end if

    iyear=CFSData4.iYear   
    imonth=CFSData4.imonth  
    iday=CFSData4.iday   
    !13.2. read CFS precipitation data   
    call ReadCFSop(CFSData4)
    !13.3. IDW interpolate CFS to 1/8 degree
    if (CFSData4.alive ) then
      CFSData4.iyear=iyear
      CFSData4.imonth=imonth
      CFSData4.iday=iday
      if (trim(CFSData4.EName).eq.trim('tmp2m')) then
        CFSData4.EName= trim('Tmax')
        call CFSinterpolate(CFSData4)!max temperature
        CFSData4.EName= trim('tmp2m')
      else
        call CFSinterpolate(CFSData4)  !precipitation prate
      end if
       
      if (trim(CFSData4.EName).eq.trim('tmp2m')) then
        CFSData4.EName= trim('Tmin')
        do l=1,CFSData4.leadT
          CFSData4.MValue(l,:,:)=CFSData4.MValue1(l,:,:) !min temperature
        end do     
        do l=1,CFSData4.leadT
          CFSData4.MValueN1(l,:,:)=CFSData4.MValueN(l,:,:)  !max temperature
        end do           
        call CFSinterpolate(CFSData4)       
        CFSData4.EName= trim('tmp2m')
      end if
    end if !file exist 

    deallocate(CFSData4.PValue)
    deallocate(CFSData4.PIValue)
    deallocate(CFSData4.MValue00)
    deallocate(CFSData4.MValue06)
    deallocate(CFSData4.MValue12)
    deallocate(CFSData4.MValue18)
    if (trim(CFSdata4.EName).eq.'tmp2m')  deallocate(CFSData4.MValue1)
    deallocate(CFSData4.w)


    !13.4. operational CFS to forecast

    !epp2.CoFile=ControlFile
    !call ReadCHpara(epp2)
    if (flagOs .eq. 2) then
      call slashchange(epp2.EpFile)
      call slashchange(epp2.ObFile)
      call slashchange(epp2.SiFile)
      call slashchange(epp2.SiFile1)
      call slashchange(epp2.PaFile)
      call slashchange(epp2.DoFile)
      call slashchange(epp2.ReFile)
    end if

    epp2.iyear=iYear   
    epp2.imonth=imonth  
    epp2.iday=iday   

    epp2.flagOs=flagOs

    !jday1=JULDAY(1,1,epp2.BeginYear)
    epp2.ndays=1
    epp2.ndays1=JULDAY(12,31,epp2.EndYearEns)-JULDAY(1,1,epp2.BeginYearEns)+1
    epp2.ncols=obw.ncols
    epp2.nrows=obw.nrows

    call NewspaceEP(epp2)
     
    rEPara1.Events=epp2.Events
    rEPara1.pthresh_obs=>epp2.pthresh_obs
    rEPara1.pthresh_fcst=>epp2.pthresh_fcst
    rEPara1.popobs=>epp2.popobs
    rEPara1.cavgobs=>epp2.cavgobs
    rEPara1.ccvobs=>epp2.ccvobs
    rEPara1.popfcst=>epp2.popfcst
    rEPara1.cavgfcst=>epp2.cavgfcst
    rEPara1.ccvfcst=>epp2.ccvfcst
    rEPara1.rho_est=>epp2.rho_est

    rData1.flagT=epp2.flagT
    rData1.flagOs=flagOs
    rData1.EpFile=epp2.EpFile
    rData1.ObFile=epp2.ObFile
    rData1.SiFile=epp2.SiFile
    rData1.BeginYear=epp2.BeginYear
    rData1.EndYear=epp2.EndYear
    rData1.BeginYearEns=epp2.BeginYearEns
    rData1.EndYearEns=epp2.EndYearEns
    rData1.leadT=epp2.leadT
    rData1.ndays=epp2.ndays
    rData1.ep=>epp2.ep
    rData1.ob=>epp2.ob
    rData1.si=>epp2.si
    rData1.nparint=1
	epp2.PensArea=-999.9
    epp2.EName=CFSdata4.EName
    write(Tdate,'(i4,2i2.2)') epp2.iYear,epp2.imonth,epp2.iday

    if (trim(CFSdata4.EName).eq.'prate') then
      write(outputfile,'(5a)')trim(epp2.ReFile),trim(epp2.EName),'\',Tdate,'.txt'     
      if (epp2.flagOs .eq. 2) then
        call slashchange(outputfile)
 	    write(inputfile,'(2a)')'rm -rf ',trim(outputfile)
	    call system(inputfile)     
      else
	    write(inputfile,'(2a)')'del ',trim(outputfile)
	    call system(inputfile)
      end if
    end if

    if (trim(CFSdata4.EName).eq.'tmp2m') then

      write(outputfile,'(4a)')trim(epp2.ReFile),'Tmax\',Tdate,'.txt' 
      if (epp2.flagOs .eq. 2) then
        call slashchange(outputfile)
 	    write(inputfile,'(2a)')'rm -rf ',trim(outputfile)
	    call system(inputfile)     
      else
	    write(inputfile,'(2a)')'del ',trim(outputfile)
	    call system(inputfile)
      end if

      write(outputfile,'(4a)')trim(epp2.ReFile),'Tmin\',Tdate,'.txt' 
      if (epp2.flagOs .eq. 2) then
        call slashchange(outputfile)
 	    write(inputfile,'(2a)')'rm -rf ',trim(outputfile)
	    call system(inputfile)     
      else
	    write(inputfile,'(2a)')'del ',trim(outputfile)
	    call system(inputfile)
      end if

    end if
    if (trim(CFSdata4.EName).ne.'prate') rData1.flagT=1
    do ii= epp2.x1,epp2.x2 ! 141,141  !464,464 !3,3  ! obw.ncols,1,-1    
      do jj= epp2.y1,epp2.y2 !55,55 !162,162 !186,186  ! obw.nrows,1,-1  
        if (abs(obw.MValue(ii,jj)-obw.NODATA_value)>0.01) then
          write(*,'(a2,i4,tr1,a2,i4)')'x=',ii,'y=',jj 
          epp2.x=ii
          epp2.y=jj
          !6.1 read observed data and forecasting data
          rData1.x=epp2.x 
          rData1.y=epp2.y
          rData1.nparint=epp2.nparint
          rData1.leadT=epp2.leadT
          if (epp2.flagT.eq.0) then
            call readDataEp(rData1)
          else     		      
            call readDataEpT(rData1)
          end if

          !read CFS
          rData1.si(1,:)=CFSdata4.MValueN(:,ii,jj)

          !13.5 read parameter of epp2
          if (trim(CFSdata4.EName).eq.'prate') then
            if (flagOs .eq. 2)   then
              write(rEPara1.PaFile,'(3a,i3.3,a1,i3.3,a4)')trim(epp2.PaFile),trim(epp2.EName),'/',rData1.x,'/',rData1.y,'.txt'        
            else
              write(rEPara1.PaFile,'(3a,i3.3,a1,i3.3,a4)')trim(epp2.PaFile),trim(epp2.EName),'\',rData1.x,'\',rData1.y,'.txt' 
            end if
          else
            if (flagOs .eq. 2)   then
              write(rEPara1.PaFile,'(2a,i3.3,a1,i3.3,a4)')trim(epp2.PaFile),'Tmin/',rData1.x,'/',rData1.y,'.txt'        
            else
              write(rEPara1.PaFile,'(2a,i3.3,a1,i3.3,a4)')trim(epp2.PaFile),'Tmin\',rData1.x,'\',rData1.y,'.txt' 
            end if           
            epp2.EName="Tmin"
          end if
          call ReadEpara(rEPara1)
          call Operation(epp2)
          call WriteOperation(epp2)
          !int(epp.Pens(i,j,epp.nems+1)*100),((char(9),int(epp.Pens(i,j,imem)*100)),imem=1,epp.nems) 
          !call WriteOperationb(epp2) PensArea(:,:,:,:) !Ensemble [Lon,Lat,lead time,members]
          epp2.PensArea(ii,jj,:,:)=epp2.Pens(1,:,:)
		  !write(*,*)epp.Pens(1,1,1) ! epp.PensArea(ii,jj,1,1)
        end if
      end do
    end do
    call WriteOperationb(epp2)

    ! caculate max   temperature 
    if (trim(CFSdata4.EName).ne.'prate') then
	  rData1.flagT=2
      do ii= epp2.x1,epp2.x2     
        do jj= epp2.y1,epp2.y2  
          if (abs(obw.MValue(ii,jj)-obw.NODATA_value)>0.01) then
            write(*,'(a2,i4,tr1,a2,i4)')'x=',ii,'y=',jj 
            epp2.x=ii
            epp2.y=jj
            !6.1 read observed data and forecasting data
            rData1.x=epp2.x 
            rData1.y=epp2.y
            rData1.nparint=epp2.nparint
            rData1.leadT=epp2.leadT     
            call readDataEpT(rData1)             
            !read CFS
            rData1.si(1,:)=CFSdata4.MValueN1(:,ii,jj)
            !13.5 read parameter of epp2
            if (flagOs .eq. 2)   then
              write(rEPara1.PaFile,'(2a,i3.3,a1,i3.3,a4)')trim(epp2.PaFile),'Tmax/',rData1.x,'/',rData1.y,'.txt'        
            else
              write(rEPara1.PaFile,'(2a,i3.3,a1,i3.3,a4)')trim(epp2.PaFile),'Tmax\',rData1.x,'\',rData1.y,'.txt' 
            end if           
            epp2.EName="Tmax"       
            call ReadEpara(rEPara1)
            call Operation(epp2)
            call WriteOperation(epp2)
            epp2.PensArea(ii,jj,:,:)=epp2.Pens(1,:,:)
          end if
        end do
      end do
      call WriteOperationb(epp2)
    end if
    call deletspaceEP(epp2)
    deallocate(obw.MValue)  
    deallocate(CFSData4.MValueN)
    if (trim(CFSdata4.EName).eq.'tmp2m')   deallocate(CFSData4.MValueN1)




!***************************************************************************************************** 
!14. get sub-domain MASK file 
!China province mask
!
!***************************************************************************************************** 

  case(14)
    epp2.CoFile=ControlFile
   
    !CFSData5.CFile E:\EPP4\DATA\CFSlonlat.txt
    !CFSData5.MName E:\EPP4\DATA\Provincexyz.txt
    !CFSData5.WFile E:\EPP4\DATA\CProvinceMask.txt

    if (flagOs .eq. 2)   call slashchange(CFSData5.CFile)  
    if (flagOs .eq. 2)   call slashchange(CFSData5.WFile)   
    call ReadCFScoord(CFSData5)
    if (flagOs .eq. 2)   call slashchange(CFSData5.MName)
    call  calMask(CFSdata5)



!***************************************************************************************************** 
!15.from historical observed data and forecast to get parameter of EPP model  for CFS in China
!                  
!***************************************************************************************************** 
  case(15)
    CanE4.CoFile=ControlFile
    call ReadCEparaCC(CanE4)
    if (flagOs .eq. 2) then
      call slashchange(CanE4.ObFile)
      call slashchange(CanE4.SiFile)
      call slashchange(CanE4.SiFile1)
      call slashchange(CanE4.ReFile)
      call slashchange(CanE4.DoFile)
    end if
    CanE4.flagOs=flagOs
    call  ReadASC(obw,CanE4.DoFile)  
    jday1=JULDAY(1,1,CanE4.BeginYear)
    CanE4.ndays=JULDAY(12,31,CanE4.EndYear)-jday1+1
    call NewspaceCE(CanE4)
     
    rData1.flagOs=flagOs
    rData1.ObFile=CanE4.ObFile
    rData1.SiFile=CanE4.SiFile
    rData1.BeginYear=CanE4.BeginYear
    rData1.EndYear=CanE4.EndYear
    rData1.leadT=CanE4.leadT
    rData1.ndays=CanE4.ndays
    rData1.ob=>CanE4.ob
    rData1.si=>CanE4.si
    rData1.nparint=1
    rData1.flagT=CanE4.flagT
    do ii=CanE4.x1,CanE4.x2  !3,3 !464,464 ! obw.ncols,1,-1    
      if (flagOs.eq.2) then
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'mkdir ',trim(CanE4.ReFile),ii             
      else
        write(outputfile,'(2a,i3.3,a1,i3.3,a4)')'md ',trim(CanE4.ReFile),ii   
      end if
      call system(outputfile)

      do jj=CanE4.y1,CanE4.y2 !186,186 !162,162 ! obw.nrows,1,-1
        if (abs(obw.MValue(ii,jj)-obw.NODATA_value)>0.01) then
          write(*,'(a2,i4,tr1,a2,i4)')'x=',ii,'y=',jj 
          CanE4.x=ii
          CanE4.y=jj
          !5.1 read observed data and forecasting data
          rData1.x=CanE4.x 
          rData1.y=CanE4.y   
          rData1.leadT=CanE4.leadT
          if (CanE4.flagT.eq.0) then
            call readDataobChina(rData1)
          else
            call readDataobT(rData1)        
          end if
          !read CFS
          rData1.si=-99.0
          rData1.nparint=CanE4.nparint
          rData1.SiFile=CanE4.SiFile1  
          rData1.leadT=CanE4.leadT
          call readDatafcst(rData1)
          !if (CanE4.flagT.ne.0)  rData1.si=rData1.si-200.0
          if (CanE4.flagCG.eq.1) then
            !read GFS
            rData1.nparint=0
            rData1.SiFile=CanE4.SiFile 
            rData1.leadT=CanE4.leadT1
            call readDatafcst(rData1)
          end if 

          !5.2 from historical observed data and forecast to get parameter of EPP model
          k=CanE4.nparint
          CanE4.nparint=5
          call EppPara(CanE4)
          CanE4.nparint=k
          !5.3 Write  paramater information of EPP
          call WriteCEpara(CanE4)

        end if
      end do
    end do

    call deletspaceCE(CanE4)




  end select
end
