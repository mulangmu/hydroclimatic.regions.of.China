implicit none 
  integer, parameter :: r4 = selected_real_kind(5)
  integer, parameter :: r8 = selected_real_kind(8)
  integer, parameter :: i8 = selected_int_kind(13)
  integer, parameter :: c50 = 50
  integer, parameter :: c100 =220
  integer, parameter :: c1000 =4000



  type ::  CFSData    !CFS data
    character(c100)  ObFile       !  observed data file  
    character(c100)  CFile       !Name of coordinate file   
    character(c100)  DFile       !Name of data file   
    character(c100)  WFile       !Name of weight file   
    character(c100)  OutFile     !Name of output data file   
    character(c100)  MName       !Name of merge file
    character(c100)  EName       !Name of Forcing 'tmp2m' 'prate'
    integer    ifile            !file no.
    logical    alive            !file exist
    integer    flagOs           !select Operating system, 1-windows  2-linux 
    integer    flagAF           !judge whether record exist, 1-Yes ,0-No
    integer    ntime            !time total
    integer    nstep            !time step /hours
    integer    ncols            !column  total
    integer    nrows            !row total     
	integer    ix              !point longitude   
	integer    iy              !point latitude
	integer    Crow            !current row in binary file
	integer    Crow1            !current row in binary file
    integer    x1,x2       !column  number
    integer    y1,y2       !row number 

    integer    leadT           !Lead time 
    integer    BeginYear
    integer    EndYear 
    integer    iYear 
    integer    imonth 
    integer    iday 
    integer    cday           !current day 1-365
    integer    ensm           !current member 1-4
    integer    ihour 
    integer    monthday(12,32)

    integer    ndays
	integer    nparint !interval of calculation in 365 days
	  
    real(r8) NODATA_value  
    real(r8) ,pointer :: x(:)    
    real(r8) ,pointer :: y(:)       
	type(Stack), pointer :: position(:,:) !(ix,iy) point position in raw CFS data 
    type(ASCFile)  ::  obw
    integer  ,pointer :: mask(:,:) !(384,190)!mask of each CFS grid
    real  ,pointer :: w(:,:,:) !(464,224,4)!weight of least 4 points (lon, lat,4) real   w(464,224,4)
    real  ,pointer :: MValue(:,:,:) !CFS  raw value (lead time, lon, lat)
    real  ,pointer :: MValue1(:,:,:) !CFS  raw value (lead time, lon, lat)
    real  ,pointer :: MValue00(:,:,:) !CFS  raw value (lead time, lon, lat)
    real  ,pointer :: MValue06(:,:,:) !CFS  raw value (lead time, lon, lat)
    real  ,pointer :: MValue12(:,:,:) !CFS  raw value (lead time, lon, lat)
    real  ,pointer :: MValue18(:,:,:) !CFS  raw value (lead time, lon, lat)
    real  ,pointer :: PValue(:)   !CFS new point value ( lead time)
    real  ,pointer :: MValueN(:,:,:)   !CFS new point value ( lead time, lon, lat)
    real  ,pointer :: MValueN1(:,:,:)   !CFS new point value ( lead time, lon, lat)

    integer  ,pointer :: PLValue(:,:,:)   !CFS file point value (time,lat,lead time) (time,4,lead time)
    character(10),pointer :: PIValue(:)   !CFS file point value ( lead time)
    real(r8),pointer ::   cor(:,:,:)       ! correlation coefficient between obs and fcst members (includes zeros)

  end type CFSData 



!!********************************!!******************************************************* 
!74 Read CFS File  
!!********************************!!******************************************************* 
subroutine ReadCFS(CFSdata1)
use PrecStru
implicit none
integer i,j,k ,ii,jj
type(CFSdata) ::  CFSdata1
type(Scaldat) ::  caldat1
character(C1000)  asFile,tempc 
integer JULDAY,jday1
logical alive
   

  jday1=JULDAY(1,1,CFSdata1.iYear) 
  caldat1.julian=jday1+CFSdata1.cday-1
  call caldat(caldat1)  
  CFSdata1.imonth=caldat1.mm
  CFSdata1.iday=caldat1.id
  alive=.false.
  k=CFSdata1.ihour
  do while (.not.alive.and.k.lt.20)
    write(asFile,'(3a,i4,i2.2,3a,i4,3i2.2,a)') trim(CFSdata1.DFile),trim(CFSdata1.EName),'\',&
      CFSdata1.iYear,CFSdata1.imonth,'\',trim(CFSdata1.EName),'.',CFSdata1.iYear,CFSdata1.imonth,&
      CFSdata1.iday,k,'.time.grb2'    
	if (CFSData1.flagOs.eq.2)    call slashchange(asFile)   
    inquire(file=asFile,exist=alive)
    if (.not.alive )k=k+6
  end do
  CFSdata1.alive=alive
 
  if (.not.alive.or.k.gt.18) then
     
    print *,trim(asFile),' is not exist.'
    return
  end if
  if (CFSData1.flagOs.eq.2) then 
    write(asFile,'(4a,i4,i2.2,3a,i4,3i2.2,5a,i2)')'./wgrib2  ',trim(CFSdata1.DFile),trim(CFSdata1.EName),'\',&
      CFSdata1.iYear,CFSdata1.imonth,'\',trim(CFSdata1.EName),'.',CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,&
	   k,'.time.grb2 -bin ',trim(CFSdata1.OutFile),trim(CFSdata1.EName),'\','fort.',CFSData1.ifile   
  else
    read(CFSdata1.DFile,'(a2,a)') tempc,asFile
    write(asFile,'(a,a1,3a,i4,i2.2,3a,i4,3i2.2,a,i2)')'wgrib2  \cygdrive\',tempc,trim(asFile),trim(CFSdata1.EName),'\',&
      CFSdata1.iYear,CFSdata1.imonth,'\',trim(CFSdata1.EName),'.',CFSdata1.iYear,CFSdata1.imonth,CFSdata1.iday,&
      k,'.time.grb2 -bin fort.',CFSData1.ifile   
  end if	       
   
  call slashchange(asFile)
  call system(asFile)
  if (CFSData1.flagOs.eq.2) then
    write(asFile,'(4a,i2)')trim(CFSdata1.OutFile),trim(CFSdata1.EName),'\','fort.',CFSData1.ifile        
    call slashchange(asFile)
  else
    write(asFile,'(a,i2)') 'fort.',CFSData1.ifile        

  end if
  open(unit=20,file=asFile,form='unformatted',access='sequential')

  j=CFSdata1.leadT*4  
  CFSdata1.MValue=>CFSdata1.MValue00
  do i=1,j     
    write(*,*)"i=",i
    read(20,end=20) CFSdata1.MValue(i,:,:) 
  end do
  20 close(20)
  ! 6 hours data to daily data
  if (trim(CFSdata1.EName).eq.'prate')  then
    do k=2,4-CFSdata1.ihour/6    
      CFSdata1.MValue(1,:,:) =CFSdata1.MValue(1,:,:) +CFSdata1.MValue(k,:,:) 
    end do
    CFSdata1.MValue(1,:,:) =CFSdata1.MValue(1,:,:) *3600*6
    do i=2,CFSdata1.leadT  
      write(*,*)"leadT=",i
      CFSdata1.MValue(i,:,:) =0
      do k=1,4
        j=(i-1)*4+k-CFSdata1.ihour/6
        CFSdata1.MValue(i,:,:) =CFSdata1.MValue(i,:,:) +CFSdata1.MValue(j,:,:) 
      end do
      CFSdata1.MValue(i,:,:) =CFSdata1.MValue(i,:,:) *3600*6
    end do 
  end if


  if (trim(CFSdata1.EName).eq.'tmp2m') then

    do k=2,4    
      CFSdata1.MValue1(1,:,:) = CFSdata1.MValue(1,:,:)  
      do ii=1,CFSData1.ncols   
        do jj=1,CFSData1.nrows
          if (CFSdata1.MValue(1,ii,jj).lt. CFSdata1.MValue(k,ii,jj)) &
            CFSdata1.MValue(1,ii,jj) = CFSdata1.MValue(k,ii,jj) 
          if (CFSdata1.MValue1(1,ii,jj).gt. CFSdata1.MValue(k,ii,jj)) &
            CFSdata1.MValue1(1,ii,jj) = CFSdata1.MValue(k,ii,jj) 
        end do
      end do
    end do
    
    do i=2,CFSdata1.leadT  
      write(*,*)"leadT=",i
      CFSdata1.MValue(i,:,:) =-999
      CFSdata1.MValue1(i,:,:) =999
      do k=1,4
        j=(i-1)*4+k
        do ii=1,CFSData1.ncols   
          do jj=1,CFSData1.nrows
            if (CFSdata1.MValue(i,ii,jj).lt. CFSdata1.MValue(j,ii,jj)) &
              CFSdata1.MValue(i,ii,jj) = CFSdata1.MValue(j,ii,jj) 
            if (CFSdata1.MValue1(i,ii,jj).gt. CFSdata1.MValue(j,ii,jj)) &
              CFSdata1.MValue1(i,ii,jj) = CFSdata1.MValue(j,ii,jj) 
            
          end do
        end do
      end do       
    end do 

    
  end if
end subroutine ReadCFS