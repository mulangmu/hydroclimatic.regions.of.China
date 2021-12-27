load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 


begin


;check to make sure netcdf file doesn't exist; if it does, remove it
;then open netcdf file to write to
system("rm basins17check.nc")
out = addfile("basins17check.nc", "c")

;create new variables to hold ascii data
;we use two because we have to do some re-arranging of the data

basin=asciiread ("chinabasin17Mask.txt", (/190,384/), "float")
basins    = basin
basins@_FillValue = 9.96921e+36

do i = 0, 189
   do j = 0, 383
      if (basins(i,j) .eq. 1 ) then
         basins(i,j)=3
      else if (basins(i,j) .eq. 2 ) then
         basins(i,j)=1
      else if (basins(i,j) .eq. 3 ) then
         basins(i,j)=2
      else if (basins(i,j) .eq. 4 ) then
         basins(i,j)=4
      else if (basins(i,j) .eq. 5 ) then
         basins(i,j)=17
      else if (basins(i,j) .eq. 6 ) then
         basins(i,j)=6
      else if (basins(i,j) .eq. 7 ) then
         basins(i,j)=12
      else if (basins(i,j) .eq. 8 ) then
         basins(i,j)=11
      else if (basins(i,j) .eq. 9 ) then
         basins(i,j)=10
      else if (basins(i,j) .eq. 10 ) then
         basins(i,j)=15
      else if (basins(i,j) .eq. 11 ) then
         basins(i,j)=13
      else if (basins(i,j) .eq. 12 ) then
         basins(i,j)=9
      else if (basins(i,j) .eq. 13 ) then
         basins(i,j)=14
      else if (basins(i,j) .eq. 14 ) then
         basins(i,j)=5
      else if (basins(i,j) .eq. 15 ) then
         basins(i,j)=8
      else if (basins(i,j) .eq. 16 ) then
         basins(i,j)=7
      else if (basins(i,j) .eq. 17 ) then
         basins(i,j)=16
      else
         basins(i,j)=basins@_FillValue
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
   end do
end do

;if you were going to plot, you should assign
;named dimensions and coordinate variables
basins!0="lat"
basins!1="lon"

lon=fspan(0,359.062,384)
lat=fspan(89.277,-89.277,190)
basins&lat=lat
basins&lon=lon

lon@units="degrees_east"
lat@units="degrees_north"
lat@La1="89.277"
lat@La2="-89.277"
lat@Di="0.9374987"
lat@grid_type="Gaussian Latitude/Longitude"
lat@long_name="latitude"
lon@Lo1="0"
lon@Lo2="359.062"
lon@Di="0.9374987"
lon@grid_type="Gaussian Latitude/Longitude"
lon@long_name="longitude"

;copy coordinate variable info from variable we created above from 
;netcdf file to our new variable so that we can read the output in grads
;this copy_VarCoords function is in contributed.edited.ncl and saves a lot of time
;make sure the dims match between the two variables
;copy_VarCoords(ts,basins)

;output vars to .nc file - you must also put in lat/lon information
out->basins = basins
out->lat = lat
out->lon = lon

end

