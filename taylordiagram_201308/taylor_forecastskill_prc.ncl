;**********************************
; taylor_4.ncl
;langyang 20130420
;**********************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"     
load "./taylor_diagram_prc.ncl"
;load "./taylor_metrics_table.ncl"
;load "./taylor_diagram_enoble.ncl"

;**********************************
begin
;**********************************
; Assume the following have already been computed:
;   _ratio are the ratio:  Case_Variance/Reference_Variance
;   _cc    are the cross correlation coef of Case to Reference
; In this example, these are derived for annual mean climatologies.
;**********************************

; Cases [Model]
  ;case      = (/ "ASO(day 200-300)", "JFM(day 350-080)" /) 
  case =(/"JFM","FMA","MAM","AMJ","MJJ","JJA","JAS","ASO","SON","OND","NDJ","DJF"/)
  nCase     = dimsizes(case )                 ; # of Cases [Cases]

; variables compared
  var       = (/ "Basin1","Basin2" ,"Basin3","Basin4","Basin5","Basin6", "Basin7", "Basin8" \ 
               , "Basin9" ,"Basin10","Basin11","Basin12","Basin13","Basin14","Basin15" , "Basin16", "Basin17" /) 
  nVar      = dimsizes(var)                   ; # of Variables

; more info to be added [all are bogus]
  source     =(/"JFM:001-090","FMA:032-120","MAM:060-151","AMJ:091-181","MJJ:121-212","JJA:152-243" \
                  ,"JAS:182-273","ASO:213-304","SON:244-334","OND:274-365","NDJ:305-031","DJF:335-059"/)
               
 ratio      = new ((/nVar, nCase/),float)  
 cor         = new ((/nVar, nCase/),float)   


 cor_data=readAsciiTable("/mnt/LVM/langyang/EA/taylordiagram_201308/basins17_prec_cor.txt",12,"float",1)
 ratio_data=readAsciiTable("/mnt/LVM/langyang/EA/taylordiagram_201308/basins17_prec_stdratios.txt",12,"float",1)

; "Basin 17"   
 do i=0, 16,1                     
    ratio(i,:)   = ratio_data(i,:)
    cor(i,:)     = cor_data(i,:)
 end do

;**********************************
; create plot
;**********************************

  wks   = gsn_open_wks("eps","taylor_prc")
  plot = new(18,graphic)
  
  ty_opt   = True                           ; taylor diagram with options
        
  ty_opt@Markers       = (/4,6,8,4,6,8,4,6,8,4,6,8/)         ; make all solid fill
  ty_opt@gsMarkerSizeF   = 0.008   ; marker size
  ty_opt@gsMarkerThicknessF = 1.8
  ty_opt@txFontHeightF    = 0.025
  ty_opt@txFontThicknessF = 5
   
  ty_opt@Colors        = (/ "blue" , "blue", "green","green" , "green", "red", "red" , "red", "brown","brown","brown", "blue" /)
 ;                          , "torquoise", "brown", "yellow","blue" , "red", "green", "cyan"/)
         
  ;ty_opt@varLabels     =(/"JFM","FMA","MAM","AMJ","MJJ","JJA","JAS","ASO","SON","OND","NDJ","DJF"/)
  ;ty_opt@caseLabels    =(/"JFM","FMA","MAM","AMJ","MJJ","JJA","JAS","ASO","SON","OND","NDJ","DJF"/)

  ;ty_opt@varLabelsYloc = 0.65                ; Move location of variable labels [default 0.45]
  ty_opt@caseLabelsFontHeightF = 0.05      ; make slight larger   [default=0.12 ]
  ;ty_opt@varLabelsFontHeightF  = 0         ; make slight smaller  [default=0.013]
 
  ;ty_opt@tiMainString  = "Two seasonal forecast skill of Precipitation for 17 Basins"           ; title
  ;ty_opt@tiMainFontHeightF = 0.014 ; tiMainString size

  ty_opt@stnRad        = (/ 0.5, 1.5 /)     ; additional standard radii
  ty_opt@ccRays        = (/ 0.1,0.2,0.3,0.4,0.5,0.6 ,0.7,0.8,0.9,0.95,0.99/)     ; correllation rays
  ty_opt@centerDiffRMS = True                ; RMS 'circles'

  
  ;plot  = taylor_diagram(wks,ratio,cc,ty_opt)
  ty_opt@taylorDraw  = False                          ; don't draw
  ty_opt@taylorFrame = False                          ; don't advance frame
  ; Generate one plot for each basin.
  ; *******************************************************************
  plots      = 17
  do i=0,plots-1
     ;res@tiMainString    = plot_title(i)              ; title
     ratio_plot=onedtond(ratio(i,:),(/12,1/))
     cor_plot=onedtond(cor(i,:),(/12,1/))
     plot(i)  = taylor_diagram(wks,ratio_plot,cor_plot,ty_opt)
  end do
     ;ty_opt1   = True
     ;ty_opt1@taylorDraw  = False                          ; don't draw
     ;ty_opt1@taylorFrame = False                          ; don't advance frame
     ;ty_opt1@Markers       = (/4,6,8,4,6,8,4,6,8,4,6,8/)         ; make all solid fill
     ;ty_opt1@gsMarkerSizeF   = 0.005   ; marker size
     ;ty_opt1@gsMarkerThicknessF = 0.8
     ;ty_opt1@Colors        = (/ "blue" , "blue", "green","green" , "green", "red", "red" , "red", "brown","brown","brown", "blue" /)

     ;plot(17) = taylor_diagram(wks,ratio_plot,cor_plot,ty_opt1)
     ;ty_opt1   = False
     ;plot(17) = gsn_blank_plot(wks,False)
     
    ;plot(17) = gsn_csm_xy(wks, x,y,xyres)
    ;newplot= gsn_attach_plots(plot(16),(/plot(16)/),True,True)
  ; *******************************************************************
  ; Create demo panels
  ; *******************************************************************
  resP                                  = True
  resP@gsnFrame         = False                  ; don't advance panel plot
  ;resP@gsnMaximize                      = False ;True
  ;resP@txString            = "Seasonal forecast skill of Precipitation for 17 Basins"
  ;resP@gsnPanelLabelBar    = True                ; add common colorbar
  resP@gsnPanelFigureStrings            = (/"1)","2)","3)","4)","5)","6)","7)","8)","9)","10)","11)","12)","13)","14)","15)","16)","17)"/)
  resP@gsnPanelFigureStringsFontHeightF = 0.01 
  resP@amJust                           = "TopLeft"


   Markers = (/4,6,8,4,6,8,4,6,8,4,6,8/)         
   Colors  = (/ "blue" , "blue", "green","green" , "green", "red", "red" , "red", "brown","brown","brown", "blue" /)
   gsMarkerSizeF      = ty_opt@gsMarkerSizeF  
   ;gsMarkerThicknessF = 2
   caseLabelsFontHeightF = 0.20
    
   
      lgres                    = True
      lgres@lgMarkerColors     = Colors        ; colors of markers
      lgres@lgMarkerIndexes    = Markers       ; Markers 
      lgres@lgMarkerSizeF      = gsMarkerSizeF ; Marker size
      lgres@lgMarkerThicknessF = 1.8
      lgres@lgItemType         = "Markers"     ; draw markers only
      lgres@lgTitleString      =  "CFSv2 Precipitation forecasts"
      lgres@lgTitleFontHeightF       = 0.03
      lgres@lgLabelFontHeightF = caseLabelsFontHeightF  ; font height of legend case labels
      lgres@vpWidthF       = 0.15           ; width of legend (NDC)
      lgres@vpHeightF      = 0.06*nCase   ; height of legend (NDC)
      lgres@caseLabels     =(/"JFM","FMA","MAM","AMJ","MJJ","JJA","JAS","ASO","SON","OND","NDJ","DJF"/)

      lgres@lgPerimOn          = False         ; turn off perimeter
      nModel                   = dimsizes( lgres@caseLabels)
      lbid = gsn_create_legend(wks,nModel,lgres@caseLabels,lgres)
	 
      amres = True
      amres@amParallelPosF     =  1;0.35           
      amres@amOrthogonalPosF   = 0.08;-0.35             
      annoid1 = gsn_add_annotation(plot(16),lbid,amres)	; add legend to plot


  gsn_panel(wks,plot,(/3,6/),resP)
  txres               = True
  txres@txFontHeightF = 0.01
  ;gsn_text_ndc(wks,"Seasonal predictive performance of CFSv2 Precipitation forecasts over 17 large basins",0.5,0.25,txres)
  frame(wks)



end 
