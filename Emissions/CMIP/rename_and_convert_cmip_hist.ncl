;
; read regridded CEDS/CMIP6 and biomass burning fire  emissions files
; convert from kg/m2/s to molecules/cm2/s
; sum all anthro sectors to one variable
;  (requires gnu compiler, not intel, on cheyenne/geyser)
;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; USER CHANGES start
; change <user> to your user name
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 resolution = "ne30np4"   ; adjust to your grid
 dir_inp = "/glade/scratch/<user>/emissions/historical_"+resolution+"/"
 dir_out = "/glade/scratch/<user>/cmip6inputs/historical_ne30np4/"
; resolution= "ne30pg3"
; resolution = "C96"
; resolution = "ne0CONUSne30x8"

 creation_date = systemfunc ("date +%Y%m%d")

 avog = 6.022e23
 earth_rad  = 6.371e8 ; in cm
 earth_area = 4. * acos(-1.) * earth_rad * earth_rad
 PI  = get_pi("float" )
 sector = (/"Anthropogenic","Biomass Burning","Aviation"/)
 sector1 = (/"anthro","bb","air"/)

; process anthro, biomass burning and aircraft emissions
;
  system("mkdir -p " + dir_out)


do nspc=0,33  ; nspc gt 30 require other emissions to be processed 
;
; loop over sectors
;
  do sec=0,1   ; anthro, bb, air

    if ( nspc.eq.0 ) then
      spc_name = "BC"
      spc_name2 = "bc_a4"
      sect_anthro = 1
      sect_bb= 1
      sect_air= 1
      mw = 12.
      sf = 1.0
    end if
    if ( nspc.eq.1 ) then
      spc_name = "CO"
      mw = 28.
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 1
    end if
    if ( nspc.eq.2 ) then
      spc_name = "NH3"
      mw = 17.
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
    end if
    if ( nspc.eq.3 ) then
      spc_name = "NOx"
      if sec.eq.1 then
       mw = 30.
      else
       mw = 46.
      end if 
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 1
    end if
    if ( nspc.eq.4 ) then
      spc_name = "OC"
      spc_name2 = "pom_a4"
      mw = 12.
      sf = 1.
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
    end if
    if ( nspc.eq.5 ) then
      spc_name = "SO2"
      spc_name2 = "so4_a1"
      mw = 64.
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 1
    end if
;   if ( nspc.eq.6 ) then
;     spc_name = "CO2"
;     mw = 44.
;   end if
;   if ( nspc.eq.6 ) then
;     spc_name = "CH4"
;     mw = 16.
;   end if
;   if ( nspc.ge.7 ) then
;     cnt_voc = nspc - 6

    if ( nspc.eq.6 ) then
     spc_name = "C2H6"
     spec_ceds = (/"VOC02-ethane"/)
     spec_bb= (/"NMVOC_C2H6"/)
     mw = 30.
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
    end if
    if ( nspc.eq.7 ) then
     spc_name = "C3H8"
     spec_ceds = (/"VOC03-propane"/)
     spec_bb= (/"NMVOC_C3H8"/)
     mw = 44.
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.8 ) then
     spc_name = "BIGALK"
      spec_ceds = (/"VOC04-butanes","VOC05-pentanes","VOC06-hexanes-pl","VOC18-esters","VOC19-ethers"/)
      spec_bb= (/"NMVOC_Higher_Alkenes"/)
      mw = 72.
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.9 ) then
     spc_name = "C2H4"
      spec_ceds = (/"VOC07-ethene"/)
      spec_bb= (/"NMVOC_C2H4"/)
      mw = 28.
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.10 ) then
     spc_name = "C3H6"
      spec_ceds  = (/"VOC08-propene"/)
      spec_bb= (/"NMVOC_C3H6"/)
      mw = 42.
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.11 ) then
     spc_name = "C2H2"
      spec_ceds  = (/"VOC09-ethyne"/)
      spec_bb= (/"NMVOC_C2H2"/)
      mw = 26.
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.12) then
     spc_name = "BIGENE"
      spec_ceds  = (/"VOC12-other-alke"/)
      spec_bb= (/"NMVOC_Higher_Alkanes"/)
      mw = 56.
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.13 ) then
     spc_name = "BENZENE"
      spec_ceds = (/"VOC13-benzene"/)
      spec_bb= (/"NMVOC_C6H6"/)
      mw = 78.
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.14 ) then
     spc_name = "TOLUENE"
      spec_ceds  = (/"VOC14-toluene"/)
      spec_bb= (/"NMVOC_C7H8"/)
      mw = 92.
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.15 ) then
     spc_name = "XYLENES"
      spec_ceds  = (/"VOC15-xylene","VOC16-trimethylb","VOC17-other-arom"/)
      spec_bb= (/"NMVOC_C8H10"/)
      mw = 106.
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.16 ) then
     spc_name = "CH2O"
      spec_ceds  = (/"VOC21-methanal"/)
      spec_bb= (/"NMVOC_CH2O"/)
      mw = 30.
      sf = 1.0
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.17 ) then
     spc_name = "CH3CHO"
      spec_ceds  = (/"VOC22-other-alka"/)
      spec_bb= (/"NMVOC_C2H4O"/)
      mw = 44.
      sf = 1.
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.18 ) then
     spc_name = "CH3OH"
      spec_ceds  = (/"VOC01-alcohols"/)
      spec_bb= (/"NMVOC_CH3OH"/)
      sf = 0.15
      mw = 32.
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.19 ) then
     spc_name = "C2H5OH"
      spec_ceds  = (/"VOC01-alcohols"/)
      spec_bb= (/"NMVOC_C2H5OH"/)
      sf = 0.85
      mw = 46.
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.20 ) then
     spc_name = "CH3COCH3"
      spec_ceds  = (/"VOC23-ketones"/)
      spec_bb= (/"NMVOC_C3H6O"/)
      sf = 0.2
      mw = 58.
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.21 ) then
     spc_name = "MEK"
      spec_ceds  = (/"VOC23-ketones"/)
      spec_bb= (/"NMVOC_MEK"/)
      sf = 0.8
      mw = 72.
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.22 ) then
     spc_name = "HCOOH"
      spec_ceds  = (/"VOC24-acids"/)
      spec_bb= (/"NMVOC_HCOOH"/)
      sf = 0.5
      mw = 46.
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.23 ) then
     spc_name = "CH3COOH"
      spec_ceds  = (/"VOC24-acids"/)
      spec_bb= (/"NMVOC_CH3COOH"/)
      sf = 0.5
      mw = 60.
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.24 ) then
     spc_name = "ISOP"
      spec_ceds= (/"NMVOC_C5H8"/)
      spec_bb= (/"NMVOC_C5H8"/)
      sf = 1.
      mw = 68.
      sect_anthro = 0
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.25 ) then
     spc_name = "MTERP"
      spec_ceds= (/"NMVOC_C10H16"/)
      spec_bb= (/"NMVOC_C10H16"/)
      sf = 1.
      mw = 136.
      sect_anthro = 0
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.26 ) then
      spc_name = "DMS"
      spec_ceds = "DMS"
      spec_bb= (/"NMVOC_C2H6S"/)
      sf = 1.
      mw = 62.
      sect_anthro = 0
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.27 ) then
     spc_name = "CH3COCHO"
      spec_ceds= (/"NMVOC_CH3COCHO"/)
      spec_bb= (/"NMVOC_CH3COCHO"/)
      sf = 1.0
      mw = 60.
      sect_anthro = 0
      sect_bb= 1
      sect_air= 0
   end if
    if ( nspc.eq.28 ) then
     spc_name = "GLYALD"
      spec_ceds= (/"NMVOC_HOCH2CHO"/)
      spec_bb= (/"NMVOC_HOCH2CHO"/)
      sf = 1.
      mw = 60.
      sect_anthro = 0
      sect_bb= 1
      sect_air= 0
    end if
    if ( nspc.eq.29 ) then
     spc_name = "HCN"
      spec_ceds  = (/"CO"/)
      spec_bb= (/"NMVOC_HCN"/)
      sf = 0.003  ; 0.5*0.006 assuming half is biofuel
      mw = 27.
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
    end if
    if ( nspc.eq.30 ) then
     spc_name = "CH3CN"
      spec_ceds  = (/"CO"/)
      spec_bb= (/"NMVOC_HCN"/)
      mw = 41.
      sf = 0.002
      sf_bb = 0.33
      sect_anthro = 1
      sect_bb= 1
      sect_air= 0
    end if
    if ( nspc.eq.31) then
     spc_name = "IVOC"
     sf = 1.
     mw = 184.
      sect_anthro = 0
      sect_bb= 0
      sect_air= 0
    end if
    if ( nspc.eq.32) then
     spc_name = "SVOC"
     sf = 1.
     mw = 310.
      sect_anthro = 0
      sect_bb= 0
      sect_air= 0
    end if
    if ( nspc.eq.33 ) then
     spc_name = "SOAGx1.5"
      spec_ceds  = (/"BIGALK","BIGENE","TOLUENE","BENZENE","XYLENES","ISOP","MTERP"/)
      yields = (/0.05, 0.05, 0.15, 0.15, 0.15, 0.04, 0.25/)
      mw = 12.
      sf = 1.0
      sect_anthro = 0
      sect_bb= 0
      sect_air= 0
   end if

   run_sec = 1
   if sec.eq.0 .and. sect_anthro.eq.0 then
      run_sec = 0
   end if
   if sec.eq.1 .and. sect_bb.eq.0 then
      run_sec = 0
   end if
   if sec.eq.2 .and. sect_air.eq.0 then
      run_sec = 0
   end if

  if run_sec.eq.1 then

  nspec_ceds = 1
  if sect_anthro.eq.1 then
    if ( nspc.le.5 ) then
      nspec_ceds = 1
    else
      nspec_ceds = dimsizes(spec_ceds)
      str1 = str_join(spec_ceds,"+")
      if sf .lt. 1 then 
         str2 = sprintf("%4.2e", sf)
      else 
        str2= " "
      end if
      hist_spec = "CEDS species: "+str2+str1
      print(hist_spec)
    end if
  end if
  do ced=0,nspec_ceds-1
   if sec.eq.0 .and. sect_anthro.eq.1 then
    if ( nspc.le.5 ) then
      files = systemfunc("ls " + dir_inp + spc_name + "-em-anthro*")
      ;files = systemfunc("ls " + dir_inp + spc_name + "-em-anthro*200001-201412*")
    else
     files = systemfunc("ls " + dir_inp + spec_ceds(ced)+ "*-em-speciated-VOC*")
     if sec.eq.0  .and. spc_name .eq. "HCN" then
      delete(files)
      files = systemfunc("ls " + dir_inp + spec_ceds(ced) + "-em-anthro*")
     end if
     if sec.eq.0  .and. spc_name .eq. "CH3CN" then
      delete(files)
      files = systemfunc("ls " + dir_inp + spec_ceds(ced) + "-em-anthro*")
     end if
    end if
   end if
   if sec.eq.1 .and. sect_bb.eq.1 then
    if ( nspc.le.5 ) then
       spec_ceds = spc_name
    else
     delete(spec_ceds)
     spec_ceds=spec_bb
    end if
    files = systemfunc("ls " + dir_inp + spec_ceds+ "-em-biomassburning*")
    hist_spec = spec_ceds
   end if ;bb
   if sec.eq.2 .and. sect_air.eq.1 then
    ; files = systemfunc("ls " + dir_inp + spc_name+ "-em-AIR-anthro*200001-201412.nc")
     files = systemfunc("ls " + dir_inp + spc_name+ "-em-AIR-anthro*")
   end if

  
; should be just 1 file per species
    nfil=0
;
      nfils = dimsizes(files)
      print(nfils)
      file_inp1 = files(nfil)
      file_inp2 = files(nfils-1)
      in = addfile(file_inp1,"r")
      fin=addfiles(files,"r")
;
      print("Processing " + files)
;
; find emissions variable name
;
      vnames = getfilevarnames(in)
      nnames = dimsizes(vnames)
     
     if sec.eq.1 then 
       rk=2
     else
       rk=3
     end if
      do n=0,nnames-1
        print(vnames(n))
        var = in->$vnames(n)$
        dimv = dimsizes(var)
        rank = dimsizes(dimv)
        if (rank.eq.rk) then 
           var_name = vnames(n)
        end if
        delete(var)
        delete(dimv)
        delete(rank)
      end do

;
; read coordinates
;
      lat= in->lat
      lon= in->lon
      time= in->time
      area = in->area
      mw_s = mw
     if ( nspc.gt.5 .and. sec.eq.0) then
      if spec_ceds(ced).ne."CO" then
        mw_spec = in@molecular_weight
      else
        mw_spec = 28
      end if
      if mw_spec.ne.mw_s then 
         print(mw_s)
         print(mw_spec)
         mw_s=mw_spec
      end if
     end if
     if spc_name .eq. "CH3CN" .and. sec.eq.1 then
      mw_s = 27.
     end if
     if spc_name .eq. "BIGALK" .and. sec.eq.1 then
      mw_s = 12.
     end if
     if spc_name .eq. "BIGENE" .and. sec.eq.1 then
      mw_s = 12.
     end if
        
      if sec.eq.2 .and. sect_air.eq.1 then
       level = in->level
      end if
;
; read emissions - all sectors in one variable
;
      wk = addfiles_GetVar(fin,files,var_name)
  if(any(ismissing(wk))) then
      print("missing data")
      wk1D      = ndtooned(wk)
      ibad = ind(ismissing(wk1D))
      wk1D(ibad)=0.
      wk = onedtond(wk1D, dimsizes(wk))
  else
      print("no missing data")
  end if
printVarSummary(wk)
print(max(wk))
print(min(wk))

  
      time1 = addfiles_GetVar(fin,files,"time")
      nwk = dimsizes(wk)
      ntim = nwk(0)
      ncol = nwk(rk-1)
;
; sum over all sectors if anthro
;
printVarSummary(wk)

   if sec.eq.0 then
      ;emiss_inp = dim_sum_n_Wrap(wk,1)
      emiss_inp = dim_sum_n_Wrap(wk(:,:,:),1)
printVarSummary(emiss_inp)
   else
      emiss_inp = wk
      sf = 1.
   end if
   if sec.eq.1  .and. spc_name .eq. "CH3CN" then
      sf = 0.33
   end if
;  emiss_inp = wk(:,2,:)

;
; convert from kg/m2/s to molec/cm2/s
;
;                            kg to g         m2 to cm2
;                               |                |


     emiss_inp = emiss_inp * 1.e3 * avog/mw_s * 1.e-4
     if ced.eq.0 then
      emiss_out = emiss_inp
      emiss_out = sf * emiss_inp
     else
      emiss_out = emiss_out + sf * emiss_inp
     end if


if sec.eq.0 then
; compute global average (in molec/cm2/s)
; annual average in Tg/year
;
     ;area_emis = conform(emiss_inp,area,1)
     area_emis = conform(emiss_out,area,1)
     emiss_out_gl = emiss_out*area_emis
     emiss_out_gl1 = dim_sum_n_Wrap(emiss_out_gl,1) / dim_sum_n_Wrap(area_emis,1)

      dims = dimsizes(emiss_out_gl1)
      nyr = dims(0)/12
      emiss_out_gl_ann = emiss_out_gl1(0:nyr-1)
      do n=0,nyr-1
        n1 = n * 12
        n2 = n1 + 11
        emiss_out_gl_ann(n) = dim_avg_Wrap(emiss_out_gl1(n1:n2))
      end do
      emiss_out_gl_ann = emiss_out_gl_ann * earth_area * 86400. * 365 * mw * 1.e-12/avog
;
     delete(emiss_out_gl)
     delete(emiss_out_gl1)
     delete(emiss_out_gl_ann)
     delete(area_emis)
end if

   delete(files)
  end do ;ced
; altidude considereation for aircarft emissions
    if sec.eq.2 .and. sect_air.eq.1 then
     altitude_int_air = (/0, 0.61, 1.22, 1.83, 2.44, 3.05, 3.66, 4.27, 4.88, 5.49, 6.1, \
     6.71, 7.32, 7.93, 8.54, 9.150001, 9.76, 10.37, 10.98, 11.59, 12.2, 12.81, \
     13.42, 14.03, 14.64, 15.25/)
     altitude_air = (/0.305, 0.915, 1.525, 2.135, 2.745, 3.355, 3.965, 4.575, 5.185, \
     5.795, 6.405, 7.015, 7.625, 8.235001, 8.845, 9.455001, 10.065, 10.675, \
     11.285, 11.895, 12.505, 13.115, 13.725, 14.335, 14.945/)
     altitude_air!0="altitude"
     altitude_int_air!0="altitude_int"
     altitude_air@units="km"
     altitude_air@long_name="Altitude"
     altitude_int_air@units="km"
     altitude_int_air@long_name="Altitude interfaces"
     nlev = dimsizes(altitude_air)
    do ilev=0,nlev-1
     dz = (altitude_int_air(ilev+1) - altitude_int_air(ilev))*1.e5   ;km -> cm
     emiss_out(:,ilev,:) = emiss_out(:,ilev,:)*dz
     emiss_out!1="altitude"
    end do
    end if
printVarSummary(emiss_out)
print(max(emiss_out))
print(min(emiss_out))

      emiss_out@molecular_weight=mw
     if sec.eq.2 then
      emiss_out@units="molecules/cm3/s"
     else
      emiss_out@units="molecules/cm2/s"
     end if
      if nspc.gt.5 then    
        emiss_out@long_name= spc_name+" "+sector1(sec)+" emissions"
        emiss_out@history= hist_spec 
      end if
printVarSummary(emiss_out)

  if spc_name .eq. "SO2"  then
      emiss_out1= emiss_out
      emis_num= emiss_out
      if sec.ne.2 then
       emiss_out= emiss_out*0.975
      end if
      emiss_out1 = emiss_out1*0.025
      diam = 0.134e-6
      rho = 1770.
      mw2 = 115.
      mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
      emis_num= emiss_out1 *mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
      emis_num@units="(particles/cm2/s)(molecules/mole)(g/kg)"
      emis_num@long_name = "number of SO4: 0.025* "+spc_name+" "+sector(sec)+" Emissions"
      emiss_out1@long_name = "SO4: 0.025* "+spc_name+" "+sector(sec)+" Emissions"
      if sec.ne.2 then
       emiss_out@long_name =  "0.975* "+spc_name+" "+sector1(sec)+" emissions"
      else
       emiss_out@long_name =  spc_name+" "+sector1(sec)+" emissions"
      end if
      ;if sec.eq.2 then
      ; ncdf->altitude          = altitude_air
      ; ncdf->altitude_int         = altitude_int_air 
      ;end if 
  end if
    

  if spc_name .eq. "SO2" .and. sec.eq.0 then
; for SO2 add sector specific output
      sect_agriculture  = emiss_out
      sect_energy = emiss_out
      sect_industry = emiss_out
      sect_transport = emiss_out
      sect_resident = emiss_out
      sect_solvents = emiss_out
      sect_waste = emiss_out
      sect_shipping = emiss_out

      sect_agriculture  = wk(:,0,:)
      sect_energy = wk(:,1,:)
      sect_industry = wk(:,2,:)
      sect_transport = wk(:,3,:)
      sect_resident = wk(:,4,:)
      sect_solvents = wk(:,5,:)
      sect_waste = wk(:,6,:)
      sect_shipping = wk(:,7,:)
 
      sect_agriculture = sect_agriculture * 1.e3 * avog/mw * 1.e-4
      sect_energy = sect_energy * 1.e3 * avog/mw * 1.e-4
      sect_industry = sect_industry * 1.e3 * avog/mw * 1.e-4
      sect_transport = sect_transport * 1.e3 * avog/mw * 1.e-4
      sect_resident = sect_resident * 1.e3 * avog/mw * 1.e-4
      sect_solvents = sect_solvents * 1.e3 * avog/mw * 1.e-4
      sect_waste = sect_waste * 1.e3 * avog/mw * 1.e-4
      sect_shipping = sect_shipping * 1.e3 * avog/mw * 1.e-4

      mw_so2 = 64.
      mw_so4_bam = 96.
      mw_so4_mam = 115.
      ag_sol_was_so2 = sect_agriculture 
      ag_sol_was_so4 = sect_agriculture 
      ag_sol_was_so2 = (sect_agriculture+sect_solvents+sect_waste)*0.975
      ag_sol_was_so4 = (sect_agriculture+sect_solvents+sect_waste)*0.025
      ship_so2 = sect_shipping
      ship_so4 = sect_shipping
      ship_so2 = sect_shipping*0.975
      ship_so4 = sect_shipping*0.025
      res_tran_so2 = sect_resident
      res_tran_so4 = sect_resident
      res_tran_so2 = (sect_resident+sect_transport)*0.975
      res_tran_so4 = (sect_resident+sect_transport)*0.025
      ene_ind_so2_surf = sect_energy
      ene_ind_so2_surf = sect_energy+sect_industry

  ; write so4_a1: sulfate accumulation mode for energy+industrial
  ; vertically distributed: 100-300m
     nalt = 8
     dz = 0.05 ;km
     dalt = 2.e4 ;cm = 200m
     altitude = fspan(0,nalt-1,nalt) * dz + 0.5*dz
     altitude_int = fspan(0,nalt,nalt+1) * dz
     altitude!0="altitude"
     altitude_int!0="altitude_int"
     altitude@units="km"
     altitude@long_name="Altitude"
     altitude_int@units="km"
     altitude_int@long_name="Altitude interfaces"
     ene_ind_so2_vert = new((/ntim,nalt,ncol/),"float")
     ene_ind_so2_vert(:,:,:) = 0.
     ene_ind_so2_surf_vol = ene_ind_so2_surf/dalt
     ene_ind_so2_surf_vol=where(ene_ind_so2_surf.eq.0.,0.,ene_ind_so2_surf_vol)
     do ialt=3,6
       ene_ind_so2_vert(:,ialt,:) = ene_ind_so2_surf_vol
     end do
     

     ene_ind_so2_vert!0="time"
     ene_ind_so2_vert!1="altitude"
     ene_ind_so2_vert!2="ncol"
     ene_ind_so2_vert@units="molecules/cm3/s"
     ene_ind_so2_vert@long_name="SO2 Anthropogenic Emissions from Industry"
     ene_ind_so2_vert&time=emiss_out&time
     ene_ind_so2_vert&altitude=altitude

     ene_ind_so2_surf!0="time"
     ene_ind_so2_surf!1="ncol"
     ene_ind_so2_surf@units="molecules/cm2/s"
     ene_ind_so2_surf@long_name="SO2 Anthropogenic Emissions from Industry"
     ene_ind_so2_surf&time=emiss_out&time
     ene_ind_so4_surf = ene_ind_so2_surf
     ene_ind_so4_surf@long_name="SO4 Anthropogenic Emissions from Industry"

     ene_ind_so4_vert = ene_ind_so2_vert
     ene_ind_so4_vert@long_name="SO4 Anthropogenic Emissions from Industry"
     
      ene_ind_so2_surf= 0.975*(ene_ind_so2_surf)
      ene_ind_so4_surf = 0.025*(ene_ind_so2_surf)
      ene_ind_so2_vert= 0.975*(ene_ind_so2_vert)
      ene_ind_so4_vert = 0.025*(ene_ind_so2_vert)
printVarSummary(ene_ind_so4_surf)
printVarSummary(ene_ind_so4_vert)

      PI  = get_pi("float" ) 
   ; calculate number emissions
  ;varname = 'emiss_ag_sol_was'
      diam = 0.134e-6
      rho = 1770.
      mw2 = 115.
      mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
      num_ag_sol_was_so4 = ag_sol_was_so4
      num_ag_sol_was_so4 = ag_sol_was_so4 *mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
      num_ag_sol_was_so4@units="(particles/cm2/s)(molecules/mole)(g/kg)"
   ; varname = 'emiss_shipping'
      diam = 0.261e-6
      rho = 1770.
      mw2 = 115.
      mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
      num_ship_so4 = ship_so4
      num_ship_so4 = ship_so4 *mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
      num_ship_so4@units="(particles/cm2/s)(molecules/mole)(g/kg)"
  ;varname = 'emiss_res_tran'
      diam = 0.0504e-6
      rho = 1770.
      mw2 = 115.
      mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
      num_res_tran_so4 = res_tran_so4
      num_res_tran_so4 = res_tran_so4 *mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
      num_res_tran_so4@units="(particles/cm2/s)(molecules/mole)(g/kg)"
  ;varname = 'emiss_ene_ind'
      diam = 0.261e-6
      rho = 1770.
      mw2 = 115.
      mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
      num_ene_ind_so4_surf = ene_ind_so4_surf
      num_ene_ind_so4_vert = ene_ind_so4_vert
      num_ene_ind_so4_surf = ene_ind_so4_surf *mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
      num_ene_ind_so4_vert = ene_ind_so4_vert *mw2 /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
      num_ene_ind_so4_vert@units="(particles/cm2/s)(molecules/mole)(g/kg)"
  end if

  if spc_name .eq. "BC" .or. spc_name .eq. "OC" then
   PI  = get_pi("float" )
   ; calculate number emissions
      diam = 0.134e-6
     if spc_name .eq. "BC" then
      rho = 1700.
     else 
      rho = 1000.
     end if
      emiss_out1 = emiss_out
      emiss_out1 = 1.4*emiss_out1
      emis_num = emiss_out1
      emiss_out1=sf*emiss_out1
      mw2 = 12.
      sector = (/"Anthropogenic","Biomass Burning","Aviation"/)
      mass_particle = rho *(PI/6.) *(diam)^3  ;mass per particle (kg/particle)
      emis_num = emiss_out*mw2  /mass_particle  ;(particles/cm2/s)(molecules/mole)(g/kg)
      emis_num@units="(particles/cm2/s)(molecules/mole)(g/kg)"
      emis_num@long_name = "number of BC "+sector(sec)+" Emissions"
      if spc_name .eq. "OC" then
       emiss_out1@long_name = "(1.4*OC) POM "+sector(sec)+" Emissions"
       emis_num@long_name = "number of (1.4*OC) POM "+sector(sec)+" Emissions"
      end if
   end if

;
  if sec.lt.1 then
;
; compute global average (in molec/cm2/s)
; annual average in Tg/year
;
     ;area_emis = conform(emiss_inp,area,1)
     area_emis = conform(emiss_out,area,1)
     emiss_out_gl = emiss_out*area_emis
     emiss_out_gl1 = dim_sum_n_Wrap(emiss_out_gl,1) / dim_sum_n_Wrap(area_emis,1)

      dims = dimsizes(emiss_out_gl1)
      nyr = dims(0)/12
      emiss_out_gl_ann = emiss_out_gl1(0:nyr-1)
      do n=0,nyr-1
        n1 = n * 12
        n2 = n1 + 11
        emiss_out_gl_ann(n) = dim_avg_Wrap(emiss_out_gl1(n1:n2))
      end do
      emiss_out_gl_ann = emiss_out_gl_ann * earth_area * 86400. * 365 * mw * 1.e-12/avog
;
     print(emiss_out_gl_ann)
     delete(emiss_out_gl)
     delete(emiss_out_gl1)
     delete(emiss_out_gl_ann)
     delete(area_emis)
   end if
;
; extract time range string from input file name
;
      length = strlen(file_inp1(0))
      char   = stringtochar(file_inp1(0))
      time_range1 = char(length-16:length-10)
      length2 = strlen(file_inp2(0))
      char2   = stringtochar(file_inp2(0))
      time_range2 = char2(length-9:length-4)
      delete(char)
      delete(char2)
;
; define creation date
;
      creation_date = systemfunc ("date +%Y%m%d")
;
; output to netcdf file
;
      ;file_out = dir_out + "emissions-cmip6"+"_" + spc_name + "_anthro_surface_mol_" + time_range + "_" + resolution + "_c" + tostring(creation_date) + ".nc"
    if sec.le.1 then
      file_out = dir_out + "emissions-cmip6"+"_" + spc_name + "_"+sector1(sec)+"_surface_" + time_range1 + time_range2 +  "_" + resolution + "_c" + tostring(creation_date) + ".nc"
    else
      file_out = dir_out + "emissions-cmip6"+"_" + spc_name + "_"+sector1(sec)+"_vertical_" + time_range1 + time_range2 +  "_" + resolution + "_c" + tostring(creation_date) + ".nc"
    end if
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp1,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
     
      sector = (/"Anthropogenic","Biomass Burning","Aviation"/)
      file_att               = True
      file_att@data_title = sector(sec)+" emissions of "+spc_name+" for CMIP6" 
      file_att@molecular_weight = mw
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_creator = USER NAME (user_email)
      file_att@cesm_contact = USER NAME (user_email)
    if sec.ne.1 then
      file_att@data_summary = "Emissions from the Community Emission Data System (CEDS) have been manipulated for use in CESM2 for CMIP6."
      file_att@data_source_url = "http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/" ;
      file_att@data_reference = "Hoesly et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)" 
    else
      file_att@data_summary = "Biomassburning Emissions  based on GFED4s, FIRE-MIP, visibility-observations and GCD data that have been manipulated for use in CESM2 for CMIP6."
      file_att@data_reference = "Van Marle et al (2016), http://www.geosci-model-dev.net/special_issue590.html" 
    end if
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_hist/rename_and_convert_cmip_hist.ncl"
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)
;
; define date array
;
      utc_date = cd_calendar(time1, 0)
      date = doubletointeger(time1)
      date@units = "YYYYMMDD"
      date@long_name = "date"
      date!0 = "time"
      date&time = time1
      date = tointeger(utc_date(:,0)) * 10000 + tointeger(utc_date(:,1)) * 100 + tointeger(utc_date(:,2))
;
      ncdf->emiss= emiss_out
      ncdf->date         = date
      ncdf->time         = time1
      ncdf->lat          = lat
      ncdf->lon          = lon
      ncdf->area = area 
      ;ncdf->time_bnds    = time_bnds
      if sec.eq.2 then
       ncdf->altitude          = altitude_air
       ncdf->altitude_int         = altitude_int_air 
      end if 
;
; clean up
      delete(ncdf)

printVarSummary(emiss_out)
print(max(emiss_out))
print(min(emiss_out))

 if spc_name .eq. "BC" .or. spc_name .eq. "OC" .or. spc_name .eq. "SO2" then
  ;
; output to netcdf file
;
    if sec.le.1 then
      file_out = dir_out + "emissions-cmip6"+"_" + spc_name2 + "_"+sector1(sec)+"_surface_" + time_range1 + time_range2 +  "_" + resolution + "_c" + tostring(creation_date) + ".nc"
    else
      file_out = dir_out + "emissions-cmip6"+"_" + spc_name2 + "_"+sector1(sec)+"_vertical_" + time_range1 + time_range2 +  "_" + resolution + "_c" + tostring(creation_date) + ".nc"
    end if
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp1,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_creator = USER NAME (user_email)
      file_att@cesm_contact = USER NAME (user_email)
      file_att@data_title = "Number of "+sector(sec)+" emissions of "+spc_name+" for CMIP6" 
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)
;
; define date array
;
      utc_date = cd_calendar(time1, 0)
      date = doubletointeger(time1)
      date@units = "YYYYMMDD"
      date@long_name = "date"
      date!0 = "time"
      date&time = time1
      date = tointeger(utc_date(:,0)) * 10000 + tointeger(utc_date(:,1)) * 100 + tointeger(utc_date(:,2))
;
      ncdf->emiss= emiss_out1
      ncdf->date         = date
      ncdf->time         = time1
      ncdf->lat          = lat
      ncdf->lon          = lon
      ncdf->area = area 
      ;ncdf->time_bnds    = time_bnds
      if sec.eq.2 then
       ncdf->altitude          = altitude_air
       ncdf->altitude_int         = altitude_int_air
      end if

;
; clean up
      delete(ncdf)
      delete(emiss_out1)


   ; output to netcdf num file
;
    if sec.le.1 then
      file_out = dir_out + "emissions-cmip6"+"_num_" + spc_name2 + "_"+sector1(sec)+"_surface_" + time_range1 + time_range2 +  "_" + resolution + "_c" + tostring(creation_date) + ".nc"
    else
      file_out = dir_out + "emissions-cmip6"+"_num_" + spc_name2 + "_"+sector1(sec)+"_vertical_" + time_range1 + time_range2 +  "_" + resolution + "_c" + tostring(creation_date) + ".nc"
    end if
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_creator = USER NAME (user_email)
      file_att@cesm_contact = USER NAME (user_email)
      file_att@data_title = "Number of "+sector(sec)+" emissions of "+spc_name2+" for CMIP6" 
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)
;
; define date array
;
;
      ncdf->emiss_num = emis_num 
      ncdf->date         = date
      ncdf->time         = time1
      ncdf->lat          = lat
      ncdf->lon          = lon
      ncdf->area = area 
      if sec.eq.2 then
       ncdf->altitude          = altitude_air
       ncdf->altitude_int         = altitude_int_air
      end if


      ;ncdf->time_bnds    = time_bnds
;
; clean up
      delete(ncdf)
      delete(emis_num)

  end if ; bc and pom and so4


  if spc_name .eq. "SO2" .and. sec.eq.0 then
 ;  print out additional files    
  spec_name_sec = (/"SO2_anthro-ag-ship-res_surface_mol_","SO2_anthro-ene_surface_mol_","so4_a1_anthro-ene_vertical_mol_","so4_a1_anthro-ag-ship_surface_mol_","so4_a2_anthro-res_surface_mol_","num_so4_a1_anthro-ene_vertical_mol_","num_so4_a1_anthro-ag-ship_surface_mol_","num_so4_a2_anthro-res_surface_mol_"/)
   
;s1var = (/2,5/)
;  do s2 = 0,1
;  s1 = s1var(s2)
  do s1 = 0,7
     spec_ceds = spc_name
     if s1.eq.0 then
       emiss_out1 =  ag_sol_was_so2
       emiss_out2 =  res_tran_so2
       emiss_out3 =  ship_so2
       emiss_out1@long_name = "SO2*0.975 anthro agriculture+waste emissions"
       emiss_out2@long_name = "SO2*0.975 anthro residential+transportation emissions"
       emiss_out3@long_name = "SO2*0.975 anthro shipping emissions"
     end if
     if s1.eq.1 then
       emiss_out1 = ene_ind_so2_surf 
       emiss_out1@long_name = "SO2*0.975 anthro energy+industrial emissions"
     end if
     if s1.eq.2 then
       emiss_out1 = ene_ind_so4_vert
       emiss_out1@long_name = "SO2*0.025 anthro energy+industry emissions"
     end if
     if s1.eq.3 then
       emiss_out1 =  ag_sol_was_so4
       emiss_out2 =  ship_so4
       emiss_out1@long_name = "SO2*0.025 anthro agriculture+waste emissions"
       emiss_out2@long_name = "SO2*0.025 anthro shipping emissions"
     end if
     if s1.eq.4 then
       emiss_out1 =  res_tran_so4
       emiss_out1@long_name = "SO2*0.975 anthro residential+transportation emissions"
     end if
     if s1.eq.5 then
       emiss_out1 = num_ene_ind_so4_vert
       emiss_out1@long_name = "particle number emissions of SO2*0.025 anthro energy+industry"
     end if
     if s1.eq.6 then
       emiss_out1 =  num_ag_sol_was_so4
       emiss_out2 =  num_ship_so4
       emiss_out1@long_name = "particle number emissions of SO2*0.025 anthro agriculture+waste"
       emiss_out2@long_name = "particle number emissions of SO2*0.025 anthro shipping"
     end if
     if s1.eq.7 then
       emiss_out1 =  num_res_tran_so4
       emiss_out1@long_name = "particle number emissions of SO2*0.975 anthro residential+transportation"
     end if

print(spec_name_sec(s1))
      spec_name11 = spec_name_sec(s1)
      file_out = dir_out + "emissions-cmip6"+"_" + spec_name11+ time_range1 + time_range2 +  "_" + resolution + "_c" + tostring(creation_date) + ".nc"
      system("/bin/rm -f " + file_out) ; remove any pre-existing file
      ncdf = addfile(file_out,"c")
;
; find file name (without path)
;
      strs = str_split(file_inp1,"/")
      indx = dimsizes(strs)
      file_name = strs(indx-1)
      delete(strs)
;
; define global attributes
;
      file_att               = True
      file_att@source_file   = file_name
      file_att@creation_date = creation_date
      file_att@data_title = "Number of "+sector(sec)+" emissions of "+spec_name11+" for CMIP6" 
      file_att@data_creator = USER NAME (user_email)
      file_att@cesm_contact = USER NAME (user_email)
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)
;
; define date array
;
      utc_date = cd_calendar(time1, 0)
      date = doubletointeger(time1)
      date@units = "YYYYMMDD"
      date@long_name = "date"
      date!0 = "time"
      date&time = time1
      date = tointeger(utc_date(:,0)) * 10000 + tointeger(utc_date(:,1)) * 100 + tointeger(utc_date(:,2))
;
    if s1.le.4 then
      ncdf->emiss= emiss_out1
      file_att@molecular_weight = mw
    else
      ncdf->emiss_num= emiss_out1
      file_att@molecular_weight = mw2
    end if
     if s1.eq.0 .or. s1.eq.3 .or. s1.eq.6 then
      ncdf->emiss2= emiss_out2
      delete(emiss_out2)
     end if
     if s1.eq.0 then
      ncdf->emiss3 = emiss_out3
      delete(emiss_out3)
     end if
      ncdf->date         = date
      ncdf->time         = time1
      ncdf->lat          = lat
      ncdf->lon          = lon
      ncdf->area = area 
      if s1.eq.2 .or. s1.eq.5 then
       ncdf->altitude = altitude
       ncdf->altitude_int = altitude_int
       printVarSummary(altitude_int)
       printVarSummary(altitude)
      end if
      ;ncdf->time_bnds    = time_bnds
      if sec.eq.0 then 
       delete(spec_ceds)
      end if
      delete(ncdf)
      delete(emiss_out1)
  end do
  delete(altitude)
  delete(altitude_int)


  end if
  


;
; clean up
;
      delete(wk)
      delete(emiss_inp)
      delete(emiss_out)
      delete(in)
      delete(file_inp1)
      delete(file_inp2)
      delete(time)
     ; delete(time_bnds)
      delete(utc_date)
      delete(date)
 ;  delete(spc_name)
    delete(vnames)
    delete(time1)
    delete(nnames)
    delete(nwk)
    delete(ntim)
    delete(ncol)
;
;
  end if ;   sectors exist
  end do    ; sector

if nspc.gt.30 then
do sec=0,1
 if spc_name.eq. "IVOC"
   specs_in = (/"C3H6", "C3H8", "C2H6", "C2H4", "BIGENE", "BIGALK", "CH3COCH3", "MEK", "CH3CHO", "CH2O", "BENZENE", "TOLUENE", "XYLENES"/)
   hist_spec = "IVOC=0.2*("+str_join(specs_in,"+")+")"
   fact = 0.2
 end if
  if spc_name.eq. "SVOC"
   specs_in = "pom_a4"
   hist_spec = "SVOC=0.6*pom_a4"
   fact = 0.6
  end if
 if spc_name.eq. "SOAGx1.5"
  if sec.eq.0 then
   specs_in  = (/"BIGALK","BIGENE","TOLUENE","BENZENE","XYLENES"/)
   yield  = (/0.05, 0.05, 0.15, 0.15, 0.15/)
   hist_spec = "SOAG=1.5* yields(0.05, 0.05, 0.15, 0.15, 0.15) from ("+str_join(specs_in,"+")+")"
  end if
  if sec.eq.1 then
   specs_in  = (/"BIGALK","BIGENE","TOLUENE","BENZENE","XYLENES","ISOP","MTERP"/)
   yield  = (/0.05, 0.05, 0.15, 0.15, 0.15, 0.04, 0.25/)
   hist_spec = "SOAG=1.5* yields(0.05, 0.05, 0.15, 0.15, 0.15, 0.04, 0.25) from ("+str_join(specs_in,"+")+")"
  end if
   fact = 1.5
 end if
 ; read required files 
   cnt_s = 0
   nsp = dimsizes(specs_in)
   do sp = 0,nsp-1
    files = systemfunc("ls " + dir_out + "emissions-cmip6"+"_" + specs_in(sp)+ "_"+sector1(sec)+"_surface_*"+ resolution + "_c*.nc")
    files1=files(0)
print(files1)
    in = addfile(files1,"r")
    if cnt_s.eq.0 then
      time1 = in->time
      date = in->date
      area = in->area
      length = strlen(files1(0))
      char   = stringtochar(files1(0))
      time_range = char(length-34:length-22)
      print(time_range)
      delete(char)
      lat = in->lat
      lon = in->lon
      wk = in->emiss
      printVarSummary(wk)
print("Max: "+max(wk))
print("Min: "+min(wk))
      mw_hc =  wk@molecular_weight
      if spc_name.eq."SOAGx1.5" then
       wk = yield(sp)*wk*mw_hc/mw
      else
       wk = wk*mw_hc/mw
     end if
    else
      emis_hc1 = in->emiss
      printVarSummary(emis_hc1)
      mw_hc =  emis_hc1@molecular_weight
      print(mw_hc)
      print(mw)
      printVarSummary(wk)
      if spc_name.eq."SOAGx1.5" then
      wk = wk + yield(sp)*emis_hc1*mw_hc/mw
     else
      wk = wk + emis_hc1*mw_hc/mw
     end if
print("Max2: "+max(wk))
print("Min2: "+min(wk))
      delete(emis_hc1)
    end if
    ;IVOC = 0.2*(HCs)
    ;SVOC = 0.6*pom_a4
    wk = fact*wk
    cnt_s = cnt_s+1 
   end do    
  emiss_out = wk
  emiss_out@history= hist_spec
  file_out = dir_out + "emissions-cmip6"+"_" + spc_name+ "_"+sector1(sec)+"_surface_" + time_range +  "_" + resolution + "_c" + tostring(creation_date) + ".nc"
  printVarSummary(emiss_out)
; emiss_out@units="molecules/cm2/s"
  emiss_out@long_name= spc_name+" "+sector1(sec)+" emissions"
  system("/bin/rm -f " + file_out) ; remove any pre-existing file
  ncdf = addfile(file_out,"c")
  
      file_att               = True
      file_att@data_title = sector1(sec)+" emissions of "+spc_name+" for CMIP6"
      file_att@molecular_weight = mw
      file_att@creation_date = creation_date
      file_att@data_creator = USER NAME (user_email)
      file_att@cesm_contact = USER NAME (user_email)
      file_att@summary = "Emissions from the Community Emission Data System (CEDS) have been manipulated for use in CESM2 for CMIP6."
      file_att@data_script = "/glade/u/home/tilmes/scripts/CMIP6/emissions_hist/rename_and_convert_cmip_hist.ncl"
      file_att@data_source_url = "http://www.globalchange.umd.edu/ceds/ceds-cmip6-data/" ;
      file_att@data_reference = "Hoesly et al., GMD, 2017 (http://www.geosci-model-dev-discuss.net/gmd-2017-43/)"
      fileattdef(ncdf,file_att)
      filedimdef(ncdf,"time",-1,True)
  ncdf->emiss= emiss_out
  ncdf->date         = date
  ncdf->time         = time1
  ncdf->lat          = lat
  ncdf->lon          = lon
  ncdf->area = area 

  delete(ncdf)
  delete(time1)
  delete(date)
  delete(emiss_out)
  delete(wk)
  delete(specs_in)
 
 if spc_name.eq."SOAGx1.5" then
  delete(yield)
 end if
 end do
 end if  ; SOAG, IVOC, SVOC

  if nspc .ge.6 .and. nspc.le.30 then
   delete(spec_ceds)
  end if


  end do    ; species  
;
end
