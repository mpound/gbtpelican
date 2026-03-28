pro vanecal,scan1,sky=sky, ifnum=ifnum,maint=maint,nfd=nfd
;
;;Computes Tsys values for Argus beams for ifnum 
;
;;Inputs:
;;scan1 = vane scan (sky scan assume to be scan1+1
;;ifnum = IFnum of spectral window
;;plnum = pol-number =0 for argus
;;nfd = number of argus beams (16 default)
;;maint=1 during maintenance, ignore for normal observing 
;;
;;Output:
;;Prints approximate effective Tsys* for each beam (Tsys* = Tsys *exp(tau)/eta_l)
;;Tsys*=Tcal[Coff]/[Con-Coff]
;;Also prints the mapping between fdnum and beam number as well as
;;weather conditions during the scan 

if (n_elements(sky) eq 0) then sky = scan1+1
if (n_elements(ifnum) eq 0) then ifnum = 0
if (n_elements(maint) eq 0) then maint=0
if (n_elements(nfd) eq 0) then nfd = 16

;;get center beam scan for ATM parameters
gettp,scan1,ifnum=ifnum,fdnum=9,/quiet
;;Compute Tcal
;;twarm in C in header
twarm=!g.s[0].twarm+273.15
if (twarm ge 370.) then begin
    print, FORMAT='(A15, F8.4)', 'twarm was = ', twarm
    twarm = !g.s[0].tambient + 1.5
    print, FORMAT='(A15, F8.4)', 'twarm set to = ', twarm
endif

;twarm=280.0
time=!g.s[0].mjd
el=!g.s[0].elevation
freq = !g.s[0].reference_frequency/1.e9
;elevation not available in data during maintenance
if (maint ne 0) then el=77.85
getatmos,el=el,mjd=time,freq=freq,outvals
;getatmos2, mjd=time, freq=freq, outvals2

if (n_elements(tau) eq 0) then tau=outvals(0)
atmtsys = outvals(1)
tatm = outvals(2)
am=1./sin(el*!pi/180.)
;;
tbg = 2.725
tcal = (tatm -tbg) + (twarm-tatm)*exp(tau*am)

print
print, FORMAT='(A7, F6.2)', 'freq = ', freq
print, FORMAT='(A7, F5.3)', 'ztau = ', tau 
;if (freq lt 100) then print,'ztau =    ', tau
;if (freq ge 100) then print,'ztau =     ', tau
print
;print, '   FDNUM        BEAM       Tsys*[K]'

minTsys = 1000
maxTsys = 0

for i=0,nfd-1 do begin
  gettp,scan1,ifnum=ifnum,fdnum=i,/quiet
  copy,0,2
  gettp,sky,ifnum=ifnum,fdnum=i,/quiet
  copy,0,1
  subtract,2,1
  divide,0,1
  ;print,'fdnum, beam, Tsys*[K]: ',i,!g.s[0].feed,tcal/median(getdata(0))    

  tsys = tcal/median(getdata(0))
  ;if (tsys gt 10000 or tsys lt 0) then tsys = !VALUES.F_NAN

  if (tsys lt minTsys) then begin
    minTsys = tsys
    minFd = i
  endif

  if (tsys gt maxTsys) then begin
    maxTsys = tsys
    maxFd = i
  endif

  if (i MOD 4) eq 0 then tsys0 = tsys
  if (i MOD 4) eq 1 then tsys1 = tsys
  if (i MOD 4) eq 2 then tsys2 = tsys
  if (i MOD 4) eq 3 then tsys3 = tsys

  if (i MOD 4) eq 3 then begin
    print, FORMAT='(F10.1, F10.1, F10.1, F10.1)', tsys0, tsys1, tsys2, tsys3 
    ;print,i,!g.s[0].feed,tcal/median(getdata(0))
  endif
endfor

print
print, FORMAT='(A15, F6.1, A8, I2, A8, I2, A1)', 'lowest Tsys  = ', minTsys, '  (Beam ', minFd+1, ', Fdnum ', minFd, ')'
print, FORMAT='(A15, F6.1, A8, I2, A8, I2, A1)', 'highest Tsys = ', maxTsys, '  (Beam ', maxFd+1, ', Fdnum ', maxFd, ')'
print
print, FORMAT='(A19, F8.2, F8.2, F8.2)', 'Tcal, Twarm, tatm: ', tcal, twarm, tatm
print
return
end



