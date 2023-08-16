pro pipeall

; dir contains the download directory
; the directory needs to contain the "rate" file and the "cal" file
; for each observation

dir='download/01250/*/'

; searach for all cal files
fs=file_search(dir+'*cal*',count=count)
sp=strpos(fs,'_cal.fits')

fsr=fs
for i=0,count-1 do fsr[i]=strmid(fs[i],0,sp[i])+'_rate.fits'
newfile=intarr(count)


; check that all cal files have rate files associated with them
missing=0
for i=0,count-1 do begin
	if not file_test(fsr[i]) then begin
		print,'missing a file:'+fsr[i]
		missing+=1
	endif
endfor

if missing ne 0 then begin
	message,'missing '+string(missing)+' files'
endif

for i=0,count-1 do begin
	h=headfits(fsr[i],/silent)
	; construct a new file name
	fn='data/'+strcompress(sxpar(h,'targprop'),/re)+'_'+strcompress(sxpar(h,'grating'),/re)+'_'+strcompress(sxpar(h,'exposure'),/re)+'_'+strcompress(sxpar(h,'detector'),/re)+'_'+strcompress(sxpar(h,'observtn'),/re)+'_rate.fits'
	
	print,'processing ',fn
	imc=readfits(fs[i],ext=1,/silent)
	im=readfits(fsr[i],ext=1,/silent)

	; use the cal file to define where the good pixels are
	w=where(finite(imc))
	imb=im
	imb[w]=!values.f_nan
	imb[*,926:1134]=!values.f_nan  ; manually mask the central portion

	; fit and remove the 1/f patter noise
	; doing a vertical smoothing too
	smsize=150 ; +/- 150 pixels
	m=fltarr(2048,2048)
	for j=0,2047 do begin
		sp=imb[*,(j-smsize)>0:(j+smsize)<2047]
		sp=median(sp,dim=2)
		m[*,j]=sp
	endfor
	im=im-m
	; replace rate file data with new 1/f removed data
	file_copy,fsr[i],fn,/overwrite
	modfits,fn,im,exten=1

	; file is now ready to be processed. 
	; through the JWST Level 2 pipeline
	spawn,'pipeone.py '+fn,a,b
	file_delete,fsr[i],fs[i]
endfor

end
