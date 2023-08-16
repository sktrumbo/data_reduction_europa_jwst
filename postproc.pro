pro postproc


; take the 3d data cubes and combine them using some light
; bad pixel suppression. 

; use appropriate G star data for each file
restore,'../gstars/gstars.sav'


files='data/EUR*'+['140*NRS1','140*NRS2','235*NRS1','235*NRS2','395*NRS1','395*NRS2']+'*internal_s3d.fits'
vrs=['s1','s2','m1','m2','l1','l2']
strs=['s140h1','s140h2','s235h1','s235h2','s395h1','s395h2']

; read in the dither offsets
readcol,'offsets.txt',file,xxs,yys,format='a,i,i'

for fi=0,n_elements(files)-1 do begin
	void=execute('star='+strs[fi])

	fs=file_search(files[fi],count=count)
	xs=fltarr(n_elements(fs))
	ys=xs
	for i=0,n_elements(fs)-1 do begin
		; mask the obvious dead pixels
		im=readfits(fs[i],ext=1,head)
		err=readfits(fs[i],ext=2)
		mask=readfits(fs[i],ext=3)
		mask2=readfits(fs[i],ext=4)
		w=where(mask ne 0 or im le 0)
		im[w]=!values.f_nan

		; mask cosmic rays
		sz=size(im)
		sz[3]=n_elements(star)<sz[3]
		print,'removing bad pixels and cosmic rays'
		for j=0,sz[1]-1 do begin
			for k=0,sz[2]-1 do begin
				s=reform(im[j,k,*]/star)
				diff=s-median(s,20)
				for l=0,sz[3]-1 do begin
					h=diff[(l-100)>0:(l+100)<(sz[3]-1)]
					h=h(sort(h))
					szh=n_elements(h)
					onesig=(-h[szh*.14]+h[szh*.84])/2.
					if abs(diff[l]) gt 5*onesig then begin
						im[j,k,l]=!values.f_nan
					endif
				endfor
			endfor
		endfor


		sz=size(im)
		sz=size(im)
		w=where(file_basename(fs[i]) eq file,c)
		xxss=xxs[w[0]]
		yyss=yys[w[0]]

		if i eq 0 then begin
			bim=median(im[*,*,200:250],dim=3)
			bim=bim[*,2:-2]
			bspec=fltarr(sz[1],sz[2],sz[3],n_elements(fs))
			bspec[*,*,*,0]=im
			wave=getwave(head)
		endif else begin
			xs=xxss
			ys=yyss
		endelse
		for k=0,sz[3]-1 do begin
			bspec[*,*,k,i]=xyshift(im[*,*,k],xs,ys)
		endfor
	endfor
	sz=size(bspec)
	cspec=fltarr(sz[1],sz[2],sz[3])
	for x=0,sz[1]-1 do begin
		for y=0,sz[2]-1 do begin
			; first median normalize the spectra
			for j=1,sz[4]-1 do begin
				m=median(bspec[x,y,*,j]/bspec[x,y,*,0])
				bspec[x,y,*,j]=bspec[x,y,*,j]/m
			endfor
			for j=0,sz[3]-1 do begin
				pix=bspec[x,y,j,*]
				w=where(finite(pix),c)
				if c eq 0 then begin
					cspec[x,y,j]=!values.f_nan
				endif else begin
					cspec[x,y,j]=total(pix[w])/c
				endelse
			endfor
		endfor
	endfor
	void=execute(vrs[fi]+'=cspec')
endfor
; the spectra are saved into short, medium, long wavelength
; with NRS1 and NRS2 saved separately.
save,s1,s2,m1,m2,l1,l2,fi='europa.sav'
end

