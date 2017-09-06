 // 先高斯平滑Iq，之后计算Fq，对Fq进行damping，实时显示G(r))

//Guassian 平滑函数
realimage smooth_Fq(object self, realimage img, number fwhm, number integ)
{
  realimage Gauss, Gauss_matrix, img_matrix, blur_img, kernal;
  number xsize,ysize, norm, alpha,n=2
  
  //缩小图像
  img.GetSize(xsize,ysize)
  number scale=1/n //只用缩略图
  Image SmallImg:= RealImage("",4,xsize*scale,1) 
  SmallImg =warp(img,icol/scale,irow/1)

  number xsize1,ysize1
  SmallImg.GetSize(xsize1,ysize1)
  
  //定义高斯函数
  fwhm=fwhm*scale
  Gauss = exprsize(1,xsize1,0)
  alpha = log(2)*4.0/(fwhm)**2;
  Gauss = exp(-alpha*((icol-ysize1/2)**2+(irow-xsize1/2)**2));
  Gauss/=sum(Gauss)
    
  kernal       = exprsize(xsize1,xsize1,0);
  blur_img     = exprsize(xsize1,0);
  Gauss_matrix = exprsize(xsize1,xsize1,0);
  img_matrix   = exprsize(xsize1,xsize1,0);

  slice2(Gauss_matrix, 0, 0, 0, 0, xsize1, 1, 1, xsize1, 1)   = Gauss[0,irow+xsize1/2-xsize1/xsize1*icol]; 
  slice2(img_matrix, 0, 0, 0,   0, xsize1, 1, 1, xsize1, 1)   = SmallImg[irow,0]; 
 
  kernal = Gauss_matrix*img_matrix;
  if((xsize1+1)%2==0) //Don't use Simpson if this condition is not satisfied.  Just becomes skyscraper summation instead
    kernal *= (4*(irow%2)+2*((irow+1)%2))/3.0*tert(irow==0,0,1)*tert(irow==xsize1-1,0,1)+tert(irow==0,1,0)+tert(irow==xsize1-1,1,0)
  blur_img[icol,0] += slice2(kernal, 0, 0, 0, 0, xsize1, 1, 1, xsize1, 1); //integration summation

  //还原
  realimage OutImg=img*0
  OutImg=warp(blur_img,icol*scale,irow*1)

  return OutImg 
}

//定义G(r)函数 
image RIF2Gr(image damped,number xscale)
{
//显示G(r)
number j, k, s, rmax, smax,  nums,ysize
number PI = pi()
number rsize=400
rmax=15

getsize(damped, nums, ysize)

number scale=1/5  //只用缩略图
Image dampImg:= RealImage("",4,nums*scale,1) 
dampImg = warp(damped, icol/scale,irow/1)

getsize(dampImg, nums, ysize)

smax=nums*xscale/scale

if((nums+1)%2!=0) nums-=1; 

image Gr, Gr_matrix, phi_matrix,f_sr
f_sr = exprsize(nums,rsize,0)
Gr_matrix = exprsize(nums,rsize,0);
Gr = exprsize(rsize,0);
phi_matrix = exprsize(nums,rsize,0)
		
f_sr  = sin(icol*smax/nums * irow*rmax/rsize)
f_sr *=2/Pi()  //这是F(Q)时的因子，若为F(S)则为8*Pi
f_sr*=(4*(icol%2)+2*((icol+1)%2))/3.0*tert(icol==0,0,1)*tert(icol==nums-1,0,1)+tert(icol==0,1,0)+tert(icol==nums-1,1,0)

slice2(phi_matrix, 0, 0, 0, 0, nums, 1, 1, rsize, 1) = dampImg[icol,0] 
Gr_matrix = f_sr*phi_matrix
Gr[irow,0] += slice2(Gr_matrix, 0, 0, 0, 0, nums, 1, 1, rsize, 1) //integration summation
Gr *= smax/nums; //integration interval

return gr
}

//定义RIFFilter函数 
image RIFFilter(number filterxsize, number filtertype, number mode, number lowerthreshold, number upperthreshold, string &filtername)
{

// The imgfilternumber specifies the filter type
// 0=Boxcar
// 1= Triangular
// 2=Trapezoidal
// 3=Happ-Genzel
// 4=3-Term Blackman Harris
// 5=4-Term Blackman Harris


image bigfilterprofile

// Filter is Boxcar
if(filtertype==0) 
	{
	bigfilterprofile=integerimage("",2,1,filterxsize*10,1)
	bigfilterprofile=1

		// Trap for an unsuitable boxcar upper threshold >1 is 1

		if (upperthreshold>filterxsize-1) upperthreshold=filterxsize-1

		bigfilterprofile=0
		bigfilterprofile[0,0,1,upperthreshold]=1
		filtername="Boxcar"
	}


// Filter is triangular
if(filtertype==1) 
	{
	bigfilterprofile=integerimage("",2,1,filterxsize*10,1)
	bigfilterprofile=1
	
		bigfilterprofile=1-(icol/upperthreshold)
		bigfilterprofile=tert(bigfilterprofile<0,0,bigfilterprofile)
		filtername="Triangular"
	}
	

// Filter is Trapezoidal
if(filtertype==2) 
	{
	bigfilterprofile=integerimage("",2,1,filterxsize*10,1)
	bigfilterprofile=1
	
		// Compute the sloping part of the trapezoid as a straight line equation
		// then use tert to set everything below the lower threshold to 1 and above it to 0	

		number slope=-1/(upperthreshold-lowerthreshold)
		number const=1-(slope*lowerthreshold)

		bigfilterprofile=slope*(icol)+const
		bigfilterprofile=tert(bigfilterprofile>1, 1, bigfilterprofile)
		bigfilterprofile=tert(bigfilterprofile<0,0,bigfilterprofile)

		filtername="Trapezoidal"
	}


// Filter is Happ-Genzel
if(filtertype==3) 
	{
	bigfilterprofile=integerimage("",2,1,filterxsize*10,1)
	bigfilterprofile=1
	
		bigfilterprofile=0.54+(0.46*cos(pi()*icol/upperthreshold))
		bigfilterprofile=tert(icol>upperthreshold, 0,bigfilterprofile)
		filtername="Happ-Genzel"
	}


// Filter is 3-Term Blackman Harris
if(filtertype==4) 
	{
	bigfilterprofile=integerimage("",2,1,filterxsize*10,1)
	bigfilterprofile=1
	
		bigfilterprofile=0.42323+(0.49755*cos(pi()*icol/upperthreshold))+(0.07922*(cos(pi()*2*icol/upperthreshold)))
		bigfilterprofile=tert(icol>upperthreshold,0,bigfilterprofile)
		filtername="3-T B-H"
	}

// Filter is 4-Term Blackman Harris
if(filtertype==5) 
	{
	bigfilterprofile=integerimage("",2,1,filterxsize*10,1)
	bigfilterprofile=1
	
		bigfilterprofile=0.35875+(0.48829*cos(1*pi()*(icol-lowerthreshold)/upperthreshold))+(0.14128*(2*cos(pi()*(icol-lowerthreshold)/upperthreshold))+(0.01168*(cos(4*pi()*(icol-lowerthreshold)/upperthreshold))))
		//bigfilterprofile=tert(icol>upperthreshold,0,bigfilterprofile)
		//bigfilterprofile=tert(icol<lowerthreshold,0,bigfilterprofile)
		//bigfilterprofile=bigfilterprofile/bigfilterprofile.max()
		filtername="4-T B-H"
	}

	/*A 7-term Blackman-Harris window has the highest dynamic range; it is ideal for signal-to-noise ratio applications. A 7-term Blackman Harris window is applied to the waveform using the following equation:
y[i] = x[i] × [a0 – a1cos(w) + a2cos(2w) – a3cos(3w) + a4cos(4w) – a5cos(5w) + a6cos(6w)]
where w = (2)i/n
           n is the waveform size
           a0 = 0.27105140069342
           a1 = 0.43329793923448
           a2 = 0.21812299954311
           a3 = 0.06592544638803
           a4 = 0.01081174209837
           a5 = 0.00077658482522
           a6 = 0.00001388721735.	
*/

// Filter is 7-Term Blackman Harris
if(filtertype==6) 
	{
	bigfilterprofile=integerimage("",2,1,filterxsize*10,1)
	bigfilterprofile=1
	
		bigfilterprofile=0.27105140069342+0.43329793923448*cos(1*pi()*(icol-lowerthreshold)/upperthreshold)+0.21812299954311*cos(2*pi()*(icol-lowerthreshold)/upperthreshold)+ 0.06592544638803*cos(3*pi()*(icol-lowerthreshold)/upperthreshold)+ 0.01081174209837*cos(4*pi()*(icol-lowerthreshold)/upperthreshold)+0.00077658482522*cos(5*pi()*(icol-lowerthreshold)/upperthreshold)+0.00001388721735*cos(6*pi()*(icol-lowerthreshold)/upperthreshold)
		
		bigfilterprofile=tert(icol>upperthreshold,0,bigfilterprofile)
		//bigfilterprofile=tert(icol<lowerthreshold,0,bigfilterprofile)
		//bigfilterprofile=bigfilterprofile/bigfilterprofile.max()
		filtername="7-T B-H"
	}

/*	
The Low Side Lobe window further reduces the size of the main lobe. The following equation defines the Low Side Lobe window: 
where N is the length of window
           w = (2n)/N
           a0 = 0.323215218
           a1 = 0.471492057
           a2 = 0.17553428
           a3 = 0.028497078
           a4 = 0.001261367
*/
		   
// Filter is Low side lobe
if(filtertype==7) 
	{
	bigfilterprofile=integerimage("",2,1,filterxsize*10,1)
	bigfilterprofile=1
	
		bigfilterprofile=0.323215218+(-1)*0.471492057*cos(1*pi()*(icol-lowerthreshold)/upperthreshold)+(-1)**2*0.17553428*cos(2*pi()*(icol-lowerthreshold)/upperthreshold)+(-1)**3*0.028497078*cos(3*pi()*(icol-lowerthreshold)/upperthreshold)+(-1)**4*0.001261367*cos(4*pi()*(icol-lowerthreshold)/upperthreshold)
		bigfilterprofile=tert(icol>upperthreshold,0,bigfilterprofile)
		//bigfilterprofile=tert(icol<lowerthreshold,0,bigfilterprofile)
		//bigfilterprofile=bigfilterprofile/bigfilterprofile.max()
		filtername="Low-Side-Lobe"
	}

	
//Gaussian damping function: B=exp(-rQ)2/2)	
if(filtertype==8) 
	{
		bigfilterprofile=integerimage("",2,1,filterxsize*1,1)
		bigfilterprofile=1
	
		number xscale=lowerthreshold
		number ratio=30/(filterxsize*xscale)  //假设Q=30 时衰减到0
		number Qdamp=upperthreshold
		bigfilterprofile=exp(-(icol*xscale*Qdamp*ratio)**2/2)

		filtername="Gaussian Filter"
	}
	
	// bigfilterprofile is a 10x expanded in x image to allow for filters which run outside the 
	// range of the normal image. For example a Happ-Genzel may reach zero at 3.5 rather than 1.
	// The filter is computed for the large image, then is cropped to the same size as the profile	
	image filterprofile=integerimage("",2,1,filterxsize,1)
	filterprofile=bigfilterprofile[0,0,1,filterxsize]

	// implement the mode 0=normal, 1=flipped left to right, 2=mirror symmetric
	
	// ensure the passed in mode value is sensible ie 0, 1 or 2
	
	mode=round(abs(mode))
	if (mode>2) mode=2
	
	
	// Create the mirror image ie the filterimage flipped about the vertical axis
	image flipfilter=imageclone(filterprofile)
	flipfilter=0
	flipfilter=filterprofile[filterxsize-icol,irow]
	
	
	// Modify the calculate profile to reflect the mode	
	if(mode==0) // normal mode
		{
			return filterprofile
		}
	if(mode==1) // flipped left to right mode
		{
			return flipfilter
		}
	
	if(mode=2) // mirror symmetric mode
		{
			// the filterprofile and its mirror image need to be subsampled by 2
			// so that the former comprises the first half of the profile and the 
			// latter the second half			
			image twohalves=imageclone(filterprofile)
			twohalves=0			
			
			// set the first half to every second pixel in the filter image
			twohalves[0,0,1, filterxsize/2]=filterprofile[filterxsize-(2*icol),irow]			
			
			// set the second half to every second pixel in the flipped image			
			twohalves[0,filterxsize/2,1,filterxsize]=flipfilter[filterxsize-(2*icol),irow]
			return twohalves
		}
}

 


number token11   //listener11  用于N listener
number token12  //listener12  用于取消N listener
number token21   //listener21  用于c listener
number token22  //listener22  用于取消c listener
number token31   //listener31  用于damp listener
number token32  //listener32  用于取消damp listener
number token41   //listener41  用于smooth listener
number token42  //listener42  用于取消smooth listener


// 定义类
class ImageDisplayEventListener : object
{
//定义变量		
ROI 	theroi, theroi2,theroi3, theroi4
number  left, right, SpecialROIID,SpecialROIID2,SpecialROIID3,SpecialROIID4,counter
image spectrum,grimg
	
	// 扑捉ROI	
	ROI GetWin(object self) 
		{
			return theroi
			return theroi2
			return theroi3
			return theroi4
		}

//N roi is changed
void ROIChanged1( Object self, Number event_flags, ImageDisplay imgdisp, Number roi_change_flags, Number roi_disp_change_flags, ROI theroi )
{	
counter=counter+1  //用于记录执行的次数

if(counter==1)
{
number xsize, ysize,raw_N_factor,raw_c_factor
spectrum.getsize(xsize,ysize)
number xscale=spectrum.ImageGetDimensionScale(0)

number thisroiid=theroi.roigetid()
if(thisroiid!=SpecialROIID)return

image Iq=spectrum
image f2qimg=imgdisp{"<f>2"}
image fq2img=imgdisp{"<f2>"}
image Fq=spectrum*0
image DampFilter=imgdisp{"Filter"}
image damped=spectrum*0
image Smooth=spectrum

number k,t,b,l,r,N,c,hi,lo,sample_l,sample_r,n0,N_factor,c_factor,SmoothWidth,SmoothRight
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:N0",n0)
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:C",c) 
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:N factor",N_factor) 
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:C factor",c_factor)
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:Smooth:FWHM",SmoothWidth) 
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:Smooth:Rightpart",SmoothRight) 

theroi.ROIGetRange(l,r)
if(0.5*(l+r)<0.1*xsize)
{
N_factor=0.5*(r+l)*N_factor/(0.89*xsize)
GetNumber("Set the value of N factor:",N_factor,N_factor)
Setnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:N factor",N_factor) 
theroi.ROIsetRange(0.89*xsize,0.89*xsize)
}

if(0.5*(l+r)>0.9*xsize)
{
N_factor=0.5*(r+l)*N_factor/(0.11*xsize)
GetNumber("Set the value of N factor:",N_factor,N_factor)
Setnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:N factor",N_factor) 
theroi.ROIsetRange(0.11*xsize,0.11*xsize)
}

N=0.5*(r+l)
N=N*N_factor

theroi.roisetlabel("N="+format(N*n0,"%4.4f"))
setnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:N",n) 

//对Iq平滑
image smoothIq=self.smooth_Fq(spectrum,SmoothWidth, xsize)
Smooth[0,SmoothRight,1,xsize]=smoothIq[0,SmoothRight,1,xsize]

Fq=(icol*xscale)*(N*Smooth+c-fq2img)/(f2qimg) //采用F(Q)，以便与G(r)统一
damped=Fq*DampFilter //过滤

//更新图层内容
imgdisp{"F(Q)"}=Fq
imgdisp{"Smooth"}=Smooth
imgdisp{"Damped"}=damped

number min1,max1
damped.minmax(min1,max1)
imgdisp.lineplotimagedisplaysetdoautosurvey(0,0)
imgdisp.lineplotimagedisplaysetcontrastlimits(1.3*min1,1.1*max1) //y轴显示范围

grImg=RIF2Gr(damped,xscale)
ImageDisplay grdisp=grImg.ImageGetImageDisplay(0)

number grsizex,grsizey
grImg.GetSize(grsizex,grsizey)
grImg[0,0.1*grsizex,1,grsizex].minmax(min1,max1)
grdisp.LinePlotImageDisplaySetDoAutoSurvey( 0, 0 ) 
grdisp.LinePlotImageDisplaySetContrastLimits(1.1*min1,1.1*max1)

}

else if(counter>1&&counter<=3) //不执行
{

}

else if(counter>3)
{
counter=0
}
}

// C roi is changed
void ROIChanged2( Object self, Number event_flags, ImageDisplay imgdisp, Number roi_change_flags, Number roi_disp_change_flags, ROI theroi2 )
{		
counter=counter+1  //用于记录执行的次数

if(counter==1)
{
number xsize, ysize,N_factor,c_factor,raw_N_factor,raw_c_factor
spectrum.getsize(xsize,ysize)
number xscale=spectrum.ImageGetDimensionScale(0)

number thisroiid=theroi2.roigetid()
if(thisroiid!=SpecialROIID2)return

image Iq=spectrum
image f2qimg=imgdisp{"<f>2"}
image fq2img=imgdisp{"<f2>"}
image Fq=spectrum*0
image DampFilter=imgdisp{"Filter"}
image damped=spectrum*0
image Smooth=spectrum

number k,t,b,l,r,N,c,hi,lo,sample_l,sample_r,n0, SmoothWidth, SmoothRight
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:N factor",N_factor) 
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:C factor",c_factor)
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:N0",n0)
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:N",n) 
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:Smooth:FWHM",SmoothWidth) 
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:Smooth:Rightpart",SmoothRight) 

theroi2.ROIGetRange(l,r)
c=(0.5*(l+r)-0.5*xsize)*c_factor   //极值的10倍

theroi2.roisetlabel("c="+format(c,"%4.3f"))
setnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:C",c) 

//对Iq平滑
image smoothIq=self.smooth_Fq(spectrum,SmoothWidth, xsize)
Smooth[0,SmoothRight,1,xsize]=smoothIq[0,SmoothRight,1,xsize]

Fq=(icol*xscale)*(N*smooth+c-fq2img)/(f2qimg) //采用F(Q)，以便与G(r)统一

//过滤
damped=Fq*DampFilter 

//更新图层内容
imgdisp{"F(Q)"}=Fq
imgdisp{"Smooth"}=smooth
imgdisp{"Damped"}=damped

number min1,max1
damped.minmax(min1,max1)
imgdisp.lineplotimagedisplaysetdoautosurvey(0,0)
imgdisp.lineplotimagedisplaysetcontrastlimits(1.3*min1,1.1*max1) //y轴显示范围

grImg=RIF2Gr(damped,xscale)
ImageDisplay grdisp=grImg.ImageGetImageDisplay(0)

number grsizex,grsizey
grImg.GetSize(grsizex,grsizey)
grImg[0,0.1*grsizex,1,grsizex].minmax(min1,max1)
grdisp.LinePlotImageDisplaySetDoAutoSurvey( 0, 0 ) 
grdisp.LinePlotImageDisplaySetContrastLimits(1.1*min1,1.1*max1)  
}

else if(counter>1&&counter<=3) //不执行
{

}

else if(counter>3)
{
counter=0
}
}

//damped roi is changed
void ROIChanged3( Object self, Number event_flags, ImageDisplay imgdisp, Number roi_change_flags, Number roi_disp_change_flags, ROI theroi3 )
{		
counter=counter+1  //用于记录执行的次数

if(counter==1)
{
number xsize, ysize,N_factor,c_factor,raw_N_factor,raw_c_factor
spectrum.getsize(xsize,ysize)
number xscale=spectrum.ImageGetDimensionScale(0)

number thisroiid=theroi3.roigetid()
if(thisroiid!=SpecialROIID3)return

image Iq=spectrum
image f2qimg=imgdisp{"<f>2"}
image fq2img=imgdisp{"<f2>"}
image Fq=spectrum*0
image DampFilter=spectrum*0
image damped=spectrum*0
image Smooth=spectrum

number k,t,b,l,r,N,c,hi,lo,sample_l,sample_r,n0, SmoothWidth, SmoothRight
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:N factor",N_factor) 
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:C factor",c_factor)
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:N0",n0)
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:N",n) 
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:C",c) 
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:Smooth:FWHM",SmoothWidth) 
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:Smooth:Rightpart",SmoothRight) 

number filtervalue
getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Damping Method", filtervalue)
number filtermode=0
string filterstring

if(filtervalue==8)
{
theroi3.ROIGetRange(l,r)

if(l<0.6*xsize)theroi3.ROIsetRange(0.6*xsize,0.6*xsize)  //只能为正

lo=xscale   
number position=0.5*(l+r)   //右侧为衰减参数
hi=0.001*(position-0.6*xsize)   //最大值为10
theroi3.roisetlabel("Damp="+format(hi,"%3.3f"))

Setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Damping Value", hi)
}

else
{
theroi3.ROIGetRange(l,r)
lo=l   //左侧为衰减起始值
number position=r   //右侧为衰减参数

hi=10*(position-l)   //最大值为10
theroi3.roisetlabel("Damp="+format(hi/(xsize-l),"%3.2f"))

Setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Damping Value", hi)
}

image filterimg=RIFFilter(xsize,filtervalue, filtermode, lo, hi, filterstring)  //只是profile中的一部分
DampFilter=filterimg   //显示filter的大小

if(filtervalue!=8)
{
filterimg[0,0,1,lo]=filterimg.GetPixel(lo,0)  //l处为最大值1
filterimg=filterimg/filterimg.max()
}

//对Iq平滑
image smoothIq=self.smooth_Fq(spectrum,SmoothWidth, xsize)
Smooth[0,SmoothRight,1,xsize]=smoothIq[0,SmoothRight,1,xsize]

Fq=(icol*xscale)*(N*smooth+c-fq2img)/(f2qimg) //采用F(Q)，以便与G(r)统一

//过滤
damped=Fq*DampFilter 

//更新图层内容
imgdisp{"F(Q)"}=Fq
imgdisp{"Smooth"}=smooth
imgdisp{"Filter"}=filterimg
imgdisp{"Damped"}=damped

number min1,max1
damped.minmax(min1,max1)
imgdisp.lineplotimagedisplaysetdoautosurvey(0,0)
imgdisp.lineplotimagedisplaysetcontrastlimits(1.3*min1,1.1*max1) //y轴显示范围

grImg=RIF2Gr(damped,xscale)
ImageDisplay grdisp=grImg.ImageGetImageDisplay(0)

number grsizex,grsizey
grImg.GetSize(grsizex,grsizey)
grImg[0,0.1*grsizex,1,grsizex].minmax(min1,max1)
grdisp.LinePlotImageDisplaySetDoAutoSurvey( 0, 0 ) 
grdisp.LinePlotImageDisplaySetContrastLimits(1.1*min1,1.1*max1)  
}

else if(counter>1&&counter<=3) //不执行
{

}

else if(counter>3)
{
counter=0
}

}

//Smooth roi is changed
void ROIChanged4( Object self, Number event_flags, ImageDisplay imgdisp, Number roi_change_flags, Number roi_disp_change_flags, ROI theroi4 )
{		
counter=counter+1  //用于记录执行的次数

if(counter==1)
{
number thisroiid=theroi4.roigetid()
if(thisroiid!=SpecialROIID4)return

image Iq=spectrum
image f2qimg=imgdisp{"<f>2"}
image fq2img=imgdisp{"<f2>"}
image Fq=spectrum*0
image DampFilter=imgdisp{"Filter"}
image damped=spectrum*0
image Smooth=spectrum

number t,b,l,r,xsize, ysize,fwhm
spectrum.getsize(xsize,ysize)
number xscale=spectrum.ImageGetDimensionScale(0)

number n,c
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:N",n) 
getnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:C",c) 

theroi4.roigetrange(l,r)
fwhm=0.10*abs(r-l)
theroi4.roisetlabel("Smooth="+fwhm)
setnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:Smooth:FWHM",fwhm) 
setnumbernote(spectrum,"ElectronDiffraction Tools:PDF:Parameters:Smooth:Rightpart",r) 

//对Iq平滑
image smoothIq=self.smooth_Fq(spectrum, fwhm, xsize)
Smooth[0,r,1,xsize]=smoothIq[0,r,1,xsize]

//计算Fq
Fq=(icol*xscale)*(N*smooth+c-fq2img)/(f2qimg)

//对其过滤
damped=Fq*DampFilter //过滤

imgdisp{"F(Q)"}=Fq
imgdisp{"Smooth"}=smooth
imgdisp{"Damped"}=damped

number min1,max1
damped.minmax(min1,max1)
imgdisp.lineplotimagedisplaysetdoautosurvey(0,0)
imgdisp.lineplotimagedisplaysetcontrastlimits(1.3*min1,1.1*max1) //y轴显示范围

grImg=RIF2Gr(damped,xscale)
ImageDisplay grdisp=grImg.ImageGetImageDisplay(0)

number grsizex,grsizey
grImg.GetSize(grsizex,grsizey)
grImg[0,0.1*grsizex,1,grsizex].minmax(min1,max1)
grdisp.LinePlotImageDisplaySetDoAutoSurvey( 0, 0 ) 
grdisp.LinePlotImageDisplaySetContrastLimits(1.1*min1,1.1*max1)  
}

else if(counter>1&&counter<=3) //不执行
{

}

else if(counter>3)
{
counter=0
}

}

// 新建roi时删除listener，删除所有roi，提取出peaks
void ROIRemoved( Object self, Number event_flags, ImageDisplay imgdisp, Number roi_change_flags, Number roi_disp_change_flags, ROI theroi )
{
// 只对特征roi操作			
number thisroiid=theroi.roigetid()
if(thisroiid!=SpecialROIID) return	

number l,r
theroi.roigetrange(l,r)			
//移除token
imgdisp.ImageDisplayRemoveEventListener(token11)
imgdisp.ImageDisplayRemoveEventListener(token12)
imgdisp.ImageDisplayRemoveEventListener(token21)
imgdisp.ImageDisplayRemoveEventListener(token22)
imgdisp.ImageDisplayRemoveEventListener(token31)
imgdisp.ImageDisplayRemoveEventListener(token32)
imgdisp.ImageDisplayRemoveEventListener(token41)
imgdisp.ImageDisplayRemoveEventListener(token42)
			
//删除临时缓存
TagGroup tg = spectrum.ImageGetTagGroup()
TagGroup tg1
if(tg.taggroupgettagastaggroup("ElectronDiffraction Tools:Background", tg1))
deletenote(spectrum, "ElectronDiffraction Tools:Background")

//删除roi
number roino= imgdisp.ImageDisplayCountROIS()
for (number i=0; i<roino; i++)
{
number index
ROI currentROI = imgdisp.ImageDisplayGetROI( index )
imgdisp.imagedisplaydeleteROI(currentROI)
}

//删除文字标记
number annots=imgdisp.componentcountchildrenoftype(13)
for (number i=0;i<annots;i++)
{
component annotid=imgdisp.componentgetnthchildoftype(13,0)
annotid.componentremovefromparent()
}

grimg.deleteimage()
}
	
		
// This function sources the ROI on the passed in image
		object init(object self, image front,image gr)
		{
			spectrum:=front
			grimg:=gr
			//grimg:=front.findnextimage()
			imagedisplay imgdisp=spectrum.ImageGetImageDisplay(0)
			
			//寻找roi
			number annots= imgdisp.ImageDisplayCountROIS()
			for (number i=0; i<annots; i++)
				{
					ROI currentROI = imgdisp.ImageDisplayGetROI( i )
					string thischar=currentROI .ROIGetLabel()
					string roilabel=left(thischar,1)					
					
					if(roilabel=="c")
					{
					theroi2 = currentROI
					SpecialROIID2=theroi2.roigetid()
					}					

					if(roilabel=="D")  //衰减参数
					{
					theroi3 = currentROI
					SpecialROIID3=theroi3.roigetid()
					}				

					if(roilabel=="S")  //平滑
					{
					theroi4 = currentROI
					SpecialROIID4=theroi4.roigetid()
					
					number l,r
					theroi4.ROIGetRange(l,r)
					theroi4.roisetrange(l-1,r-1)
					}					
					
					if(roilabel=="N")
					{
					theroi = currentROI
					SpecialROIID=theroi.roigetid()					
					}
				}

			// 特征roi的响应			
			number dummy=0			
			self.roichanged1(dummy, imgdisp, dummy, dummy, theroi)
			self.roichanged2(dummy, imgdisp, dummy, dummy, theroi2)
			self.roichanged3(dummy, imgdisp, dummy, dummy, theroi3)
			self.roichanged4(dummy, imgdisp, dummy, dummy, theroi4)
			
			return self
		}
		
	// Constructor
	ImageDisplayEventListener(object self)
		{

		}
		
	// Destructor
	~ImageDisplayEventListener(object self)
		{
		}
}



// function to check the displayed image is appropriate and create and apply
// the event listeners

void main()
{
image img:=getfrontimage()
imagedisplay imgdisp=imagegetimagedisplay(img, 0)
string imgname=img.getname()
image Fq=imgdisp{"F(Q)"}

number rmax=15  //缩略图的大小 ?
number rsize=400
image gr
gr=RealImage("",4,rsize,1)
gr.showimage()
gr.SetName("Live G(r)-"+imgname)

gr.ImageSetDimensionScale(0,rmax/rsize)
gr.ImageSetDimensionUnitString(0,"r (A)")
imagesetIntensityUnitString(gr,"Reduced G(r)")
ImageDisplay grdisp=gr.ImageGetImageDisplay(0)
lineplotimagedisplaysetslicedrawingstyle(grdisp,0,2)   //图层0：G(r), 填充
lineplotimagedisplaysetslicecomponentcolor(grdisp, 0, 1,0,0.6,0.6)

object sliceid=grdisp.imagedisplaygetsliceidbyindex(0)
grdisp.imagedisplaysetslicelabelbyid(sliceid, "G(r)")

imagedocument imgdoc=getfrontimagedocument()
documentwindow tempwin=imgdoc.imagedocumentgetwindow()
windowsetframeposition(tempwin, 185,0)
tempwin.windowsetframesize(550,350)
		
img.ShowImage()

//提示语
result("--------------------------------------------------------------------------\nMove ROI: to change parameter N, c and damping factor Hi\nCtrl+ ROI: to calculate G(r)\nAlt+'S': to set parameters N, c, Hi\n"+"--------------------------------------------------------------------------\n")

//删除已有图层
number noslices=imgdisp.imagedisplaycountslices()
for(number i=noslices-1; i>-1; i--)
{
object sliceid=imgdisp.imagedisplaygetsliceidbyindex(i)
string slicename=imgdisp.imagedisplaygetslicelabelbyid(sliceid)
if(slicename=="Filter"||slicename=="Damped"||slicename=="G(r)"||slicename=="Smooth") imgdisp.imagedisplaydeleteslicewithid(sliceid)
}

//添加一些图层
lineplotimagedisplaysetslicedrawingstyle(imgdisp,3,2)   //图层3：存放RIF或F(Q)图层，填充
lineplotimagedisplaysetslicecomponentcolor(imgdisp, 3, 1,0,0.6,0.6)  //填充颜色
//lineplotimagedisplaysetslicecomponentcolor(imgdisp, 2, 0,0.5,0.5,0.5)  //边缘颜色
imgdisp.lineplotimagedisplaysetgridon(0)  //gridon=0
imgdisp.lineplotimagedisplaysetbackgroundon(0)  //bgon=0

image slice6=imgdisp{"I(Q)"}
imgdisp.imagedisplayaddimage(slice6, "Smooth")  //图层4：存放Smooth图层
lineplotimagedisplaysetslicedrawingstyle(imgdisp,4,1)
lineplotimagedisplaysetslicecomponentcolor(imgdisp, 4, 0,1,0,1)

image slice4=img*0
slice4=1
imgdisp.imagedisplayaddimage(slice4, "Filter")  //图层4：存放Filter
lineplotimagedisplaysetslicedrawingstyle(imgdisp,5,1)
lineplotimagedisplaysetslicecomponentcolor(imgdisp, 5, 0,256,0,256)

image slice5=img*0
imgdisp.imagedisplayaddimage(slice5, "Damped")  //图层5：存放damped F(Q)
imgdisp.LinePlotImageDisplaySetLegendShown(1)
lineplotimagedisplaysetslicedrawingstyle(imgdisp,6,1)
lineplotimagedisplaysetslicecomponentcolor(imgdisp, 6, 0,1,0,0)

//获取N值，用于设置N-ROI的相对位置
number xsize,ysize,n,n_factor,c_factor,c,xscale,Ratio
img.getsize(xsize,ysize)
xscale=img.ImageGetDimensionScale(0)

Getnumbernote(img,"ElectronDiffraction Tools:PDF:Parameters:N",N)   
Getnumbernote(img,"ElectronDiffraction Tools:PDF:Parameters:C",c) 
c_factor=c/(0.2*xsize)
n_factor=N/(0.4*xsize)

number filtervalue,dampVal=0
getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Damping Method", filtervalue)
getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Damping Value", dampVal)

roi roi1=createroi()
ROI1.ROISetrange(n/n_factor-4,n/n_factor+4)   //常数N，0.1*xsize处为拟合值
roi1.roisetcolor(1,0,0)
roi1.roisetvolatile(0)
roi1.roisetlabel("N")	
imgdisp.ImageDisplayAddROI( ROI1)

roi1=createroi()
ROI1.ROISetrange(0.5*xsize-10,0.5*xsize+10)   //Smooth
roi1.roisetcolor(0,1,0)
roi1.roisetvolatile(0)
roi1.roisetlabel("Smooth")	
imgdisp.ImageDisplayAddROI( ROI1)

setnumbernote(img,"ElectronDiffraction Tools:PDF:Parameters:Smooth:FWHM",0.10*20) 
setnumbernote(img,"ElectronDiffraction Tools:PDF:Parameters:Smooth:Rightpart",0.5*xsize+10) 


Setnumbernote(img,"ElectronDiffraction Tools:PDF:Parameters:N factor",N_factor)   //N的比例系数
Setnumbernote(img,"ElectronDiffraction Tools:PDF:Parameters:C factor",c_factor)   //c的比例系数

//result("N range: ( 0,  "+xsize*n_factor+" );      C range: ( -"+0.5*xsize*c_factor+",   "+0.5*xsize*c_factor+" )\n")

roi1=createroi()
ROI1.ROISetrange(c/c_factor+0.5*xsize-4,c/c_factor+0.5*xsize+4)   //常数c，0.5*xsize处为0
roi1.roisetvolatile(0)
roi1.roisetcolor(0,0,1)
roi1.roisetlabel("c")	
imgdisp.ImageDisplayAddROI( ROI1)

if(filtervalue==8)  //Gaussian damping
{
roi1=createroi()
ROI1.ROISetrange(0.6*xsize+1000*dampval,0.6*xsize+1000*dampval) //默认值为缓存值
roi1.roisetvolatile(0)
roi1.roisetcolor(1,0,1)
roi1.roisetlabel("Damp")	
imgdisp.ImageDisplayAddROI( ROI1)
}

else
{
number x,y,maxval
maxval=Fq.max(x,y)
roi1=createroi()
ROI1.ROISetrange(10,10+dampVal/10)   //参数hi，最大值为10，默认1
roi1.roisetvolatile(0)
roi1.roisetcolor(1,0,1)
roi1.roisetlabel("Damp")	
imgdisp.ImageDisplayAddROI( ROI1)
}

// listener-image关联
// Listener for ROI removal
string messagemap11="roi_removed:ROIRemoved"
object ROIRemovalListener11=alloc(ImageDisplayEventListener).init(img,gr)
token11 =imgdisp.ImageDisplayAddEventListener( RoiRemovalListener11, messagemap11)
string messagemap12="roi_changed:ROIChanged1"
object ROIChangeListener12=alloc(ImageDisplayEventListener).init(img,gr)
token12 =imgdisp.ImageDisplayAddEventListener( RoiChangeListener12, messagemap12)	

string messagemap21="roi_removed:ROIRemoved"
object ROIRemovalListener21=alloc(ImageDisplayEventListener).init(img,gr)
token21 =imgdisp.ImageDisplayAddEventListener( RoiRemovalListener21, messagemap21)
string messagemap22="roi_changed:ROIChanged2"
object ROIChangeListener22=alloc(ImageDisplayEventListener).init(img,gr)
token22 =imgdisp.ImageDisplayAddEventListener( RoiChangeListener22, messagemap22)

string messagemap31="roi_removed:ROIRemoved"
object ROIRemovalListener31=alloc(ImageDisplayEventListener).init(img,gr)
token31 =imgdisp.ImageDisplayAddEventListener( RoiRemovalListener31, messagemap31)
string messagemap32="roi_changed:ROIChanged3"
object ROIChangeListener32=alloc(ImageDisplayEventListener).init(img,gr)
token32 =imgdisp.ImageDisplayAddEventListener( RoiChangeListener32, messagemap32)	

string messagemap41="roi_removed:ROIRemoved"
object ROIRemovalListener41=alloc(ImageDisplayEventListener).init(img,gr)
token41 =imgdisp.ImageDisplayAddEventListener( RoiRemovalListener41, messagemap41)
string messagemap42="roi_changed:ROIChanged4"
object ROIChangeListener42=alloc(ImageDisplayEventListener).init(img,gr)
token42 =imgdisp.ImageDisplayAddEventListener( RoiChangeListener42, messagemap42)	
}

// Main script
main()
