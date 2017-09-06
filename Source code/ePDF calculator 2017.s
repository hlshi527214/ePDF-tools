// PDF Calculator是基于D.R.G. Mitchell and T. Petersen的RDF calculator V1.0, 2010版本构建的，添加了Image preparation(HRTEM to SAED, SAED to profile, profile calibration, ROI value)和data save(能保存RIF、G以及PASD中的相关数据)

/*************************************************************************
可通过以下代码设置特征峰位和密度
number defaultpeak=2.494
number defaultaveragedensity=2.711
setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Average Density",defaultaveragedensity)
*************************************************************************/
 number widthfitting=15  // the width of the window for fitting (3-99) // must be odd
 if(!getpersistentnumbernote("ElectronDiffraction Tools:Fit:Width for fitting", widthfitting))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:Fit:Width for fitting", widthfitting)
	}
	
 number width=100  // the FWHM of Gaussian
 if(!getpersistentnumbernote("ElectronDiffraction Tools:Fit:Gaussian width",width))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:Fit:Gaussian width",width)
	}
	
 number polyorder=2  // the polynomial order 2 (for 1D plots) to 3 (for derivatives)
 if(!getpersistentnumbernote("ElectronDiffraction Tools:Fit:Polynomial order", polyorder))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:Fit:Polynomial order", polyorder)
	}	
	
	
//-------------------------自定义函数-----------------------------------
//函数：2D高斯函数
image Create2DGaussian(number xsize, number ysize, number centrex, number centrey, number intensity, number sigmax, number sigmay)
	{
		image gaussian=realimage("", 4,xsize, ysize)
		if(centrex>xsize || centrex<0 || centrey>ysize || centrey<0) return gaussian
		gaussian=intensity*exp(-(icol-centrex)**2/(2*sigmax**2))*exp(-(irow-centrey)**2/(2*sigmay**2))

		return gaussian
	}


// 函数：技术x2，拟合误差
number computechisqrd(image original, image fit, number centrex, number centrey, number centralpcnt)
	{
		number xsize, ysize
		getsize(original, xsize, ysize)
		
		
		// Compare only subarea of the image - centred on the target centre.
		
		if(centralpcnt<25) centralpcnt=25
		if(centralpcnt>75) centralpcnt=75
		
		// Keep the size of the central region sensible >=10 pixels		
		number xdim=xsize*(centralpcnt/100)
		if(xdim<10) xdim=10
		
		number ydim=ysize*(centralpcnt/100)
		if(ydim<10) ydim=10
		
		
		// Define the search area centred on the target centre
		
		number t=round(centrey-(ydim/2))
		number l=round(centrex-(xdim/2))
		number b=round(centrey+(ydim/2))
		number r=round(centrex+(xdim/2))
		
		if(t<0) t=0
		if(l<0) l=0
		if(b>ysize-1) b=ysize-1
		if(r>xsize-1) r=xsize-1
		

		// measure the difference between the original and the fit (gaussian) images within their subareal regions
		// and compute chisqrd		
		image difference=(original[t,l,b,r]-fit[t,l,b,r])**2/(original[t,l,b,r]+0.000001)
		number chisqrd=sum(difference)
		return chisqrd
	}


// 函数：用MonteCarlo方法自动寻找σ参数（高斯函数）
number FindSigmaMonteCarlo(image front, number centrex, number centrey, number subareapcnt, number intensity, number initialsigma, number iterations, number &minchisqrd, number &sigma)
	{
		// Source the image size

		number xsize, ysize, chisqrd
		getsize(front, xsize, ysize)
		number i
		
		//if(centrex==0) centrex=xsize/2
		//if(centrey==0) centrey=ysize/2


		// random walk - vary the initial sigma with a random gaussian function and compute the fit for the resulting sigma		
		for(i=0; i<iterations; i++)
			{
				number sigmatest=initialsigma+gaussianrandom()
				if(sigmatest<1) sigmatest=1
				image testgaussian=Create2DGaussian(xsize, ysize, centrex, centrey, intensity, sigmatest, sigmatest)			
				chisqrd=computechisqrd(front, testgaussian, centrex, centrey, subareapcnt)

				if(chisqrd<minchisqrd) 
					{
						minchisqrd=chisqrd
						sigma=sigmatest
						initialsigma=sigmatest
					}	

			}
		return minchisqrd
	}


// 函数：Monte Carlo自动找中心
number FindCentreMonteCarlo(image front, number sigma, number intensity, number iterations, number centrex, number centrey, number subareapcnt, number &minchisqrd, number &fitcentrex, number &fitcentrey)
	{
		number xsize, ysize, chisqrd, newx, newy
		getsize(front, xsize, ysize)
		image testgaussian=front*0
	
	
		// If no centre information is provided  estimate the centre as the geometric centre
		
		if(centrex==0) centrex=xsize/2
		if(centrey==0) centrey=ysize/2
		number i, randcentrex, randcentrey
		fitcentrex=centrex
		fitcentrey=centrey
	
	
		// Do a random walk looking for an improved fit
		
		for(i=0; i<iterations; i++)
			{
				// vary the centres with a random walk
				
				randcentrex=gaussianrandom()
				randcentrey=gaussianrandom()
				newx=randcentrex+fitcentrex
				newy=randcentrey+fitcentrey
				
				
				// Ignore any values which take the centre outside the bounds of the image
				
				if(newx<0) randcentrex=0
				if(newx>xsize-1) randcentrex=0
				if(newy<0) randcentrey=0
				if(newy>ysize-1) randcentrey=0
			
			
				// test the fit with the new centre
				
				testgaussian=Create2DGaussian(xsize, ysize, newx,newy, intensity, sigma, sigma)
				chisqrd=computechisqrd(front, testgaussian, newx, newy, subareapcnt)

				if(chisqrd<minchisqrd) 
					{
						minchisqrd=chisqrd
						fitcentrex=newx
						fitcentrey=newy
					}	
			}
		
		return minchisqrd
	}
	

// ------------------------定义class ------------------------
// 初始值
number noofelements=3 // 化学元素的数目
number beamvoltage=200// the default beam voltage in kV
number defaultuseatorwtpercent=0 // the default value for the at or wt% radio - 0=at%, 1=wt%
number defaultposition=1 // the variable which defines the position of the profile 1=North, 2=NE, 3=E, 4=SE  . . 8=NW - values are 1-8 incl.

number defaultprofileaspectratio=1.6 // the default aspect ratio of the displayed profiles
number defaultprofileheight=550 // the default height of the displayed profiles

number defaultsorqunits=1 // sets whether profiles are displayed in units of s (0) or Q (1)
number defaultdampingmethod=8 // set the default value of the radio button setting the damping method
number defaultfiltermode=0 // set the filter to normal, 1=flipper about a vertical axis, 2= mirror symmetric

number defaultaveragedensity=2.5 // A value for the average density of the material
number defaultminimumdensity=0.01 // The lowest permissible density value for the density field
number defaultfittotopchannelpcnt=40 // the scattering factor curve is initially fitted to the experimental curve by minimising the difference between
// it and the experimental curve over the top x% of channels - this sets that percentage.
number noiterations=50 // the number of iterations with which to fit a scattering factor plot to an intensity profile
number defaultresolution=2000 // the resolution of the G(r) and related plots in pixels
number rmax=20   //size of the HRTEM image
number defaultminimumradius=1.8 // This is the minimum value the radius field can take - set it lower
// and it will default to this value
number defaultminimumresolution=100 // This is the smallest resolution - the Resolution field value
// can take - set it below this and it defaults to 100

number defaultpeak=2.7741 // 特征峰的d值

// 把以上数据添加到缓存
 number Positiveg=1  // 标定时g(r)是否非负
 if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Positive g(r)", Positiveg))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Positive g(r)",Positiveg)
	}
	
 number autoContrast=1  // 自动衬度
 if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Auto Contrast", autoContrast))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Auto Contrast",autoContrast)
	}

number flagBG=0   // auto BG or not?
if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Auto BG", flagBG))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Auto BG", flagBG)
	}
	
	
number flagrGr=0   // to display rG(r) image
if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:rG(r)", flagrGr))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:rG(r)", flagrGr)
	}
	

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:R_max", rmax))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:R_max", rmax)
	}
	
if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Number of Elements", noofelements))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Number of Elements", noofelements)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Beam Voltage (kV)", beamvoltage))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Beam Voltage (kV)", beamvoltage)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Use At% (0) or Wt% (1)", defaultuseatorwtpercent))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Use At% (0) or Wt% (1)", defaultuseatorwtpercent)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Damping Method", defaultdampingmethod))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Damping Method", defaultdampingmethod)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Filter Mode (0-2)", defaultfiltermode))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Filter Mode (0-2)", defaultfiltermode)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Profile Position (1-8)", defaultposition))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Profile Position (1-8)", defaultposition)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Profile Aspect Ratio", defaultprofileaspectratio))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Profile Aspect Ratio", defaultprofileaspectratio)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Profile Height", defaultprofileheight))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Profile Height", defaultprofileheight)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Use s (0) or Q (1)", defaultsorqunits))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Use s (0) or Q (1)", defaultsorqunits)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Average Density",defaultaveragedensity))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Average Density",defaultaveragedensity)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Minimum Density",defaultminimumdensity))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Minimum Density",defaultminimumdensity)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Fit to Top channel %",defaultfittotopchannelpcnt))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Fit to Top channel %",defaultfittotopchannelpcnt)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:No of Iterations to Fit",noiterations))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:No of Iterations to Fit",noiterations)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Resolution",defaultresolution))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Resolution",defaultresolution)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Minimum Resolution",defaultminimumresolution))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Minimum Resolution",defaultminimumresolution)
	}

if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Peak",defaultpeak))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Peak",defaultpeak)
	}

	number defaultLeft=190,defaultRight=24  //PDF对话框的位置
	if(!getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Dialog Position:Left",defaultLeft))
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Dialog Position:Left",defaultLeft)
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Dialog Position:Right",defaultRight)
	}
	
/*************************************************************************
Class PDFcalculator - // This section responds to changes made in the dialog box, such as pushing a button.
*************************************************************************/
class PDFcalculator : uiframe
{
// Variables
number hushflag// a flag used to temporarily halt chnaged method responses
image sourceimg // the original sadp from which the profile was created
image imgf2Q // the image which holds the the weighted scattering factor squared ie sum(C x F**2)
image imgfQ2 // the image which holds the sum of the weighted scattering factor squared ie sum(C x F)**2
image errplot // the fraction error plot derived from the RIF plot - note this is global as it has a listener
image scatteringfactorcurve // the scattering factor curve

/*************************************************************************
Add  ROIs
*************************************************************************/
void AddROI(object  self)
{
image img:=getfrontimage()
ImageDisplay imgdisp=img.ImageGetImageDisplay(0)

Result("-----------------------------------------------------------------------------------------\nSingle click: to Add the ROIs\nAlt+click: to read and save ROIs\n-----------------------------------------------------------------------------------------\n")

if(OptionDown()) //读取并保持roi
{
//删除临时缓存
TagGroup tg =GetPersistentTagGroup( )
TagGroup tg1
if(tg.taggroupgettagastaggroup("ElectronDiffraction Tools:PDF:ROIs", tg1))deletepersistentnote("ElectronDiffraction Tools:PDF:ROIs")

number annots= imgdisp.ImageDisplayCountROIS()
for (number i=0; i<annots; i++)
{
ROI currentROI = imgdisp.ImageDisplayGetROI(i)
number l,r,x
currentROI.roigetrange(l,r)
x=r/1000 //为方便排序
setpersistentnumbernote("ElectronDiffraction Tools:PDF:ROIs:"+x+":X",r)
} 

result("ROIs are saved.")
}

else    //添加rois
{
TagGroup tg =GetPersistentTagGroup( )
TagGroup tg1
if(!tg.taggroupgettagastaggroup("ElectronDiffraction Tools:PDF:ROIs", tg1))
{
ShowAlert("No ROIs will be added!",1)
exit(0)
}

number mcpno=tg1.taggroupcounttags()
for (number i=0; i<mcpno-1; i++) //不添加尾roi
{
string mcplabel = tg1.TagGroupGetTaglabel( i )
number value=val(mcplabel)

number x
getpersistentnumbernote("ElectronDiffraction Tools:PDF:ROIs:"+value+":X",x)

roi roi1=createroi()
ROI1.ROISetrange(x,x)
roi1.roisetvolatile(0)
roi1.roisetcolor(1,0,1)
imgdisp.ImageDisplayAddROI( ROI1)
} 

result("ROIs are added.")
}

}

/*************************************************************************
profile to RIF
*************************************************************************/
void Profile2RIF(object  self)
{
image img:=getfrontimage()
ImageDisplay imgdisp=img.ImageGetImageDisplay(0)
img.showimage()

Result("-----------------------------------------------------------------------------------------\nTo get a simple RIF profile.\nAlt+click: to subtract Gaussian background\n-----------------------------------------------------------------------------------------\n")

if(OptionDown())	
{
ChooseMenuItem("ePDF Tools","","BG-Removed-Gauss")
}

else
{
ChooseMenuItem("ePDF Tools","","Profile-RIF")
}

}

/*************************************************************************
从profile提取各图层
*************************************************************************/
void Profile2sliceimage(object  self)
{
Result("-----------------------------------------------------------------------------------------\nGet slices from the profile image\n-----------------------------------------------------------------------------------------\n")

image front:=getfrontimage()
imagedisplay imgdisp=front.imagegetimagedisplay(0)
string imgname=getname(front)
number xsize, ysize
getsize(front, xsize, ysize)

//获取各图层，在result window中显示
number n=1
number noslices=imgdisp.imagedisplaycountslices()
for(number i=0;i<noslices;i++)
{
object sliceid=imgdisp.imagedisplaygetsliceidbyindex(i)
string slicename=imgdisp.imagedisplaygetslicelabelbyid(sliceid)

Result("Slice "+n+": "+slicename+"\n")
n=n+1
}

number i, found
number no=1
GetNumber("No of the slice:\nAll(0)",no,no)
noslices=imgdisp.imagedisplaycountslices()

if(no==0)
{
for(i=0; i<noslices; i++)
	{
		// Get the slice		
		object sliceid=imgdisp.imagedisplaygetsliceidbyindex(i)
		string slicename=imgdisp.imagedisplaygetslicelabelbyid(sliceid)

		// Extract it as an image and display it		
		image slice:=imageclone(front{i})// clone the actual slice to pick up all the image calibrations, tags etc
		showimage(slice)
		
		// label and display the component spectrum		
		setname(slice, imgname+" ("+slicename+")")
		documentwindow docwin=getdocumentwindow(0)
		docwin.windowsetframeposition(found*30, found*30)
		imagecopycalibrationfrom(slice, front)
		found=found+1
	}
}

else
{
for(i=0; i<noslices; i++)
	{
	
	if(i==no-1)
	{
		// Get the slice		
		object sliceid=imgdisp.imagedisplaygetsliceidbyindex(i)
		string slicename=imgdisp.imagedisplaygetslicelabelbyid(sliceid)

		// Extract it as an image and display it		
		image slice:=imageclone(front{i})// clone the actual slice to pick up all the image calibrations, tags etc
		showimage(slice)
		
		// label and display the component spectrum		
		setname(slice, imgname+" ("+slicename+")")
		documentwindow docwin=getdocumentwindow(0)
		docwin.windowsetframeposition(found*30, found*30)
		imagecopycalibrationfrom(slice, front)
		
		break
	}
		found=found+1
	}

}
}

/*************************************************************************
设置PDF参数
*************************************************************************/
void SetPreferences(object  self)
{
Result("-----------------------------------------------------------------------------------------\nSet default parameters of PDF calculator\n-----------------------------------------------------------------------------------------\n")

TagGroup tg = GetPersistentTagGroup( )  
TagGroup mcptags
tg.taggroupgettagastaggroup("ElectronDiffraction Tools:PDF:Default Values", mcptags)
mcptags.TagGroupOpenBrowserWindow( 0 )
}

/*************************************************************************
对Profile进行平滑
*************************************************************************/
void smoothIq(object  self)
{
image img:=GetFrontImage()
img.showimage()
	
ChooseMenuItem("ePDF Tools","","I(Q)-Smooth")
}

/*************************************************************************
图层替换
*************************************************************************/
void ReplaceIq(object  self)
{
Result("-----------------------------------------------------------------------------------------\nReplace button: \n1) Copy (Ctrl+C) and paster (Ctrl+V) the desired profile to the target profile;\n2) Run this button and set the desired slice number.\n-----------------------------------------------------------------------------------------\n")

image img:=getfrontimage()
imagedisplay imgdisp=imagegetimagedisplay(img, 0)

//获取各图层
number noslices=imgdisp.imagedisplaycountslices()
for(number i=0;i<noslices;i++)
{
object sliceid=imgdisp.imagedisplaygetsliceidbyindex(i)
string slicename=imgdisp.imagedisplaygetslicelabelbyid(sliceid)

Result("Slice "+i+": "+slicename+"\n")
}

number sliceno=0
GetNumber("Which slice will be replaced?",sliceno,sliceno)
img{sliceno}=img{noslices-1}  //最后图层替换目标图层

//删除粘贴的图层
object sliceid=imgdisp.imagedisplaygetsliceidbyindex(noslices-1)
imgdisp.imagedisplaydeleteslicewithid(sliceid)
}

/*************************************************************************
F(Q)->I(Q)和SAED
*************************************************************************/
void Fq2Iqbutton(object  self)
{
image img:=getfrontimage()
imagedisplay imgdisp=imagegetimagedisplay(img, 0)
string imgname=img.getname()

image Fq=imgdisp{"Damped"}
image fq2=imgdisp{"<f2>"}
image F2q=imgdisp{"<f>2"}

number N,c,xscale
xscale=img.ImageGetDimensionScale(0)
Getnumbernote(img,"ElectronDiffraction Tools:PDF:Parameters:N",N) 
Getnumbernote(img,"ElectronDiffraction Tools:PDF:Parameters:C",c) 
//image Iq=2*Pi()*Fq*N*F2q/(icol*xscale)+N*Fq2-c   //转为I(s)
image Iq=Fq*N*F2q/(icol*xscale)+N*Fq2-c //转为I(Q)
Iq.SetPixel(0,0,0)   //第一像素为0

Iq.ShowImage()
Iq.SetName("I(Q) - "+imgname)
Iq.ImageCopyCalibrationFrom(img)
SetStringNote(Iq,"ElectronDiffraction Tools:PDF:Image name",imgname)  
setstringnote(Iq,"Radial Distribution Function","Calculated Intensity profile")
imagesetIntensityUnitString(Iq,"I(Q)")

//转为SAED，注意单位为1/nm
number xsize,ysize
img.GetSize(xsize,ysize)
image SAED=RealImage("",4,2*xsize,2*xsize)
SAED=0

SAED = Iq[iradius,0] 
SAED.ShowImage()

xscale=10*xscale/(2*Pi())
SAED.setorigin(xsize+0.5,xsize+0.5)
SAED.SetScale(xscale,xscale)
SAED.SetUnitString("1/nm")
SAED.SetName("Calculated SAED - "+imgname)
SetStringNote(SAED,"ElectronDiffraction Tools:PDF:Image name",imgname)  

SAED.SetNumberNote("ElectronDiffraction Tools:Center:X",xsize+0.5)
SAED.SetNumberNote("ElectronDiffraction Tools:Center:Y",xsize+0.5)
}

/*************************************************************************
Live Gaussian Fit
*************************************************************************/
void LiveGaussian(object  self)
{
image img:=getfrontimage()
ImageDisplay imgdisp=img.ImageGetImageDisplay(0)
img.showimage()
	
ChooseMenuItem("ePDF Tools","","Live-Gaussian")
}

/*************************************************************************
Close button
*************************************************************************/
void CloseImages(object  self)
{
Result("-----------------------------------------------------------------------------------------\nClose button: \n1) Close all images without saving;\n2) ALT+ Close images with saving dialog. \n-----------------------------------------------------------------------------------------\n")

Number kWINDOWTYPE_IMAGEWINDOW = 5
number numberDocs = CountDocumentWindowsOfType(kWINDOWTYPE_IMAGEWINDOW)

for(number  i = 0; i < numberDocs; ++ i )
{
ImageDocument imgDoc = GetImageDocument( 0 )
image img:=getfrontimage() 

if(optiondown()) img.closeimage()   //ALT+，保存关闭
else img.deleteimage()   //不保存，直接关闭
}
}

/*************************************************************************
BG button
*************************************************************************/
void BGRemoved(object  self)
{
image img:=getfrontimage()
ImageDisplay imgdisp=img.ImageGetImageDisplay(0)
img.showimage()

number norois=imgdisp.imagedisplaycountrois()
	if(norois==0)
	{
	ShowAlert("Drag an ROI first!",2)
	break	
	}
	
ChooseMenuItem("ePDF Tools","","BG-Removed")
}


/*************************************************************************
I(Q) to F(Q) button
*************************************************************************/
void Iq2Fq(object  self)
{
Result("-----------------------------------------------------------------------------------------\nSinlge click: to get F(Q)\nCtrl+ to set the value of N.\n\n")

if(controlDown())
{
for(number k=1;k<3;k++)
{
image img:=GetFrontImage()
ImageDisplay imgdisp=img.ImageGetImageDisplay(0)
object sliceid=imgdisp.ImageDisplayGetSliceIDByIndex(0)

number N
Getnumbernote(img,"ElectronDiffraction Tools:PDF:Parameters:N",N)

number x,cursorstate,min1,max1
lineplotimagedisplaygetcursorstate(imgdisp,sliceid,x,cursorstate)

img[0,x-5,1,x+5].minmax(min1,max1)

if(abs(min1)>abs(max1))
{
SetPersistentNumberNote("ElectronDiffraction Tools:Temp:"+k+":N",N)
SetPersistentNumberNote("ElectronDiffraction Tools:Temp:"+k+":Y",min1)
//Result(n+"\t"+min1+"\n")
}

else
{
SetPersistentNumberNote("ElectronDiffraction Tools:Temp:"+k+":N",N)
SetPersistentNumberNote("ElectronDiffraction Tools:Temp:"+k+":Y",max1)
//Result(N+"\t"+max1+"\n")
}
 
img.deleteimage()		
}

//计算该峰正好为0时的x值
number x,x1,y1,x2,y2
getPersistentNumberNote("ElectronDiffraction Tools:Temp:1:N",x1)
getPersistentNumberNote("ElectronDiffraction Tools:Temp:1:Y",y1)
getPersistentNumberNote("ElectronDiffraction Tools:Temp:2:N",x2)
getPersistentNumberNote("ElectronDiffraction Tools:Temp:2:Y",y2)

result("("+x1+","+y1+")\t("+x2+","+y2+")\n")
x=(0-y1)/(y2-y1)*(x2-x1)+x1
Result("N roi will be set to: "+x+"\n")

//找到intensity profile
number nodocs=countdocumentwindowsoftype(5)
image profile
imagedocument imgdoc
string imgname

number N
for(number i=0; i<nodocs; i++)
	{
		imagedocument imgdoc=getimagedocument(i)
		image tempimg:=imgdoc.imagedocumentgetimage(0)
		ImageDisplay profiledisp=tempimg.ImageGetImageDisplay(0)

		TagGroup tg =tempimg.ImageGetTagGroup( )  //缓存中是否有文件名，
		if(tg.taggroupdoestagexist("Radial Distribution Function"))
			{
			getStringNote(tempimg,"Radial Distribution Function",imgname)
						
			if(imgname=="Intensity Profile")
			{
			number N_factor,n0,l,r
			getnumbernote(tempimg,"ElectronDiffraction Tools:PDF:Parameters:N factor",N_factor) 
			getnumbernote(tempimg,"ElectronDiffraction Tools:PDF:Parameters:N0",n0)
			
			//寻找N roi
			number annots=profiledisp.ImageDisplayCountROIS()
			for(number j=0;j<annots;j++)
			{
			ROI currentROI =profiledisp.ImageDisplayGetROI( j )
			string thischar=currentROI .ROIGetLabel()
			string roilabel=left(thischar,1)
			
			if(roilabel=="N")
				{
				x=x/(N_factor*n0)
				currentROI .roisetrange(x,x)
				}
			}

			}
			}
}
			
deletepersistentnote("ElectronDiffraction Tools:Temp")	
}

else
{
image img:=getfrontimage()
ImageDisplay imgdisp=img.ImageGetImageDisplay(0)
img.showimage()

taggroup dspacingtags
taggroup ptags=img.ImageGetTagGroup()
if(ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:PDF:Parameters", dspacingtags)) 
{
ChooseMenuItem("ePDF Tools","","Iq-Fq")
}
	
else break
}

}


//手动输入定义中心
void SetCenter(object  self)
{
Result("Set the new center of the image!\n")

number sizex, sizey, centery, centerx, scalex, scaley,x,y,centerx0,centery0
image img:=getfrontimage()
imagedisplay imgdisp = img.ImageGetImageDisplay(0)

//读取已有中心
img.getorigin(centerx0,centery0)

if (!GetNumber( "Set New Center X (pix)\nX0= "+centerx0, centerx0, centerx0 ))exit(0)
if (!GetNumber( "Set New Center Y (pix)\nY0= "+centery0, centery0, centery0 ))exit(0)

setorigin(img,centerx0,centery0)
setNumberNote(img,"ElectronDiffraction Tools:Center:X", centerx0 )
setNumberNote(img,"ElectronDiffraction Tools:Center:Y", centery0)

Result("New center is: ("+centerx0+", "+centery0+")\n")

ROI roi_1 = NewROI()
ROISetRectangle(roi_1, centery0-5, centerx0-5, centery0+5, centerx0+5 ) 
roi_1.roisetvolatile(0)
roi_1.roisetcolor(1,0,0)
imgdisp.ImageDisplayAddROI(roi_1)
}

//用圆定义中心
void SetCenterFromOval(object  self)
{
Result("Set center using the oval annotation!\n")
image img:=getfrontimage()
ImageDisplay imgdisp=img.imagegetimagedisplay(0)

number t,l,b,r,x0,y0,Radius
number annots=imgdisp.componentcountchildrenoftype(6)  
if(annots!=0)  //如果是圆,圆心为中心
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.ComponentGetBoundingRect(t,l,b,r)

x0=l+(r-l)/2
y0=t+(b-t)/2

//设置中心
if (!GetNumber( "Set New Center X (pix)\nX0= "+x0, x0, x0 ))exit(0)
if (!GetNumber( "Set New Center Y (pix)\nY0= "+y0, y0, y0 ))exit(0)

setorigin(img,x0,y0)
setNumberNote(img,"ElectronDiffraction Tools:Center:X", x0 )
setNumberNote(img,"ElectronDiffraction Tools:Center:Y", y0)

Result("New center is: ("+x0+", "+y0+")\n")

ROI roi_1 = NewROI()
ROISetRectangle(roi_1, y0-5, x0-5, y0+5, x0+5 ) 
roi_1.roisetvolatile(0)
roi_1.roisetcolor(1,0,0)
imgdisp.ImageDisplayAddROI(roi_1)
}


//如果是roi
number roiannots=imgdisp.ImageDisplayCountROIS()
if(roiannots!=0)  //如果是圆,圆心为中心
{
ROI currentROI = imgdisp.ImageDisplayGetROI( 0 )
currentROI.roigetrectangle(t,l,b,r)

x0=l+(r-l)/2
y0=t+(b-t)/2

imgdisp.imagedisplaydeleteROI(currentROI)

//设置中心
if (!GetNumber( "Set New Center X (pix)\nX0= "+x0, x0, x0 ))exit(0)
if (!GetNumber( "Set New Center Y (pix)\nY0= "+y0, y0, y0 ))exit(0)

setorigin(img,x0,y0)
setNumberNote(img,"ElectronDiffraction Tools:Center:X", x0 )
setNumberNote(img,"ElectronDiffraction Tools:Center:Y", y0)

Result("New center is: ("+x0+", "+y0+")\n")

ROI roi_1 = NewROI()
ROISetRectangle(roi_1, y0-5, x0-5, y0+5, x0+5 ) 
roi_1.roisetvolatile(0)
roi_1.roisetcolor(1,0,0)
imgdisp.ImageDisplayAddROI(roi_1)
}

}

//对图像进行标定
void ImgCalibration(object  self)
{
Result("Calibrate the image using Line ROI, Line, or user-defined two circles!\n")
image img:=getfrontimage()
ImageDisplay imgdisp=img.imagegetimagedisplay(0)

string units=img.ImageGetDimensionUnitString(0)
number x1,y1,x2,y2,xscale,xorigin,yorigin
img.getorigin(xorigin,yorigin)

number annots=imgdisp.componentcountchildrenoftype(2)  
if(annots!=0)  //如果是线
{
component annotid=imgdisp.componentgetnthchildoftype(2,0)
annotid.componentgetcontrolpoint(4,x1,y1)  //线的起点
annotid.componentgetcontrolpoint(5,x2,y2)  //线的终点
number R=sqrt((y2-y1)**2+(x2-x1)**2)

number stringno=1,dspacing=1
GetNumber("Set Unit string\nnm (1), 1/nm (2)",stringno,stringno)
GetNumber("Set d-spacing (A)",dspacing,dspacing)

if(stringno==1)  //图像
{
xscale=0.1*dspacing/R
ImageSetDimensionCalibration(img,0,xorigin,xscale,"nm",1)
ImageSetDimensionCalibration(img,1,yorigin,xscale,"nm",1)

result("Two points of the line: ("+x1+", "+y1+"), ("+x2+", "+y2+")\nR= "+R+",  d= "+dspacing+" A,  Scale= "+xscale+" nm\n")
}

else if(stringno==2) //衍射
{
xscale=10/(dspacing*R)
ImageSetDimensionCalibration(img,0,xorigin,xscale,"1/nm",1)
ImageSetDimensionCalibration(img,1,yorigin,xscale,"1/nm",1)

result("Two points of the line: ("+x1+", "+y1+"), ("+x2+", "+y2+")\nR= "+R+",  d= "+dspacing+" A,  Scale= "+xscale+" 1/nm\n")
}
}

annots=imgdisp.ImageDisplayCountROIs()
if(annots!=0)  //如果是ROI
{
ROI currentROI = imgdisp.ImageDisplayGetROI( 0 )
if(ROIIsLine(currentROI))
{
currentROI.ROIGetLine(x1,y1,x2,y2)
number R=sqrt((y2-y1)**2+(x2-x1)**2)

number stringno=1,dspacing=1
GetNumber("Set Unit string\nnm (1), 1/nm (2)",stringno,stringno)
GetNumber("Set d-spacing (A)",dspacing,dspacing)

if(stringno==1)  //图像
{
xscale=0.1*dspacing/R
ImageSetDimensionCalibration(img,0,xorigin,xscale,"nm",1)
ImageSetDimensionCalibration(img,1,yorigin,xscale,"nm",1)

result("Two points of the line: ("+x1+", "+y1+"), ("+x2+", "+y2+")\nR= "+R+",  d= "+dspacing+" A,  Scale= "+xscale+" nm\n")
}

else if(stringno==2) //衍射
{
xscale=10/(dspacing*R)
ImageSetDimensionCalibration(img,0,xorigin,xscale,"1/nm",1)
ImageSetDimensionCalibration(img,1,yorigin,xscale,"1/nm",1)

result("Two points of the line: ("+x1+", "+y1+"), ("+x2+", "+y2+")\nR= "+R+",  d= "+dspacing+" A,  Scale= "+xscale+" 1/nm\n")
}
}
}

annots=imgdisp.componentcountchildrenoftype(6)  //如果是圆
if(annots!=0)  
{
for (number i=1; i<annots+1; i++)  //读取圆心
{
number t,l,b,r	
component annotid=imgdisp.componentgetnthchildoftype(6,i-1)
annotid.ComponentGetBoundingRect(t,l,b,r)

number x0=l+(r-l)/2
number y0=t+(b-t)/2
SetNumberNote(img,"ElectronDiffraction Tools:Temp:X"+i+"", x0 )
SetNumberNote(img,"ElectronDiffraction Tools:Temp:Y"+i+"", y0 )
}

GetNumberNote(img,"ElectronDiffraction Tools:Temp:X1", x1 )
GetNumberNote(img,"ElectronDiffraction Tools:Temp:Y1", y1 )
GetNumberNote(img,"ElectronDiffraction Tools:Temp:X2", x2 )
GetNumberNote(img,"ElectronDiffraction Tools:Temp:Y2", y2 )
number R=sqrt((y2-y1)**2+(x2-x1)**2)

number stringno=1,dspacing=1
GetNumber("Set Unit string\nnm (1), 1/nm (2)",stringno,stringno)
GetNumber("Set d-spacing (A)",dspacing,dspacing)

if(stringno==1)  //图像
{
xscale=0.1*dspacing/R
ImageSetDimensionCalibration(img,0,xorigin,xscale,"nm",1)
ImageSetDimensionCalibration(img,1,yorigin,xscale,"nm",1)

result("Two points of the line: ("+x1+", "+y1+"), ("+x2+", "+y2+")\nR= "+R+",  d= "+dspacing+" A,  Scale= "+xscale+" nm\n")
}

else if(stringno==2) //衍射
{
xscale=10/(dspacing*R)
ImageSetDimensionCalibration(img,0,xorigin,xscale,"1/nm",1)
ImageSetDimensionCalibration(img,1,yorigin,xscale,"1/nm",1)

result("Two points of the line: ("+x1+", "+y1+"), ("+x2+", "+y2+")\nR= "+R+",  d= "+dspacing+" A,  Scale= "+xscale+" 1/nm\n")
}
}
}


//函数：获取图像中的标定参数，以便origin中画图
void GetCalibration(object  self)
{
Result("-----------------------------------------------------------------------------------------\nSinlge click: to get the image scale\nAlt+ to set the image scale.\n\n")
if(OptionDown())   //s-Q互换
{
image img:=GetFrontImage()
number xorigin,xscale,yorigin,yscale,xsize,ysize,s_xscale,Q_xscale,s_yscale,Q_yscale
string xunitstring,yunitstring

img.GetSize(xsize,ysize)
If(ysize<=1)
{
img.ImageGetDimensionCalibration(0,xorigin,xscale,xunitstring,1)
Result("\nOrigin: "+xorigin+";    Scale: "+xscale+";    Units: "+xunitstring+"\n")

s_xscale=xscale/(2*Pi()*0.1)
Q_xscale=(2*Pi()*0.1)*xscale

number flag=1
if(xunitstring=="1/nm")flag=1
if(xunitstring=="Q (1/A)")flag=2
GetNumber("0: xscale="+xscale+", "+xunitstring+"\n------------------------------------------------\n1: xscale="+Q_xscale+", Q (1/A)"+"\n2: xscale="+s_xscale+", 1/nm",flag,flag)
if(flag==1)
{
img.ImageSetDimensionCalibration(0,xorigin,Q_xscale,"Q (1/A)",1)
Result("Origin: "+xorigin+";    Scale: "+Q_xscale+";    Units: Q (1/A)\n")
}

if(flag==2)
{
img.ImageSetDimensionCalibration(0,xorigin,s_xscale,"1/nm",1)
Result("Origin: "+xorigin+";    Scale: "+s_xscale+";    Units: 1/nm\n")
}
}


If(ysize>1)
{
img.ImageGetDimensionCalibration(0,xorigin,xscale,xunitstring,1)
img.ImageGetDimensionCalibration(1,yorigin,yscale,yunitstring,1)
Result("\nOrigin: ("+xorigin+", "+yorigin+");    Scale: ("+xscale+", "+yscale+");    Units: ("+xunitstring+", "+yunitstring+")\n")

s_xscale=xscale/(2*Pi()*0.1)
Q_xscale=(2*Pi()*0.1)*xscale

number flag=1
if(xunitstring=="1/nm")flag=1
if(xunitstring=="Q (1/A)")flag=2
GetNumber("0: xscale="+xscale+", "+xunitstring+"\n------------------------------------------------\n1: xscale="+Q_xscale+", Q (1/A)"+"\n2: xscale="+s_xscale+", 1/nm",flag,flag)
if(flag==1)
{
img.ImageSetDimensionCalibration(0,xorigin,Q_xscale,"Q (1/A)",1)
img.ImageSetDimensionCalibration(1,yorigin,Q_xscale,"Q (1/A)",1)
Result("Origin: ("+xorigin+", "+yorigin+");    Scale: ("+Q_xscale+", "+Q_xscale+");    Units: (Q (1/A), Q (1/A))\n")
}

if(flag==2)
{
img.ImageSetDimensionCalibration(0,xorigin,s_xscale,"1/nm",1)
img.ImageSetDimensionCalibration(1,yorigin,s_xscale,"1/nm",1)
Result("Origin: ("+xorigin+", "+yorigin+");    Scale: ("+s_xscale+", "+s_xscale+");    Units: (1/nm, 1/nm)\n")
}
}

}


//获取单位信息
else
{
image img:=GetFrontImage()
number xorigin,xscale,yorigin,yscale,xsize,ysize
string xunitstring,yunitstring

img.GetSize(xsize,ysize)
If(ysize<=1)
{
img.ImageGetDimensionCalibration(0,xorigin,xscale,xunitstring,1)
Result("\nOrigin: "+xorigin+";    Scale: "+xscale+";    Units: "+xunitstring+"\n")

Result("colint x0:=-"+xorigin*xscale+" inc:="+xscale+"\n")   //输出origin中的interval数据
}

If(ysize>1)
{
img.ImageGetDimensionCalibration(0,xorigin,xscale,xunitstring,1)
img.ImageGetDimensionCalibration(1,yorigin,yscale,yunitstring,1)
Result("\nOrigin: ("+xorigin+", "+yorigin+");    Scale: ("+xscale+", "+yscale+");    Units: ("+xunitstring+", "+yunitstring+")\n")
}
}

}




//保存G(r)文件
void saveGr(object  self)
{
Image profile:=GetFrontImage()
number xsize,ysize,xscale,yscale,xorigin,yorigin
xsize=profile.imagegetdimensionsize(0)
xscale=profile.imagegetdimensionscale(0)
xorigin=profile.imagegetdimensionorigin(0)

number startvector=xorigin,endvector=xsize*xscale+xorigin
number xend=endvector-startvector

string filename, othername
string imagename=profile.getname()

If (!SaveAsDialog("Save text file as", GetApplicationDirectory(2,0) +imagename+ ".txt", filename)) Exit(0)

number fileID = CreateFileForWriting(filename)
Result("Exported file:  "+imagename+"\n")
number yval, Q

for(number i=0; i<xsize; i++)
	{
		yval=getpixel(profile, i, 0)
		number xcal=xorigin*xscale+(i*xscale)

		if(xcal>=0&&xcal<=xend)
		{
        Q=1*(xcal+startvector)
		WriteFile(fileID, Q+"        "+yval+"        "+0+"        "+0+"\n")
		}
	}
CloseFile(fileID)

profile.closeimage()
}



/*************************************************************************
Zoom button
*************************************************************************/
void ZoomImage(object self)
{
Result("-----------------------------------------------------------------------------------------\nZoom button: \n1) Zoom image;\n2) ALT+  Scale image.\n-----------------------------------------------------------------------------------------\n")

if(OptionDown())
{
image img:=getfrontimage()
string imgname=img.GetName()

img.ShowImage()
ChooseMenuItem("Process","","Scale...")

image scaledImg:=getfrontimage()
scaledImg.SetName(imgname)
}

else
{
Image myImage := GetFrontImage()
imagedisplay imgdisp=myimage.imagegetimagedisplay(0)

number screenWidth, screenHeight, imagewidth, imageheight, zoomfactor
number xPos, yPos, zoom

GetSize( myImage, imagewidth, imageheight)
GetScreenSize( screenWidth, screenHeight )

number workingareawidth = screenwidth-250
number workingareaheight = screenheight-250
number scalewidth = workingareawidth/imagewidth
number scaleheight = workingareaheight/imageheight

 If(scalewidth>scaleheight)
{
zoomfactor = scaleheight
}

else
{
zoomfactor = scalewidth
}

SetWindowPosition( myImage, 190,30 )

if(imageheight<=1)
{
SetWindowSize( myImage, screenwidth-800, screenheight-250)

number min1,max1
myImage.minmax(min1,max1)
imgdisp.LinePlotImageDisplaySetDoAutoSurvey( 0, 0 ) 
imgdisp.LinePlotImageDisplaySetContrastLimits(1.1*min1,1.1*max1)  
}
else
{
SetWindowSize( myImage, imagewidth*zoomfactor, imageheight*zoomfactor)
setzoom(myimage, zoomfactor)
}

SetImagePositionWithinWindow( myImage, 0, 0 )
}
}


/*************************************************************************
ROI image button
*************************************************************************/
void ROIImage(object self)
{
Result("-----------------------------------------------------------------------------------------\nROIImg button: \n1) Get ROI image;\n2) ALT+  Get the largest 2^n image around center.\n-----------------------------------------------------------------------------------------\n")

//得到最大的2n image
if(OptionDown())
{
image RIF:=getfrontimage()
number minx,miny,minval,x0,y0,xsize,ysize,xscale,yscale
RIF.getsize(xsize,ysize)
rif.getscale(xscale,yscale)

string imgname=rif.GetName()
GetNumberNote(rif,"ElectronDiffraction Tools:Center:X", x0 )
GetNumberNote(rif,"ElectronDiffraction Tools:Center:Y", y0 )

minx=min(xsize-x0,x0)
miny=min(ysize-y0,y0)
minval=min(minx,miny)

image img=RIF[y0-minval,x0-minval,y0+minval,x0+minval]
img.ImageCopyCalibrationFrom(rif)
img.getsize(xsize,ysize)

number no=0
for(number i=0;i<20;i++)
{
if(round(xsize)>=2)
{
xsize=xsize/2
no=no+1
}
}

img.getsize(xsize,ysize)
number scale=2**no/xsize

image new=realimage("the rescaling image",4,xsize*scale,ysize*scale)
new=warp(img,icol/scale,irow/scale)
new.showimage()
new.ImageCopyCalibrationFrom(rif)
new.setscale(xscale/scale,xscale/scale)

new.GetSize(xsize,ysize)
new.SetOrigin(xsize/2+0.5,ysize/2+0.5)
setNumberNote(new,"ElectronDiffraction Tools:Center:X",xsize/2+0.5 )
setNumberNote(new,"ElectronDiffraction Tools:Center:Y",ysize/2+0.5)
new.SetName("ROI image-"+imgname)
}

//得到roi image
else
{
image front:=getfrontimage()
imagedisplay imgdisp=front.imagegetimagedisplay(0)
number annots= imgdisp.ImageDisplayCountROIS()

string imgname
getname(front, imgname)

number t,l,b,r,xsize,ysize,xscale
front.getsize(xsize,ysize)
front.getselection(t,l,b,r)

if(annots==0||annots>1)
{
if(ysize<=1)
{
xscale=front.imagegetdimensionscale(0)

r=19  //右侧的长度
getnumber("Right side of analysis region",r,r)
r=r/xscale

image cropped=realimage("", 4,r, ysize)
cropped=front[0,0,1,r]
//cropped=cropped-cropped.min()

showimage(cropped)
imagecopycalibrationfrom(cropped, front)
setname(cropped, imgname)

TagGroup tag11,tag2,tag22
TagGroup tag1=front.imagegettaggroup()
if(tag1.taggroupdoestagexist("ElectronDiffraction Tools"))
	{
		tag1.taggroupgettagastaggroup("ElectronDiffraction Tools",tag11)
	}

setnumbernote(cropped,"ElectronDiffraction Tools:X",1)
tag2=cropped.imagegettaggroup()
tag2.taggroupgettagastaggroup("ElectronDiffraction Tools",tag22)
tag22.TagGroupCopyTagsFrom(tag11)
deletenote(cropped,"ElectronDiffraction Tools:X")

tag1=front.imagegettaggroup()
if(tag1.taggroupdoestagexist("Radial Distribution Function"))
{
string imglabel=getstringnote(front,"Radial Distribution Function")
setstringnote(cropped,"Radial Distribution Function",imglabel)
}

imagedisplay newdisp=cropped.ImageGetImageDisplay(0)
number noslices=imgdisp.imagedisplaycountslices()
for(number i=0;i<noslices;i++)
{
object sliceid=imgdisp.imagedisplaygetsliceidbyindex(i)
string slicename=imgdisp.imagedisplaygetslicelabelbyid(sliceid)

if(i==0)
{
newdisp.imagedisplaysetslicelabelbyID(sliceid, slicename)
}

else
{
image new=front{i}
image new1=new[0,0,1,r]
newdisp.imagedisplayaddimage(new1, slicename)
newdisp.LinePlotImageDisplaySetLegendShown(1)
}
}

}
else
{
image cropped=front
showimage(cropped)
imagecopycalibrationfrom(cropped, front)
setname(cropped, imgname)
}
}

if(annots==1)
{
if(ysize<=1)
{
image cropped=realimage("", 4,r, ysize)
front.setselection(0,0,1,r)
cropped=front[]
showimage(cropped)
imagecopycalibrationfrom(cropped, front)
setname(cropped, imgname)

TagGroup tag11,tag2,tag22
TagGroup tag1=front.imagegettaggroup()
if(tag1.taggroupdoestagexist("ElectronDiffraction Tools"))
	{
		tag1.taggroupgettagastaggroup("ElectronDiffraction Tools",tag11)
	}

setnumbernote(cropped,"ElectronDiffraction Tools:X",1)
tag2=cropped.imagegettaggroup()
tag2.taggroupgettagastaggroup("ElectronDiffraction Tools",tag22)
tag22.TagGroupCopyTagsFrom(tag11)
deletenote(cropped,"ElectronDiffraction Tools:X")

tag1=front.imagegettaggroup()
if(tag1.taggroupdoestagexist("Radial Distribution Function"))
{
string imglabel=getstringnote(front,"Radial Distribution Function")
setstringnote(cropped,"Radial Distribution Function",imglabel)
}

imagedisplay newdisp=cropped.ImageGetImageDisplay(0)
number noslices=imgdisp.imagedisplaycountslices()
for(number i=0;i<noslices;i++)
{
object sliceid=imgdisp.imagedisplaygetsliceidbyindex(i)
string slicename=imgdisp.imagedisplaygetslicelabelbyid(sliceid)

if(i==0)
{
newdisp.imagedisplaysetslicelabelbyID(sliceid, slicename)
}

else
{
image new=front{i}
image new1=new[0,0,1,r]
newdisp.imagedisplayaddimage(new1, slicename)
newdisp.LinePlotImageDisplaySetLegendShown(1)
}
}
}

else
{
image cropped=front[]
showimage(cropped)
imagecopycalibrationfrom(cropped, front)
setname(cropped, imgname)
}
}
}

}

/*************************************************************************
Center button
*************************************************************************/
void Center(object  self)
{
Result("-----------------------------------------------------------------------------------------\nSet center button: \n1) SPACE bar to find the center of the 2D pattern by Gaussian fitting;\n2) ALT+ to set the center by the oval component or the origin.\n-----------------------------------------------------------------------------------------\n")

if(OptionDown())
{
image img:=getfrontimage()
ImageDisplay imgdisp=img.imagegetimagedisplay(0)

number t,l,b,r,x0,y0,Radius
number annots=imgdisp.componentcountchildrenoftype(6)  
if(annots!=0)  //如果是圆,圆心为中心
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.ComponentGetBoundingRect(t,l,b,r)

x0=l+(r-l)/2
y0=t+(b-t)/2

//设置中心
if (!GetNumber( "Set New Center X (pix)\nX0= "+x0, x0, x0 ))exit(0)
if (!GetNumber( "Set New Center Y (pix)\nY0= "+y0, y0, y0 ))exit(0)

setorigin(img,x0,y0)
setNumberNote(img,"ElectronDiffraction Tools:Center:X", x0 )
setNumberNote(img,"ElectronDiffraction Tools:Center:Y", y0)

Result("New center is: ("+x0+", "+y0+")\n")

ROI roi_1 = NewROI()
ROISetRectangle(roi_1, y0-5, x0-5, y0+5, x0+5 ) 
roi_1.roisetvolatile(0)
roi_1.roisetcolor(1,0,0)
imgdisp.ImageDisplayAddROI(roi_1)
}

if(annots==0)  //如果没有圆
{
number sizex, sizey, centery, centerx, scalex, scaley,x,y,centerx0,centery0

//读取已有中心
img.getorigin(centerx0,centery0)
if(centerx0==0||centery0==0)
{
img.GetSize(sizex,sizey)
centerx0 =0.5*sizex+0.5
centery0 =0.5*sizey+0.5
}

if (!GetNumber( "Set New Center X (pix)\nX0= "+centerx0, centerx0, centerx0 ))exit(0)
if (!GetNumber( "Set New Center Y (pix)\nY0= "+centery0, centery0, centery0 ))exit(0)

setorigin(img,centerx0,centery0)
setNumberNote(img,"ElectronDiffraction Tools:Center:X", centerx0 )
setNumberNote(img,"ElectronDiffraction Tools:Center:Y", centery0)

Result("New center is: ("+centerx0+", "+centery0+")\n")

ROI roi_1 = NewROI()
ROISetRectangle(roi_1, centery0-5, centerx0-5, centery0+5, centerx0+5 ) 
roi_1.roisetvolatile(0)
roi_1.roisetcolor(1,0,0)
imgdisp.ImageDisplayAddROI(roi_1)
}
}
	
else
{	
image  img:=getfrontimage()
ImageDisplay imgdisp = img.ImageGetImageDisplay(0)
ImageDocument  imgdoc=getfrontimagedocument()
documentwindow imgwin=ImageDocumentGetWindow(imgdoc)
string imgname=img.getname()

number centrex, centrey, x1,y1,x2,y2,x,y,xsize,ysize,xscale,yscale, top,left,bottom,right

img.getsize(xsize,ysize)

number width, times=1
getpersistentnumbernote("ElectronDiffraction Tools:Fit:Gaussian width",width)

number k=1,m=1,noclick=1,n=1 //m用于while循环，noclick用于if取点，n用于对称点for循环
while(2>m)
{
number keypress=getkey()
		
if(keypress==32) //鼠标取点
{
 getwindowsize(img, xsize, ysize)
 getscale(img,xscale,yscale)
 
number mouse_win_x, mouse_win_y
Number img_view_t, img_view_l, img_view_b, img_view_r
Number v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y
Number img_win_t, img_win_l, img_win_b, img_win_r
Number i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y
Number i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y
Number mouse_img_x, mouse_img_y

windowgetmouseposition(imgwin, mouse_win_x, mouse_win_y)
imgdisp.ImageDisplayGetImageRectInView( img_view_t, img_view_l, img_view_b, img_view_r )
imgdoc.ImageDocumentGetViewToWindowTransform( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y )
ObjectTransformTransformRect( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                            , img_view_t, img_view_l, img_view_b, img_view_r\
                            , img_win_t, img_win_l, img_win_b, img_win_r );

imgdisp.ComponentGetChildToViewTransform( i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y )
		
ObjectTransformCompose( v2w_off_x, v2w_off_y, v2w_scale_x, v2w_scale_y\
                      , i2v_off_x, i2v_off_y, i2v_scale_x, i2v_scale_y\
                      , i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y )
		
ObjectTransformUntransformPoint( i2w_off_x, i2w_off_y, i2w_scale_x, i2w_scale_y\
                               , mouse_win_x, mouse_win_y, mouse_img_x, mouse_img_y );

centrex=mouse_img_x
centrey=mouse_img_y

 //-------------开始用高斯拟合确定该点中心------------
for(n=1;n<2;n++) 
{
//获取选区图
img := GetFrontImage()
img.setselection(centrey-width,centrex-width,centrey+width,centrex+width)
img.getselection(top, left, bottom, right)

//平滑
image temp=img[]
temp=medianfilter(temp,3,0)

temp.getsize(xsize,ysize)

times=1/times
image new := RealImage("target",4,xsize/times,ysize/times)
new = temp[].warp(icol*times,irow*times)
temp.deleteimage()
centrex=centrex/times
centrey=centrey/times

//归一化
number minval, maxval
minmax(new, minval, maxval)
number range=maxval-minval
new=(new-minval)/range
new.getsize(xsize,ysize)
centrex=xsize/2
centrey=ysize/2

number subareapcnt=50 // central percent of area tested for goodness of fit (25-75%)
number fitcentrex, fitcentrey
number i, chisqrd
number intensity=1 // set to 1 for normalised images
number initialsigma
number iterations=50 // The number of random steps taken to refine the value

// Note the processing time (for a ca 50 x 50 blob image) is a function of iterations
// multiplied by refinementloops (below). Values of 50 iterations x 5 refinements takes about 0.4s.
number sigma=xsize/4 // Sigma of the Gaussian - estimate an initial value before refinement
initialsigma=sigma
number minchisqrd=1e99

// Loop to refine sigma, then Centre then Sigma etc
number refinementloops=30
number firstchi, secondchi
number counter
for(i=0; i<refinementloops; i++)
	{
		// Refine sigma
		minchisqrd=FindSigmaMonteCarlo(new, centrex, centrey, subareapcnt, intensity, initialsigma, iterations, minchisqrd, sigma)		
		firstchi=minchisqrd		
		
		// Refine the centre
		image gaussian=Create2DGaussian(xsize, ysize, centrex, centrey, intensity, sigma,sigma)
		minchisqrd=FindCentreMonteCarlo(new, sigma, intensity, iterations, centrex, centrey, subareapcnt, minchisqrd, fitcentrex, fitcentrey)
		centrex=fitcentrex
		centrey=fitcentrey
		
		// If the chisqrd value for both sigma and centre optimisation are the same, increment a counter.
		// If that counter reaches 3, then assume the iteration is optimised and curtail further iteration.				
		if(firstchi==secondchi) counter=counter+1
		else counter=0
		
		if(counter==3) i=refinementloops
		
		if(i==refinementloops-1)  //最后一次保存到缓存
		{
		setpersistentnumbernote("ElectronDiffraction Tools:maxpeak:"+minchisqrd+":X",centrex)
		setpersistentnumbernote("ElectronDiffraction Tools:maxpeak:"+minchisqrd+":Y",centrey)
		}
	}

//删除 roi
img.clearselection()
new.deleteimage()

noclick=noclick+1

if(noclick==2)m=2
}
 if(keypress>0 && keypress!=32) // Any key except space pressed - cancel
{
m=2
}
} 
}

img:=getfrontimage()
imgdisp = img.ImageGetImageDisplay(0)

//找最小x2的参数
taggroup mcptags
taggroup ptags=getpersistenttaggroup()
ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:maxpeak", mcptags)

String mcplabel = mcptags.TagGroupGetTaglabel(0)
number value=val(mcplabel)

getpersistentnumbernote("ElectronDiffraction Tools:maxpeak:"+value+":X",x1) 
getpersistentnumbernote("ElectronDiffraction Tools:maxpeak:"+value+":Y",y1) 

x1=x1/times+left
y1=y1/times+top
component dot=NewBoxAnnotation(y1-5, x1-5, y1+5, x1+5)
dot.componentsetfillmode(2)
dot.componentsetdrawingmode(2)
dot.componentsetforegroundcolor(255,0,255)
imgdisp.componentaddchildatend(dot)

 // 设置中心
setNumberNote(img,"ElectronDiffraction Tools:Center:X", x1 )
setNumberNote(img,"ElectronDiffraction Tools:Center:Y", y1)
setorigin(img, x1, y1)
result("Center is: "+"("+x1+", "+y1+")"+"\n")

//删除缓存标记
deletepersistentnote("ElectronDiffraction Tools:maxpeak")
}
}



/*************************************************************************
Profile button
*************************************************************************/
void profile_intensity(object  self)
{
number xsize, ysize, centrex, centrey, minivalue,xscale,yscale

image frontimage:=getfrontimage()
GetSize( frontimage, xsize, ysize )
frontimage.getscale(xscale,yscale)
string imgname=frontimage.getname()
string unitstring=frontimage.GetUnitString()

Result("-----------------------------------------------------------------------------------------\nSet Profile button: \n1) to get intensity profile;\n2) ALT+ to get linear image overlaid live profile.\n-----------------------------------------------------------------------------------------\n")

//得到linear image
GetNumberNote(frontimage,"ElectronDiffraction Tools:Center:X", centrex )
GetNumberNote(frontimage,"ElectronDiffraction Tools:Center:Y", centrey )

if(OptionDown())  //ALT+得到linear image, intensity profile
{
number samples=ysize
number k=2*Pi()/samples

number minor = min( xsize, ysize )
image linearimg := RealImage( "Linear image-"+imgname, 4, minor, samples )
linearimg = warp( frontimage,	(icol * sin(irow * k))+ centrex, (icol * cos(irow * k))+ centrey )

getsize(linearimg,xsize,ysize)

if(unitstring=="1/nm")xscale=frontimage.ImageGetDimensionScale(0)/10*2*pi()
if(unitstring=="Q (1/A)")xscale=xscale

linearimg.imagesetdimensioncalibration(0,0,xscale,"Q (1/A)",0)
imagesetdimensioncalibration( linearimg, 1, 0, 360/ysize, "deg.", 0)
linearimg.ShowImage()

ImageDisplay lindisp=linearimg.ImageGetImageDisplay(0)
lindisp.ImageDisplaySetCaptionOn(1)

//添加roi
NewLiveProfile( lindisp,0,0.25*minor,minor,0.25*minor, 0.1*minor )
}

else  //默认直接得到intensity profile
{
number xwidth,ywidth,minwidth
xwidth=tert(centrex-0.5*xsize>0,xsize-centrex,centrex) 
ywidth=tert(centrey-0.5*ysize>0,ysize-centrey,centrey) 
minwidth=min(xwidth,ywidth)

number samples=2*minwidth
number k=2*Pi()/samples

image linearimg := RealImage( "Linear image-"+imgname, 4, samples, samples )
linearimg = warp( frontimage[centrey-minwidth,centrex-minwidth,centrey+minwidth,centrex+minwidth],	(icol * sin(irow * k))/2+ minwidth+0.5, (icol * cos(irow * k))/2+ minwidth+0.5 )

//投影到1D
image profile:= RealImage("Intensity profile of "+imgname,4,samples,1) 
profile[icol,0] += linearimg
profile/=samples

profile.ShowImage()

xscale=frontimage.ImageGetDimensionScale(0)
unitstring=frontimage.ImageGetDimensionUnitString(0)
if(unitstring=="1/nm")xscale=frontimage.ImageGetDimensionScale(0)/10*pi()
if(unitstring=="Q (1/A)")xscale=frontimage.ImageGetDimensionScale(0)

profile.imagesetdimensionscale(0,xscale)
profile.imagesetdimensionunitstring(0,"Q (1/A)") //标定为Q
				
profile.imagesetdimensionorigin(0,0)
profile.ImagesetDimensionCalibration(1,0,1,"Counts",0)
setwindowposition(profile,150, 320)

imagedisplay imgdisp=profile.imagegetimagedisplay(0)
number profsizex, profsizey, maxx, maxy
getsize(profile, profsizex, profsizey)

profile=profile-profile[0,profsizex-8,1,profsizex].min()

number maxval=max(profile[0,(0.5/xscale),1, profsizex])
imgdisp.LinePlotImageDisplaySetContrastLimits( 0, maxval)  
imgdisp.LinePlotImageDisplaySetDoAutoSurvey( 0, 0 ) 
imgdisp.lineplotimagedisplaysetgridon(0)  //gridon=0
imgdisp.lineplotimagedisplaysetbackgroundon(0)  //bgon=0
}
}


/*************************************************************************
Calibration button：用已知d值标定profile
*************************************************************************/
void Calibrate(object  self)
{
image img:=getfrontimage()
Number xsize, ysize,xscale,yscale,t,l,b,r
img.getsize(xsize,ysize)
img.getselection(t,l,b,r)
xscale=img.imagegetdimensionscale(0)

Result("-----------------------------------------------------------------------------------------\nProfile calibration using the known peak: \n1) Drag a ROI at the desired peak;\n2) Set the know d-spacing (A), the middle location of the ROI will be used to calibrate.\n-----------------------------------------------------------------------------------------\n")

result("xscale (before) is: "+xscale+"\n")

number d
getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Peak",d)
getnumber("d-spacing (A) of the peak",d,d)
xscale=2*Pi()*1/(0.5*(l+r)*d)
img.imagesetdimensionscale(0,xscale)
result("xscale (after) is: "+xscale+"\n")
img.clearselection()
}


/*************************************************************************
Save XY button
*************************************************************************/
void XY(object  self)
{
Image profile:=GetFrontImage()
number xsize,ysize,xscale,yscale,xorigin,yorigin
xsize=profile.imagegetdimensionsize(0)
xscale=profile.imagegetdimensionscale(0)
xorigin=profile.imagegetdimensionorigin(0)

number startvector=xorigin,endvector=xsize*xscale+xorigin
number xend=endvector-startvector

string filename, othername
string imagename=profile.getname()

If (!SaveAsDialog("Save text file as", GetApplicationDirectory(2,0) +imagename+ ".txt", filename)) Exit(0)

number fileID = CreateFileForWriting(filename)
Result("Exported file:  "+imagename+"\n")
number yval, Q

for(number i=0; i<xsize; i++)
	{
		yval=getpixel(profile, i, 0)
		number xcal=xorigin*xscale+(i*xscale)

		if(xcal>=0&&xcal<=xend)
		{
        Q=1*(xcal+startvector)
		WriteFile(fileID, Q+"        "+yval+"\n")
		}
	}
CloseFile(fileID)

profile.closeimage()
}

/*************************************************************************
Clear button
*************************************************************************/
void Clear(object  self)
{
image front:=getfrontimage()
component imgdisp=imagegetimagedisplay(front, 0)
number i

//删除pdf标记
deletenote(front,"Radial Distribution Function")

// 文字
number annots=imgdisp.componentcountchildrenoftype(13)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(13,0)
annotid.componentremovefromparent()
}

//线
annots=imgdisp.componentcountchildrenoftype(2)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(2,0)
annotid.componentremovefromparent()
}

// 箭头
annots=imgdisp.componentcountchildrenoftype(3)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(3,0)
annotid.componentremovefromparent()
}

// 框
annots=imgdisp.componentcountchildrenoftype(5)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(5,0)
annotid.componentremovefromparent()
}

//椭圆
annots=imgdisp.componentcountchildrenoftype(6)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(6,0)
annotid.componentremovefromparent()
}

// this bit removes the spot mask annotations
annots=imgdisp.componentcountchildrenoftype(8)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(8,0)
annotid.componentremovefromparent()
}

// this bit removes the arry mask annotations
annots=imgdisp.componentcountchildrenoftype(9)
for (i=0; i<annots; i++)
{
component annotid=imgdisp.componentgetnthchildoftype(9,0)
annotid.componentremovefromparent()
}


// this bit removes the roi annotations
annots= imgdisp.ImageDisplayCountROIS()
for (i=0; i<annots; i++)
{
number index
ROI currentROI = imgdisp.ImageDisplayGetROI( index )
imgdisp.imagedisplaydeleteROI(currentROI)

}

//删除图层标志
deletenote(front,"Radial Distribution Function")
}

/*************************************************************************
Get the atom number
*************************************************************************/
number GetAtomNumber(object self, string symb, number &result)
{
	result = 1
if(symb=="")     return 0
if(symb=="H"||symb== "h")          return 1
if(symb=="He"||symb=="he")         return 2
if(symb=="Li"||symb=="li")         return 3
if(symb=="Be"||symb=="be")         return 4
if(symb=="B"||symb== "b")          return 5
if(symb=="C"||symb== "c")          return 6
if(symb=="N"||symb== "n")          return 7
if(symb=="O"||symb== "o")          return 8
if(symb=="F"||symb== "f")          return 9
if(symb=="Ne"||symb=="ne")         return 10
if(symb=="Na"||symb=="na")         return 11
if(symb=="Mg"||symb=="mg")         return 12
if(symb=="Al"||symb=="al")         return 13
if(symb=="Si"||symb=="si")         return 14
if(symb=="P"||symb== "p")          return 15
if(symb=="S"||symb== "s")          return 16
if(symb=="Cl"||symb=="cl")         return 17
if(symb=="Ar"||symb=="ar")         return 18
if(symb=="K"||symb== "k")          return 19
if(symb=="Ca"||symb=="ca")         return 20
if(symb=="Sc"||symb=="sc")         return 21
if(symb=="Ti"||symb=="ti")         return 22
if(symb=="V"||symb== "v")          return 23
if(symb=="Cr"||symb=="cr")         return 24
if(symb=="Mn"||symb=="mn")         return 25
if(symb=="Fe"||symb=="fe")         return 26
if(symb=="Co"||symb=="co")         return 27
if(symb=="Ni"||symb=="ni")         return 28
if(symb=="Cu"||symb=="cu")         return 29
if(symb=="Zn"||symb=="zn")         return 30
if(symb=="Ga"||symb=="ga")         return 31
if(symb=="Ge"||symb=="ge")         return 32
if(symb=="As"||symb=="as")         return 33
if(symb=="Se"||symb=="se")         return 34
if(symb=="Br"||symb=="br")         return 35
if(symb=="Kr"||symb=="kr")         return 36
if(symb=="Rb"||symb=="rb")         return 37
if(symb=="Sr"||symb=="sr")         return 38
if(symb=="Y"||symb== "y")          return 39
if(symb=="Zr"||symb=="zr")         return 40
if(symb=="Nb"||symb=="nb")       return 41
if(symb=="Mo"||symb=="mo")      return 42
if(symb=="Tc"||symb=="tc")         return 43
if(symb=="Ru"||symb=="ru")         return 44
if(symb=="Rh"||symb=="rh")         return 45
if(symb=="Pd"||symb=="pd")         return 46
if(symb=="Ag"||symb=="ag")         return 47
if(symb=="Cd"||symb=="cd")         return 48
if(symb=="In"||symb=="in")       	  return 49
if(symb=="Sn"||symb=="sn")         return 50
if(symb=="Sb"||symb=="sb")         return 51
if(symb=="Te"||symb=="te")         return 52
if(symb=="I"||symb== "i")         	 return 53
if(symb=="Xe"||symb=="xe")         return 54
if(symb=="Cs"||symb=="cs")         return 55
if(symb=="Ba"||symb=="ba")         return 56
if(symb=="La"||symb=="la")         return 57
if(symb=="Ce"||symb=="ce")         return 58
if(symb=="Pr"||symb=="pr")         return 59
if(symb=="Nd"||symb=="nd")         return 60
if(symb=="Pm"||symb=="pm")         return 61
if(symb=="Sm"||symb=="sm")         return 62
if(symb=="Eu"||symb=="eu")         return 63
if(symb=="Gd"||symb=="gd")         return 64
if(symb=="Tb"||symb=="tb")         return 65
if(symb=="Dy"||symb=="dy")         return 66
if(symb=="Ho"||symb=="ho")         return 67
if(symb=="Er"||symb=="er")         return 68
if(symb=="Tm"||symb=="tm")         return 69
if(symb=="Yb"||symb=="yb")         return 70
if(symb=="Lu"||symb=="lu")         return 71
if(symb=="Hf"||symb=="hf")         return 72
if(symb=="Ta"||symb=="ta")         return 73
if(symb=="W"||symb== "w")          return 74
if(symb=="Re"||symb=="re")         return 75
if(symb=="Os"||symb=="os")         return 76
if(symb=="Ir"||symb=="ir")         return 77
if(symb=="Pt"||symb=="pt")         return 78
if(symb=="Au"||symb=="au")         return 79
if(symb=="Hg"||symb=="hg")         return 80
if(symb=="Tl"||symb=="tl")         return 81
if(symb=="Pb"||symb=="pb")         return 82
if(symb=="Bi"||symb=="bi")         return 83
if(symb=="Po"||symb=="po")         return 84
if(symb=="At"||symb=="at")         return 85
if(symb=="Rn"||symb=="rn")         return 86
if(symb=="Fr"||symb=="fr")         return 87
if(symb=="Ra"||symb=="ra")         return 88
if(symb=="Ac"||symb=="ac")         return 89
if(symb=="Th"||symb=="th")         return 90
if(symb=="Pa"||symb=="pa")         return 91
if(symb=="U"||symb== "u")          return 92
if(symb=="Np"||symb=="np")         return 93
if(symb=="Pu"||symb=="pu")         return 94
if(symb=="Am"||symb=="am")         return 95
if(symb=="Cm"||symb=="cm")         return 96
if(symb=="Bk"||symb=="bk")         return 97
if(symb=="Cf"||symb=="cf")         return 98
if(symb=="Es"||symb=="es")         return 99
if(symb=="Fm"||symb=="fm")         return 100
if(symb=="Md"||symb=="md")         return 101
if(symb=="No"||symb=="no")         return 102
if(symb=="Lr"||symb=="lr")          return 103	
	
	result = 0
	return 0
	}

/*************************************************************************
Get the element name
*************************************************************************/
string GetElementName(object self, number atNum, number &result)
	{
	result = 1
	if(atnum==0) return ""
	if(atNum==1) return "Hydrogen"
	if(atNum==2) return "Helium"
	if(atNum==3) return "Lithium"
	if(atNum==4) return "Beryllium"
	if(atNum==5) return "Boron"
	if(atNum==6) return "Carbon"
	if(atNum==7) return "Nitrogen"
	if(atNum==8) return "Oxygen"
	if(atNum==9) return "Fluorine"
	if(atNum==10) return "Neon"
	if(atNum==11) return "Sodium"
	if(atNum==12) return "Magnesium"
	if(atNum==13) return "Aluminum"
	if(atNum==14) return "Silicon"
	if(atNum==15) return "Phosphorus"
	if(atNum==16) return "Sulfur"
	if(atNum==17) return "Chlorine"
	if(atNum==18) return "Argon"
	if(atNum==19) return "Potassium"
	if(atNum==20) return "Calcium"
	if(atNum==21) return "Scandium"
	if(atNum==22) return "Titanium"
	if(atNum==23) return "Vanadium"
	if(atNum==24) return "Chromium"
	if(atNum==25) return "Manganese"
	if(atNum==26) return "Iron"
	if(atNum==27) return "Cobalt"
	if(atNum==28) return "Nickel"
	if(atNum==29) return "Copper"
	if(atNum==30) return "Zinc"
	if(atNum==31) return "Gallium"
	if(atNum==32) return "Germanium"
	if(atNum==33) return "Arsenic"
	if(atNum==34) return "Selenium"
	if(atNum==35) return "Bromine"
	if(atNum==36) return "Krypton"
	if(atNum==37) return "Rubidium"
	if(atNum==38) return "Strontium"
	if(atNum==39) return "Yttrium"
	if(atNum==40) return "Zirconium"
	if(atNum==41) return "Niobium"
	if(atNum==42) return "Molybdenum"
	if(atNum==43) return "Technetium"
	if(atNum==44) return "Ruthenium"
	if(atNum==45) return "Rhodium"
	if(atNum==46) return "Palladium"
	if(atNum==47) return "Silver"
	if(atNum==48) return "Cadmium"
	if(atNum==49) return "Indium"
	if(atNum==50) return "Tin"
	if(atNum==51) return "Antimony"
	if(atNum==52) return "Tellurium"
	if(atNum==53) return "Iodine"
	if(atNum==54) return "Xenon"
	if(atNum==55) return "Cesium"
	if(atNum==56) return "Barium"
	if(atNum==57) return "Lanthanum"
	if(atNum==58) return "Cerium"
	if(atNum==59) return "Praseodymium"
	if(atNum==60) return "Neodymium"
	if(atNum==61) return "Promethium"
	if(atNum==62) return "Samarium"
	if(atNum==63) return "Europium"
	if(atNum==64) return "Gadolinium"
	if(atNum==65) return "Terbium"
	if(atNum==66) return "Dysprosium"
	if(atNum==67) return "Holmium"
	if(atNum==68) return "Erbium"
	if(atNum==69) return "Thulium"
	if(atNum==70) return "Ytterbium"
	if(atNum==71) return "Lutetium"
	if(atNum==72) return "Hafnium"
	if(atNum==73) return "Tantalum"
	if(atNum==74) return "Tungsten"
	if(atNum==75) return "Rhenium"
	if(atNum==76) return "Osmium"
	if(atNum==77) return "Iridium"
	if(atNum==78) return "Platinum"
	if(atNum==79) return "Gold"
	if(atNum==80) return "Mercury"
	if(atNum==81) return "Thallium"
	if(atNum==82) return "Lead"
	if(atNum==83) return "Bismuth"
	if(atNum==84) return "Polonium"
	if(atNum==85) return "Astatine"
	if(atNum==86) return "Radon"
	if(atNum==87) return "Francium"
	if(atNum==88) return "Radium"
	if(atNum==89) return "Actinium"
	if(atNum==90) return "Thorium"
	if(atNum==91) return "Protactinium"
	if(atNum==92) return "Uranium"
	if(atNum==93) return "Neptunium"
	if(atNum==94) return "Plutonium"
	if(atNum==95) return "Americium"
	if(atNum==96) return "Curium"
	if(atNum==97) return "Berkelium"
	if(atNum==98) return "Californium"
	if(atNum==99) return "Einsteinium"
	if(atNum==100) return "Fermium"
	if(atNum==101) return "Mendelevium"
	if(atNum==102) return "Nobelium"
	if(atNum==103) return "Lawrencium"
	
	result = 0
	return ""
}

/*************************************************************************
Get the atom symbol
*************************************************************************/
string GetElementSymbol(object self, number atNum, number &result)
	{
	result = 1
	if(atnum==0) return ""
	if(atNum==1) return "H"
	if(atNum==2) return "He"
	if(atNum==3) return "Li"
	if(atNum==4) return "Be"
	if(atNum==5) return "B"
	if(atNum==6) return "C"
	if(atNum==7) return "N"
	if(atNum==8) return "O"
	if(atNum==9) return "F"
	if(atNum==10) return "Ne"
	if(atNum==11) return "Na"
	if(atNum==12) return "Mg"
	if(atNum==13) return "Al"
	if(atNum==14) return "Si"
	if(atNum==15) return "P"
	if(atNum==16) return "S"
	if(atNum==17) return "Cl"
	if(atNum==18) return "Ar"
	if(atNum==19) return "K"
	if(atNum==20) return "Ca"
	if(atNum==21) return "Sc"
	if(atNum==22) return "Ti"
	if(atNum==23) return "V"
	if(atNum==24) return "Cr"
	if(atNum==25) return "Mn"
	if(atNum==26) return "Fe"
	if(atNum==27) return "Co"
	if(atNum==28) return "Ni"
	if(atNum==29) return "Cu"
	if(atNum==30) return "Zn"
	if(atNum==31) return "Ga"
	if(atNum==32) return "Ge"
	if(atNum==33) return "As"
	if(atNum==34) return "Se"
	if(atNum==35) return "Br"
	if(atNum==36) return "Kr"
	if(atNum==37) return "Rb"
	if(atNum==38) return "Sr"
	if(atNum==39) return "Y"
	if(atNum==40) return "Zr"
	if(atNum==41) return "Nb"
	if(atNum==42) return "Mo"
	if(atNum==43) return "Tc"
	if(atNum==44) return "Ru"
	if(atNum==45) return "Rh"
	if(atNum==46) return "Pd"
	if(atNum==47) return "Ag"
	if(atNum==48) return "Cd"
	if(atNum==49) return "In"
	if(atNum==50) return "Sn"
	if(atNum==51) return "Sb"
	if(atNum==52) return "Te"
	if(atNum==53) return "I"
	if(atNum==54) return "Xe"
	if(atNum==55) return "Cs"
	if(atNum==56) return "Ba"
	if(atNum==57) return "La"
	if(atNum==58) return "Ce"
	if(atNum==59) return "Pr"
	if(atNum==60) return "Nd"
	if(atNum==61) return "Pm"
	if(atNum==62) return "Sm"
	if(atNum==63) return "Eu"
	if(atNum==64) return "Gd"
	if(atNum==65) return "Tb"
	if(atNum==66) return "Dy"
	if(atNum==67) return "Ho"
	if(atNum==68) return "Er"
	if(atNum==69) return "Tm"
	if(atNum==70) return "Yb"
	if(atNum==71) return "Lu"
	if(atNum==72) return "Hf"
	if(atNum==73) return "Ta"
	if(atNum==74) return "W"
	if(atNum==75) return "Re"
	if(atNum==76) return "Os"
	if(atNum==77) return "Ir"
	if(atNum==78) return "Pt"
	if(atNum==79) return "Au"
	if(atNum==80) return "Hg"
	if(atNum==81) return "Tl"
	if(atNum==82) return "Pb"
	if(atNum==83) return "Bi"
	if(atNum==84) return "Po"
	if(atNum==85) return "At"
	if(atNum==86) return "Rn"
	if(atNum==87) return "Fr"
	if(atNum==88) return "Ra"
	if(atNum==89) return "Ac"
	if(atNum==90) return "Th"
	if(atNum==91) return "Pa"
	if(atNum==92) return "U"
	if(atNum==93) return "Np"
	if(atNum==94) return "Pu"
	if(atNum==95) return "Am"
	if(atNum==96) return "Cm"
	if(atNum==97) return "Bk"
	if(atNum==98) return "Cf"
	if(atNum==99) return "Es"
	if(atNum==100) return "Fm"
	if(atNum==101) return "Md"
	if(atNum==102) return "No"
	if(atNum==103) return "Lr"

	result = 0
	return ""
}

/*************************************************************************
Get the atom mass
*************************************************************************/
number GetAtomicWt(object self, number atNum)
{
	if(atnum==0) return 0
	if(atNum==1) return 1.008
	if(atNum==2) return 4.0026
	if(atNum==3) return 6.941
	if(atNum==4) return 9.01218
	if(atNum==5) return 10.81
	if(atNum==6) return 12.011
	if(atNum==7) return 14.0067
	if(atNum==8) return 15.9994
	if(atNum==9) return 18.9984
	if(atNum==10) return 20.179
	if(atNum==11) return 22.9898
	if(atNum==12) return 24.305
	if(atNum==13) return 26.9815
	if(atNum==14) return 28.086
	if(atNum==15) return 30.9738
	if(atNum==16) return 32.06
	if(atNum==17) return 35.453
	if(atNum==18) return 39.948
	if(atNum==19) return 39.102
	if(atNum==20) return 40.08
	if(atNum==21) return 44.9559
	if(atNum==22) return 47.9
	if(atNum==23) return 50.9414
	if(atNum==24) return 51.996
	if(atNum==25) return 54.938
	if(atNum==26) return 55.847
	if(atNum==27) return 58.9332
	if(atNum==28) return 58.71
	if(atNum==29) return 63.546
	if(atNum==30) return 65.37
	if(atNum==31) return 69.72
	if(atNum==32) return 72.59
	if(atNum==33) return 74.9216
	if(atNum==34) return 78.96
	if(atNum==35) return 79.904
	if(atNum==36) return 83.8
	if(atNum==37) return 85.4678
	if(atNum==38) return 87.62
	if(atNum==39) return 88.9059
	if(atNum==40) return 91.22
	if(atNum==41) return 92.9064
	if(atNum==42) return 95.94
	if(atNum==43) return 99
	if(atNum==44) return 101.07
	if(atNum==45) return 102.9055
	if(atNum==46) return 106.4
	if(atNum==47) return 107.868
	if(atNum==48) return 112.4
	if(atNum==49) return 114.82
	if(atNum==50) return 118.69
	if(atNum==51) return 121.75
	if(atNum==52) return 127.6
	if(atNum==53) return 126.9045
	if(atNum==54) return 131.3
	if(atNum==55) return 132.905
	if(atNum==56) return 137.34
	if(atNum==57) return 138.905
	if(atNum==58) return 140.12
	if(atNum==59) return 140.9077
	if(atNum==60) return 144.24
	if(atNum==61) return 147
	if(atNum==62) return 150.4
	if(atNum==63) return 151.96
	if(atNum==64) return 157.25
	if(atNum==65) return 158.9254
	if(atNum==66) return 162.5
	if(atNum==67) return 164.9303
	if(atNum==68) return 167.26
	if(atNum==69) return 168.9342
	if(atNum==70) return 173.04
	if(atNum==71) return 174.97
	if(atNum==72) return 178.49
	if(atNum==73) return 180.9479
	if(atNum==74) return 183.85
	if(atNum==75) return 186.2
	if(atNum==76) return 190.2
	if(atNum==77) return 192.22
	if(atNum==78) return 195.09
	if(atNum==79) return 196.9665
	if(atNum==80) return 200.59
	if(atNum==81) return 204.37
	if(atNum==82) return 207.2
	if(atNum==83) return 208.9806
	if(atNum==84) return 210
	if(atNum==85) return 210
	if(atNum==86) return 222
	if(atNum==87) return 223
	if(atNum==88) return 226
	if(atNum==89) return 227
	if(atNum==90) return 232.0381
	if(atNum==91) return 231
	if(atNum==92) return 238.029
	if(atNum==93) return 237
	if(atNum==94) return 242
	if(atNum==95) return 243
	if(atNum==96) return 247
	if(atNum==97) return 247
	if(atNum==98) return 251
	if(atNum==99) return 254
	if(atNum==100) return 253
	if(atNum==101) return 256
	if(atNum==102) return 254
	if(atNum==103) return 257
}

/*************************************************************************
The atom symbol is changed
*************************************************************************/
void ChangedSymbol(object self, taggroup tg) 
	{
//删除缓存
deletepersistentnote("ElectronDiffraction Tools:PDF:Temporary Data")
deletepersistentnote("ElectronDiffraction Tools:PDF:Temp")
	
string atomID
tg.dlggetidentifier(atomID)   //elementfield1,elementfield2...

number strlen=len(atomID)
atomID=right(atomID,strlen-12)   //elementfield共12字节
number menuid=val(atomID)   //第几个元素,1,2,3...

//找相关元素
taggroup atomnofield=self.lookupelement("atomnofield"+menuid)
taggroup atomwtfield=self.lookupelement("atomwtfield"+menuid)
taggroup Elementfield=self.lookupelement("elementfield"+menuid)

string symbol=tg.dlggetstringvalue()  //输入的元素符号，可小写
number dummy,i,noofelements

number radiovalue   //百分比的类型
getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Use At% (0) or Wt% (1)", radiovalue)

getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Number of Elements",noofelements)

if(symbol=="")
{
dlgvalue(atomnofield,0)
dlgvalue(atomwtfield,0)

self.setelementisenabled("atomnofield"+menuid,0)
dlgvalue(self.lookupelement("atomnofield"+menuid),0)

self.setelementisenabled("atomwtfield"+menuid,0)
dlgvalue(self.lookupelement("atomwtfield"+menuid),0)

self.setelementisenabled("wtpcntfield"+menuid,0)
dlgvalue(self.lookupelement("wtpcntfield"+menuid),0)

self.setelementisenabled("atompcntfield"+menuid,0)
dlgvalue(self.lookupelement("atompcntfield"+menuid),0)
self.validateview()
}
			
//原子序数大于0
if(symbol!="")
{
//更新原子序数
number atomNo=self.GetAtomNumber(symbol,dummy)
dlgvalue(atomnofield,atomNo)

//更新元素符号
dlgvalue(Elementfield,self.GetElementSymbol(atomno,dummy))

//更新原子量
dlgvalue(atomwtfield,self.getatomicwt(atomNo))

//设置某些按钮可用、可输入
if(radiovalue==0)  //原子比
{
self.setelementisenabled("atompcntfield"+menuid,1)	
}

if(radiovalue==1)// 质量比
{
self.setelementisenabled("wtpcntfield"+menuid,1)	
}
			
self.setelementisenabled("setsfsbutton",1)
self.setelementisenabled("calculatebutton",1)		
}
}	

/*************************************************************************
Make atom number field
*************************************************************************/
taggroup makeatomnofield(object self, number fieldid)
	{
		taggroup atomnofield=dlgcreateintegerfield(0,4).dlgidentifier("atomnofield"+fieldid).dlgenabled(0)

		return atomnofield
	}

/*************************************************************************
Make element field
*************************************************************************/
taggroup makeElementfield(object self, number fieldid)
	{
	string changemethod="ChangedSymbol"
		taggroup elementfield=dlgcreatestringfield("",6,changemethod).dlgidentifier("elementfield"+fieldid)

		return elementfield
	}
	
/*************************************************************************
Atom weight
*************************************************************************/
taggroup makeatomweightfield(object self, number fieldid)
	{
		taggroup atomweightfield=dlgcreaterealfield(0,6,4).dlgidentifier("atomwtfield"+fieldid).dlgenabled(0)
		return atomweightfield
	}

/*************************************************************************
Atom%
*************************************************************************/
taggroup makeatompcntfield(object self, number fieldid)
	{
		taggroup atompcntfield=dlgcreaterealfield(0,6,4).dlgidentifier("atompcntfield"+fieldid).dlgenabled(0).dlgchangedmethod("atfieldcalc")
		return atompcntfield
	}

/*************************************************************************
Wt. %
*************************************************************************/
taggroup makeweightpcntfield(object self, number fieldid)
	{
		taggroup weightpcntfield=dlgcreaterealfield(0,6,4).dlgidentifier("wtpcntfield"+fieldid).dlgenabled(0).dlgchangedmethod("wtfieldcalc")
		return weightpcntfield
	}

/*************************************************************************
Calculate the scattering factor
*************************************************************************/
number calculatescatteringfactor(object self, number element, number kvalue)
{
// Some traps
if(element<1) element=1
if(element>103) element=103
if(kvalue<=0) kvalue=0

// Variables
number i, j, q, k
number h, aval, bval, cval, dval
number  scatf = 0.0

// Step through the calculation of the three lorenztians and three gaussians which make up the approximation
// the aval, bval, cval and dval parameters are stored in the Global Info
    for( i=0; i<3; i++)
		{    
        //Lorenztians        
        getpersistentnumbernote("ElectronDiffraction Tools:PDF:Scattering Factors:Element "+element+":a"+(i+1),aval)
		getpersistentnumbernote("ElectronDiffraction Tools:PDF:Scattering Factors:Element "+element+":b"+(i+1),bval)
        scatf=scatf+aval/( kvalue**2 + bval )
		
		// Gaussians		
		getpersistentnumbernote("ElectronDiffraction Tools:PDF:Scattering Factors:Element "+element+":c"+(i+1),cval)
		getpersistentnumbernote("ElectronDiffraction Tools:PDF:Scattering Factors:Element "+element+":d"+(i+1),dval)
        scatf=scatf+cval*exp( -kvalue**2 * dval)
   	}
	
return scatf
}

/*************************************************************************
Calculate button
*************************************************************************/
void Calculate(object self) // 计算化学组分
{
	number i, j,noofelementsselected, noofelements, sumatno,sumwtpcnt, sumatpcnt
	number thisatwt, massaccumulator, thiswtpcnt, thisatpcnt, summoles

	// 是at% 还是wt%
	number atorwtpcntstatus
	getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Use At% (0) or Wt% (1)", atorwtpcntstatus)

	//元素个数
	getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Number of Elements",noofelements)
	
	//由元素%来识别有多少个元素
	number statusaccumulator	
	for(i=1;i<noofelements+1;i++)
		{
			taggroup atno=self.lookupelement("atomnofield"+i)
			sumatno=sumatno+atno.dlggetvalue()

			taggroup wtpcnt=self.lookupelement("wtpcntfield"+i)
			sumwtpcnt=sumwtpcnt+wtpcnt.dlggetvalue()

			taggroup atpcnt=self.lookupelement("atompcntfield"+i)
			sumatpcnt=sumatpcnt+atpcnt.dlggetvalue()
		}

		Setpersistentnumbernote("ElectronDiffraction Tools:PDF:Temp:Sum of atoms", sumatpcnt)
		
		// The atom % radio button is checked but no values have been entred in the atom % column
		if(atorwtpcntstatus==0 && sumatpcnt==0)
			{
				showalert("You must enter atomic percentages.",2)
				self.setelementisenabled("calculatebutton",1)
				return
			}

		// The weight % radio button is checked but no values have been entred in the weight % column
		if(atorwtpcntstatus==1 && sumwtpcnt==0)
			{
				showalert("You must enter weight percentages.",2)
				self.setelementisenabled("calculatebutton",1)
				return
			}							

	// Update the progress bar
	self.dlgsetprogress("progbar",0.33)
	self.validateview()

	//number radiovalue=dlggetvalue(self.lookupelement("atwtradiosetting"))
	number radiovalue
	getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Use At% (0) or Wt% (1)", radiovalue)

if(radiovalue==0) // 原子比
{
	//元素个数
	getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Number of Elements",noofelements)
	
	// 计算总质量
	sumatpcnt=0
	for(i=1; i<noofelements+1;i++)
		{
			thisatwt=dlggetvalue(self.lookupelement("atomwtfield"+i))
			thisatpcnt=dlggetvalue(self.lookupelement("atompcntfield"+i))
			if(thisatpcnt>0) sumatpcnt=sumatpcnt+thisatpcnt
			if(thisatwt>0) massaccumulator=massaccumulator+(thisatpcnt*thisatwt)
		}

		sumwtpcnt=0
	for(i=1; i<noofelements+1;i++)
		{
			thisatwt=dlggetvalue(self.lookupelement("atomwtfield"+i))
			thisatpcnt=dlggetvalue(self.lookupelement("atompcntfield"+i))

			// Trap for divide by zero error when 0 is entered into the atom % field
			if(thisatwt>0 && massaccumulator>0) thiswtpcnt=((thisatpcnt*thisatwt)/massaccumulator)*100
			else thiswtpcnt=0
	
			if(thisatwt>0) dlgvalue(self.lookupelement("wtpcntfield"+i),thiswtpcnt)	
			if(thisatwt>0) sumwtpcnt=sumwtpcnt+thiswtpcnt
		}

		Setpersistentnumbernote("ElectronDiffraction Tools:PDF:Temp:Sum of weight", sumwtpcnt)			
}

if(radiovalue==1) // 质量比
{
	// Get the number of rows of elemental data to check
	getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Number of Elements",noofelements)

	sumatpcnt=0
	sumwtpcnt=0
	summoles=0

	// Calculate the total number of moles in the specimen
	for(i=1; i<noofelements+1;i++)
		{
			thisatwt=dlggetvalue(self.lookupelement("atomwtfield"+i))
			thiswtpcnt=dlggetvalue(self.lookupelement("wtpcntfield"+i))
			if(thisatwt>0) summoles=summoles+(thiswtpcnt/thisatwt)
			if(thisatwt>0) sumwtpcnt=sumwtpcnt+thiswtpcnt
		}

	for(i=1; i<noofelements+1;i++)
		{
			thisatwt=dlggetvalue(self.lookupelement("atomwtfield"+i))
			thiswtpcnt=dlggetvalue(self.lookupelement("wtpcntfield"+i))

			// Note need to trap for a value of 0 entered into the wt % field  divide by zero errors occur
			if(thisatwt>0 && summoles>0) thisatpcnt=((thiswtpcnt/thisatwt)/summoles)*100
			else thisatpcnt=0

			if(thisatwt>0) dlgvalue(self.lookupelement("atompcntfield"+i),thisatpcnt)	
			if(thisatwt>0) sumatpcnt=sumatpcnt+thisatpcnt
		}

	Setpersistentnumbernote("ElectronDiffraction Tools:PDF:Temp:Sum of atoms", sumatpcnt)	
	Setpersistentnumbernote("ElectronDiffraction Tools:PDF:Temp:Sum of weight", sumwtpcnt)
	result(sumatpcnt+", "+sumwtpcnt+"\n")	

	self.validateview()
}


// Update the progress bar
self.dlgsetprogress("progbar",0.66)
self.validateview()							
							
	// Prepare some results for outputting
	string date, time
	date=getdate(0)
	time=gettime(1)

	result("\n\nThe composition are:\n-----------------------------------------------------------------------------------------\n")
	result("Symbol\t\tAt. No.\t\tAt. Wt.\t\t\tAt. %\t\t\tWt. %")

	//输出元素组分
	for(i=1; i<noofelements+1;i++)
		{
			number atno=dlggetvalue(self.lookupelement("atomnofield"+i))
			number atwt=dlggetvalue(self.lookupelement("atomwtfield"+i))
			number wtpcnt=dlggetvalue(self.lookupelement("wtpcntfield"+i))
			number atpcnt=dlggetvalue(self.lookupelement("atompcntfield"+i))
		
			// Get the element symbol
			number dummy
			string symbol=self.GetElementSymbol(atno,dummy)

			if(atno>0&&(wtpcnt>0||atpcnt>0)) result("\n  "+symbol+"\t\t"+format(atno,"%3.0f")+"\t\t\t"+format(atwt,"%6.3f")+"\t\t\t"+format(atpcnt,"%6.3f")+"\t\t\t"+format(wtpcnt,"%6.3f"))
		}

	// update the progrss bar
	self.dlgsetprogress("progbar",1)
	self.validateview()							
							
// Look up the summed atomic percentages - so that the individual atom% values can be scaled accordingly for calculating the mean scattering factor.
number sumatpcntvalue
getpersistentnumbernote("ElectronDiffraction Tools:PDF:Temp:Sum of atoms", sumatpcntvalue)

number noelements, noscatteringfactors, k
string pulldownelement
getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Number of Elements",noelements)

		for(j=1; j<noelements+1; j++)
			{
				//获取元素名称
				number result
				number atno=dlggetvalue(self.lookupelement("atomnofield"+j))
				string pulldownelement=self.GetElementName(atno,result)
				
				//获取元素含量
				number atompcntfieldvalue=dlggetvalue(self.lookupelement("atompcntfield"+j))
				number atomfraction=atompcntfieldvalue/sumatpcntvalue

				// 写入缓存
				if(pulldownelement!="")
					{
						setpersistentnumbernote("ElectronDiffraction Tools:PDF:Temporary Data:"+pulldownelement+":Atomic Fraction",atomfraction)
						setpersistentnumbernote("ElectronDiffraction Tools:PDF:Temporary Data:"+pulldownelement+":Atomic Number",atno)					
					}
				
				self.validateview()
			}		

		// Get the density and compute the no of atoms per unit volume (A3)		
		number densityvalue=dlggetvalue(self.lookupelement("densityfield"))
		number atwt, wtpercent, wtfraction,molespercc, molesperA3,atomsperA3, totalatomsperA3
		getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Number of Elements", noofelements)

		for (i=0; i<noofelements;i++)
			{
				// Get the atomic weight and wt percentage values for each of the elements selected
				atwt=dlggetvalue(self.lookupelement("atomwtfield"+(i+1)))
				wtpercent=dlggetvalue(self.lookupelement("wtpcntfield"+(i+1)))
				if(atwt==0) break
				
				// Calculate the number of moles of each element present in 1cm3
				wtfraction=wtpercent/100
				molespercc=(densityvalue*wtfraction)/atwt
				
				// Convert that into moles (and then atoms) per cubic Angstrom
				molesperA3=molespercc/1e24
				atomsperA3=molesperA3*6.02214E23
				totalatomsperA3=totalatomsperA3+atomsperA3
			}

			// Output the density calculations to the results window
			number atomNumber
			Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Temp:Sum of atoms", atomNumber)
			//totalatomsperA3=totalatomsperA3/atomNumber    //应除以总原子数
			
			densityvalue=dlggetvalue(self.lookupelement("densityfield"))  //由密度、分子量直接计算数密度
			totalatomsperA3=densityvalue/massaccumulator*6.02214/10  //输出时应*e24，其中10为1cm^3=1e24A^3的差别

			
			result("\n\nDensity Calculations:\n-----------------------------------------------------------------------------------------\n")
			result("Material Density = "+densityvalue+" g/cm3,\tAtomic Density = "+totalatomsperA3+" E24 atoms/cm3\n")
			Result("Molecular weight = "+massaccumulator+" g/mol\n")
			
			Setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Atomic Density",totalatomsperA3)
			Setpersistentnumbernote("ElectronDiffraction Tools:PDF:Temp:Molecular weight", massaccumulator)	
			
	self.setelementisenabled("Scatterbutton",1)

	self.dlgsetprogress("progbar",0.0)
	self.validateview()	
}

/*************************************************************************
Scatter button: Add scattering curves on the profile
*************************************************************************/
void Scatter(object self)
{
	// Check to see if any images are displayed
	number nodocs=countdocumentwindowsoftype(5)

	if(nodocs==0)
		{
			showalert("There are no images or profiles displayed.",2)
			return
		}

	// Get the frontmost image and check if it is an image or a profile
	sourceimg:=getfrontimage()
	number xsize, ysize
	getsize(sourceimg, xsize, ysize)
	string imgname=getname(sourceimg)
	
	if(ysize==1) // It's a profile
		{
		self.setelementisenabled("Scatterbutton",0)
			
		// Check to ensure that the profile has been calibrated - get the units string and
		// ensure that the leftmost character is 1 ie either 1/A or 1/nm
		number origin, scale
		string unitsstring, unitsshort
		imagegetdimensioncalibration(sourceimg, 0, origin, scale, unitsstring, 0)
		if(len(unitsstring)>0) unitsshort=left(unitsstring,1)

		if(unitsshort!="1"&&unitsshort!="Q")  //单位不以1开头，比如1/nm, 1/A和Q开头
			{
				if(!twobuttondialog("This profile should be calibrated in reciprocal distance (1/A or 1/nm).","Proceed","Cancel"))
					{
						return
					}
			}

		// Do the image calibration - based on the radio settings		
		image profileimg=imageclone(sourceimg)
		imagecopycalibrationfrom(profileimg, sourceimg)  //复制-显示profile	

		number angle
		if(getnumbernote(sourceimg,"ElectronDiffraction Tools:PDF:Angle",angle))   //提取角度数据，以便调用
			{
				setnumbernote(profileimg,"ElectronDiffraction Tools:PDF:Angle",angle)
			}
			
		// 数据标定
		string currentunits=imagegetdimensionunitstring(profileimg,0)
		string idstring

		if(currentunits=="1/?")   //单位为1/A
			{
				number xscale=profileimg.imagegetdimensionscale(0)
				profileimg.imagesetdimensionscale(0,xscale*2*pi())
				profileimg.imagesetdimensionunitstring(0,"Q (1/A)")
				idstring="Q"
			}

		if(currentunits=="1/nm")  //单位为1/nm
			{
				number xscale=profileimg.imagegetdimensionscale(0)
				profileimg.imagesetdimensionscale(0,xscale/10*2*pi())
				profileimg.imagesetdimensionunitstring(0,"Q (1/A)")
				idstring="Q"
			}			

		// Set the name of the profile
		string sourcename=getname(sourceimg)
		setname(profileimg, sourcename)
		SetStringNote(profileimg,"ElectronDiffraction Tools:PDF:Image name",imgname)  //以备调用名字

		// add a marker tag to identify it for future deletion
		setstringnote(profileimg, "Radial Distribution Function","Intensity Profile")
		
		//-----------接下来计算结构因子曲线--------------
		//计算电子速度
		number voltage=dlggetvalue(self.lookupelement("voltagefield"))
		number relvelocity, c
		c=2.998e8
		relvelocity=c*sqrt(1-(1/(1+(voltage/511))**2))

		// Use that to calculate the relative relativistic mass of the electron ie mrel/m0
		// scattering factors are multiplied by the mrel/m0 factor
		number emassrest, emassrel
		emassrest=9.109e-31
		emassrel=1/sqrt((1-(relvelocity**2/c**2)))		

		//检查标定
		string unitstring=profileimg.imagegetdimensionunitstring(0)
		unitstring=left(unitstring,1)

		number xscale=profileimg.imagegetdimensionscale(0)
		getsize(profileimg, xsize, ysize)
		number maxxvalue=xsize*xscale // the maximum value of K (s) or q in the experimental profile		
	
		// This is the image to which the scaled (by the atom fraction) scattering factor plots are added.
		// Note if we are working in units of Q we have to divide the max x value by 2 x Pi
		if(unitstring=="Q") maxxvalue=maxxvalue/(2*pi())

		// 元素的数目
		number noelements, noscatteringfactors, k
		string pulldownelement

		getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Number of Elements",noelements)

		// 开始计算结构因子: 
		image temp:=realimage("",4,xsize,1) // a temporary image for storing data
		imgf2q:=realimage("",4,xsize,1) // RIF component image of N*f**2 as a function of q
		imgfq2:=realimage("",4,xsize,1) // RIF component image of f as a function of (Nfq) **2
		number plotstep=maxxvalue/xsize
		number sfvalue

		// During the calculation function information on the elements selected and the atomic fraction was stored in some temporary tags - souce this info for the calculation
		taggroup ptags=getpersistenttaggroup()
		taggroup temptags
		ptags.taggroupgettagastaggroup("ElectronDiffraction Tools:PDF:Temporary Data", temptags)
		number counter=0
		number noofchosenelements=temptags.taggroupcounttags()	

if(noofchosenelements==1) // 单一元素
	{
		// Get some information on the selected element from Temporary data tags (created by the 'Calculate' button				
		string elementname=temptags.taggroupgettaglabel(0)
		number thisatomnumber		
		
		// Note we do not need to source the atomic fraction, since for a single element system
		// this is, by definition, equal to 1				
		getpersistentnumbernote("ElectronDiffraction Tools:PDF:Temporary Data:"+elementname+":Atomic Number",thisatomnumber)

		// Calculate the scattering factor plot from Kirkland's data
		for(number i=0;i<xsize;i++)
			{
				sfvalue=self.calculatescatteringfactor(thisatomnumber, ((i/xsize)*maxxvalue))
				
						if(mod(counter, round(xsize/5))==0)
						{
							self.dlgsetprogress("progbar",i/xsize)
							self.validateview()
						}
				
				setpixel(temp,i,0,sfvalue)
			}
			
		temp=temp*emassrel

		// Since this is a single component system (atomic fraction=1) both the f2q and fq2 images have the same value 
		imgf2q=temp**2
		imgfq2=temp**2
		deleteimage(temp)
		
		self.dlgsetprogress("progbar",0)
		self.validateview()
}

else // 多元素
{
number counter=1
	for(number i=0; i<noofchosenelements; i++)
		{
			string elementname=temptags.taggroupgettaglabel(i)
			number thisatomfraction, thisatomnumber
			getpersistentnumbernote("ElectronDiffraction Tools:PDF:Temporary Data:"+elementname+":Atomic Fraction",thisatomfraction)
			getpersistentnumbernote("ElectronDiffraction Tools:PDF:Temporary Data:"+elementname+":Atomic Number",thisatomnumber)

			number j
			for(j=0;j<xsize;j++)
				{
					sfvalue=self.calculatescatteringfactor(thisatomnumber, ((j/xsize)*maxxvalue))
					
					if(mod(counter, round((noofchosenelements*xsize/10)))==0)
						{
							self.dlgsetprogress("progbar",(counter/(xsize*noofchosenelements)))
							self.validateview()
						}
						
					setpixel(temp,j,0,sfvalue)
					counter=counter+1
				}

			// Scale the scattering factor plot by the relativistic factor calculated earlier			
			temp=temp*emassrel	
	
			imgf2q=imgf2q+temp*thisatomfraction
			imgfq2=imgfq2+temp**2*thisatomfraction
		}
		
		imgf2q*=imgf2q
		
		
		// clean up		
		deleteimage(temp)
		self.dlgsetprogress("progbar",0)
		self.validateview()
		
}		//结构因子计算完毕
					
		// Get the default settings for height and aspect ratio of the profile from th global info
		// and display the profile accordingly
		number profileaspectratio, profileheight
		getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Profile Aspect Ratio",profileaspectratio)
		getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Default Profile Height", profileheight)
		
		profileimg.ShowImage()
		profileimg.setwindowposition(142,24)
	
		imagedocument imgdoc=getfrontimagedocument()
		documentwindow tempwin=imgdoc.imagedocumentgetwindow()
		profileimg.imagesetintensityscale(1)		
		
		tempwin.windowsetframesize(profileheight*profileaspectratio, profileheight)

		//第一图层为I(Q)
		lineplotimagedisplay lindisp=profileimg.imagegetimagedisplay(0)	
		object sliceid=lindisp.imagedisplaygetsliceidbyindex(0)
		lindisp.imagedisplaysetslicelabelbyid(sliceid, "Raw")
		lindisp.lineplotimagedisplaysetlegendshown(1)
		
//第二图层为f(q)，第三图层为fq2，第四图层为rif
lindisp.imagedisplayaddimage(imgf2q, "<f>2")
lindisp.imagedisplayaddimage(imgfq2, "<f2>")

lindisp.lineplotimagedisplaysetgridon(0)  //gridon=0
lindisp.lineplotimagedisplaysetbackgroundon(0)  //bgon=0

imagesetIntensityUnitString(profileimg,"Intensity")
				
		// Enable/disable the dialog buttons appropriately
		self.setelementisenabled("Scatterbutton",0)
		self.setelementisenabled("fitprofilebutton",1)
		self.setelementisenabled("nupbutton",1)
		self.setelementisenabled("ndownbutton",1)
		self.setelementisenabled("rdfbutton",1)

		//This is the end of the section which deals with a profile being foremost				
		return
	}


	// Display the selected image

	showimage(sourceimg)
	setwindowposition(sourceimg, 142,24)
	updateimage(sourceimg)

	// Attach a flag to the selected image (for subsequent cleaning up operations)

	setstringnote(sourceimg, "Radial Distribution Function","Source Image")
}


/*************************************************************************
Reduced G(r)
*************************************************************************/
void grbutton(object self)
{
image img:=GetFrontImage()
ImageDisplay imgdisp=img.ImageGetImageDisplay(0)

//-------------获取参数c和N---------------
number N,N_factor,n0,c,c_factor,l,r
getnumbernote(img,"ElectronDiffraction Tools:PDF:Parameters:N factor",N_factor) 
getnumbernote(img,"ElectronDiffraction Tools:PDF:Parameters:N0",n0)
getnumbernote(img,"ElectronDiffraction Tools:PDF:Parameters:C factor",c_factor)
			
//寻找N roi
number annots=imgdisp.ImageDisplayCountROIS()
for(number j=0;j<annots;j++)
	{
		ROI currentROI =imgdisp.ImageDisplayGetROI( j )
		string thischar=currentROI .ROIGetLabel()
		string roilabel=left(thischar,1)
			
		if(roilabel=="N")
			{
				currentROI.ROIGetRange(l,r)
				N=(r+l)/(2)
				N=N*N_factor*n0
			}
			
		if(roilabel=="c")	
		{
		currentROI.ROIGetRange(l,r)
		c=(0.5*(l+r)-0.5*img.ImageGetDimensionSize(0))*c_factor  
		}
	}
	
number atomicdensity,maxradius,resolution
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Atomic Density",atomicdensity)
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:R_max",maxradius)
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Resolution",resolution)
		
image rifimage
string imgname

TagGroup tg =img.ImageGetTagGroup( )  //缓存中是否有文件名，
TagGroup tg1
if(tg.taggroupgettagastaggroup("ElectronDiffraction Tools:PDF", tg1))
{
getStringNote(img,"ElectronDiffraction Tools:PDF:Image name",imgname)
}

else
{
imgname=img.getname()
}	
		
		number noslices=imgdisp.imagedisplaycountslices() //寻找rif图层
		for(number i=noslices-1; i>-1; i--)
			{
				object sliceid=imgdisp.imagedisplaygetsliceidbyindex(i)
				string slicename=imgdisp.imagedisplaygetslicelabelbyid(sliceid)
				if(slicename=="Damped") 
					{
						rifimage=imgdisp{"Damped"}
						break
					}
					
				if(slicename=="F(Q)")
					{
						rifimage=imgdisp{"F(Q)"}
						break					
					}
					
				if(noslices==1)
					{
						rifimage=imgdisp{0}
						break					
					}
			}
		
		//phi is the reduced intensity function (RIF).  For monatomic solids: phi(s) = (I(s)/Nf(s)^2-1)*s
		number i,j,  rmax, smax, numr, nums,qscale
		rmax=maxradius
		numr=resolution

		//The RDF function expects the data to be in terms of s - so test to see if it is in terms of Q
		// and if so, provide an appropriate scaling factor
		string isqstring=img.imagegetdimensionunitstring(0)

		if(isqstring=="1/nm"||left(isqstring,1)=="s")   //针对intensity profile, xscale=0.1*xscale
		{
		qscale=0.1*img.imagegetdimensionscale(0)  //1/nm到1/A
		qscale=2*Pi()*qscale   //再到Q
		}
		
		else if(left(isqstring,1)=="Q")
		{
		qscale=img.imagegetdimensionscale(0)
		}

		number ysize
		getsize(img, nums, ysize) 
		
		smax=qscale*nums
		
		// decrement the number of integration points if not odd for Simpson's rule
	    if((nums+1)%2!=0) nums-=1; 

		image Gr, Gr_matrix, phi_matrix,f_sr
		f_sr = exprsize(nums,numr,0)
		Gr_matrix = exprsize(nums,numr,0);
		Gr = exprsize(numr,0);
		phi_matrix = exprsize(nums,numr,0)
		
		f_sr  = sin(icol*smax/nums * irow*rmax/numr) //sin(Q*r)
		f_sr *=2/Pi(); //F(Q)时的系数，若为F(s)，则为8*pi
        // below, the weightings of each icol are adjusted to implement Simpson's rule
        f_sr*=(4*(icol%2)+2*((icol+1)%2))/3.0*tert(icol==0,0,1)*tert(icol==nums-1,0,1)+tert(icol==0,1,0)+tert(icol==nums-1,1,0)

		slice2(phi_matrix, 0, 0, 0, 0, nums, 1, 1, numr, 1) = rifimage[icol,0] 
		Gr_matrix = f_sr*phi_matrix
		Gr[irow,0] += slice2(Gr_matrix, 0, 0, 0, 0, nums, 1, 1, numr, 1) //integration summation
		Gr *= smax/nums; //integration interval		

//Gr to RDF
image grimage=Gr
grimage.showimage()
ImageDisplay grdisp=grimage.ImageGetImageDisplay(0)		
setstringnote(grimage,"Radial Distribution Function","Reduced Density Function")	
setwindowposition(grimage,195,40)
		
		// Label the G(r) image appropriately
		setname(grimage,"Reduced G(r) - "+imgname)
		imagesetIntensityUnitString(grimage,"Reduced G(r)")
		ImageSetDimensionCalibration(grimage,0,0,rmax/numr,"r (A)",1)
		SetStringNote(grimage,"ElectronDiffraction Tools:PDF:Image name",imgname)  //保存名字，以备调用		
		
		number x,y,xsize
		grimage.getsize(xsize,ysize)
		grimage[0,0.05*xsize,1,xsize].minmax(x,y)
		grdisp.lineplotimagedisplaysetdoautosurvey(0,0)
		grdisp.lineplotimagedisplaysetcontrastlimits(1.1*x,1.1*y)
		
		setnumbernote(grimage,"ElectronDiffraction Tools:PDF:Parameters:N",N)
		setnumbernote(grimage,"ElectronDiffraction Tools:PDF:Parameters:c",c)

		number angle
		if(getnumbernote(img,"ElectronDiffraction Tools:PDF:Angle",angle))   //提取角度数据，以便调用
			{
				setnumbernote(grimage,"ElectronDiffraction Tools:PDF:Angle",angle)
			}	
}

/*************************************************************************
G(r) calibration
*************************************************************************/
void grcalibration(object self)
{
ChooseMenuItem("ePDF Tools","","G(r)-Calibration")
}

/*************************************************************************
g(r) to G(r)
*************************************************************************/
void gr2ReducedGr(object self)
{
image gr:=getfrontimage()
image Reducedgr=gr*0
image rGr=gr*0
image RDF=gr*0

string imgname
gr.GetStringNote("ElectronDiffraction Tools:PDF:Image name",imgname)

number atomicdensity,maxradius,resolution
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Atomic Density",atomicdensity)
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:R_max",maxradius)
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Resolution",resolution)

Reducedgr=(gr-1)*4*pi()*atomicdensity*(icol*maxradius/resolution)
Reducedgr.showimage()
Reducedgr.SetName("Reduced G(r) - "+imgname)
Reducedgr.ImageCopyCalibrationFrom(gr)
SetStringNote(Reducedgr,"ElectronDiffraction Tools:PDF:Image name",imgname)  
setstringnote(Reducedgr,"Radial Distribution Function","Reduced Density Function")
imagesetIntensityUnitString(Reducedgr,"Reduced G(r)")
setwindowposition(Reducedgr, 150,60)

number flagrGr
getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:rG(r)", flagrGr)
if(flagrGr==1)
{
rGr=Reducedgr*(icol*maxradius/resolution)
rGr.showimage()
rGr.SetName("rG(r) - "+imgname)
rGr.ImageCopyCalibrationFrom(gr)
SetStringNote(rGr,"ElectronDiffraction Tools:PDF:Image name",imgname)  
imagesetIntensityUnitString(rGr,"rG(r)")
setwindowposition(rGr, 142,24)
}

RDF=gr*4*pi()*atomicdensity*(icol*maxradius/resolution)**2
RDF.showimage()
RDF.SetName("RDF - "+imgname)
RDF.ImageCopyCalibrationFrom(gr)
SetStringNote(RDF,"ElectronDiffraction Tools:PDF:Image name",imgname)  
setstringnote(RDF,"Radial Distribution Function","Radial Distribution Function")
imagesetIntensityUnitString(RDF,"RDF")
setwindowposition(RDF, 150,60)

Reducedgr.ShowImage()
}

/*************************************************************************
G(r) to F(Q)
*************************************************************************/
void Gr2FQ(object self)
{
image grimg:=GetFrontImage()
image f2q,fq2,Iq

string imgname
grimg.GetStringNote("ElectronDiffraction Tools:PDF:Image name",imgname)

string idstring,unitstring
getstringnote(grimg, "Radial Distribution Function",idstring)  //仅对Gr操作
if(idstring!="Reduced Density Function")
	{
		break
	}
	
//找F(q)
number i, grfoundflag, markedimgs,Qsize,ysize,rscale,Qscale
image tempimg

number t,b,l,r,N,c,N_factor,c_factor,hi,lo,sample_l,sample_r


rscale=grimg.ImageGetDimensionScale(0)

number nodocs=countdocumentwindowsoftype(5)
for(i=0; i<nodocs; i++)
	{
		imagedocument imgdoc=getimagedocument(i)
		image tempimg:=imgdoc.imagedocumentgetimage(0)
		getstringnote(tempimg, "Radial Distribution Function",idstring)

		if(idstring=="Intensity Profile")
			{
			Iq=tempImg*0
			tempImg.GetSize(Qsize,ysize)
			Qscale=tempImg.ImageGetDimensionScale(0)			
			unitstring=tempimg.imagegetdimensionunitstring(0)
			
			tempImg.ShowImage()
			ImageDisplay imgdisp=tempImg.ImageGetImageDisplay(0)
			f2q=imgdisp{"<f>2"}
			fq2=imgdisp{"<f2>"}
			
			//读取N和c值
			getnumbernote(tempimg,"ElectronDiffraction Tools:PDF:Parameters:N factor",N_factor) 
			getnumbernote(tempimg,"ElectronDiffraction Tools:PDF:Parameters:C factor",c_factor)

			number annots= imgdisp.ImageDisplayCountROIS()
			for (number i=0; i<annots; i++)
				{
					ROI currentROI = imgdisp.ImageDisplayGetROI( i )
					string thischar=currentROI .ROIGetLabel()
					string roilabel=left(thischar,1)

					if(roilabel=="c")
						{
							currentROI.ROIGetRange(l,r)
							c=(l-0.5*Qsize)*c_factor   //极值的10倍
						}

					if(roilabel=="N")
						{
							currentROI.ROIGetRange(l,r)
							N=(r+l)/(2)
							N=N*N_factor
						}
				}
				
			break
			}
	}
	
	// numr is the resolution of the calculated plot
	number nums=Qsize
	number  numr
	getsize(grimg,numr,r);

	image phi, phi_matrix, Gr_matrix;
	image f_sr

    // decrement the number of integration points if not odd for Simpson's rule
	if((numr+1)%2!=0) numr-=1; 

	f_sr = exprsize(numr,nums,0);
	phi_matrix = exprsize(numr,nums,0);
	Gr_matrix = exprsize(numr,nums,0);
	phi = exprsize(nums,0);
	f_sr = sin(icol*rscale*irow*Qscale);
	f_sr /= (2*pi()); //constants in accordance with Cockayne & McKenzie Acta Cryst. 1988
    // below, the weightings of each icol are adjusted to implement Simpson's rule
    f_sr*=(4*(icol%2)+2*((icol+1)%2))/3.0*tert(icol==0,0,1)*tert(icol==numr-1,0,1)+tert(icol==0,1,0)+tert(icol==numr-1,1,0)

	slice2(Gr_matrix, 0, 0, 0, 0, numr, 1, 1, nums, 1) = grimg[icol,0] 
	phi_matrix = f_sr*Gr_matrix
	phi[irow,0] += slice2(phi_matrix, 0, 0, 0, 0, numr, 1, 1, nums, 1) //integration summation
	phi *=rscale; //integration interval
	
	// If data are calibrated in Q then modify the scaling - the above calculation was done in terms
	// of s, so now change the scaling so that it is in terms of Q	
	if(unitstring=="Q (1/A)") ImageSetDimensionCalibration(phi,0,0,Qscale,unitstring,1)
	
	// Set the labelling on the image	
	imagesetintensityunitinfo(phi, "F(Q)",1)	
	
	setname(phi,"F(Q) from G(r) - "+imgname)
	setstringnote(phi,"Radial Distribution Function","Calculated Reduced Interference Function")
	phi.SetStringNote("ElectronDiffraction Tools:PDF:Image name",imgname)

	showimage(phi)
	setwindowposition(phi,202,84)
	
	//开始计算Iq
	Iq=phi*N*f2q+N*fq2-c
	Iq.showimage()
	
	Qscale=10*Qscale/(2*Pi())
	unitstring="Diffraction vector (1/nm)"
	Iq.ImageSetDimensionCalibration(0,0,Qscale,unitstring,1)
	
	// Set the labelling on the image	
	imagesetintensityunitinfo(Iq, "Calculated I(Q)",1)	
	
	setname(Iq,"I(Q) from G(r) - "+imgname)
	setstringnote(Iq,"Radial Distribution Function","Calculated diffraction profile")
	Iq.SetStringNote("ElectronDiffraction Tools:PDF:Image name",imgname)	
	
	setwindowposition(Iq, 142,24)
	setwindowposition(phi, 150,60)
	phi.ShowImage()
}

/*************************************************************************
g(r) to RDF
*************************************************************************/
void gr2RDF(object self)
{
image gr:=getfrontimage()
image RDF=gr*0

string imgname
gr.GetStringNote("ElectronDiffraction Tools:PDF:Image name",imgname)

number atomicdensity,maxradius,resolution
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Atomic Density",atomicdensity)
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:R_max",maxradius)
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Resolution",resolution)

RDF=gr*4*pi()*atomicdensity*(icol*maxradius/resolution)**2
RDF.showimage()
RDF.SetName("RDF - "+imgname)
RDF.ImageCopyCalibrationFrom(gr)
SetStringNote(RDF,"ElectronDiffraction Tools:PDF:Image name",imgname)  
setstringnote(RDF,"Radial Distribution Function","Radial Distribution Function")
imagesetIntensityUnitString(RDF,"RDF")
}


//G>R按钮：function grtorifbutton()根据gr计算rif曲线
void grtorifbutton(object self)
{
  //This function computes the RIF = 8*PI*integral(G(r)*sin(2*PI*sr))dr
  //using a simple discrete summation over the finite s sampling.

  //RIF is the reduced intensity function.  For monatomic solids: 
  //RIF(s) = (I(s)/Nf(s)^2-1)*s.  Gr is the reduced density function.
  //s here is 2*sin(theta)/lambda. The pair separation r is in the
  //same units as the wavelength lambda.

// Sort through the open images looking for a G(r) image
// ensure that the appropriate G(r) image is selected
image grimg
number nodocs=countdocumentwindowsoftype(5)
string imgname
getpersistentstringnote("ElectronDiffraction Tools:PDF:Temporary Data:Image name",imgname)

// Count the number of shown images to make sure more at least 1 is shown
if(nodocs<1) 
	{
		showalert("Ensure a G(r) image is present.",2)
		return
	}

number i, grfoundflag, markedimgs
string idstring
image tempimg

// Go through the shown images to identify the G(r) image and to flag any 
// CalcRIF images for deletion - if the ALT key is NOT held down ie the default is
// to always replace an existing CalcRIF image - if ALT is held down keep the old one as well
for(i=0; i<nodocs; i++)
	{
		imagedocument imgdoc=getimagedocument(i)
		image tempimg:=imgdoc.imagedocumentgetimage(0)
		getstringnote(tempimg, "Radial Distribution Function",idstring)

		if(idstring=="Calculated Reduced Interference Function" && !optiondown())
			{
				number imgid=tempimg.imagegetid()
				setpersistentnumbernote("ElectronDiffraction Tools:PDF:Temporary Data:"+markedimgs,imgid)
				markedimgs=markedimgs+1
			}
	}


// Delete the Calculated RIF image(s) if the ALT key is NOT held down
for(i=0; i<markedimgs;i++)
	{
		number imgdeleteid
		getpersistentnumbernote("ElectronDiffraction Tools:PDF:Temporary Data:"+i,imgdeleteid)
		
		if(imgdeleteid!=0)
			{
				image temp:=findimagebyid(imgdeleteid)
				//deleteimage(temp)
			}
	}

// Sort through the images again to find the G(r) image
number grflag
nodocs=countdocumentwindowsoftype(5)

for(i=0; i<nodocs; i++)
	{
		imagedocument imgdoc=getimagedocument(i)
		image tempimg:=imgdoc.imagedocumentgetimage(0)

		getstringnote(tempimg, "Radial Distribution Function",idstring)

		if(idstring=="Reduced Density Function")
			{
				grflag=1 // found the RDF image so set a flag
				grimg:=tempimg
				break
			}
	}

		// If the G(r) image has not been found the flag value is zero - bail out		
		if(grflag==0)
			{
				showalert("A G(r) image must be displayed.",2)
				return
			}

	// Sort through the open images to find the RIF image and from it obtain the xscale
	// value (smax)with which to calibrate the calculated RIF	
	number smax,sxsize, sysize
	string unitstring
	
	for(i=0; i<nodocs; i++)
		{
			imagedocument imgdoc=getimagedocument(i)
			image tempimg:=imgdoc.imagedocumentgetimage(0)
			getstringnote(tempimg, "Radial Distribution Function",idstring)

			if(idstring=="Reduced Interference Function"||idstring=="Intensity Profile")
				{
					smax=tempimg.imagegetdimensionscale(0)
					getsize(tempimg, sxsize, sysize)
					unitstring=tempimg.imagegetdimensionunitstring(0)
					break
				}
		}		
	
	// The maximum value of s in the original RIF image is the pixel calibration 
	// multiplied by the number of pixels in the image. If this value is zero for
	// any reason, the user is prompted to manually enter a value	
	smax=smax*sxsize
	
	if(smax==0) // Valid calibration data for s could not be sourced from the original RIF image
		{
			showalert("Ensure a RIF image is displayed.",2)
			return
		}	
	
	// numr is the resolution of the calculated plot
	number nums=dlggetvalue(self.lookupelement("resolutionfield"))

	// rmax is the maximum radius
	number rmax=dlggetvalue(self.lookupelement("radiusfield"))

	// If data are calibrated in Q then modify the scaling - this calculation expects 
	// data to be in terms of s	
	if(unitstring=="Q (1/A)") smax=smax/(2*pi())

	number r, numr
	getsize(grimg,numr,r);

	image phi, phi_matrix, Gr_matrix;
	image f_sr

    // decrement the number of integration points if not odd for Simpson's rule
	if((numr+1)%2!=0) numr-=1; 

	f_sr = exprsize(numr,nums,0);
	phi_matrix = exprsize(numr,nums,0);
	Gr_matrix = exprsize(numr,nums,0);
	phi = exprsize(nums,0);
	f_sr = sin(icol*1.0/numr*rmax*irow*2*pi()*1.0/nums*smax);
	f_sr /= (2*pi()); //constants in accordance with Cockayne & McKenzie Acta Cryst. 1988
    // below, the weightings of each icol are adjusted to implement Simpson's rule
    f_sr*=(4*(icol%2)+2*((icol+1)%2))/3.0*tert(icol==0,0,1)*tert(icol==numr-1,0,1)+tert(icol==0,1,0)+tert(icol==numr-1,1,0)

	slice2(Gr_matrix, 0, 0, 0, 0, numr, 1, 1, nums, 1) = grimg[icol,0] 
	phi_matrix = f_sr*Gr_matrix
	phi[irow,0] += slice2(phi_matrix, 0, 0, 0, 0, numr, 1, 1, nums, 1) //integration summation
	phi *= rmax/numr; //integration interval
	
	// If data are calibrated in Q then modify the scaling - the above calculation was done in terms
	// of s, so now change the scaling so that it is in terms of Q	
	if(unitstring=="Q (1/A)") ImageSetDimensionCalibration(phi,0,0,(smax/nums)*2*pi(),unitstring,1)
	else ImageSetDimensionCalibration(phi,0,0,(smax/nums),unitstring,1)	
	
	// Set the labelling on the image	
	imagesetintensityunitinfo(phi, "Q[S(Q)-1]",1)
	setname(phi,"RIF from G(r) - "+imgname)
	setstringnote(phi,"Radial Distribution Function","Calculated Reduced Interference Function")

	showimage(phi)
	setwindowposition(phi,202,84)
}



/*************************************************************************
Data button
*************************************************************************/
void databutton(object self)
{
image img:=GetFrontImage()
img.showimage()
	
ChooseMenuItem("ePDF Tools","","Data")
}

/*************************************************************************
Version information
*************************************************************************/
string createversionstring(object self)
	{		
		string todaysdate, microscopestring="1", temp
		getdate(0,todaysdate)

		number length=len(todaysdate)
		if(length<10) todaysdate="0"+todaysdate

		string year = right(todaysdate,4)
		string month = mid(todaysdate,3,2)
		string day = left(todaysdate,2)

		string currentdatestring=year+month+day
		return currentdatestring
	}

	/*************************************************************************
r_max is changed
*************************************************************************/
	void radiuschanged(object self, taggroup tg)
{
	number currentradius=tg.dlggetvalue()
	setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:R_max",currentradius)
	
	number minradius,defaultradiuslimit
	getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:R_max",defaultradiuslimit)
	if(currentradius<minradius)
		{
			tg.dlgvalue(defaultradiuslimit)
		}
}

/*************************************************************************
Density is changed
*************************************************************************/
void densitychanged(object self, taggroup tg)
{
	number currentdensity=tg.dlggetvalue()
	number minimumdensity
	getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Minimum Density",minimumdensity)
	
	setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Average Density",currentdensity)
	
	if(currentdensity<minimumdensity)
		{
			tg.dlgvalue(minimumdensity)
			setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Average Density",minimumdensity)
		}
		
		// disable the J(r), G>R and Data buttons to force a recalculation of G(r)		
		self.setelementisenabled("jrbutton",0)
		self.setelementisenabled("grtorifbutton",0)
		self.setelementisenabled("databutton",0)
}

/*************************************************************************
Resolution is changed
*************************************************************************/
void resolutionchanged(object self, taggroup tg)
{
	number currentresolution=tg.dlggetvalue()
	setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Resolution",currentresolution)
	
	number minresolution
	getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Minimum Resolution",minresolution)
	if(currentresolution<minresolution)
		{
			tg.dlgvalue(minresolution)
		}
		
		// disable the J(r), G>R and Data buttons to force a recalculation of G(r)		
		self.setelementisenabled("jrbutton",0)
		self.setelementisenabled("grtorifbutton",0)
		self.setelementisenabled("databutton",0)
}


/*************************************************************************
The total atom percent
*************************************************************************/
void atfieldcalc(object self, taggroup tg)
{
	number thisatwt, massaccumulator, thiswtpcnt, thisatpcnt, sumwtpcnt, sumatpcnt, noofelements, summoles, i

	// Get the number of rows of elemental data to check
	getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Number of Elements",noofelements)
	
	sumatpcnt=0
	// Calculate the total mass in the specimen
	for(i=1; i<noofelements+1;i++)
		{
			thisatpcnt=dlggetvalue(self.lookupelement("atompcntfield"+i))
			sumatpcnt=sumatpcnt+thisatpcnt
		}

	self.validateview()
	return
}


/*************************************************************************
Atom weight calculation
*************************************************************************/
void wtfieldcalc(object self, taggroup tg)
{
	number thisatwt, massaccumulator, thiswtpcnt, thisatpcnt, sumwtpcnt, sumatpcnt, noofelements, summoles, i

	// Get the number of rows of elemental data to check
	getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Number of Elements",noofelements)

	// Calculate the total wt% based on the current values
	sumwtpcnt=0
	for(i=1; i<noofelements+1;i++)
		{

		thiswtpcnt=dlggetvalue(self.lookupelement("wtpcntfield"+i))
		sumwtpcnt=sumwtpcnt+thiswtpcnt

		}

	self.validateview()
	return
}


/*************************************************************************
Set flag state
*************************************************************************/
void setflagstatus(object self, number elementno, number statusvalue)
	{
		taggroup container=self.lookupelement("colourflag"+elementno)
		taggroup bitmap=container.dlggetelement(0)
		image temp:=rgbimage("",4,18,18)
		rgbnumber colour

		if(statusvalue==0) // this menas that the pulldown had been set to the null value ("Package Name" - set the flag to grey
			{
				colour=rgb(120,120,120)
			}

		if(statusvalue==1) // this means that the package status is current - set the flag to green
			{
				colour=rgb(0,255,0)
			}
		
		if(statusvalue==2) // this menas that the package status is not current - set the flag to red
			{
				colour=rgb(255,0,0)
			}

		temp=colour
		bitmap.DLGBitmapData(temp)
	}

/*************************************************************************
Close ePDF tools window
*************************************************************************/
void AboutToCloseDocument( object self, number test)
{
	number xpos, ypos
	documentwindow dialogwindow=getframewindow(self)
	windowgetframeposition(dialogwindow, xpos, ypos)

	//软件窗口在可视范围内
	if(OptionDown())
	{
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Dialog Position:Left",xpos)
		setpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Dialog Position:Right",ypos)
	}
}

/*************************************************************************
function supplywtoratpercent()
*************************************************************************/
void SupplyWtorAtPercent(object self, taggroup tg) // When the 'Calculate Using' radio button is changed the atorwtercent variable is set
	{
		number atorwtpercent = tg.DLGGetValue()
	}
	
	
//---------------------------------------生成  UI  界面---------------------------------------
/*************************************************************************
生成Composition box
*************************************************************************/
TagGroup CompositionBox(object self)
{
//-----------化学组分-------------
number i
TagGroup box1_items
TagGroup box1 = DLGCreateBox("Composition", box1_items).dlgexternalpadding(0,10)

TagGroup element_items //用于盛放化学组分
TagGroup elementbox = DLGCreateBox("", element_items)

// Create the top row of pulldown, fields and labels
taggroup elementlabel=dlgcreatelabel("Symb.").dlgexternalpadding(6,0).dlginternalpadding(0,0).dlgexternalpadding(4,9)
taggroup atompcntlabel=dlgcreatelabel("At. %").dlgexternalpadding(5,0).dlginternalpadding(0,0)
taggroup wtpcntlabel=dlgcreatelabel("Wt. %").dlgexternalpadding(6,0).dlginternalpadding(0,0)
taggroup atomnolabel=dlgcreatelabel("No.").dlgexternalpadding(6,0).dlginternalpadding(0,0)

taggroup atomwtlabel=dlgcreatelabel("Mass").dlgexternalpadding(8,0).dlginternalpadding(0,0)

taggroup group1=dlggroupitems(elementlabel,atompcntlabel, wtpcntlabel, atomnolabel).dlgtablelayout(4,1,0).dlganchor("West")
taggroup group2=dlggroupitems(group1, atomwtlabel).dlgtablelayout(2,1,0).dlganchor("West")
element_items.DLGAddElement(group2).dlganchor("West")

//多元素，自动生成列表
taggroup elementgroup
for(i=1; i<noofelements+1; i++)
	{
		taggroup elementmenu=self.makeElementfield(i).dlgexternalpadding(0,0)
		taggroup atompcntfield=self.makeatompcntfield(i).dlgexternalpadding(4,0)//原子数比
		taggroup wtpcntfield=self.makeweightpcntfield(i)  //质量比		
		taggroup atomnofield=self.makeatomnofield(i).dlgexternalpadding(4,0)   //原子序数
		
		taggroup atomwtfield=self.makeatomweightfield(i).dlgexternalpadding(0,0) //原子量		
		
		taggroup fieldgroup1=dlggroupitems(elementmenu, atompcntfield, wtpcntfield, atomnofield).dlgtablelayout(4,1,0)
		elementgroup=dlggroupitems(fieldgroup1, atomwtfield).dlgtablelayout(2,1,0).dlganchor("West")
		
		element_items.dlgaddelement(elementgroup)		
	}

//-----------基本参数-----------------	
	
	TagGroup Parameters_items //用于加速电压、calc
	TagGroup Parametersbox = DLGCreateBox("", Parameters_items)

	//加速电压
	number defaultvoltage
	getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Beam Voltage (kV)",defaultvoltage)
	taggroup voltslabel=DLGCreatelabel("   kV  ")
	taggroup voltagefield = DLGCreateRealField(defaultvoltage, 6, 0).DLGAnchor("Center").dlgidentifier("voltagefield").dlgexternalpadding(0,0)
	taggroup voltagegroup=dlggroupitems(voltslabel, voltagefield).dlgtablelayout(2,1,0)
	
	//密度
		number defaultdensityvalue
		getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Average Density",defaultdensityvalue)
		taggroup densitylabel=DLGCreatelabel("Dens.")
		taggroup densityfield = DLGCreaterealField(defaultdensityvalue,6,3).dlgchangedmethod("densitychanged").DLGAnchor("Center").dlgidentifier("densityfield").dlgexternalpadding(0,0)
		taggroup densitygroup=dlggroupitems(densitylabel, densityfield).dlgtablelayout(2,1,0).dlgexternalpadding(0,0)
	
	//R-max
		number defaultradiuslimit
		getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:R_max",defaultradiuslimit)
		taggroup radiuslabel=DLGCreatelabel("r_max")
		taggroup radiusfield = DLGCreaterealField(defaultradiuslimit,6,1).dlgchangedmethod("radiuschanged").DLGAnchor("Center").dlgidentifier("radiusfield").dlgexternalpadding(0,0)
		taggroup radiusgroup=dlggroupitems(radiuslabel, radiusfield).dlgtablelayout(2,1,0).dlgexternalpadding(0,0)

	// The resolution of the G(r) plot in pixels		
		number defaultresolution
		getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Resolution",defaultresolution)
		taggroup resolutionlabel=DLGCreatelabel("  Size ")
		taggroup resolutionfield = DLGCreateIntegerField(defaultresolution,6).dlgchangedmethod("resolutionchanged").DLGAnchor("Center").dlgidentifier("resolutionfield").dlgexternalpadding(0,0)
		taggroup resolutiongroup=dlggroupitems(resolutionlabel, resolutionfield).dlgtablelayout(2,1,0).dlgexternalpadding(0,0)
		
	taggroup ParametersGroup=dlggroupitems(voltagegroup, densitygroup,radiusgroup,resolutiongroup).dlgtablelayout(1,4,0).dlgexternalpadding(0,6)
	Parameters_items.dlgaddelement(ParametersGroup)
	
//	组合后添加入box1中
	taggroup totgroup=dlggroupitems(elementbox,Parametersbox).dlgtablelayout(2,1,0).dlganchor("North")
	box1_items.dlgaddelement(totgroup)

//-------------其它----------------
	// 质量比和原子比的选择
	number radiovalue, invradiovalue
	getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Use At% (0) or Wt% (1)", radiovalue)

	TagGroup radio_items
	taggroup radioList = DLGCreateRadioList( radio_items, radiovalue, "SupplyWtorAtpercent" ).dlgidentifier("atwtradiosetting").dlgexternalpadding(0,0)

	if(radiovalue==1) invradiovalue=0
	else invradiovalue=1

	radio_items.DLGAddElement( DLGCreateRadioItem( "At.%", 0) ).DLGSide("Left")
	radio_items.DLGAddElement( DLGCreateRadioItem( "Wt.%", 1 ) ).DLGSide("Left")	

	// 总比例
	taggroup sumlabel=dlgcreatelabel("% Sums")
	taggroup sumatpcntfield=dlgcreaterealfield(0,6,4).dlgidentifier("sumatpcnt")
	taggroup sumwtpcntfield=dlgcreaterealfield(0,6,4).dlgidentifier("sumwtpcnt").dlgexternalpadding(4,0)
	taggroup sumfieldgroup=dlggroupitems(sumlabel, sumatpcntfield,sumwtpcntfield,radioList).dlgtablelayout(4,1,0).dlgexternalpadding(0,0).dlginternalpadding(6,0).DLGAnchor("West")

	return box1;
}

/*************************************************************************
生成PDF tools box
*************************************************************************/
taggroup CreatePDFbox(object self)
{
		taggroup PDFbox_items
		taggroup PDFbox=dlgcreatebox("PDF tools", PDFbox_items).dlginternalpadding(10,8).dlgexternalpadding(0,0)
		
		TagGroup CalculateButton = PDFbox_items.dlgaddelement(DLGCreatePushButton("Calc.", "Calculate")).dlgenabled(0).dlgidentifier("calculatebutton").dlgexternalpadding(2,2).dlginternalpadding(3,0)
		taggroup Scatterbutton=PDFbox_items.dlgaddelement(DLGCreatePushButton("Scatt.","Scatter")).dlginternalpadding(2,0).dlgexternalpadding(2,2).dlgidentifier("Scatterbutton").dlgenabled(0)
		taggroup BGbutton=PDFbox_items.dlgaddelement(dlgcreatepushbutton("BG", "BGRemoved")).dlgexternalpadding(2,2).dlginternalpadding(8,0)
		taggroup Fqbutton=PDFbox_items.dlgaddelement(dlgcreatepushbutton("F(Q)", "Iq2Fq")).dlgexternalpadding(2,2).dlginternalpadding(6,0)
		taggroup grbutton=PDFbox_items.dlgaddelement(dlgcreatepushbutton("G(r)", "grbutton")).dlgexternalpadding(2,2).dlginternalpadding(9,0)
		taggroup Calibrationbutton=PDFbox_items.dlgaddelement(dlgcreatepushbutton("Calib.", "grcalibration")).dlgexternalpadding(2,2).dlginternalpadding(3,0)
		
		taggroup gr2Grbutton=PDFbox_items.dlgaddelement(dlgcreatepushbutton("->G(r)", "gr2reducedGr")).dlgexternalpadding(2,2).dlginternalpadding(1,0)
		taggroup F2Ibutton=PDFbox_items.dlgaddelement(dlgcreatepushbutton("F->I", "Fq2Iqbutton")).dlgexternalpadding(2,2).dlginternalpadding(8,0)
		taggroup gr2RDFbutton=PDFbox_items.dlgaddelement(dlgcreatepushbutton("Iq > rif", "Profile2RIF")).dlgexternalpadding(2,2).dlginternalpadding(1,0)
		taggroup Iq2SAEDbutton=PDFbox_items.dlgaddelement(dlgcreatepushbutton("+ ROI", "AddROI")).dlgexternalpadding(2,2).dlginternalpadding(3,0)		
		taggroup databutton=PDFbox_items.dlgaddelement(dlgcreatepushbutton("Data", "databutton")).dlgexternalpadding(2,2).dlginternalpadding(5,0)		
		
		taggroup closebutton=PDFbox_items.dlgaddelement(dlgcreatepushbutton("Prefer.", "SetPreferences")).dlgexternalpadding(2,2).dlginternalpadding(0,0)	

		PDFbox.DLGTableLayOut(6,2,0)

		return PDFbox
}

/*************************************************************************
生成Image processing box
*************************************************************************/
taggroup ImageProcessingBox(object self)
{
taggroup ImageProcessingBox_items
taggroup ImageProcessingBox=dlgcreatebox("Image processing", ImageProcessingBox_items).dlginternalpadding(10,10).dlgexternalpadding(0,10)

ImageProcessingBox_items.DLGAddElement(DLGCreatePushButton("Center","center").dlgexternalpadding(0,2).dlginternalpadding(3,0))
ImageProcessingBox_items.DLGAddElement(DLGCreatePushButton("Profile", "Profile_Intensity").dlgexternalpadding(0,2).dlginternalpadding(0,0))
ImageProcessingBox_items.DLGAddElement(DLGCreatePushButton("Calib.", "Calibrate").dlgexternalpadding(0,2).dlginternalpadding(4,0))
ImageProcessingBox_items.DLGAddElement(DLGCreatePushButton("Slices","Profile2sliceimage").dlgexternalpadding(0,2).dlginternalpadding(4,0))
ImageProcessingBox_items.DLGAddElement(DLGCreatePushButton("Smooth","SmoothIq").dlgexternalpadding(0,2).dlginternalpadding(1,0))
ImageProcessingBox_items.DLGAddElement(DLGCreatePushButton("Repl.","ReplaceIq").dlgexternalpadding(0,2).dlginternalpadding(1,0))

ImageProcessingBox_items.DLGAddElement(DLGCreatePushButton("ROIImg","ROIImage").dlgexternalpadding(0,2).dlginternalpadding(0,0))
ImageProcessingBox_items.DLGAddElement(DLGCreatePushButton("Zoom", "ZoomImage").dlgexternalpadding(0,2).dlginternalpadding(1,0))
ImageProcessingBox_items.DLGAddElement(DLGCreatePushButton(">Calib.","GetCalibration").dlgexternalpadding(0,2).dlginternalpadding(0,0))
ImageProcessingBox_items.DLGAddElement(DLGCreatePushButton("SaveGr","SaveGr").dlgexternalpadding(0,2).dlginternalpadding(0,0))
ImageProcessingBox_items.DLGAddElement(DLGCreatePushButton("SaveXY","XY").dlgexternalpadding(0,2).dlginternalpadding(0,0))
ImageProcessingBox_items.DLGAddElement(DLGCreatePushButton("Close","CloseImages").dlgexternalpadding(0,2).dlginternalpadding(0,0))

ImageProcessingBox.DLGTableLayOut(6,2,0)

return ImageProcessingBox
}

/*************************************************************************
版权申明
*************************************************************************/
Taggroup addfooter(object self)
	{
		TagGroup VersionLabel = DLGCreateLabel("H.L.Shi (honglongshi@outlook.com), v1.0, Jan. 2017.");
		VersionLabel.dlgexternalpadding(0,5)
		return VersionLabel
	}	

/*************************************************************************
// Puts all the dialog components together as a single taggroup
*************************************************************************/
taggroup MakePDFdialog(object self)
	{
		TagGroup dialog_items;	
		TagGroup PDFdialog = DLGCreateDialog("", dialog_items)
	
		dialog_items.dlgaddelement(self.CompositionBox()) //参数设置box
		dialog_items.dlgaddelement(self.CreatePDFbox())  //RDF box
		dialog_items.dlgaddelement(self.ImageProcessingBox())
		dialog_items.dlgaddelement(self.addfooter())  //脚注 box	
		
		return PDFdialog
	}	
	
/*************************************************************************
// Dialog constructor
*************************************************************************/	
PDFcalculator(object self)
	{
		self.init(self.MakePDFdialog())
		self.display("ePDF Tools")		
			
		// Position the dialog
		documentwindow dialogwin=getdocumentwindow(0)
		number xpos, ypos
		getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Dialog Position:Left",xpos)
		getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Dialog Position:Right",ypos)

		// if the dialog position persistant tags have been set position the dialog
		if(xpos!=0 && ypos!=0) 
			{
				windowsetframeposition(dialogwin, xpos, ypos)
			}
	}
	
/*************************************************************************
// Dialog destructor
*************************************************************************/	
~PDFcalculator(object self)
	{

	}

}

	

/*************************************************************************
// Main program
*************************************************************************/

alloc(PDFcalculator)