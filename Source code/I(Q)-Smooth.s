//实时拟合eels背底：指数函数

number token1, token2

class ImageDisplayEventListener : object
{
// Some variables
ROI 	theroi
number 	left, right, SpecialROIID,counter
image spectrum

//Guassian 平滑函数
realimage smooth_Is(object self, realimage img, number fwhm, number integ)
{
  realimage Gauss, Gauss_matrix, img_matrix, blur_img, kernal;
  number xsize,ysize, norm, alpha;
  
  //缩小图像
  img.GetSize(xsize,ysize)
  number scale=1/5 //只用缩略图
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



// A function which returns the ROI
ROI GetWin(object self)
{
return theroi
}


// ROIChanged
void ROIChanged( Object self, Number event_flags, ImageDisplay imgdisp, Number roi_change_flags, Number roi_disp_change_flags, ROI theroi )
{
counter=counter+1  //用于记录执行的次数

//只对特征roi操作
number roiflag=1

number thisroiid=theroi.roigetid()
if(thisroiid!=SpecialROIID) //不是特征，显示roi的值
{
number left, right,xscale,flagBG

getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Auto BG", flagBG)

xscale=spectrum.ImageGetDimensionScale(0)
theroi.roigetrange(left, right)
} //不是特征roi操作结束


if(thisroiid==SpecialROIID) //特征roi，指数背底
{
// Source the position of the ROI
number left, right,xscale,xsize,ysize
spectrum.GetSize(xsize,ysize)
xscale=spectrum.ImageGetDimensionScale(0)
theroi.roigetrange(left, right)
number width=right-left
theroi.roisetlabel("Width="+width)

// 拟合背底
if(counter==1)
{
number fwhm=0.1*abs(right-left)
image background=self.smooth_Is(spectrum, fwhm, xsize)
imgdisp.lineplotimagedisplaysetlegendshown(1)
imgdisp{"Background"}=background

if(OptionDown())
{
image new=spectrum
number flag=1
GetNumber("Which part of the ROI will be replaced, left, right, or all (0, 1, 2)",flag,flag)

if(flag==0)  //left
{
new[0,0,1,left]=background[0,0,1,left]
}

if(flag==1)//right
{
new[0,right,1,xsize]=background[0,right,1,xsize]
}

if(flag==2)
{
new[0,0,1,xsize]=background[0,0,1,xsize]
}

imgdisp{"Smoothed"}=new
}

}

else if(counter>1&&counter<=3) //不执行
{

}

else if(counter>3)
{
counter=0
}

}
}

// ROI is removed
void ROIRemoved( Object self, Number event_flags, ImageDisplay imgdisp, Number roi_change_flags, Number roi_disp_change_flags, ROI theroi )
{
// 只对特征roi操作
number thisroiid=theroi.roigetid()
if(thisroiid!=SpecialROIID) return

image img=imgdisp{"Smoothed"}
img.imagecopycalibrationfrom(spectrum)
img.setname("Gaussian smoothed")
//setstringnote(img,"Radial Distribution Function","Reduced Density Function")	
img.ShowImage()


//删除两个listener，删除图层
imgdisp.ImageDisplayRemoveEventListener(token1)
imgdisp.ImageDisplayRemoveEventListener(token2)
number noslices=imgdisp.imagedisplaycountslices()
number i
for(i=noslices-1; i>-1; i--)
{
object sliceid=imgdisp.imagedisplaygetsliceidbyindex(i)
string slicename=imgdisp.imagedisplaygetslicelabelbyid(sliceid)
if(slicename=="Background" ) imgdisp.imagedisplaydeleteslicewithid(sliceid)
}

spectrum.deleteimage()
}

// 初始化
object init(object self, image front)
{
spectrum:=front
imagedisplay imgdisp=spectrum.ImageGetImageDisplay(0)
theroi = imgdisp.ImageDisplayGetROI(0)

// 设置特征roi
SpecialROIID=theroi.roigetid()
theroi.roisetvolatile(0)
theroi.roisetcolor(1,0,0)
theroi.roisetlabel("Width")

// 特征roi的响应

number dummy=0
self.roichanged(dummy, imgdisp, dummy, dummy, theroi)
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


// 主函数
void main()
{
// source the front-most image and check that it is a 1D profile
image img:=GetFrontImage()

image front=img
front.SetName("Gaussian BG")
front.ShowImage()
imagedisplay imgdisp=front.imagegetimagedisplay(0)

front.ImageCopyCalibrationFrom(img)
number xscale=img.ImageGetDimensionScale(0)

string isqstring=front.imagegetdimensionunitstring(0)
if(isqstring=="1/?")   //单位为1/A
{
xscale=front.imagegetdimensionscale(0)
front.imagesetdimensionscale(0,xscale*2*pi())
front.imagesetdimensionunitstring(0,"Q (1/A)")
}

if(isqstring=="1/nm")  //单位为1/nm
{
xscale=front.imagegetdimensionscale(0)
front.imagesetdimensionscale(0,xscale/10*2*pi())
front.imagesetdimensionunitstring(0,"Q (1/A)")
}

if(isqstring=="Q (1/A)")  //单位为Q (1/A)
{
xscale=front.imagegetdimensionscale(0)
front.imagesetdimensionscale(0,xscale)
front.imagesetdimensionunitstring(0,"Q (1/A)")
}

number xsize, ysize
getsize(front, xsize, ysize)

Result("-----------------------------------------------------------------------------------------\nAlt+ROI: Add the right of ROI to Intensity profile\nCtrl+ROI: to set auto background.\n-----------------------------------------------------------------------------------------\n")

//删除已有图层
number noslices=imgdisp.imagedisplaycountslices()
for(number i=noslices-1; i>-1; i--)
{
object sliceid=imgdisp.imagedisplaygetsliceidbyindex(i)
string slicename=imgdisp.imagedisplaygetslicelabelbyid(sliceid)
if(slicename=="Background") imgdisp.imagedisplaydeleteslicewithid(sliceid)
}

object sliceid=imgdisp.imagedisplaygetsliceidbyindex(0)
imgdisp.imagedisplaysetslicelabelbyid(sliceid, "Raw")
imgdisp.LinePlotImageDisplaySetLegendShown(1)

//添加新图层
image slice1=img*0
imgdisp.imagedisplayaddimage(slice1, "Background")  //图层1：background
lineplotimagedisplaysetslicedrawingstyle(imgdisp,1,1)
lineplotimagedisplaysetslicecomponentcolor(imgdisp, 1, 0,1,0,0)

image slice2=img*0
imgdisp.imagedisplayaddimage(slice2, "Smoothed")  //图层2：smoothed
lineplotimagedisplaysetslicedrawingstyle(imgdisp,2,1)
lineplotimagedisplaysetslicecomponentcolor(imgdisp, 2, 0,1,0,1)

// 添加roi
number t,l,b,r
img.GetSelection(t,l,b,r)
if(l!=0&&r!=xsize)front.setselection(t,l,b,r)

else front.setselection(0,0.1*xsize,1,0.1*xsize+50)

// listener-image关联
// Listener for ROI removal
string messagemap1="roi_removed:ROIRemoved"
object ROIRemovalListener=alloc(ImageDisplayEventListener).init(front)
token1 = imgdisp.ImageDisplayAddEventListener( RoiRemovalListener, messagemap1)

// Listener for ROI change
string messagemap2="roi_changed,roi_Added:ROIChanged"
object ROIChangeListener=alloc(ImageDisplayEventListener).init(front)
token2 = imgdisp.ImageDisplayAddEventListener( RoiChangeListener, messagemap2)
}

// Main script
main()