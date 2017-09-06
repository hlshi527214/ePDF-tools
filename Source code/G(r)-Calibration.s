 // 以自动寻找的数据点为基础，进行手动实时拟合背底
 //确保存在ElectronDiffraction Tools:Fit:Smooth factor的缓存，否则运行以下两列
 //number s=2,start=120,end=1200
//setpersistentnumbernote("ElectronDiffraction Tools:Fit:Smooth factor",s)
//setpersistentnumbernote("ElectronDiffraction Tools:Fit:Start",start)
//setpersistentnumbernote("ElectronDiffraction Tools:Fit:End",end)
 


number token1   //listener1  用于roichanged
number token2  //listener2  用于取消listener
number token3  //用于添加roi

// 定义类
class ImageDisplayEventListener : object
{
//定义变量		
ROI 	theroi, theroi2
number  left, right, SpecialROIID
image spectrum
	
	// 扑捉ROI	
	ROI GetWin(object self) 
		{
			return theroi
		}


void ROIChanged( Object self, Number event_flags, ImageDisplay imgdisp, Number roi_change_flags, Number roi_disp_change_flags, ROI theroi )
{		
number xsize, ysize,atomicdensity,maxradius,resolution,xscale,l,r,slope,x,y,ming,roi_r,k,Molweight,Positiveg
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Atomic Density",atomicdensity)
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:R_max",maxradius)
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Resolution",resolution)
getpersistentnumbernote("ElectronDiffraction Tools:PDF:Temp:Molecular weight", Molweight)	
getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Positive g(r)", Positiveg)

spectrum.getsize(xsize,ysize)
xscale=spectrum.ImageGetDimensionScale(0)

//先转为g(r)
image g=spectrum*0
image gr=spectrum*0
image new=abs(spectrum)
g=1+spectrum/(4*pi()*1*(icol*xscale)) //假设密度为1

//再归一化
number thisroiid=theroi.roigetid()
if(thisroiid==SpecialROIID)
{
theroi.roigetrange(l,r)

k=g.GetPixel(r,0)
g=g-k
g=g/abs(1-k)
g.SetPixel(0,0,0)

spectrum.SetNumberNote("ElectronDiffraction Tools:PDF:Parameters:Density",1-k)

if(Positiveg==1)  //roi左侧为0
{
g[0,0,1,l]=0
}

number numbDens=abs(spectrum.GetPixel(r,0))/(4*pi()*xscale*r)
number density=10*numbDens*Molweight/6.02214
theROI.roisetlabel("Dens.="+format(density,"%6.2f")+" \nNumber= "+format(numbDens,"%6.4f"))  //显示密度值

gr=(4*pi()*atomicdensity*(icol*xscale))*(g-1)
gr.SetPixel(0,0,0)

if(OptionDown())  //自动设置roi的位置
{
theroi.roigetrange(l,r)

//寻找左侧的位置
number xmin,ymin,xmax,ymax
g[0,l,1,r].max(xmax,ymax)
xmax=xmax+l

g[0,l,1,xmax].min(xmin,ymin)
xmin=l+xmin

g[0,r-10,1,r+10].min(x,y)
r=r-10+x

number Ii,Ii1
For(number i=xmin;i<xmax;i++)
{
Ii=g.GetPixel(i,0)
Ii1=g.GetPixel(i+1,0)

if(Ii<0&&Ii1>=0) theroi.roisetrange(i,r) //设置roi的位置
}
}

//更新图层g(r)
imgdisp{"g(r)"}=g
imgdisp{"G(r)"}=gr

//imgdisp.lineplotimagedisplaysetdoautosurvey(0,0)
//imgdisp.lineplotimagedisplaysetcontrastlimits(1.1*min1,-1.1*min1)
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

number atomicdensity,maxradius,resolution,Positiveg
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Atomic Density",atomicdensity)
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:R_max",maxradius)
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Resolution",resolution)
getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Positive g(r)", Positiveg)

string imgname
spectrum.GetStringNote("ElectronDiffraction Tools:PDF:Image name",imgname)

//移除token
imgdisp.ImageDisplayRemoveEventListener(token1)
imgdisp.ImageDisplayRemoveEventListener(token2)

number dens
spectrum.getNumberNote("ElectronDiffraction Tools:PDF:Parameters:Density",dens)
image normGr=spectrum/dens
normGr.ShowImage()

normGr.SetName("Normalized G(r) - "+imgname)
normGr.ImageCopyCalibrationFrom(spectrum)
SetStringNote(normGr,"ElectronDiffraction Tools:PDF:Image name",imgname)  
setstringnote(normGr,"Radial Distribution Function","Normalized Reduced Density Function")
imagesetIntensityUnitString(normGr,"Normalized G(r)")
SetWindowPosition( normGr, 190,30 )

if(Positiveg==1) //衰减至零
{
number ir=normGr.getpixel(r,0)   // ROI的左侧趋于零
number k=ir/r

for(number i=0;i<l;i++)  //G(r)线性变化
{
normGr.setpixel(i,0,k*i)
}
}


//删除roi
number roino= imgdisp.ImageDisplayCountROIS()
for (number i=0; i<roino; i++)
{
number index
ROI currentROI = imgdisp.ImageDisplayGetROI( index )
imgdisp.imagedisplaydeleteROI(currentROI)
}

image gr=imgdisp{"g(r)"}
gr.ShowImage()
ImageDisplay grdisp=gr.ImageGetImageDisplay(0)

number x,y,xsize,ysize,xscale
gr.GetSize(xsize,ysize)
gr.minmax(x,y)
grdisp.LinePlotImageDisplaySetDoAutoSurvey( 0, 0 ) 
grdisp.LinePlotImageDisplaySetContrastLimits(-0.1*y,1.1*y)

gr.SetName("g(r) - "+imgname)
gr.ImageCopyCalibrationFrom(spectrum)
SetStringNote(gr,"ElectronDiffraction Tools:PDF:Image name",imgname)  
setstringnote(gr,"Radial Distribution Function","Density Function")
imagesetIntensityUnitString(gr,"g(r)")
SetWindowPosition( gr, 190,30 )

image RDF=gr*0
xscale=spectrum.ImageGetDimensionScale(0)

RDF=4*Pi()*atomicdensity*(icol*xscale)**2*gr
RDF.ShowImage()
RDF.SetName("RDF - "+imgname)
RDF.ImageCopyCalibrationFrom(spectrum)
SetStringNote(RDF,"ElectronDiffraction Tools:PDF:Image name",imgname)  
setstringnote(RDF,"Radial Distribution Function","Radial Distribution Function")
imagesetIntensityUnitString(RDF,"RDF(r)")
SetWindowPosition(RDF, 190,30 )

image ReducedGr=gr*0
ReducedGr=4*Pi()*atomicdensity*(icol*xscale)*(gr-1)
ReducedGr.ShowImage()
ReducedGr.SetName("Reduced G(r) - "+imgname)
ReducedGr.ImageCopyCalibrationFrom(spectrum)
SetStringNote(ReducedGr,"ElectronDiffraction Tools:PDF:Image name",imgname)  
setstringnote(ReducedGr,"Radial Distribution Function","Reduced Density Function")
imagesetIntensityUnitString(ReducedGr,"Reduced G(r)")
SetWindowPosition(ReducedGr, 190,30 )

normGr.ShowImage()
gr.ShowImage()

//spectrum.deleteimage()
}
	
		
// This function sources the ROI on the passed in image
		object init(object self, image front)
		{
			spectrum:=front
			imagedisplay imgdisp=spectrum.ImageGetImageDisplay(0)
			theroi = imgdisp.ImageDisplayGetROI(0)			
			
			// 把第一个roi设置特征roi			
			SpecialROIID=theroi.roigetid()
			//theroi.roisetrange(20,70)
			theroi.roisetvolatile(0)
			theroi.roisetcolor(1,0,0)
			theroi.roisetlabel("Range")			
			
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



// function to check the displayed image is appropriate and create and apply
// the event listeners

void main()
{
image img:=getfrontimage()
imagedisplay imgdisp=imagegetimagedisplay(img, 0)
image splineimg=img*0
string imgname=img.getname()
number xsize, ysize
getsize(img, xsize, ysize)

//提示语
Result("---------------------------------------------------------------------------\nG(r) calibration:\n\n1) Mark the first physical peak using ROI.\n2) ALT + ROI: Set ROI to the optimum position.\n---------------------------------------------------------------------------\n")

//先添加图层
//删除已有图层
number noslices=imgdisp.imagedisplaycountslices()
for(number i=noslices-1; i>-1; i--)
{
object sliceid=imgdisp.imagedisplaygetsliceidbyindex(i)
string slicename=imgdisp.imagedisplaygetslicelabelbyid(sliceid)
if(slicename=="G(r)"||slicename=="g(r)") imgdisp.imagedisplaydeleteslicewithid(sliceid)
}

object sliceid=imgdisp.imagedisplaygetsliceidbyindex(0)
imgdisp.imagedisplaysetslicelabelbyid(sliceid, "Raw")

image slice1=img*0
imgdisp.imagedisplayaddimage(slice1, "G(r)")  //图层1：存放New
lineplotimagedisplaysetslicedrawingstyle(imgdisp,1,1)
lineplotimagedisplaysetslicecomponentcolor(imgdisp, 1, 0,256,0,0)

image slice2=img*0
imgdisp.imagedisplayaddimage(slice2, "g(r)")  //图层2：存放damped F(Q)
imgdisp.LinePlotImageDisplaySetLegendShown(1)
lineplotimagedisplaysetslicedrawingstyle(imgdisp,2,1)
lineplotimagedisplaysetslicecomponentcolor(imgdisp,2, 0,0,0,1)

number t,l,b,r
img.GetSelection(t,l,b,r)
img.ClearSelection()

number x,y
img[0,r-10,1,r+10].min(x,y)
r=r-10+x

//添加roi
roi roi1=createroi()
ROI1.ROISetrange(l,r)
roi1.roisetvolatile(0)
imgdisp.ImageDisplayAddROI( ROI1)

// listener-image关联
// Listener for ROI removal
string messagemap1="roi_removed:ROIRemoved"
object ROIRemovalListener=alloc(ImageDisplayEventListener).init(img)
token1 =imgdisp.ImageDisplayAddEventListener( RoiRemovalListener, messagemap1)

// Listener for ROI change
string messagemap2="roi_changed:ROIChanged"
object ROIChangeListener=alloc(ImageDisplayEventListener).init(img)
token2 =imgdisp.ImageDisplayAddEventListener( RoiChangeListener, messagemap2)		
}


// Main script
main()
