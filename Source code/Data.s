 // 根据曲线属性自动计算相应参数
 


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
image RDF=spectrum*0

number xsize, ysize,atomicdensity,maxradius,resolution,xscale,l,r,slope
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Atomic Density",atomicdensity)
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:R_max",maxradius)
Getpersistentnumbernote("ElectronDiffraction Tools:PDF:Default Values:Resolution",resolution)

number flag,r1,r2,angle,coord
string CurveLabel
Getstringnote(spectrum,"Radial Distribution Function",CurveLabel)
if(CurveLabel=="Density Function")flag=0
if(CurveLabel=="Reduced Density Function")flag=1
if(CurveLabel=="Radial Distribution Function")flag=2

spectrum.getsize(xsize,ysize)
xscale=spectrum.ImageGetDimensionScale(0)

number thisroiid=theroi.roigetid()
if(thisroiid==SpecialROIID)
{
theroi.roigetrange(l,r)
r1=l*xscale
r2=r*xscale
angle=2*(asin(r2/r1/2))*(180/pi())

if(flag==0) //density function
{
//计算RDF
RDF=4*Pi()*atomicdensity*(icol*xscale)**2*spectrum

coord=RDF[0,l,1,r].sum()*xscale
theroi.roisetlabel("r1= "+format(r1,"%4.3f")+" A \nr2= "+format(r2,"%4.3f")+" A \nCoord.= "+format(coord,"%4.2f")+"\nAngle= "+format(angle,"%4.1f")+" deg." )
}

if(flag==1) //Reduced density function
{
//计算RDF
RDF=4*Pi()*atomicdensity*(icol*xscale)**2*(1+spectrum/(4*Pi()*atomicdensity*(icol*xscale)))

coord=RDF[0,l,1,r].sum()*xscale
theroi.roisetlabel("r1= "+format(r1,"%4.3f")+" A \nr2= "+format(r2,"%4.3f")+" A \nCoord.= "+format(coord,"%4.2f")+"\nAngle= "+format(angle,"%4.1f")+" deg." )
}

if(flag==2)  //RDF
{
coord=spectrum[0,l,1,r].sum()*xscale
theroi.roisetlabel("r1= "+format(r1,"%4.3f")+" A \nr2= "+format(r2,"%4.3f")+" A \nCoord.= "+format(coord,"%4.2f")+"\nAngle= "+format(angle,"%4.1f")+" deg." )
}
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
imgdisp.ImageDisplayRemoveEventListener(token1)
imgdisp.ImageDisplayRemoveEventListener(token2)
			
//删除roi
number roino= imgdisp.ImageDisplayCountROIS()
for (number i=0; i<roino; i++)
{
number index
ROI currentROI = imgdisp.ImageDisplayGetROI( index )
imgdisp.imagedisplaydeleteROI(currentROI)
}

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

number t,l,b,r,xsize,ysize
img.GetSize(xsize,ysize)

number roino= imgdisp.ImageDisplayCountROIS()
if(roino==0)
{
roi roi1=createroi()
ROI1.ROISetrange(0.1*xsize,0.2*xsize)
roi1.roisetvolatile(0)
imgdisp.ImageDisplayAddROI( ROI1)
}

if(roino==1)
{
img.GetSelection(t,l,b,r)
img.ClearSelection()

roi roi1=createroi()
ROI1.ROISetrange(l,r)
roi1.roisetvolatile(0)
imgdisp.ImageDisplayAddROI( ROI1)
}

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
