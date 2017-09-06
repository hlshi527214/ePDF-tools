/*-----------------------------------读取Scattering factor 文件到PDF缓存中-----------------------------*/
//Scattering factor 文件: element 1	a1	a2	a3	b1	b2	b3	c1	c2	c3	d1	d2	d3......

number i, fileref,xorigin,width, tablast,xsize,xstart,xend,no=0
string  pathname,path, extension, filename, thisline,Delimiter=" "


//选择文件，得到路径和路径名
documentwindow win
if(!opendialog(win, "Scattering factor file",path, pathname))exit(0)

// 只读 (.txt)
extension=pathextractextension(pathname,0)
filename=pathextractfilename(pathname,0)

	if(extension!="txt") 
		{
			if(twobuttondialog("This is not a text file","Cancel","Proceed"))
			{
			showalert("Cancelled by User.",2)
			exit(0)
			}
		}
		
// 开始读取文件
try
{
fileref=openfileforreading(pathname)
object file_stream= NewStreamFromFileReference(fileref, 1) 
number readresult=1
number counter=0
number yval,xval

// 总行数
while(readresult==1)
{
readresult=ReadFileLine(fileref,thisline)
width=len(thisline)	

//在行内识别读取分隔符，为方便排序，除以100
for(i=0;i<width;i++)
{
string thischar=mid(thisline,i,1)

number Bytes
if(asc(thischar)==9)  //找分隔符tab
{
Bytes=len(left(thisline,i))
setpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:txt file: line "+counter+":Tab "+no,Bytes)
//result(bytes*100+",")
no=no+1
}
}

no=0
//result("\n")

counter=counter+1
}
closefile(fileref)

//逐行读取
fileref=openfileforreading(pathname)
number step=0,tab0, tab1,tab2,tab3,tab4,tab5,tab6,tab7,tab8,tab9,tab10,tab11
while(step<counter)  
{
readresult=ReadFileLine(fileref,thisline)  
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:txt file: line "+step+":Tab 0",tab0)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:txt file: line "+step+":Tab 1",tab1)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:txt file: line "+step+":Tab 2",tab2)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:txt file: line "+step+":Tab 3",tab3)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:txt file: line "+step+":Tab 4",tab4)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:txt file: line "+step+":Tab 5",tab5)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:txt file: line "+step+":Tab 6",tab6)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:txt file: line "+step+":Tab 7",tab7)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:txt file: line "+step+":Tab 8",tab8)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:txt file: line "+step+":Tab 9",tab9)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:txt file: line "+step+":Tab 10",tab10)
getpersistentnumbernote("ElectronDiffraction Tools:Crystal Tools:Read Files:txt file: line "+step+":Tab 11",tab11)

string element
number a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3
//行内读取
number width=len(thisline)		
for( i=0;i<width;i++)
{	
// 第0个tab左侧为元素，tab1和tab0之间为a1
if(i==tab0)
{
element=left(thisline,i)

a1=val(mid(thisline,i,tab1-tab0))
SetPersistentNumberNote("ElectronDiffraction Tools:PDF:Scattering Factors:"+element+":a1",a1)
}

//a2: between tab1-tab2
if(i==tab1)
{
 a2=val(mid(thisline,i,tab2-tab1))
SetPersistentNumberNote("ElectronDiffraction Tools:PDF:Scattering Factors:"+element+":a2",a2)
}

//a3: between tab3-tab2
if(i==tab2)
{
 a3=val(mid(thisline,i,tab3-tab2))
SetPersistentNumberNote("ElectronDiffraction Tools:PDF:Scattering Factors:"+element+":a3",a3)
}

//b1: between tab1-tab2
if(i==tab3)
{
b1=val(mid(thisline,i,tab4-tab3))
SetPersistentNumberNote("ElectronDiffraction Tools:PDF:Scattering Factors:"+element+":b1",b1)
}

//b2: between tab1-tab2
if(i==tab4)
{
b2=val(mid(thisline,i,tab5-tab4))
SetPersistentNumberNote("ElectronDiffraction Tools:PDF:Scattering Factors:"+element+":b2",b2)
}

//b3: between tab1-tab2
if(i==tab5)
{
b3=val(mid(thisline,i,tab6-tab5))
SetPersistentNumberNote("ElectronDiffraction Tools:PDF:Scattering Factors:"+element+":b3",b3)
}

//b1: between tab1-tab2
if(i==tab6)
{
c1=val(mid(thisline,i,tab7-tab6))
SetPersistentNumberNote("ElectronDiffraction Tools:PDF:Scattering Factors:"+element+":c1",c1)
}

//b2: between tab1-tab2
if(i==tab7)
{
c2=val(mid(thisline,i,tab8-tab7))
SetPersistentNumberNote("ElectronDiffraction Tools:PDF:Scattering Factors:"+element+":c2",c2)
}

//b3: between tab1-tab2
if(i==tab8)
{
c3=val(mid(thisline,i,tab9-tab8))
SetPersistentNumberNote("ElectronDiffraction Tools:PDF:Scattering Factors:"+element+":c3",c3)
}

//b1: between tab1-tab2
if(i==tab9)
{
d1=val(mid(thisline,i,tab10-tab9))
SetPersistentNumberNote("ElectronDiffraction Tools:PDF:Scattering Factors:"+element+":d1",d1)
}

//b2: between tab1-tab2
if(i==tab10)
{
d2=val(mid(thisline,i,tab11-tab10))
SetPersistentNumberNote("ElectronDiffraction Tools:PDF:Scattering Factors:"+element+":d2",d2)
}

//b3: between tab1-tab2
if(i==tab11)
{
d3=val(mid(thisline,i,width-tab11))
SetPersistentNumberNote("ElectronDiffraction Tools:PDF:Scattering Factors:"+element+":d3",d3)
}

}
step=step+1
}
}

catch // If an error occurs close the file.
	{
		closefile(fileref)
		showalert("Sorry, an error occurred",2)
	}
		