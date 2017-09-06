// Script to create and install a package which contains the ePDF Tools scripts.
// H. L. Shi (honglongshi@outlook.com), Minzu University of China, Beijing, China
// version 1.0, Jan 2017.

// To install ePDF Tools:
// 1) Open this installer script in DigitalMicrograph then hit CONTROL + ENTER. Follow the prompt and open the folder called 'ePDF Tools'. 2) Click 'SAVE' and the installer will create and install the ePDF Tools package in a new menu called 'ePDF Tools'. It is then ready to use.
number linecountglobal
number level=3 // specifies that the package will be gtk - highest priority

// Function which adds a script to the package
number AddS2P(string PathName, string packagename, string CommandName, string Menu, string SubMenu)
	{
	// This function installs a script as a menu command
	string filename=pathextractbasename(pathname,0)
	result(" Installing <"+FileName+">\n\t")
	try 
		{
		AddScriptFileToPackage(PathName,PackageName,level,CommandName,Menu,SubMenu,0)
		}
	catch
		{
		beep()
		result(" => Error while installing as a menu command ("+Menu+"-"+SubMenu+"-"+CommandName+").\n")	
		return -1
		}
	result(" => Installed as a menu command ("+Menu+"-"+SubMenu+"-"+CommandName+").\n")	
	return 0
	}

// Adds a spacer to the package

number AddLine2P(string packagename, string Menu, string SubMenu)
	{
	// This function intalls a seperation line in a menu
	// Note, that it is necessary to install each line under a different 
	// "position" with regard to menu/submenu/command-name. That is the 
	// reason, why the line-count variable is needed.

	LineCountGlobal++
	result(" Installing a line\n\t")
	AddScriptToPackage("",PackageName,level,"-"+menu+"."+LineCountGlobal, menu, submenu, 0)
	result(" => Installed.\n")
	return 0
	}

// Reports progress with the package creation
void InstallationLog( String Package, String ScriptFile , number thiserror)
	{
		if(thiserror==0) // the installation went OK
			{
				Result( "\n" + "Package [" + package + "], script file [" + scriptfile + "] installed." )
				return
			}
		else // there was a problem
			{
				Result( "\n" + "Package [" + package + "], script file [" + scriptfile + "] could not be installed." )
				return
			}
	}

if(!continuecanceldialog("Install ePDF Tools package in DigitalMicrograph?")) exit(0)

//----------------------- Main script starts here---------------------
	String scriptFile, pkgName, cmd, menu, submenu
	pkgName	= "ePDF Tools"
	cmd="ePDF Tools"

//	Get the directory of scripts to be installed
	String ParentPath, Path
	If( !SaveAsDialog("", "Open ePDF Tools Scripts directory then SAVE", ParentPath) )
		{
			showalert("Installation cancelled by user.",2)
			exit(0)
		}

	ParentPath = ParentPath.PathExtractDirectory(0)
	//documentwindow reswin=getresultswindow(1)

	number erroraccumulator=0
	number thiserror=0

// Install the various ePDF Tools scripts as menu items
	menu 		= "ePDF Tools"
	submenu 	= ""
	
	scriptFile 	= "I(Q)-Smooth.s"	
	cmd 		= "I(Q)-Smooth"
	thiserror=AddS2P(parentpath+scriptfile, pkgname, cmd, Menu, SubMenu)
	erroraccumulator=erroraccumulator+thiserror
	InstallationLog( pkgname, ScriptFile , thiserror)	
	
	scriptFile 	= "BG-Removed.s"	
	cmd 		= "BG-Removed"
	thiserror=AddS2P(parentpath+scriptfile, pkgname, cmd, Menu, SubMenu)
	erroraccumulator=erroraccumulator+thiserror
	InstallationLog( pkgname, ScriptFile , thiserror)
	
	scriptFile 	= "BG-Removed-Gauss smooth.s"	
	cmd 		= "BG-Removed-Gauss"
	thiserror=AddS2P(parentpath+scriptfile, pkgname, cmd, Menu, SubMenu)
	erroraccumulator=erroraccumulator+thiserror
	InstallationLog( pkgname, ScriptFile , thiserror)	

	scriptFile 	= "Iq-Fq.s"	
	cmd 		= "Iq-Fq"
	thiserror=AddS2P(parentpath+scriptfile, pkgname, cmd, Menu, SubMenu)
	erroraccumulator=erroraccumulator+thiserror
	InstallationLog( pkgname, ScriptFile , thiserror)

	scriptFile 	= "G(r)-Calibration.s"	
	cmd 		= "G(r)-Calibration"
	thiserror=AddS2P(parentpath+scriptfile, pkgname, cmd, Menu, SubMenu)
	erroraccumulator=erroraccumulator+thiserror
	InstallationLog( pkgname, ScriptFile , thiserror)

	// add a separate line
	AddLine2P(pkgname, Menu, SubMenu)

	scriptFile 	= "Data.s"	
	cmd 		= "Data"
	thiserror=AddS2P(parentpath+scriptfile, pkgname, cmd, Menu, SubMenu)
	erroraccumulator=erroraccumulator+thiserror
	InstallationLog( pkgname, ScriptFile , thiserror)
	
	scriptFile 	= "Profile-RIF.s"	
	cmd 		= "Profile-RIF"
	thiserror=AddS2P(parentpath+scriptfile, pkgname, cmd, Menu, SubMenu)
	erroraccumulator=erroraccumulator+thiserror
	InstallationLog( pkgname, ScriptFile , thiserror)

	scriptFile 	= "Scattering Factor.s"	
	cmd 		= "Scattering Factor"
	thiserror=AddS2P(parentpath+scriptfile, pkgname, cmd, Menu, SubMenu)
	erroraccumulator=erroraccumulator+thiserror
	InstallationLog( pkgname, ScriptFile , thiserror)
	
// add a separate line
	AddLine2P(pkgname, Menu, SubMenu)

// and then add ePDF Tools menu	
	scriptFile 	= "ePDF calculator 2017.s"	
	cmd 		= "ePDF calculator 2017"
	thiserror=AddS2P(parentpath+scriptfile, pkgname, cmd, Menu, SubMenu)
	erroraccumulator=erroraccumulator+thiserror
	InstallationLog( pkgname, ScriptFile , thiserror)	
	
// Do a final error reporting
	if(erroraccumulator==0)
		{
			showalert("Installation of ePDF Tools completed successfully!\n\nIt is recommended that you  restart DigitalMicrograph.",2)
		}
	else
		{
			showalert("There were "+erroraccumulator+" errors creating the '"+pkgname+"' package. See Results window for details.",2)
		}
			