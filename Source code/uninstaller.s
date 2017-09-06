// Script to create and install a package which contains the ePDF Tools scripts.
// H. L. Shi (honglongshi@outlook.com), Minzu University of China, Beijing, China
// version 1.0, Jan 2017.

// Note this script causes DM to shut down, so all open work should be saved before running this.
beep()
if(!continuecanceldialog("This uninstaller will shut down DigitalMicrograph.\n\n  Continue only if all open work has been saved."))
	{
		showalert("Unistallation cancelled by User!",2)
		exit(0)
	}

beep()
if(!continuecanceldialog("Remove the ePDF Tools Package and all its components from Digital Micrograph?"))
	{
		showalert("Uninstallation cancelled by User!",2)
		exit(0)
	}

// Delete the RDFToolPackage file from the plugins folder
string packagepath=getapplicationdirectory(1008,0)
packagepath=packagepath+"ePDF Tools.gtk"
number packagedeleted, settingsdeleted

if(doesfileexist(packagepath)) // delete it only if it exists
	{
		deletefile(packagepath)
		result("\nePDF Tools Package deleted from "+packagepath+"\n")
		packagedeleted=1
	}
else // It doesn't exist so show an alert
	{
		showalert("ePDF Tools Package not found in Plugins folder!",1)
		result("\nePDF Tools Package not found in Plugins folder!")
		packagedeleted=0
	}

// Delete the ePDF Tools settings from the Global Info
taggroup ptags=getpersistenttaggroup()
number tagexists=ptags.taggroupdoestagexist("ElectronDiffraction Tools:PDF")

if(tagexists) // the settings file exists
	{
		//deletepersistentnote("ElectronDiffraction Tools:PDF")
		result("\nPDF Settings file deleted\n")
		settingsdeleted=1
	}
else // the settings file does not exist
	{
		showalert("ePDF Tools Settings file not found!",1)
		result("\nePDF Tools Settings file not found!\n")
		settingsdeleted=0
	}

// Show alerts that reflect which files were found/deleted
// Neither the package nor the settings file was found/deleted
if(packagedeleted==0 && settingsdeleted==0) // Neither the package nor the settings file was found/deleted
	{
		showalert("No ePDF Tools files were found!",1)
		exit(0) // no need to shut down DM.
	}
		
// Package was found settings file was not	
if(packagedeleted==1 && settingsdeleted==0)
	{
		showalert("ePDF Tools Package was deleted, Settings file was not found!",2)
	}
	
// Package was not found, settings file was
if(packagedeleted==0 && settingsdeleted==1)
	{
		showalert("ePDF Tools Package was not found, Settings file was deleted!",2)
	}

showalert("Uninstallation of ePDF Tools is complete.\n\n DigitalMicrograph will now shut down.",2)

// Close DM
choosemenuitem("File","","Exit")