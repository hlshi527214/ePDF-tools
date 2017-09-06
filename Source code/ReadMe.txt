ePDF tools, a processing and analysis package of the atomic pair distribution function for electron diffraction
-------------------------------------------


The folder 'source code' contains 11 source codes in *.s, and two txt files (Scattering factor.txt, and ReadMe.txt). 
-------------------------------------------
BG-Removed.s is a source code to perform background subtraction based on the cubic spline fitting of the sampling points. 

BG-Removed-Gauss smooth.s is a source code to perform background subtraction based on the algorithm of the Gaussian smooth. 

Data.s is a source code to perform quantitative PDF analysis in real time. 

ePDF calculator 2017.s is a graphical user interface of ePDF tools. 

G(r)-Calibration.s is a source code to perform PDF calibration. 

I(Q)-Smooth.s is a source code to filter the sparks of the profile. 

Installer.s is a source code to install ePDF tools into DigitalMicrograph. 

Iq-Fq.s is a source code to get the reduced structure function by data normalization and other corrections. 

Profile-RIF.s is a source code to add sampling points for background subtraction.

Scattering Factor.s is a source code to import the scattering factors into DigitalMicrograph. 

uninstaller.s is a source code to uninstall ePDF tools from DigitalMicrograph.

ReadMe.txt is a text file to briefly describe ePDF tools. 

Scattering factor.txt is a text file containing scattering factors of elements from 1 to 99 for electrons.
-------------------------------------------


How to install ePDF tools in DigitalMicrograph.
-------------------------------------------
keep all of 13 source codes in one folder, e.g. Source code!

1) Launch DigitalMicrograph (If your computer has no Gatan Microscopy Suite, you should install it firstly: http://www.gatan.com/products/tem-analysis/gatan-microscopy-suite-software. Software licenses are required to install this software.)

2) Drag Installer.s into DigitalMicrograph and run it by holding down the CONTROL key and hitting 'ENTER'.

3) A window is prompted to select the folder 'Source code', and choose Save button. a new menu labeled 'ePDF Tools' will be added on the menu bar.

4) In order to import the scattering data, Click the item 'ePDF Tools/ Scattering factors' and then to select the file 'Scattering factor.txt'.
-------------------------------------------


How to uninstall ePDF Tools.
-------------------------------------------
1) Launch DigitalMicrograph.

2) Drag Unstaller.s into DigitalMicrograph and run it by holding down the CONTROL key and hitting 'ENTER' to uninstall it.
-------------------------------------------

How to use ePDF Tools (See the demo video for details). 
-------------------------------------------
Launch DigitalMicrograph and open ePDF Tools (ePDF Tools/ePDF calculator 2017), a GUI named 'ePDF Tools' will be displayed. 

1) Open a SAED data, e.g., Example pattern.dm3

2) Set parameters: Input 'Au' and '1' in the fields of 'Symb.' and 'At. %' for polycrystal gold, respectively. 
							   Input 200 in 'kV' for the acceleration voltage,
							   Input 19.3 in 'Dens.' for the material density (g/cm^3),
							   Input 20 in 'r_max' for the max dimension of PDFs (Å)
							   Input 2000 in 'Size' for the Size of PDFs (pixels)

3) Profile preparation: keep SAED pattern as the foremost pattern, and then 
									a) Find the center. Click 'Center' button, move the cursor to the desired position and then press the space bar to find the Center of the SAED pattern.
									b) Get an intensity Profile by clicking 'Profile' button.
									c) Calibrate the profile(Optional).  Carefully drag an ROI (region of interest) across the desired peak and then hit the ‘Calib.’ button in the image processing box to calibrate it, Input 2.355 In the prompted window. 

4) Calculate and add the scattering curves by subsequently clicking 'Calc.' and 'Scatt.' in 'PDF tools' box.

5) Background subtraction: keep the intensity Profile as the foremost image.
											a) Drag an ROI to define the first sampling point.
											b) Click 'BG' button to entry the mode of Background subtraction. Add, change, move or delete ROI to define the sampling point.
											c) Delete 'BG' ROI for next step when the Background curve looks good.

6) Data normalization: keep the Background subtracted intensity Profile as the foremost image.	
									a) Click 'F(Q)' button to entry the normalization mode.
									b) Add, change or move ROI labeled 'N', 'c', 'Damp' and 'Smooth' to set the parameters of N, c, Gaussian damping D, and the strength of smoothing.

7) PDF calibration: keep the intensity Profile containing 'Damped' slice as the foremost image.
							   a) Click G(r) button to get the reduced G(r).
							   b) Drag an ROI spanning the first physical peak and then click 'Calib.' in 'PDF tools' box to calibrate it.
							   c) Delete the ROI to get the calibrated PDFs.

8) PDF analysis: keep one of calibrated PDF as the foremost image.
							a) Drag an ROI to define the limits of the integral.
							b) Click 'Data' button to perform quantitative PDF analysis in real time.





