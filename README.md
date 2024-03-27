Welcome to the POLAR.m MATLAB script used to process data taken in the UAH Mueller Matrix Spectropolarimeter. This is the code used in the work described in the paper E.E. James, I.M. Anderson, and D.A. Gregory,
"Spectropolarimetric characterization of vortex optics" (submitted for peer review but not published yet--plan to leave an Optica Open link here).

Short version: if you download POLAR.m, POLARfuncs, and Data into the same level folder and run the code as is, it will calculate and plot the Mueller matrix and associated polarization properties from the sample
data file included under the Data folder. The sample data is a measurement of WPV10L-532 taken on February 7, 2023. Before you do that though, here is how the POLAR.m code and associated files are structured:

Under Data/SamplesTRANSMISSION/VortexRetarderWPV10L532/WPV10L532_07Feb2023, you will see a file called PARAMETERS.txt. This contains important information about the specific measurement taken on February 7 used
by the MATLAB scripts and by the user. You will also see a raw data file with a .lvm extension. The first column of this data lists the wavelengths, and the rest of the columns are irradiance measurements for
the 61 different retarder positions discussed in the paper. Finally, you will see a .mat file called MMideal--this is an "ideal" Mueller matrix with which to compare the measured one.

In lines 80-81 of POLAR.m, you will see FL.Calibration.FLDR and .NAME. This is the location and name of the calibration data required for the data reduction process. The calibration FLDR shoudl be called 
'Data\DiagnosticRunsTRANSMISSION\24Jan2023', and NAME should be 'CAL1_24Jan2023.' Do not change it, or else the Mueller matrix cannot be calculated. Fortunately, if you accidentally delete or change it, the
information is available in PARAMETERS.txt.

In lines 83-84, there are also FL.Sample.FLDR and FL.Sample.txt. These provide the location and name of the raw data file. 

After this, you will see a few toggles. TOG.CAL=1 will actually compute the calibration parameters for the file listed under FL.Calibration.FLDR and FL.Calibration.NAME. You can do this if you want, but leave it as 
0 (or any other number) for the vortex retarder. TOG.COMP=1 compares the vortex retarder Mueller matrix to the "ideal Mueller matrix"--however, MMideal is currently set as a half-wave retarder with horizontal 
fast axis orientation, so it will not be very useful for vortex retarders. TOG.PLOT should always be set to 1 if you want to see plots of the Mueller matrix and polarization properties (among other parameters).

The rest of the scripts are well-commented and contain further details. Hopefully this was enough to get started. If there are any questions, reach out to me at eej0004@uah.edu.






