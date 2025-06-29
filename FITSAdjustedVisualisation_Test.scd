// FITSAdjustedVisualisation_Test.scd
// Test script for visualizing adjusted FITS data using FITSVisualisation class.

s.boot; // Ensure the server is booted

// Define the path to your sample FITS file
// Make sure you have a small FITS file named `sample.fits` in your FITS_SC_Project directory.
~fitsFilePath = thisProcess.platform.userAppSupportDir +/+ "SuperCollider" +/+ "Extensions" +/+ "FITS_SC" +/+ "sample.fits";

// Load the FITSFile and FITSVisualisation classes
// (These are loaded automatically if SuperCollider's class library is recompiled
// and the symbolic link is set up correctly, but explicit loading is good for testing.)
// this.executeFile(thisProcess.platform.userAppSupportDir +/+ "SuperCollider" +/+ "Extensions" +/+ "FITS_SC" +/+ "FITSFile.sc");
// this.executeFile(thisProcess.platform.userAppSupportDir +/+ "SuperCollider" +/+ "Extensions" +/+ "FITS_SC" +/+ "FITSVisualisation.sc");

( // Start of test block
	var fitsFile, visualiser, filteredFits, clippedFits, radialProfileData;

	"==================================================".postln;
	" FITS ADJUSTED VISUALISATION TEST SCRIPT ".postln;
	"==================================================".postln;

	// --- 1. Load the FITS file ---
	"\n--- 1. Loading FITS File ---".postln;
	fitsFile = FITSFile.new(~fitsFilePath);
	fitsFile.printSummary;

	// --- 2. Visualizing Gaussian Filtered Data ---
	"\n--- 2. Visualizing Gaussian Filtered Data ---".postln;
	"Applying Gaussian filter (sigma: 1.5)...".postln;
	filteredFits = fitsFile.gaussianFilter(sigma: 1.5);

	visualiser = filteredFits.visualiser; // Get visualiser for the filtered data

	"Opening interactive visualization of filtered data...".postln;
	visualiser.visualizeInteractive(downsample: 2); // Downsample for faster rendering

	"Opening histogram of filtered data...".postln;
	visualiser.visualizeHistogram(numBins: 128);

	// --- 3. Visualizing Sigma Clipped Data ---
	"\n--- 3. Visualizing Sigma Clipped Data ---".postln;
	"Applying Sigma Clipping (sigma: 2, maxIterations: 3)...".postln;
	// Note: sigmaClip returns a FloatArray, not a FITSFile instance.
	// To visualize it with FITSVisualisation, we'd need to create a new FITSFile
	// or a custom visualization method that accepts a raw FloatArray and dimensions.
	// For simplicity, we'll demonstrate by creating a temporary FITSFile instance.
	clippedFits = FITSFile.new(fitsFile.filePath); // Re-load original to get dimensions
	clippedFits.data = fitsFile.sigmaClip(sigma: 2, maxIterations: 3); // Replace data with clipped data

	visualiser = clippedFits.visualiser; // Get visualiser for the clipped data

	"Opening interactive visualization of sigma clipped data...".postln;
	visualiser.visualizeInteractive(downsample: 2);

	"Opening histogram of sigma clipped data...".postln;
	visualiser.visualizeHistogram(numBins: 128);

	// --- 4. Plotting Radial Profile ---
	"\n--- 4. Plotting Radial Profile ---".postln;
	"Calculating radial profile...".postln;
	radialProfileData = fitsFile.radialProfile(fitsFile.width / 2, fitsFile.height / 2, 50);

	visualiser = fitsFile.visualiser; // Use visualiser from original fitsFile
	"Plotting radial profile...".postln;
	visualiser.plotMultipleRadialProfiles(
		[radialProfileData],
		["Original Image Radial Profile"],
		[Color.green]
	);

	"\nNOTE: Several visualization windows should have opened. Close them manually when you are finished.".postln;
	"==================================================".postln;
) // End of test block
