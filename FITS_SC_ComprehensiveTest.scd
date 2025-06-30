// FITSAdjustedVisualisation_Test.scd
// Test script for visualizing adjusted FITS data using FITSVisualisation class.

// Define the path to your sample FITS file
~fitsFilePath = "/Users/marcinpietruszewski/Downloads/nuPG_2024_ScreenShot.fits";

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

	// sigmaClip returns a new FITSFile object, so we can use it directly
	clippedFits = fitsFile.sigmaClip(sigma: 2, maxIterations: 3);

	visualiser = clippedFits.visualiser; // Get visualiser for the clipped data

	"Opening interactive visualization of sigma clipped data...".postln;
	visualiser.visualizeInteractive(downsample: 2);

	"Opening histogram of sigma clipped data...".postln;
	visualiser.visualizeHistogram(numBins: 128);

// Replace your radial profile section (section 4) with this corrected version:
"\n--- 4. Comparing Radial Profiles: Original vs Filtered vs Clipped ---".postln;

if (fitsFile.validateData.not) {
    "Data validation failed. Cannot proceed with radial profile.".postln;
} {
    // ALL var declarations must come first
    var centerX = (fitsFile.width / 2).asInteger;
    var centerY = (fitsFile.height / 2).asInteger;
    var maxRadius = min(50, min(centerX, centerY));
    var profileDataArray = [];
    var labels = [];
    var colors = [];
    var filteredProfileData, clippedProfileData; // Declare these here

    "Calculating radial profiles for comparison...".postln;

    // Original image profile
    "Calculating original image profile...".postln;
    radialProfileData = fitsFile.radialProfile(centerX, centerY, maxRadius);
    if (radialProfileData[\radii].size > 0) {
        profileDataArray = profileDataArray.add(radialProfileData);
        labels = labels.add("Original Image");
        colors = colors.add(Color.green);
    };

    // Filtered image profile
    "Calculating filtered image profile...".postln;
    filteredProfileData = filteredFits.radialProfile(centerX, centerY, maxRadius);
    if (filteredProfileData[\radii].size > 0) {
        profileDataArray = profileDataArray.add(filteredProfileData);
        labels = labels.add("Gaussian Filtered");
        colors = colors.add(Color.blue);
    };

    // Clipped image profile
    "Calculating sigma-clipped image profile...".postln;
    clippedProfileData = clippedFits.radialProfile(centerX, centerY, maxRadius);
    if (clippedProfileData[\radii].size > 0) {
        profileDataArray = profileDataArray.add(clippedProfileData);
        labels = labels.add("Sigma Clipped");
        colors = colors.add(Color.red);
    };

    if (profileDataArray.size > 0) {
        "Calculated % radial profiles for comparison.".format(profileDataArray.size).postln;

        visualiser = fitsFile.visualiser;
        "Plotting comparison of radial profiles...".postln;

        try {
            visualiser.plotMultipleRadialProfiles(
                profileDataArray,
                labels,
                colors
            );
        } { |error|
            "Plotting failed: %".format(error).postln;
            "Displaying radial profile data instead:".postln;
            visualiser.printRadialProfileData(profileDataArray, labels);
        };
    } {
        "No valid radial profiles were calculated.".postln;
    };
};

	"\nNOTE: Several visualization windows should have opened. Close them manually when you are finished.".postln;
	"==================================================".postln;
) // End of test block
