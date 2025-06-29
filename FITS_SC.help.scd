SuperCollider Help File for FITS_SC Classes

// This help file describes the FITSFile and FITSVisualisation classes,
// which provide tools for reading, processing, and visualizing FITS
// (Flexible Image Transport System) astronomical data within SuperCollider.

// Installation:
// 1. Clone the FITS_SC repository:
//    git clone https://github.com/marcinpiet/FITS_SC.git /path/to/your/FITS_SC_Project
// 2. Create a symbolic link from the cloned directory to your SuperCollider Extensions folder:
//    ln -s /path/to/your/FITS_SC_Project /Users/yourusername/Library/Application\ Support/SuperCollider/Extensions/FITS_SC
//    (Replace /path/to/your/FITS_SC_Project with the actual path where you cloned this repository.)
// 3. Recompile Class Library in SuperCollider: Language -> Recompile Class Library.

// Usage:
// To use the classes, you will need a FITS file. For testing, you can place a small
// FITS file (e.g., a 100x100 pixel image) named `sample.fits` in the root of your
// `FITS_SC_Project` directory.

// Example FITS file download:
// You can find sample FITS files on the NASA FITS website:
// https://fits.gsfc.nasa.gov/fits_samples.html
// Look for small files like the "HST WFPC II (100 x 100 pixels)" example.

// --------------------------------------------------------------------------------
// Class: FITSFile
// Description:
// A class for loading, accessing, and processing data from FITS files.
// It provides methods for reading header information, pixel data, statistical
// analysis, image processing, source finding, WCS (World Coordinate System)
// conversions, and data export.

// Constructor:
// FITSFile.new(filePath)
//   filePath: Absolute path to the FITS file.

// Instance Methods:

// .printSummary
//   Prints a summary of the FITS file, including dimensions and header information.

// .width
//   Returns the width of the primary image data.

// .height
//   Returns the height of the primary image data.

// .data
//   Returns the raw image data as a FloatArray.

// .getPixel(x, y)
//   Returns the pixel value at the specified (x, y) coordinates.

// .getRow(y)
//   Returns the entire row of pixel data at the specified y-coordinate as a FloatArray.

// .calculateStatistics
//   Calculates and returns a Dictionary of basic statistics (mean, median, std, min, max)
//   for the image data.

// .sigmaClip(sigma: 3, maxIterations: 5)
//   Performs sigma clipping on the data to remove outliers.
//   sigma: Number of standard deviations to clip beyond.
//   maxIterations: Maximum number of clipping iterations.
//   Returns a new FloatArray with clipped data.

// .gaussianFilter(sigma: 1.0)
//   Applies a Gaussian blur filter to the image data.
//   sigma: Standard deviation of the Gaussian kernel.
//   Returns a new FITSFile instance with the filtered data.

// .findSources(threshold: 100, minSeparation: 10)
//   Identifies sources (e.g., stars, galaxies) in the image data.
//   threshold: Minimum pixel value to consider as part of a source.
//   minSeparation: Minimum pixel distance between detected sources.
//   Returns an Array of Dictionaries, each describing a source (e.g., x, y, peak value).

// .aperturePhotometry(centerX, centerY, radius, innerRadius: 0, outerRadius: 0)
//   Performs aperture photometry around a given center.
//   centerX, centerY: Coordinates of the aperture center.
//   radius: Radius of the photometric aperture.
//   innerRadius, outerRadius: Optional. For annulus photometry (sky background estimation).
//   Returns a Dictionary with net flux and sky level.

// .pixelToWorld(pixelX, pixelY)
//   Converts pixel coordinates to world coordinates (e.g., RA, Dec) using WCS information.
//   Returns an Array of world coordinates.

// .worldToPixel(worldX, worldY)
//   Converts world coordinates to pixel coordinates.
//   Returns an Array of pixel coordinates.

// .radialProfile(centerX, centerY, maxRadius, numBins: 50)
//   Calculates the radial profile of the image data around a given center.
//   centerX, centerY: Center of the radial profile.
//   maxRadius: Maximum radius to consider.
//   numBins: Number of bins for the profile.
//   Returns a FloatArray representing the radial profile.

// .getHilbertCurveDataArbitrary
//   Transforms the 2D image data into a 1D array following a Hilbert curve path.
//   Returns a FloatArray.

// .exportToCSV(filePath)
//   Exports the image data to a CSV file.

// .exportRegion(x, y, width, height, filePath)
//   Exports a specified rectangular region of the image data to a CSV file.

// .visualiser
//   Returns an instance of FITSVisualisation associated with this FITSFile.

// --------------------------------------------------------------------------------
// Class: FITSVisualisation
// Description:
// A class for visualizing FITS image data within SuperCollider, providing
// interactive plotting and data representation tools.

// Constructor:
// FITSVisualisation.new(fitsFile)
//   fitsFile: An instance of FITSFile. (Typically accessed via `fitsFile.visualiser`)

// Instance Methods:

// .visualizeInteractive(downsample: 1)
//   Opens an interactive window to visualize the FITS image data.
//   Allows for zooming, panning, and displaying pixel values on hover.
//   downsample: Factor to downsample the image for faster rendering (e.g., 2 for half resolution).

// .visualizeHistogram(numBins: 256)
//   Displays a histogram of the FITS image data.
//   numBins: Number of bins for the histogram.

// .plotMultipleRadialProfiles(profiles, labels, colors)
//   Plots multiple radial profiles on a single graph for comparison.
//   profiles: An Array of FloatArrays, where each FloatArray is a radial profile.
//   labels: An Array of Strings, labels for each profile.
//   colors: An Array of Colors, colors for each profile.

// --------------------------------------------------------------------------------
// Examples:

s.boot; // Ensure the server is booted

// Define the path to your sample FITS file
// Make sure you have a small FITS file named `sample.fits` in your FITS_SC_Project directory.
~fitsFilePath = thisProcess.platform.userAppSupportDir +/+ "SuperCollider" +/+ "Extensions" +/+ "FITS_SC" +/+ "sample.fits";

// --- Example 1: Basic FITSFile Usage ---
(
"--- Example 1: Basic FITSFile Usage ---".postln;
var fitsFile = FITSFile.new(~fitsFilePath);
fitsFile.printSummary;

// Get a pixel value
var pixelValue = fitsFile.getPixel(fitsFile.width / 2, fitsFile.height / 2);
"Pixel at center: %".format(pixelValue).postln;

// Calculate statistics
var stats = fitsFile.calculateStatistics;
"Statistics: Mean=%, Median=%, Std=%".format(stats[\mean], stats[\median], stats[\std]).postln;
)

// --- Example 2: Image Processing and Source Finding ---
(
"--- Example 2: Image Processing and Source Finding ---".postln;
var fitsFile = FITSFile.new(~fitsFilePath);

// Apply a Gaussian filter
var filteredFits = fitsFile.gaussianFilter(sigma: 1.5);
"Filtered image range: % to %".format(filteredFits.data.minItem, filteredFits.data.maxItem).postln;

// Find sources in the filtered image
var sources = filteredFits.findSources(threshold: filteredFits.data.mean + (filteredFits.data.std * 2), minSeparation: 10);
"Found % sources.".format(sources.size).postln;

if(sources.size > 0) {
    var firstSource = sources[0];
    var photometry = filteredFits.aperturePhotometry(firstSource[\x], firstSource[\y], 5);
    "Photometry for first source: Net Flux=%, Sky=%".format(photometry[\netFlux].round(0.1), photometry[\skyLevel].round(0.1)).postln;
};
)

// --- Example 3: WCS Conversion ---
(
"--- Example 3: WCS Conversion ---".postln;
var fitsFile = FITSFile.new(~fitsFilePath);
var pixelCoord = [10, 10];
var worldCoord = fitsFile.pixelToWorld(pixelCoord[0], pixelCoord[1]);
"Pixel %, % converts to World Coordinates: %".format(pixelCoord[0], pixelCoord[1], worldCoord).postln;

var convertedPixel = fitsFile.worldToPixel(worldCoord[0], worldCoord[1]);
"World Coordinates % converts back to Pixel: %".format(worldCoord, convertedPixel).postln;
)

// --- Example 4: Data Transformation and Export ---
(
"--- Example 4: Data Transformation and Export ---".postln;
var fitsFile = FITSFile.new(~fitsFilePath);

// Radial Profile
var profile = fitsFile.radialProfile(fitsFile.width/2, fitsFile.height/2, 50);
"Radial profile generated with % bins.".format(profile.size).postln;

// Hilbert Curve Data
var hilbertData = fitsFile.getHilbertCurveDataArbitrary;
"Hilbert curve data generated with size %".format(hilbertData.size).postln;

// Export to CSV
var csvPath = Platform.userAppSupportDir +/+ "SuperCollider" +/+ "Extensions" +/+ "FITS_SC" +/+ "exported_fits_data.csv";
fitsFile.exportToCSV(csvPath);
"Data exported to: %".format(csvPath).postln;
)

// --- Example 5: FITSVisualisation - Interactive and Histogram ---
(
"--- Example 5: FITSVisualisation - Interactive and Histogram ---".postln;
var fitsFile = FITSFile.new(~fitsFilePath);
var visualiser = fitsFile.visualiser;

// Interactive visualization (a new window will open)
"Opening interactive visualization...".postln;
visualiser.visualizeInteractive(downsample: 2); // Downsample for larger files

// Histogram visualization (another new window will open)
"Opening histogram visualization...".postln;
visualiser.visualizeHistogram(numBins: 128);

"Close the visualization windows manually when done.".postln;
)

// --- Example 6: FITSVisualisation - Plotting Multiple Radial Profiles ---
(
"--- Example 6: FITSVisualisation - Plotting Multiple Radial Profiles ---".postln;
var fitsFile = FITSFile.new(~fitsFilePath);
var visualiser = fitsFile.visualiser;

// Generate two different radial profiles
var radialProf1 = fitsFile.radialProfile(fitsFile.width / 2, fitsFile.height / 2, 50);
var filteredFits = fitsFile.gaussianFilter(sigma: 0.8);
var radialProf2 = filteredFits.radialProfile(filteredFits.width / 2, filteredFits.height / 2, 50);

// Plot them on the same graph
"Plotting multiple radial profiles...".postln;
visualiser.plotMultipleRadialProfiles(
    [radialProf1, radialProf2],
    ["Original Image Profile", "Filtered Image Profile"],
    [Color.red, Color.blue]
);
"A window showing both radial profiles should appear.".postln;
)