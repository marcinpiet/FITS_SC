FITSFile {
    // Instance Variables
    var <file,       // The File object for the FITS file
        <headers,    // A Dictionary to store the FITS headers
        <data,       // The raw pixel data as a flat Array
        <width,      // The width of the image in pixels
        <height,     // The height of the image in pixels
        <bitpix,     // The BITPIX value from the header (bits per pixel)
        <datatype;   // The SuperCollider data type (e.g., \float32)

    // WCS (World Coordinate System) variables
    var <crpix1, <crpix2, <crval1, <crval2, <cdelt1, <cdelt2, <ctype1, <ctype2, <cunit1, <cunit2;
    var <pc1_1, <pc1_2, <pc2_1, <pc2_2; // PC matrix for rotation

    // Factory method to create a new FITSFile object from a file path.
    *new { |path|
        ^super.new.init(path);
    }

    // Factory method to create a FITSFile object from an existing data array.
    // This is useful for testing or creating FITS objects from generated data.
    *newFromData { |dataArray, w, h, headersDict|
        ^super.new.initFromData(dataArray, w, h, headersDict);
    }

    // Initializes the FITSFile by reading and parsing the file at the given path.
    init { |path|
        this.readFITS(path);
    }

    // Initializes the FITSFile from an existing data array and metadata.
    initFromData { |dataArray, w, h, headersDict|
        data = dataArray;
        width = w;
        height = h;
        headers = headersDict ? Dictionary.new;
    }

    // Reads and parses the entire FITS file, populating headers and data.
    readFITS { |path|
        var headerSize = 2880; // Standard FITS header block size
        var headerBlocks, dataSize, rawData, bytesRead;
        var bzero = 0, bscale = 1; // Default scaling values
        var endFound = false;
        var headerBlock, headerString, numCards;
        var startIdx, endIdx, card, key, value, parts, num;

        file = File(path, "rb"); // Open the file in binary read mode
        if(file.isNil) {
            Error("Could not open FITS file: %".format(path)).throw;
        };

        // --- Read and parse headers ---
        headers = Dictionary.new;
        headerBlocks = 0;

        while { endFound.not } {
            headerBlock = Int8Array.newClear(headerSize);
            file.read(headerBlock); // Read one header block
            headerString = String.newFrom(headerBlock.collect { |b|
                if(b < 0) { (b + 256).asAscii } { b.asAscii }
            });

            headerBlocks = headerBlocks + 1;
            numCards = (headerSize / 80).asInteger; // Each card is 80 characters

            // Parse each 80-character card
            numCards.do { |i|
                startIdx = i * 80;
                endIdx = ((i + 1) * 80) - 1;
                card = headerString[startIdx..endIdx];

                if(card.beginsWith("END")) {
                    endFound = true;
                };

                if(card.contains("=")) {
                    parts = card.split($=);
                    if(parts.size >= 2) {
                        key = parts[0].stripWhiteSpace;
                        value = parts[1].stripWhiteSpace;

                        // Remove comments
                        if(value.contains("/")) {
                            value = value.split($/)[0].stripWhiteSpace;
                        };

                        // Handle string vs. numeric values
                        if(value.beginsWith("'") and: value.endsWith("'")) {
                            value = value[1..(value.size-2)];
                        } {
                            num = value.asFloat;
                            if(num.notNil and: { value != "" }) { value = num };
                        };

                        headers[key] = value;
                    };
                };
            };
        };

        // --- Extract key parameters ---
        width = headers["NAXIS1"].asInteger;
        height = headers["NAXIS2"].asInteger;
        bitpix = headers["BITPIX"].asInteger;
        bzero = headers["BZERO"] ? 0;
        bscale = headers["BSCALE"] ? 1;

        // Extract WCS parameters (if they exist)
        crpix1 = headers["CRPIX1"];
        crpix2 = headers["CRPIX2"];
        crval1 = headers["CRVAL1"];
        crval2 = headers["CRVAL2"];
        cdelt1 = headers["CDELT1"];
        cdelt2 = headers["CDEL2"];
        ctype1 = headers["CTYPE1"];
        ctype2 = headers["CTYPE2"];
        cunit1 = headers["CUNIT1"];
        cunit2 = headers["CUNIT2"];
        pc1_1 = headers["PC1_1"];
        pc1_2 = headers["PC1_2"];
        pc2_1 = headers["PC2_1"];
        pc2_2 = headers["PC2_2"];

        "FITS Headers parsed:".postln;
        "Width: %, Height: %, BITPIX: %".format(width, height, bitpix).postln;
        "BZERO: %, BSCALE: %".format(bzero, bscale).postln;

        // --- Determine data type and read data ---
        datatype = switch(bitpix)
            { 8 } { \uint8 }
            { 16 } { \int16 }
            { 32 } { \int32 }
            { -32 } { \float32 }
            { -64 } { \float64 }
            { \uint8 }; // Default

        dataSize = width * height * (bitpix.abs / 8);

        "Reading FITS: % x % pixels, BITPIX=%, datatype=%, dataSize=%"
            .format(width, height, bitpix, datatype, dataSize).postln;

        rawData = Int8Array.newClear(dataSize.asInteger);
        bytesRead = file.read(rawData);
        "Bytes read: %".format(bytesRead).postln;

        data = this.convertData(rawData, datatype, bzero, bscale);

        file.close;
        file = nil;

        "FITS loaded: % values, range % to %"
            .format(data.size, data.minItem, data.maxItem).postln;
    }

    // Converts the raw byte array from the FITS file into a properly typed and scaled array of numbers.
    convertData { |rawData, dtype, bzero, bscale|
        var converted;
        var numPixels, byte1, byte2, val;
        var b1, b2, b3, b4;
        var floatBits, sign, exponent, mantissa, floatVal;

        "Converting data: dtype=%, bzero=%, bscale=%".format(dtype, bzero, bscale).postln;

        switch(dtype)
        {
            \uint8 } {
                converted = rawData.collect { |byte|
                    val = byte;
                    if(val < 0) { val = val + 256 };
                    (val * bscale) + bzero
                };
            }
            {
            \int16 } {
                numPixels = (rawData.size / 2).asInteger;
                converted = Array.newClear(numPixels);
                numPixels.do { |i|
                    byte1 = rawData[i*2];
                    byte2 = rawData[i*2 + 1];
                    if(byte1 < 0) { byte1 = byte1 + 256 };
                    if(byte2 < 0) { byte2 = byte2 + 256 };
                    val = (byte1 << 8) | byte2; // Big-endian
                    if(val > 32767) { val = val - 65536 };
                    converted[i] = (val * bscale) + bzero;
                };
            }
            {
            \int32 } {
                numPixels = (rawData.size / 4).asInteger;
                converted = Array.newClear(numPixels);
                numPixels.do { |i|
                    b1 = rawData[i*4]; b2 = rawData[i*4+1]; b3 = rawData[i*4+2]; b4 = rawData[i*4+3];
                    if(b1 < 0) { b1 = b1 + 256 }; if(b2 < 0) { b2 = b2 + 256 };
                    if(b3 < 0) { b3 = b3 + 256 }; if(b4 < 0) { b4 = b4 + 256 };
                    val = (b1 << 24) | (b2 << 16) | (b3 << 8) | b4;
                    if(val > 2147483647) { val = val - 4294967296 };
                    converted[i] = (val * bscale) + bzero;
                };
            }
            {
            \float32 } {
                numPixels = (rawData.size / 4).asInteger;
                converted = Array.newClear(numPixels);
                numPixels.do { |i|
                    b1 = rawData[i*4]; b2 = rawData[i*4+1]; b3 = rawData[i*4+2]; b4 = rawData[i*4+3];
                    if(b1 < 0) { b1 = b1 + 256 }; if(b2 < 0) { b2 = b2 + 256 };
                    if(b3 < 0) { b3 = b3 + 256 }; if(b4 < 0) { b4 = b4 + 256 };
                    floatBits = (b1 << 24) | (b2 << 16) | (b3 << 8) | b4;
                    sign = (floatBits >> 31) & 1;
                    exponent = (floatBits >> 23) & 0xFF;
                    mantissa = floatBits & 0x7FFFFF;
                    if(exponent == 0) {
                        if(mantissa == 0) { floatVal = 0.0; } { floatVal = (mantissa / 8388608.0) * (2 ** -126); };
                    } {
                        if(exponent == 255) { floatVal = inf; } { floatVal = (1.0 + (mantissa / 8388608.0)) * (2 ** (exponent - 127)); };
                    };
                    if(sign == 1) { floatVal = floatVal.neg };
                    converted[i] = (floatVal * bscale) + bzero;
                };
            }
            {
                converted = rawData.collect { |byte|
                    val = byte;
                    if(val < 0) { val = val + 256 };
                    (val * bscale) + bzero
                };
            };
        ^converted;
    }

    // --- Data Access Methods ---

    // Retrieves the pixel value at a specific (x, y) coordinate.
    // Returns 0 if the coordinates are out of bounds.
    getPixel { |x, y|
        if(x < 0 or: { x >= width } or: { y < 0 } or: { y >= height }) { ^0; };
        ^data[y * width + x];
    }

    // Extracts a single row of pixel data.
    // Returns an empty array if the row index is out of bounds.
    getRow { |y|
        var startIdx, endIdx, rowData;
        if(y < 0 or: { y >= height }) { ^[] };
        startIdx = y * width;
        endIdx = startIdx + width - 1;
        rowData = Array.newClear(width);
        width.do { |i| rowData[i] = data[startIdx + i]; };
        ^rowData;
    }

    // Extracts a rectangular region of data and returns it as a new Array.
    getRegionData { |startX, startY, regionWidth, regionHeight|
        var regionData = Array.newClear(regionWidth * regionHeight);
        regionHeight.do { |y|
            regionWidth.do { |x|
                var pixel = this.getPixel(startX + x, startY + y);
                regionData[y * regionWidth + x] = pixel;
            };
        };
        ^regionData;
    }

    // ===== STATISTICAL ANALYSIS METHODS =====

    // Calculates fundamental statistics for the entire image.
    // Returns a Dictionary with mean, median, standard deviation, etc.
    calculateStatistics {
        var sortedData, n, median, mad, std, mean;
        "Calculating statistics...".postln;
        if(data.isNil) { "No data available".postln; ^nil; };

        sortedData = data.copy.sort;
        n = data.size;
        mean = data.sum / n;

        median = if(n % 2 == 0) {
            (sortedData[n/2-1] + sortedData[n/2]) / 2
        } {
            sortedData[n/2]
        };

        mad = sortedData.collect { |x| (x - median).abs }.sort[n/2]; // Median Absolute Deviation
        std = sqrt(data.collect { |x| (x - mean).squared }.sum / (n - 1));

        "Statistics calculated.".postln;
        ^Dictionary[
            \mean -> mean,
            \median -> median,
            \std -> std,
            \mad -> mad,
            \min -> data.minItem,
            \max -> data.maxItem,
            \range -> (data.maxItem - data.minItem),
            \count -> n
        ];
    }

    // Performs sigma clipping on the data to remove outliers.
    // Returns a new array containing only the values within the specified sigma range.
    sigmaClip { |sigma=3, maxIterations=5|
        var clippedData = data.copy;
        var mean, std;

        maxIterations.do {
            mean = clippedData.sum / clippedData.size;
            std = sqrt(clippedData.collect { |x| (x - mean).squared }.sum / (clippedData.size - 1));
            clippedData = clippedData.select { |x|
                (x - mean).abs <= (sigma * std)
            };
        };
        ^clippedData;
    }

    // ===== WCS METHODS =====

    // Converts pixel coordinates (x, y) to world coordinates (e.g., RA, Dec).
    // A basic linear transformation is used if WCS headers are present.
    pixelToWorld { |x, y|
        var ra, dec;
        if(crval1.isNil or: crval2.isNil) {
            "No WCS information available".postln;
            ^[x, y];
        };

        // Simple linear transformation (basic WCS)
        ra = crval1 + ((x - (crpix1 ? 0)) * (cdelt1 ? 1));
        dec = crval2 + ((y - (crpix2 ? 0)) * (cdelt2 ? 1));

        ^[ra, dec];
    }

    // Converts world coordinates (e.g., RA, Dec) back to pixel coordinates (x, y).
    worldToPixel { |ra, dec|
        var x, y;
        if(crval1.isNil or: crval2.isNil) {
            "No WCS information available".postln;
            ^[ra, dec];
        };

        x = (crpix1 ? 0) + ((ra - crval1) / (cdelt1 ? 1));
        y = (crpix2 ? 0) + ((dec - crval2) / (cdelt2 ? 1));

        ^[x, y];
    }

    // ===== IMAGE PROCESSING METHODS =====

    // Applies a Gaussian filter to the image data to smooth it.
    // Returns a new FITSFile object with the filtered data.
    gaussianFilter { |sigma=1.0|
        var kernelSize = (6 * sigma).ceil.asInteger;
        var halfKernel = (kernelSize / 2).floor;
        var kernel1D = Array.newClear(kernelSize);
        var normFactor = 0;
        var intermediateData = Array.newClear(data.size); // For horizontal pass results
        var filteredData = Array.newClear(data.size); // Final results

        "Applying Gaussian filter (sigma=%)".format(sigma).postln;
        "Creating 1D Gaussian kernel (size: %)".format(kernelSize).postln;

        // Create 1D Gaussian kernel
        kernelSize.do { |i|
            var x = i - halfKernel;
            var value = (x.squared.neg / (2 * sigma.squared)).exp;
            kernel1D.put(i, value);
            normFactor = normFactor + value;
        };

        // Normalize 1D kernel
        kernelSize.do { |i|
            kernel1D.put(i, kernel1D.at(i) / normFactor);
        };

        "Applying horizontal convolution pass... (This may take a while)".postln;
        // Apply horizontal convolution
        height.do { |row|
            width.do { |col|
                var sum = 0;
                kernelSize.do { |k_idx|
                    var pixelCol = col + k_idx - halfKernel;
                    if((pixelCol >= 0) and: { pixelCol < width }) {
                        var weight = kernel1D.at(k_idx);
                        sum = sum + (this.getPixel(pixelCol, row) * weight);
                    };
                };
                intermediateData[row * width + col] = sum;
            };
        };

        "Applying vertical convolution pass... (This may take a while)".postln;
        // Apply vertical convolution
        height.do { |row|
            width.do { |col|
                var sum = 0;
                kernelSize.do { |k_idx|
                    var pixelRow = row + k_idx - halfKernel;
                    if((pixelRow >= 0) and: { pixelRow < height }) {
                        var weight = kernel1D.at(k_idx);
                        sum = sum + (intermediateData[pixelRow * width + col] * weight);
                    };
                };
                filteredData[row * width + col] = sum;
            };
        };

        "Gaussian filter complete.".postln;
        ^FITSFile.newFromData(filteredData, width, height, headers);
    }

    // Estimates and subtracts the background level from the image data.
    // The background is estimated from the median of corner samples.
    subtractBackground { |boxSize=50|
        var backgroundLevel = 0;
        var samples = [];
        var safeBoxSize = boxSize.min(width/4).min(height/4);

        // Sample background from corners and edges
        safeBoxSize.do { |i|
            safeBoxSize.do { |j|
                samples = samples.add(this.getPixel(j, i)); // Top-left
                samples = samples.add(this.getPixel(width-1-j, i)); // Top-right
                samples = samples.add(this.getPixel(j, height-1-i)); // Bottom-left
                samples = samples.add(this.getPixel(width-1-j, height-1-i)); // Bottom-right
            };
        };

        samples = samples.sort;
        backgroundLevel = samples[samples.size/2]; // median
        data = data.collect { |pixel| pixel - backgroundLevel };
        "Background level % subtracted".format(backgroundLevel.round(0.1)).postln;
        ^this;
    }

    // ===== PHOTOMETRY AND MEASUREMENT METHODS =====

    // Performs aperture photometry on a specified circular region.
    // Calculates total flux, sky-subtracted net flux, and the sky background level.
    aperturePhotometry { |centerX, centerY, radius, skyInner=nil, skyOuter=nil|
        var totalFlux = 0, pixelCount = 0;
        var skyFlux = 0, skyPixels = 0;
        var netFlux, skyLevel;

        // Default sky annulus
        skyInner = skyInner ? (radius * 1.5);
        skyOuter = skyOuter ? (radius * 2.0);

        // Calculate sky background
        height.do { |y|
            width.do { |x|
                var distance = ((x - centerX).squared + (y - centerY).squared).sqrt;
                if((distance >= skyInner) and: { distance <= skyOuter }) {
                    skyFlux = skyFlux + this.getPixel(x, y);
                    skyPixels = skyPixels + 1;
                };
            };
        };

        skyLevel = if(skyPixels > 0) { skyFlux / skyPixels } { 0 };

        // Calculate aperture flux
        height.do { |y|
            width.do { |x|
                var distance = ((x - centerX).squared + (y - centerY).squared).sqrt;
                if(distance <= radius) {
                    totalFlux = totalFlux + this.getPixel(x, y);
                    pixelCount = pixelCount + 1;
                };
            };
        };

        netFlux = totalFlux - (skyLevel * pixelCount);

        ^Dictionary[
            \totalFlux -> totalFlux,
            \netFlux -> netFlux,
            \skyLevel -> skyLevel,
            \pixelCount -> pixelCount,
            \centerX -> centerX,
            \centerY -> centerY,
            \radius -> radius
        ];
    }

    // Finds point sources (local maxima) in the image above a given threshold.
    // Returns an array of Dictionaries, each describing a found source.
    findSources { |threshold=nil, minSeparation=5|
        var sources = [];
        var stats;
        "Finding sources...".postln;

        // Auto-threshold if not provided
        if(threshold.isNil) {
            stats = this.calculateStatistics();
            threshold = stats[\mean] + (3 * stats[\std]);
        };

        "Finding sources with threshold: %".format(threshold.round(0.1)).postln;

        (minSeparation..(height-minSeparation)).do { |y|
            (minSeparation..(width-minSeparation)).do { |x|
                var pixel = this.getPixel(x, y);
                var isMax = true;

                if(pixel > threshold) {
                    // Check if it's a local maximum
                    (-1..1).do { |dy|
                        (-1..1).do { |dx|
                            if(dx != 0 or: dy != 0) {
                                if(this.getPixel(x + dx, y + dy) >= pixel) {
                                    isMax = false;
                                };
                            };
                        };
                    };

                    if(isMax) {
                        // Check minimum separation from existing sources
                        var tooClose = sources.any({ |aSource|
                            var dist = ((aSource[\x] - x).squared + (aSource[\y] - y).squared).sqrt;
                            dist < minSeparation;
                        });

                        if(tooClose.not) {
                            sources = sources.add(Dictionary[\x -> x, \y -> y, \flux -> pixel]);
                        };
                    };
                };
            };
        };

        sources = sources.sort { |a, b| a[\flux] > b[\flux] };
        "Found % sources".format(sources.size).postln;
        ^sources;
    }

    // ===== DATA EXPORT METHODS =====

    // Extracts a rectangular region of the image and saves it to a new file.
    exportRegion { |startX, startY, regionWidth, regionHeight, outputPath|
        var regionData = Array.newClear(regionWidth * regionHeight);
        var file;

        regionHeight.do { |y|
            regionWidth.do { |x|
                var pixel = this.getPixel(startX + x, startY + y);
                regionData[y * regionWidth + x] = pixel;
            };
        };

        file = File(outputPath, "w");
        if(file.isNil) { "Could not create output file: %".format(outputPath).postln; ^this; };
        regionData.do { |val| file.write(val.asString ++ "\n"); };
        file.close;
        "Region exported to: %".format(outputPath).postln;
        ^this;
    }

    // Calculates the average pixel value in concentric annuli around a center point.
    // Returns a Dictionary mapping radius to average flux.
    radialProfile { |centerX, centerY, maxRadius=nil|
        var profile = Dictionary.new;
        var radiusBins;

        maxRadius = maxRadius ? (width.min(height) / 2);
        radiusBins = Array.series(maxRadius.asInteger, 0, 1);

        radiusBins.do { |radius|
            var fluxSum = 0, pixelCount = 0;

            height.do { |y|
                width.do { |x|
                    var distance = ((x - centerX).squared + (y - centerY).squared).sqrt;
                    if((distance >= radius) and: { distance < (radius + 1) }) {
                        fluxSum = fluxSum + this.getPixel(x, y);
                        pixelCount = pixelCount + 1;
                    };
                };
            };

            if(pixelCount > 0) {
                profile[radius] = fluxSum / pixelCount;
            };
        };

        ^profile;
    }

    // Helper function to map a 1D Hilbert curve index to 2D coordinates.
    hilbert_d2xy { |d, n|
        var x = 0, y = 0;
        var rx, ry, s = 1;
        var temp_x, temp_y, temp; // Declare all vars at the top
        n.do { |i|
            rx = (d >> 1) bitAnd: 1;
            ry = (d bitXor: rx) bitAnd: 1;
            // rotate and add
            temp_x = x;
            temp_y = y;
            if (ry == 0) {
                if (rx == 1) {
                    temp_x = s - 1 - x;
                    temp_y = s - 1 - y;
                };
                // swap x and y
                temp = temp_x;
                temp_x = temp_y;
                temp_y = temp;
            };
            x = temp_x + (s * rx);
            y = temp_y + (s * ry);
            d = d >> 2;
            s = s * 2;
        };
        ^[x, y];
    }

    // Remaps the 2D pixel data into a 1D sequence using a Hilbert curve.
    // This can be useful for sonification or alternative visualizations.
    // Requires the image to be square with power-of-two dimensions.
    getHilbertCurveData {
        var hilbertData, n, totalPixels;
        "Generating Hilbert curve data...".postln;

        if (width != height) {
            "Warning: Image is not square. Returning linear data.".postln;
            ^data.copy;
        };

        n = width.log2;
        if (n.isInteger.not) {
            "Warning: Image dimensions are not a power of 2. Returning linear data.".postln;
            ^data.copy;
        };

        totalPixels = width * height;
        hilbertData = Array.newClear(totalPixels);

        totalPixels.do { |i|
            var coords = this.hilbert_d2xy(i, n);
            var x = coords[0];
            var y = coords[1];
            hilbertData[i] = this.getPixel(x, y);
        };
        "Hilbert curve data generated.".postln;
        ^hilbertData;
    }

    // Remaps pixel data into a 1D sequence using a Hilbert curve, handling arbitrary dimensions
    // by padding to the next power-of-two square.
    getHilbertCurveDataArbitrary {
        var hilbertData, n, effectiveSide, totalEffectivePixels;
        "Generating Hilbert curve data for arbitrary dimensions...".postln;

        effectiveSide = (width.max(height)).nextPowerOfTwo;
        n = effectiveSide.log2;
        totalEffectivePixels = effectiveSide * effectiveSide;
        hilbertData = Array.newClear(totalEffectivePixels);

        "Using effective Hilbert grid size: % x %".format(effectiveSide, effectiveSide).postln;

        totalEffectivePixels.do { |i|
            var coords = this.hilbert_d2xy(i, n);
            var x = coords[0];
            var y = coords[1];

            if ((x >= 0) and: { x < width } and: { (y >= 0) and: { y < height } }) {
                hilbertData[i] = this.getPixel(x, y);
            } {
                hilbertData[i] = 0; // Pad with 0
            };
        };
        "Hilbert curve data generated for arbitrary dimensions.".postln;
        ^hilbertData;
    }

    // Exports the entire dataset to a simple text file, one value per line.
    exportToFile { |outputPath|
        var file;
        if(data.isNil) { "No data to export".postln; ^this; };
        file = File(outputPath, "w");
        if(file.isNil) { "Could not create output file: %".format(outputPath).postln; ^this; };
        "Exporting % values to %".format(data.size, outputPath).postln;
        data.do { |val| file.write(val.asString ++ "\n"); };
        file.close;
        "Export complete".postln;
    }

    // Exports the image data to a CSV file.
    exportToCSV { |outputPath|
        var file, rowData, rowString;
        if(data.isNil) { "No data to export".postln; ^this; };
        file = File(outputPath, "w");
        if(file.isNil) { "Could not create output file: %".format(outputPath).postln; ^this; };
        "Exporting % x % image to CSV: %".format(height, width, outputPath).postln;
        height.do { |row|
            rowData = this.getRow(row);
            rowString = rowData.collect(_.asString).join(",");
            file.write(rowString ++ "\n");
        };
        file.close;
        "CSV export complete".postln;
    }

    // Prints a detailed summary of the FITS file's properties and headers.
    printSummary {
        "=== FITS File Summary ===".postln;
        "Dimensions: % x % pixels".format(width, height).postln;
        "Data type: % (BITPIX=%)".format(datatype, bitpix).postln;
        "Total pixels: %".format(data.size).postln;
        if(data.notNil) {
            "Value range: % to %".format(data.minItem, data.maxItem).postln;
            "Mean: %".format(data.sum / data.size).postln;
        };
        "Key Headers:".postln;
        ["NAXIS1", "NAXIS2", "BITPIX", "BZERO", "BSCALE"].do { |key|
            if(headers[key].notNil) { "  %: %".format(key, headers[key]).postln; };
        };
        "WCS Headers (if present):".postln;
        ["CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2",
         "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2",
         "PC1_1", "PC1_2", "PC2_1", "PC2_2"].do { |key|
            if(headers[key].notNil) { "  %: %".format(key, headers[key]).postln; };
        };
        "".postln;
    }

    // Compares the first few values of the data array with an expected array.
    // Useful for testing data integrity.
    compareWithExpected { |expectedValues|
        var matches = 0;
        if(data.isNil) { "No data loaded".postln; ^this; };
        "=== Data Comparison ===".postln;
        "Expected first 10: %".format(expectedValues[0..9]).postln;
        "Actual first 10:   %".format(data[0..9]).postln;
        10.do { |i| if(data[i] == expectedValues[i]) { matches = matches + 1; }; };
        "Matches: %/10".format(matches).postln;
        "".postln;
    }

    // Returns a FITSVisualisation instance for this FITSFile.
    visualiser {
        ^FITSVisualisation.new(this);
    }
}