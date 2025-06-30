FITSVisualisation {
    var <fitsFile; // Reference to the FITSFile instance

    // --- Instance variables for GUI state and windows ---
    var <roiWindow, <roiStartX, <roiEndX, <roiStartY, <roiEndY, <isSelectingROI=false;
    var <currentRoiData, <currentRoiCenterX, <currentRoiCenterY, <currentRoiMaxRadius;
    var <interactiveWindow, <rowsWindow, <histogramWindow, <imageWindow, <regionWindow;
    var <colorWindow, <contourWindow, <colormapWindow, <radialProfileWindow;
    var <displayMin, <displayMax; // For interactive visualization level adjustment


    *new { |fitsFileInstance|
        ^super.new.init(fitsFileInstance);
    }

    init { |fitsFileInstance|
        fitsFile = fitsFileInstance;
    }



// Add this new method to compare multiple histograms:
visualizeMultipleHistograms { |fitsFiles, titles, colors, numBins=50, width=800, height=600|
    var windowWidth = width, windowHeight = height;
    var allData = [], allTitles = [], allColors = [];
    var globalMin = inf, globalMax = -inf;
    var histograms = [], maxCounts = [];
    var win, userView;
    var globalRange, globalMaxCount; // Declare these at the top

    // Collect all data and find global range
    fitsFiles.do { |fits, i|
        var data = fits.data;
        var title = if (titles.notNil and: (i < titles.size)) { titles[i] } { "Dataset %" };
        var color = if (colors.notNil and: (i < colors.size)) { colors[i] } { Color.hsv(i * 0.3 % 1.0, 0.8, 0.9) };

        if (data.notNil and: data.size > 0) {
            allData = allData.add(data);
            allTitles = allTitles.add(title);
            allColors = allColors.add(color);
            globalMin = min(globalMin, data.minItem);
            globalMax = max(globalMax, data.maxItem);
        };
    };

    if (allData.size == 0) {
        "No valid data for histogram comparison.".postln;
        ^this;
    };

    // Calculate histograms for each dataset
    globalRange = globalMax - globalMin;
    if (globalRange == 0) { globalRange = 1 };

    allData.do { |data, i|
        var histogram = Array.fill(numBins, 0);

        data.do { |val|
            var bin = ((val - globalMin) / globalRange * (numBins - 1)).floor.clip(0, numBins - 1);
            histogram[bin] = histogram[bin] + 1;
        };

        histograms = histograms.add(histogram);
        maxCounts = maxCounts.add(histogram.maxItem);
    };

    globalMaxCount = maxCounts.maxItem;
    if (globalMaxCount == 0) { globalMaxCount = 1 };

    {
        win = Window("FITS Histogram Comparison", Rect(100, 100, windowWidth, windowHeight + 100));
        win.background = Color.white;

        userView = UserView(win, Rect(0, 0, windowWidth, windowHeight + 100));
        userView.background = Color.white;

        userView.drawFunc = {
            var plotBounds = Rect(60, 40, windowWidth - 120, windowHeight - 80);
            var barWidth = plotBounds.width / numBins;

            // Draw axes
            Pen.color = Color.black;
            Pen.width = 2;

            // X-axis
            Pen.line(
                Point(plotBounds.left, plotBounds.bottom),
                Point(plotBounds.right, plotBounds.bottom)
            );

            // Y-axis
            Pen.line(
                Point(plotBounds.left, plotBounds.top),
                Point(plotBounds.left, plotBounds.bottom)
            );
            Pen.stroke;

            // Draw histograms (overlaid with transparency)
            histograms.do { |histogram, datasetIndex|
                var color = allColors[datasetIndex];

                Pen.color = Color(color.red, color.green, color.blue, 0.6); // Add transparency
                Pen.width = 1;

                histogram.do { |count, binIndex|
                    var barHeight = (count / globalMaxCount) * plotBounds.height;
                    var x = plotBounds.left + (binIndex * barWidth);
                    var y = plotBounds.bottom - barHeight;

                    Pen.fillRect(Rect(x, y, barWidth - 1, barHeight));
                };
            };

            // Draw axis labels
            Pen.color = Color.black;
            Pen.font = Font.monospace(10);

            // X-axis labels
            5.do { |i|
                var value = globalMin + (globalRange * i / 4);
                var x = plotBounds.left + (plotBounds.width * i / 4);
                Pen.stringAtPoint(value.round(0.1).asString, Point(x - 20, plotBounds.bottom + 10));
            };

            // Y-axis labels
            5.do { |i|
                var count = globalMaxCount * i / 4;
                var y = plotBounds.bottom - (plotBounds.height * i / 4);
                Pen.stringAtPoint(count.round(1).asString, Point(plotBounds.left - 40, y - 5));
            };

            // Title
            Pen.font = Font.monospace(14);
            Pen.stringAtPoint("FITS Data Histogram Comparison",
                Point(plotBounds.left + plotBounds.width/2 - 100, plotBounds.top - 25));

            // Legend
            allTitles.do { |title, i|
                var color = allColors[i];
                var legendY = plotBounds.bottom + 40 + (i * 20);
                var legendX = plotBounds.left;

                Pen.color = color;
                Pen.fillRect(Rect(legendX, legendY, 15, 15));

                Pen.color = Color.black;
                Pen.font = Font.monospace(10);
                Pen.stringAtPoint(title, Point(legendX + 20, legendY + 2));
            };

            // Statistics
            Pen.color = Color.black;
            Pen.font = Font.monospace(9);
            Pen.stringAtPoint("Range: % to %".format(globalMin.round(0.01), globalMax.round(0.01)),
                Point(plotBounds.right - 150, plotBounds.top - 10));
            Pen.stringAtPoint("Bins: %".format(numBins),
                Point(plotBounds.right - 80, plotBounds.top - 10));
        };

        win.front;
        userView.refresh;
        histogramWindow = win;
    }.defer;

    ^this;
}




    plotRadialProfile { |profileDict, title|
        var radii, values, plotter;

        radii = profileDict.keys.asArray.sort;
        if (radii.isEmpty) {
            "No data to plot for radial profile.".postln; ^this;
        };

        values = radii.collect { |r| profileDict[r] };

        {
            // Close previous window if it exists
            if(radialProfileWindow.notNil) {
                if(radialProfileWindow.isClosed.not) {
                    radialProfileWindow.close;
                };
            };

            plotter = Plotter(title, Rect(200, 200, 500, 350));
            plotter.plotMode = \linear;
            plotter.value = radii.collect({|r, i| [r, values[i]] });
            radialProfileWindow = plotter.parent;
            radialProfileWindow.front;
        }.defer;

        ^this;
    }

    // Plots multiple radial profiles on a single Plotter window.
    // profileDicts: An Array of Dictionaries, where each Dictionary is a radial profile.
    // titles: An Array of Strings, corresponding to the titles for each profile.
    // colors: Optional. An Array of Colors, one for each profile. If nil or size mismatch, random colors are used.
    // Add this method to your FITSVisualisation class
// Replace the plotMultipleRadialProfiles method with this corrected version
plotMultipleRadialProfiles { |profileDataArray, labels, colors|
    var plotWindow, plotView;
    var allRadii = [], allValues = [];
    var minVal = inf, maxVal = -inf;
    var maxRadius = 0;
    var valueRange, padding;

    "=== PLOTTING RADIAL PROFILES ===".postln;

    // Collect all data for scaling
    profileDataArray.do { |profileData, i|
        var radii = profileData[\radii];
        var values = profileData[\values];

        if (radii.notNil and: values.notNil and: radii.size > 0) {
            allRadii = allRadii.add(radii);
            allValues = allValues.add(values);
            minVal = min(minVal, values.minItem);
            maxVal = max(maxVal, values.maxItem);
            maxRadius = max(maxRadius, radii.maxItem);
            "Profile %: % points, range: % to %".format(
                i, radii.size, values.minItem.round(0.01), values.maxItem.round(0.01)
            ).postln;
        };
    };

    if (allRadii.size == 0) {
        "No valid profile data to plot!".postln;
        ^this;
    };

    // Add some padding to the value range for better visualization
    valueRange = maxVal - minVal;
    padding = valueRange * 0.1; // 10% padding
    minVal = minVal - padding;
    maxVal = maxVal + padding;

    // Create a larger plotting window
    try {
        plotWindow = Window("Radial Profile Plot", Rect(100, 100, 1000, 700));
        plotWindow.background = Color.white;

        // Create plot view that fills most of the window
        plotView = UserView(plotWindow, Rect(0, 0, plotWindow.bounds.width, plotWindow.bounds.height));
        plotView.background = Color.white;

        plotView.drawFunc = {
            // Define plot area with proper margins from the edges
            var plotBounds = Rect(80, 60, 700, 500); // Fixed position: left=80, top=60, width=700, height=500
            var xScale = plotBounds.width / maxRadius;
            var yScale = plotBounds.height / (maxVal - minVal);
            var numVLines, numHLines;

            // Draw axes
            Pen.color = Color.black;
            Pen.width = 2;

            // X-axis
            Pen.line(
                Point(plotBounds.left, plotBounds.bottom),
                Point(plotBounds.right, plotBounds.bottom)
            );

            // Y-axis
            Pen.line(
                Point(plotBounds.left, plotBounds.top),
                Point(plotBounds.left, plotBounds.bottom)
            );

            Pen.stroke;

            // Draw grid lines
            Pen.color = Color.gray(0.9);
            Pen.width = 1;

            // Vertical grid lines (every 5 radius units)
            numVLines = (maxRadius / 5).ceil.asInteger;
            numVLines.do { |i|
                var radius = i * 5;
                if (radius <= maxRadius) {
                    var x = plotBounds.left + (radius * xScale);
                    Pen.line(Point(x, plotBounds.top), Point(x, plotBounds.bottom));
                };
            };

            // Horizontal grid lines
            numHLines = 8;
            numHLines.do { |i|
                var value = minVal + ((maxVal - minVal) * i / (numHLines - 1));
                var y = plotBounds.bottom - ((value - minVal) * yScale);
                Pen.line(Point(plotBounds.left, y), Point(plotBounds.right, y));
            };

            Pen.stroke;

            // Draw axis labels and tick marks
            Pen.color = Color.black;
            Pen.width = 1;

            // X-axis tick marks and labels
            numVLines.do { |i|
                var radius = i * 5;
                if (radius <= maxRadius) {
                    var x = plotBounds.left + (radius * xScale);
                    // Tick mark
                    Pen.line(Point(x, plotBounds.bottom), Point(x, plotBounds.bottom + 5));
                    // Label
                    Pen.stringAtPoint(radius.asString, Point(x - 10, plotBounds.bottom + 15));
                };
            };

            // Y-axis tick marks and labels
            numHLines.do { |i|
                var value = minVal + ((maxVal - minVal) * i / (numHLines - 1));
                var y = plotBounds.bottom - ((value - minVal) * yScale);
                // Tick mark
                Pen.line(Point(plotBounds.left - 5, y), Point(plotBounds.left, y));
                // Label
                Pen.stringAtPoint(value.round(0.1).asString, Point(plotBounds.left - 50, y - 8));
            };

            Pen.stroke;

            // Draw profiles
            allRadii.do { |radii, profileIndex|
                var values = allValues[profileIndex];
                var color = if (colors.notNil and: colors.size > profileIndex) {
                    colors[profileIndex];
                } {
                    Color.hsv(profileIndex * 0.3 % 1.0, 0.8, 0.9);
                };

                Pen.color = color;
                Pen.width = 3;

                // Draw the profile line
                radii.do { |radius, i|
                    if (i < values.size) {
                        var x = plotBounds.left + (radius * xScale);
                        var y = plotBounds.bottom - ((values[i] - minVal) * yScale);

                        if (i == 0) {
                            Pen.moveTo(Point(x, y));
                        } {
                            Pen.lineTo(Point(x, y));
                        };
                    };
                };

                Pen.stroke;

                // Draw points
                Pen.color = color;
                radii.do { |radius, i|
                    if (i < values.size) {
                        var x = plotBounds.left + (radius * xScale);
                        var y = plotBounds.bottom - ((values[i] - minVal) * yScale);
                        Pen.fillOval(Rect(x-3, y-3, 6, 6));
                    };
                };
            };

            // Draw labels positioned relative to plot area
            Pen.color = Color.black;

            // Title - positioned above the plot
            Pen.stringAtPoint("Radial Profile Analysis",
                Point(plotBounds.left + plotBounds.width/2 - 70, plotBounds.top - 40));

            // Data range info - positioned below title
            Pen.stringAtPoint("Range: % to % (% points)".format(
                (minVal + padding).round(0.1),
                (maxVal - padding).round(0.1),
                allRadii[0].size
            ), Point(plotBounds.left + plotBounds.width/2 - 80, plotBounds.top - 20));

            // X-axis label - positioned below the plot
            Pen.stringAtPoint("Radius (pixels)",
                Point(plotBounds.left + plotBounds.width/2 - 50, plotBounds.bottom + 40));

            // Y-axis label - positioned to the left of the plot
            Pen.stringAtPoint("Intensity Value",
                Point(plotBounds.left - 70, plotBounds.top + plotBounds.height/2 - 10));

            // Legend - positioned to the right of the plot
            if (labels.notNil) {
                labels.do { |label, i|
                    if (i < allRadii.size) {
                        var color = if (colors.notNil and: colors.size > i) {
                            colors[i];
                        } {
                            Color.hsv(i * 0.3 % 1.0, 0.8, 0.9);
                        };

                        var legendY = plotBounds.top + (i * 25) + 20;
                        var legendX = plotBounds.right + 20;

                        Pen.color = color;
                        Pen.fillOval(Rect(legendX, legendY, 12, 12));

                        Pen.color = Color.black;
                        Pen.stringAtPoint(label, Point(legendX + 20, legendY - 2));
                    };
                };
            } {
                // Default legend if no labels provided
                allRadii.do { |radii, i|
                    var color = Color.hsv(i * 0.3 % 1.0, 0.8, 0.9);
                    var legendY = plotBounds.top + (i * 25) + 20;
                    var legendX = plotBounds.right + 20;

                    Pen.color = color;
                    Pen.fillOval(Rect(legendX, legendY, 12, 12));

                    Pen.color = Color.black;
                    Pen.stringAtPoint("Profile %".format(i + 1),
                        Point(legendX + 20, legendY - 2));
                };
            };
        };

        plotWindow.front;
        plotView.refresh;

        "Plot window created successfully with dimensions: % x %".format(
            plotWindow.bounds.width, plotWindow.bounds.height
        ).postln;
        "Plot positioned at: left=80, top=60, size=700x500".postln;

    } { |error|
        "Error creating plot window: %".format(error).postln;
        "Falling back to simple data display...".postln;
        this.printRadialProfileData(profileDataArray, labels);
    };
}





// Add this fallback method for when plotting fails
printRadialProfileData { |profileDicts, titles|
    "=== RADIAL PROFILE DATA ===".postln;

    profileDicts.do { |profileDict, i|
        var radii = profileDict[\radii];
        var values = profileDict[\values];
        var title = if (titles.notNil and: (i < titles.size)) { titles[i] } { "Profile %" };

        "Profile: %".format(title).postln;
        "Radius\tValue".postln;

        if (radii.notNil and: values.notNil) {
            radii.do { |radius, j|
                if (j < values.size) {
                    "%\t%".format(radius.round(0.1), values[j].round(0.01)).postln;
                };
            };
        };
        "".postln;
    };
}


    showROIStats { |stats, x, y, w, h|
        var statsWin, textView, statsText;

        statsText = "=== ROI ANALYSIS RESULTS ===\n\n";
        statsText = statsText ++ "Region: (%, %) to (%, %)\n".format(x, y, x+w-1, y+h-1);
        statsText = statsText ++ "Size: % x % pixels (% total)\n\n".format(w, h, w*h);
        statsText = statsText ++ "STATISTICS:\n";
        statsText = statsText ++ "Mean Value: %\n".format(stats[\mean].round(0.01));
        statsText = statsText ++ "Minimum: %\n".format(stats[\min].round(0.01));
        statsText = statsText ++ "Maximum: %\n".format(stats[\max].round(0.01));
        statsText = statsText ++ "Std Dev: %\n".format(stats[\std].round(0.01));
        statsText = statsText ++ "Range: %\n".format((stats[\max] - stats[\min]).round(0.01));
        statsText = statsText ++ "\nSNR (approx): %\n".format((stats[\mean] / stats[\std]).round(0.1));

        {
            statsWin = Window("ROI Statistics", Rect(100, 100, 400, 300));
            textView = TextView(statsWin, statsWin.bounds.insetBy(10, 10));
            textView.string = statsText;
            textView.font = Font.monospace(12);
            textView.editable = false;
            statsWin.front;
            roiWindow = statsWin;
        }.defer;

        ^this;
    }

    visualizeInteractive { |downsample=2, windowWidth=800, windowHeight=600|
    var newWidth = (fitsFile.width / downsample).floor, newHeight = (fitsFile.height / downsample).floor;
    var scaleX = windowWidth / newWidth, scaleY = windowHeight / newHeight;
    var minVal = fitsFile.data.minItem, maxVal = fitsFile.data.maxItem;
    var win, userView, infoText;
    var mouseRow = 0, mouseCol = 0, mousePixel = 0;
    var plotLightCurveButton, plotEnergySpectrumButton, plotRadialProfileButton;
    var minSlider, maxSlider, minLabel, maxLabel;
    var minControlSpec, maxControlSpec; // ControlSpec objects for the sliders

    // Initialize displayMin and displayMax to the full range of the data
    displayMin = minVal;
    displayMax = maxVal;

    // Create ControlSpec objects for the sliders
    minControlSpec = ControlSpec(minVal, maxVal, \linear, (maxVal - minVal) / 1000);
    maxControlSpec = ControlSpec(minVal, maxVal, \linear, (maxVal - minVal) / 1000);

    {
        win = Window("FITS Image - Interactive (ROI Select)", Rect(400, 400, windowWidth, windowHeight + 120)); // Increased height for sliders
        infoText = StaticText(win, Rect(10, windowHeight + 5, windowWidth - 20, 20)).string_("Click and drag to select ROI").font_(Font.monospace(12));

        // Min Value Slider - Changed to mouseUpAction
        minLabel = StaticText(win, Rect(10, windowHeight + 30, 80, 20)).string_("Min Value:");
        minSlider = Slider(win, Rect(90, windowHeight + 30, windowWidth - 100, 20))
            .value_(minControlSpec.unmap(displayMin))
            .mouseUpAction_({ |sl|  // Changed from action_ to mouseUpAction_
                displayMin = minControlSpec.map(sl.value);
                if (displayMin > displayMax) {
                    displayMax = displayMin;
                    maxSlider.value = maxControlSpec.unmap(displayMax);
                };
                userView.refresh;
                "Min value updated to: %".format(displayMin.round(0.01)).postln;
            });

        // Max Value Slider - Changed to mouseUpAction
        maxLabel = StaticText(win, Rect(10, windowHeight + 60, 80, 20)).string_("Max Value:");
        maxSlider = Slider(win, Rect(90, windowHeight + 60, windowWidth - 100, 20))
            .value_(maxControlSpec.unmap(displayMax))
            .mouseUpAction_({ |sl|  // Changed from action_ to mouseUpAction_
                displayMax = maxControlSpec.map(sl.value);
                if (displayMax < displayMin) {
                    displayMin = displayMax;
                    minSlider.value = minControlSpec.unmap(displayMin);
                };
                userView.refresh;
                "Max value updated to: %".format(displayMax.round(0.01)).postln;
            });

        plotLightCurveButton = Button(win, Rect(10, windowHeight + 90, 180, 25))
            .states_([["Plot Light Curve (32 bins)"]])
            .action_({ |btn|
                if (currentRoiData.notNil) {
                    this.plotHistogram(currentRoiData, "ROI Light Curve (Histogram, 32 bins)", 32, 512);
                } { "No ROI data available. Select an ROI first.".postln; };
            });

        plotEnergySpectrumButton = Button(win, Rect(200, windowHeight + 90, 200, 25))
            .states_([["Plot Energy Spectrum (256 bins)"]])
            .action_({ |btn|
                if (currentRoiData.notNil) {
                    this.plotHistogram(currentRoiData, "ROI Energy Spectrum (Histogram, 256 bins)", 256, 512);
                } { "No ROI data available. Select an ROI first.".postln; };
            });

        plotRadialProfileButton = Button(win, Rect(410, windowHeight + 90, 150, 25))
            .states_([["Plot Radial Profile"]])
            .action_({ |btn|
                var radialProf;
                if (currentRoiData.notNil and: currentRoiCenterX.notNil) {
                    radialProf = fitsFile.radialProfile(currentRoiCenterX, currentRoiCenterY, currentRoiMaxRadius);
                    this.plotRadialProfile(radialProf, "ROI Radial Profile");
                } { "No ROI data available. Select an ROI first.".postln; };
            });

        userView = UserView(win, Rect(0, 0, windowWidth, windowHeight)).drawFunc_({ Pen.use {
                var currentRange = displayMax - displayMin;
                if (currentRange == 0) { currentRange = 1.0; }; // Avoid division by zero

                newHeight.do({ |r|
                    newWidth.do({ |c|
                        var pixel = fitsFile.getPixel(c * downsample, r * downsample);
                        var normalized = (pixel - displayMin) / currentRange;
                        Pen.fillColor = Color.gray(normalized.clip(0, 1));
                        Pen.fillRect(Rect(c * scaleX, r * scaleY, scaleX + 1, scaleY + 1));
                    });
                });
                Pen.strokeColor = Color.red;
                Pen.width = 1;
                Pen.line(Point(mouseCol * scaleX, 0), Point(mouseCol * scaleX, windowHeight));
                Pen.line(Point(0, mouseRow * scaleY), Point(windowWidth, mouseRow * scaleY));
                Pen.stroke;

                if (isSelectingROI or: { roiStartX.notNil and: roiEndX.notNil }) {
                    Pen.strokeColor = Color.green;
                    Pen.width = 2;
                    Pen.strokeRect(Rect(
                        roiStartX.min(roiEndX) * scaleX,
                        roiStartY.min(roiEndY) * scaleY,
                        (roiEndX - roiStartX).abs * scaleX,
                        (roiEndY - roiStartY).abs * scaleY
                    ));
                    Pen.stroke;
                };
            };
        }).mouseDownAction_({ |view, x, y|
            isSelectingROI = true;
            roiStartX = (x / scaleX).floor.clip(0, newWidth - 1);
            roiStartY = (y / scaleY).floor.clip(0, newHeight - 1);
            roiEndX = roiStartX;
            roiEndY = roiStartY;
            view.refresh;
        }).mouseMoveAction_({ |view, x, y|
            var origRow, origCol;
            mouseCol = (x / scaleX).floor.clip(0, newWidth - 1);
            mouseRow = (y / scaleY).floor.clip(0, newHeight - 1);
            origRow = mouseRow * downsample;
            origCol = mouseCol * downsample;
            mousePixel = fitsFile.getPixel(origCol, origRow);
            infoText.string = "Row: %, Col: %, Value: %".format(origRow, origCol, mousePixel.round(0.1));

            if (isSelectingROI) {
                roiEndX = (x / scaleX).floor.clip(0, newWidth - 1);
                roiEndY = (y / scaleY).floor.clip(0, newHeight - 1);
            };
            view.refresh;
        }).mouseUpAction_({ |view, x, y|
            var x1, y1, x2, y2, roiWidth, roiHeight;
            isSelectingROI = false;
            roiEndX = (x / scaleX).floor.clip(0, newWidth - 1);
            roiEndY = (y / scaleY).floor.clip(0, newHeight - 1);

            x1 = roiStartX.min(roiEndX);
            y1 = roiStartY.min(roiEndY);
            x2 = roiStartX.max(roiEndX);
            y2 = roiEndY.max(roiEndY);

            roiWidth = (x2 - x1).abs + 1;
            roiHeight = (y2 - y1).abs + 1;

            currentRoiData = fitsFile.getRegionData(x1, y1, roiWidth, roiHeight);
            currentRoiCenterX = x1 + (roiWidth / 2);
            currentRoiCenterY = y1 + (roiHeight / 2);
            currentRoiMaxRadius = (roiWidth.max(roiHeight) / 2).ceil;

            infoText.string = "ROI selected. Ready for analysis.";
            "ROI Selected: x=%, y=%, w=%, h=%".format(x1, y1, roiWidth, roiHeight).postln;
            view.refresh;
        });
        win.front;
        interactiveWindow = win;
    }.defer;
}


    visualizeRows { |startRow=0, numRows=50, rowHeight=8|
        var windowWidth = fitsFile.width, windowHeight = numRows * rowHeight;
        var minVal = fitsFile.data.minItem, maxVal = fitsFile.data.maxItem;
        var win, userView;

        {
            win = Window("FITS Image - Rows", Rect(500, 500, windowWidth, windowHeight));
            userView = UserView(win, Rect(0, 0, windowWidth, windowHeight)).drawFunc_({
                Pen.use {
                    numRows.do({ |r|
                        var rowData = fitsFile.getRow(startRow + r);
                        fitsFile.width.do({ |c|
                            var pixel = rowData[c];
                            var normalized = (pixel - minVal) / (maxVal - minVal);
                            Pen.fillColor = Color.gray(normalized.clip(0, 1));
                            Pen.fillRect(Rect(c, r * rowHeight, 1, rowHeight));
                        });
                    });
                };
            });
            win.front;
        }.defer;
        rowsWindow = win;
    }

      // Replace the visualizeHistogram method in your FITSVisualisation class with this corrected version:

visualizeHistogram { |numBins=50, title="FITS Data Histogram", width=600, height=400|
    var data, histogram, minVal, maxVal, range;
    var win, userView;
    var histogramTitle = title;

    if (fitsFile.isNil or: fitsFile.data.isNil) {
        "No FITS data available for histogram.".postln;
        ^this;
    };

    data = fitsFile.data;
    if (data.size == 0) {
        "FITS data is empty.".postln;
        ^this;
    };

    // Calculate histogram
    minVal = data.minItem;
    maxVal = data.maxItem;
    range = maxVal - minVal;

    if (range == 0) {
        "Data has no variation (all values are the same).".postln;
        ^this;
    };

    histogram = Array.fill(numBins, 0);

    data.do { |val|
        var bin = ((val - minVal) / range * (numBins - 1)).floor.clip(0, numBins - 1);
        histogram[bin] = histogram[bin] + 1;
    };

    {
        win = Window(histogramTitle, Rect(100, 100, width, height));
        win.background = Color.white;

        userView = UserView(win, Rect(0, 0, width, height));
        userView.background = Color.white;

        userView.drawFunc = {
            var plotBounds = Rect(60, 40, width - 120, height - 80);
            var barWidth = plotBounds.width / numBins;
            var maxCount = histogram.maxItem;

            if (maxCount == 0) { maxCount = 1 };

            // Draw axes
            Pen.color = Color.black;
            Pen.width = 2;

            // X-axis
            Pen.line(
                Point(plotBounds.left, plotBounds.bottom),
                Point(plotBounds.right, plotBounds.bottom)
            );

            // Y-axis
            Pen.line(
                Point(plotBounds.left, plotBounds.top),
                Point(plotBounds.left, plotBounds.bottom)
            );
            Pen.stroke;

            // Draw histogram bars
            Pen.color = Color.blue;
            Pen.width = 1;

            histogram.do { |count, binIndex|
                var barHeight = (count / maxCount) * plotBounds.height;
                var x = plotBounds.left + (binIndex * barWidth);
                var y = plotBounds.bottom - barHeight;

                Pen.fillRect(Rect(x, y, barWidth - 1, barHeight));
            };

            // Draw axis labels
            Pen.color = Color.black;
            Pen.font = Font.monospace(10);

            // X-axis labels (5 tick marks)
            5.do { |i|
                var value = minVal + (range * i / 4);
                var x = plotBounds.left + (plotBounds.width * i / 4);
                Pen.stringAtPoint(value.round(0.1).asString, Point(x - 20, plotBounds.bottom + 10));
            };

            // Y-axis labels (5 tick marks)
            5.do { |i|
                var count = maxCount * i / 4;
                var y = plotBounds.bottom - (plotBounds.height * i / 4);
                Pen.stringAtPoint(count.round(1).asString, Point(plotBounds.left - 40, y - 5));
            };

            // Title
            Pen.font = Font.monospace(14);
            Pen.stringAtPoint(histogramTitle,
                Point(plotBounds.left + plotBounds.width/2 - (histogramTitle.size * 4), plotBounds.top - 25));

            // Statistics
            Pen.color = Color.black;
            Pen.font = Font.monospace(9);
            Pen.stringAtPoint("Min: %".format(minVal.round(0.01)),
                Point(plotBounds.right - 150, plotBounds.top - 10));
            Pen.stringAtPoint("Max: %".format(maxVal.round(0.01)),
                Point(plotBounds.right - 80, plotBounds.top - 10));
            Pen.stringAtPoint("Bins: %".format(numBins),
                Point(plotBounds.right - 150, plotBounds.top + 5));
            Pen.stringAtPoint("Count: %".format(data.size),
                Point(plotBounds.right - 80, plotBounds.top + 5));
        };

        win.front;
        userView.refresh;
        histogramWindow = win;
    }.defer;

    ^this;
}


    closeAllWindows {
        [imageWindow, regionWindow, colorWindow, interactiveWindow, rowsWindow, histogramWindow, contourWindow, colormapWindow, roiWindow, radialProfileWindow].do({
            |win|
            if (win.notNil) { win.close; };
        });
        "All visualization windows closed".postln;
    }

    testPlotter {
        var testData;
        testData = [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 4, 6, 3, 7, 2, 8, 1, 9, 0]; // Simple line data
        testData.plot("Test Plotter");
        "Test Plotter window opened.".postln;
    }
}
