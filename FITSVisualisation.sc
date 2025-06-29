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

    // A general-purpose histogram plotting method for any data array.
    plotHistogram { |plotData, title, numBins=50, width=600, height=400|
        var windowWidth = width, windowHeight = height;
        var minVal, maxVal, range;
        var histogram, maxCount, barWidth;
        var win, userView;

        if(plotData.isNil or: {plotData.isEmpty}) {
            "No data provided for histogram.".postln;
            ^this;
        };

        minVal = plotData.minItem;
        maxVal = plotData.maxItem;
        range = maxVal - minVal;
        if(range == 0) { range = 1 }; // Avoid division by zero

        histogram = Array.fill(numBins, 0);
        plotData.do({
            |val|
            var bin = ((val - minVal) / range * (numBins - 1)).floor.clip(0, numBins - 1);
            histogram[bin] = histogram[bin] + 1;
        });
        maxCount = histogram.maxItem;
        if(maxCount == 0) { maxCount = 1 }; // Avoid division by zero

        barWidth = windowWidth / numBins;

        {
            win = Window(title ? "Histogram", Rect(700, 200, windowWidth, windowHeight));
            if (win.notNil) {
                userView = UserView(win, Rect(0, 0, windowWidth, windowHeight)).drawFunc_({
                    Pen.use {
                        histogram.do({
                            |count, i|
                            var barHeight = (count / maxCount) * (windowHeight - 40);
                            var x = i * barWidth;
                            Pen.fillColor = Color.gray(i / (numBins - 1));
                            Pen.fillRect(Rect(x, windowHeight - barHeight - 20, barWidth - 1, barHeight));
                        });
                        Pen.strokeColor = Color.black;
                        Pen.font = Font.monospace(10);
                        Pen.stringAtPoint("Min: %".format(minVal.round(1)), Point(10, 5));
                        Pen.stringAtPoint("Max: %".format(maxVal.round(1)), Point(windowWidth - 80, 5));
                        Pen.stringAtPoint("Count: %".format(maxCount), Point(windowWidth / 2 - 30, 5));
                    };
                });
                win.front;
            } {
                "Error: Window could not be created for histogram visualization.".postln;
            };
        }.defer;
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
    plotMultipleRadialProfiles { |profileDicts, titles, colors=nil|
        var plotter, plots = Array.newClear(profileDicts.size);
        var allRadii = Array.new, allValues = Array.new;
        var globalMinY, globalMaxY, globalMinX, globalMaxX;
        var yRange, paddedMinY, paddedMaxY; 
        var currentPlot; // Moved declaration to the top

        // Collect all data and find global min/max for scaling
        profileDicts.do { |profileDict, i|
            var radii = profileDict.keys.asArray.sort;
            var values = radii.collect { |r| profileDict[r] };

            if (radii.notEmpty) {
                allRadii = allRadii.addAll(radii);
                allValues = allValues.addAll(values);
                // Create a Plot object for each profile and set its properties individually
                currentPlot = Plot.new(nil, nil, Skin.default); // Pass a default skin
                currentPlot.plotMode = \linear;
                currentPlot.resolution = 1;
                currentPlot.value = radii.collect({|r, j| [r, values[j]] });
                currentPlot.domain = [radii.minItem, radii.maxItem];
                currentPlot.range = [values.minItem, values.maxItem];
                plots.put(i, currentPlot);
            } {
                plots.put(i, Plot.new.value_([])); // Add empty plot if no data
            };
        };

        if (allRadii.isEmpty) {
            "No data to plot for multiple radial profiles.".postln;
            ^this;
        };

        globalMinX = allRadii.minItem;
        globalMaxX = allRadii.maxItem;
        globalMinY = allValues.minItem;
        globalMaxY = allValues.maxItem;

        // Add some padding to the y-range
        yRange = globalMaxY - globalMinY;
        if (yRange == 0) { yRange = globalMaxY.abs * 0.1 + 0.01 }; // Avoid zero range
        paddedMinY = globalMinY - (yRange * 0.1);
        paddedMaxY = globalMaxY + (yRange * 0.1);

        // Defer GUI creation
        {
            // Close previous window if it exists
            if(radialProfileWindow.notNil) {
                if(radialProfileWindow.isClosed.not) {
                    radialProfileWindow.close;
                };
            };

            plotter = Plotter("Multiple Radial Profiles", Rect(200, 200, 600, 400));
            plotter.plotMode = \linear; // Still use linear for each profile

            // Assign the array of Plot objects to plotter.plots
            plotter.plots = plots;

            if (colors.notNil and: { colors.size == profileDicts.size }) {
                plotter.plotColor = colors;
            } {
                // Default colors if not provided or size mismatch
                plotter.plotColor = Array.fill(profileDicts.size, { Color.rand });
            };

            // Add legends for each profile
            plotter.legend = titles;

            radialProfileWindow = plotter.parent; // Store the window reference
            radialProfileWindow.front;
        }.defer;

        ^this;
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

        // Initialize displayMin and displayMax to the full range of the data
        displayMin = minVal;
        displayMax = maxVal;

        {
            win = Window("FITS Image - Interactive (ROI Select)", Rect(400, 400, windowWidth, windowHeight + 120)); // Increased height for sliders
            infoText = StaticText(win, Rect(10, windowHeight + 5, windowWidth - 20, 20)).string_("Click and drag to select ROI").font_(Font.monospace(12));

            // Min Value Slider
            minLabel = StaticText(win, Rect(10, windowHeight + 30, 80, 20)).string_("Min Value:");
            minSlider = Slider(win, Rect(90, windowHeight + 30, windowWidth - 100, 20))
                .minVal_(minVal).maxVal_(maxVal).value_(displayMin)
                .action_({ |sl|
                    displayMin = sl.value;
                    if (displayMin > displayMax) { displayMax = displayMin; maxSlider.value = displayMax; };
                    userView.refresh;
                });

            // Max Value Slider
            maxLabel = StaticText(win, Rect(10, windowHeight + 60, 80, 20)).string_("Max Value:");
            maxSlider = Slider(win, Rect(90, windowHeight + 60, windowWidth - 100, 20))
                .minVal_(minVal).maxVal_(maxVal).value_(displayMax)
                .action_({ |sl|
                    displayMax = sl.value;
                    if (displayMax < displayMin) { displayMin = displayMax; minSlider.value = displayMin; };
                    userView.refresh;
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

    visualizeHistogram { |numBins=50|
        this.plotHistogram(fitsFile.data, "FITS Image - Histogram", numBins);
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