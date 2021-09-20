path1="/Users/masahiro/Desktop/Experiments/IF/dapi/210822/image_analysis/ESC/";
for(i=1; i<19;i++){
	fname=""+i+"_SR_high.czi";
	run("Bio-Formats Importer", "open="+path1+fname+" autoscale color_mode=Default rois_import=[ROI manager] split_timepoints split_channels view=Hyperstack stack_order=XYCZT");
	title=getTitle();
	run("Duplicate...",  "title=image duplicate");
	saveAs("Tiff", path1+i+".tif");
	close();
	selectWindow(title);
	run("Gaussian Blur...", "sigma=2 stack");
	setAutoThreshold("Default dark");
	//run("Threshold...");
	setAutoThreshold("Otsu dark");
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Analyze Particles...", "size=2000-Infinity pixel show=Masks display exclude clear include stack");
	run("Invert LUT");
	run("Clear Results");
	run("Duplicate...",  "title=nucleus_mask duplicate");
	selectWindow("nucleus_mask");
	saveAs("Tiff", path1+i+"_nucleus_mask.tif");
	selectWindow("Mask of "+title);
	run("Duplicate...", "duplicate");
	selectWindow("Mask of "+title+"-1");
	run("Duplicate...", "duplicate");
	selectWindow("Mask of "+title+"-2");
	run("Options...", "iterations=5 count=1 black do=Erode stack");
	imageCalculator("Subtract create stack", "Mask of "+title+"-1","Mask of "+title+"-2");
	selectWindow("Result of Mask of "+title+"-1");
	saveAs("Tiff", path1+i+"_nucleus_rim_mask.tif");
	selectWindow("Mask of "+title+"-1");
	run("Analyze Particles...", "size=2000-Infinity pixel show=Outlines display exclude clear include add stack");
	roiManager("Save", path1+i+"_nuclear_periphery_roiset.zip");
	selectWindow("Drawing of Mask of "+title+"-1");
	close();
	selectWindow("Mask of "+title+"-2");
	close();
	selectWindow("Mask of "+title+"-1");
	close();
	selectWindow("Mask of "+title);
	close();
	selectWindow(i+"_nucleus_mask.tif");
	close();
	selectWindow(i+"_nucleus_rim_mask.tif");
	close();
	selectWindow(title);
	close();
	run("Bio-Formats Importer", "open="+path1+fname+" autoscale color_mode=Default rois_import=[ROI manager] split_timepoints split_channels view=Hyperstack stack_order=XYCZT");
	selectWindow(title);
	run("Gaussian Blur...", "sigma=2 stack");
	run("Clear Results");
	run("Auto Threshold", "method=MaxEntropy white stack");
	run("Open", "stack");
	run("Analyze Particles...", "size=50-Infinity pixel show=Masks display exclude clear include add stack");
	saveAs("Results", path1+i+"_DAPI_dense_regions.csv");
	selectWindow("Mask of "+title);
	run("Invert LUT");
	saveAs("Tiff", path1+i+"_DAPI_dense_regions.tif");
	close();
	selectWindow(title);
	close();
	run("Clear Results");
}

path1="/Users/masahiro/Desktop/Experiments/IF/dapi/210822/image_analysis/EpiLC/";
for(i=1; i<24;i++){
	fname=""+i+"_SR_high.czi";
	run("Bio-Formats Importer", "open="+path1+fname+" autoscale color_mode=Default rois_import=[ROI manager] split_timepoints split_channels view=Hyperstack stack_order=XYCZT");
	title=getTitle();
	run("Duplicate...",  "title=image duplicate");
	saveAs("Tiff", path1+i+".tif");
	close();
	selectWindow(title);
	run("Gaussian Blur...", "sigma=2 stack");
	setAutoThreshold("Default dark");
	//run("Threshold...");
	setAutoThreshold("Otsu dark");
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Analyze Particles...", "size=2000-Infinity pixel show=Masks display exclude clear include stack");
	run("Invert LUT");
	run("Clear Results");
	run("Duplicate...",  "title=nucleus_mask duplicate");
	selectWindow("nucleus_mask");
	saveAs("Tiff", path1+i+"_nucleus_mask.tif");
	selectWindow("Mask of "+title);
	run("Duplicate...", "duplicate");
	selectWindow("Mask of "+title+"-1");
	run("Duplicate...", "duplicate");
	selectWindow("Mask of "+title+"-2");
	run("Options...", "iterations=5 count=1 black do=Erode stack");
	imageCalculator("Subtract create stack", "Mask of "+title+"-1","Mask of "+title+"-2");
	selectWindow("Result of Mask of "+title+"-1");
	saveAs("Tiff", path1+i+"_nucleus_rim_mask.tif");
	selectWindow("Mask of "+title+"-1");
	run("Analyze Particles...", "size=2000-Infinity pixel show=Outlines display exclude clear include add stack");
	roiManager("Save", path1+i+"_nuclear_periphery_roiset.zip");
	selectWindow("Drawing of Mask of "+title+"-1");
	close();
	selectWindow("Mask of "+title+"-2");
	close();
	selectWindow("Mask of "+title+"-1");
	close();
	selectWindow("Mask of "+title);
	close();
	selectWindow(i+"_nucleus_mask.tif");
	close();
	selectWindow(i+"_nucleus_rim_mask.tif");
	close();
	selectWindow(title);
	close();
	run("Bio-Formats Importer", "open="+path1+fname+" autoscale color_mode=Default rois_import=[ROI manager] split_timepoints split_channels view=Hyperstack stack_order=XYCZT");
	selectWindow(title);
	run("Gaussian Blur...", "sigma=2 stack");
	run("Clear Results");
	run("Auto Threshold", "method=MaxEntropy white stack");
	run("Open", "stack");
	run("Analyze Particles...", "size=50-Infinity pixel show=Masks display exclude clear include add stack");
	saveAs("Results", path1+i+"_DAPI_dense_regions.csv");
	selectWindow("Mask of "+title);
	run("Invert LUT");
	saveAs("Tiff", path1+i+"_DAPI_dense_regions.tif");
	close();
	selectWindow(title);
	close();
	run("Clear Results");
}


path1="/Users/masahiro/Desktop/Experiments/IF/dapi/210822/image_analysis/d2PGCLC/";
for(i=1; i<20;i++){
	fname=""+i+"_SR_high.czi";
	run("Bio-Formats Importer", "open="+path1+fname+" autoscale color_mode=Default rois_import=[ROI manager] split_timepoints split_channels view=Hyperstack stack_order=XYCZT");
	title=getTitle();
	run("Duplicate...",  "title=image duplicate");
	saveAs("Tiff", path1+i+".tif");
	close();
	selectWindow(title);
	run("Gaussian Blur...", "sigma=2 stack");
	setAutoThreshold("Default dark");
	//run("Threshold...");
	setAutoThreshold("Otsu dark");
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Analyze Particles...", "size=2000-Infinity pixel show=Masks display exclude clear include stack");
	run("Invert LUT");
	run("Clear Results");
	run("Duplicate...",  "title=nucleus_mask duplicate");
	selectWindow("nucleus_mask");
	saveAs("Tiff", path1+i+"_nucleus_mask.tif");
	selectWindow("Mask of "+title);
	run("Duplicate...", "duplicate");
	selectWindow("Mask of "+title+"-1");
	run("Duplicate...", "duplicate");
	selectWindow("Mask of "+title+"-2");
	run("Options...", "iterations=5 count=1 black do=Erode stack");
	imageCalculator("Subtract create stack", "Mask of "+title+"-1","Mask of "+title+"-2");
	selectWindow("Result of Mask of "+title+"-1");
	saveAs("Tiff", path1+i+"_nucleus_rim_mask.tif");
	selectWindow("Mask of "+title+"-1");
	run("Analyze Particles...", "size=2000-Infinity pixel show=Outlines display exclude clear include add stack");
	roiManager("Save", path1+i+"_nuclear_periphery_roiset.zip");
	selectWindow("Drawing of Mask of "+title+"-1");
	close();
	selectWindow("Mask of "+title+"-2");
	close();
	selectWindow("Mask of "+title+"-1");
	close();
	selectWindow("Mask of "+title);
	close();
	selectWindow(i+"_nucleus_mask.tif");
	close();
	selectWindow(i+"_nucleus_rim_mask.tif");
	close();
	selectWindow(title);
	close();
	run("Bio-Formats Importer", "open="+path1+fname+" autoscale color_mode=Default rois_import=[ROI manager] split_timepoints split_channels view=Hyperstack stack_order=XYCZT");
	selectWindow(title);
	run("Gaussian Blur...", "sigma=2 stack");
	run("Clear Results");
	run("Auto Threshold", "method=MaxEntropy white stack");
	run("Open", "stack");
	run("Analyze Particles...", "size=50-Infinity pixel show=Masks display exclude clear include add stack");
	saveAs("Results", path1+i+"_DAPI_dense_regions.csv");
	selectWindow("Mask of "+title);
	run("Invert LUT");
	saveAs("Tiff", path1+i+"_DAPI_dense_regions.tif");
	close();
	selectWindow(title);
	close();
	run("Clear Results");
}

path1="/Users/masahiro/Desktop/Experiments/IF/dapi/210822/image_analysis/d4c7PGCLC/";
for(i=1; i<28;i++){
	fname=""+i+"_SR_high.czi";
	run("Bio-Formats Importer", "open="+path1+fname+" autoscale color_mode=Default rois_import=[ROI manager] split_timepoints split_channels view=Hyperstack stack_order=XYCZT");
	title=getTitle();
	run("Duplicate...",  "title=image duplicate");
	saveAs("Tiff", path1+i+".tif");
	close();
	selectWindow(title);
	run("Gaussian Blur...", "sigma=2 stack");
	setAutoThreshold("Default dark");
	//run("Threshold...");
	setAutoThreshold("Otsu dark");
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Analyze Particles...", "size=2000-Infinity pixel show=Masks display exclude clear include stack");
	run("Invert LUT");
	run("Clear Results");
	run("Duplicate...",  "title=nucleus_mask duplicate");
	selectWindow("nucleus_mask");
	saveAs("Tiff", path1+i+"_nucleus_mask.tif");
	selectWindow("Mask of "+title);
	run("Duplicate...", "duplicate");
	selectWindow("Mask of "+title+"-1");
	run("Duplicate...", "duplicate");
	selectWindow("Mask of "+title+"-2");
	run("Options...", "iterations=5 count=1 black do=Erode stack");
	imageCalculator("Subtract create stack", "Mask of "+title+"-1","Mask of "+title+"-2");
	selectWindow("Result of Mask of "+title+"-1");
	saveAs("Tiff", path1+i+"_nucleus_rim_mask.tif");
	selectWindow("Mask of "+title+"-1");
	run("Analyze Particles...", "size=2000-Infinity pixel show=Outlines display exclude clear include add stack");
	roiManager("Save", path1+i+"_nuclear_periphery_roiset.zip");
	selectWindow("Drawing of Mask of "+title+"-1");
	close();
	selectWindow("Mask of "+title+"-2");
	close();
	selectWindow("Mask of "+title+"-1");
	close();
	selectWindow("Mask of "+title);
	close();
	selectWindow(i+"_nucleus_mask.tif");
	close();
	selectWindow(i+"_nucleus_rim_mask.tif");
	close();
	selectWindow(title);
	close();
	run("Bio-Formats Importer", "open="+path1+fname+" autoscale color_mode=Default rois_import=[ROI manager] split_timepoints split_channels view=Hyperstack stack_order=XYCZT");
	selectWindow(title);
	run("Gaussian Blur...", "sigma=2 stack");
	run("Clear Results");
	run("Auto Threshold", "method=MaxEntropy white stack");
	run("Open", "stack");
	run("Analyze Particles...", "size=50-Infinity pixel show=Masks display exclude clear include add stack");
	saveAs("Results", path1+i+"_DAPI_dense_regions.csv");
	selectWindow("Mask of "+title);
	run("Invert LUT");
	saveAs("Tiff", path1+i+"_DAPI_dense_regions.tif");
	close();
	selectWindow(title);
	close();
	run("Clear Results");
}

path1="/Users/masahiro/Desktop/Experiments/IF/dapi/210822/image_analysis/GSC/";
for(i=1; i<23;i++){
	fname=""+i+"_SR_high.czi";
	run("Bio-Formats Importer", "open="+path1+fname+" autoscale color_mode=Default rois_import=[ROI manager] split_timepoints split_channels view=Hyperstack stack_order=XYCZT");
	title=getTitle();
	run("Duplicate...",  "title=image duplicate");
	saveAs("Tiff", path1+i+".tif");
	close();
	selectWindow(title);
	run("Gaussian Blur...", "sigma=2 stack");
	setAutoThreshold("Default dark");
	//run("Threshold...");
	setAutoThreshold("Otsu dark");
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Analyze Particles...", "size=2000-Infinity pixel show=Masks display exclude clear include stack");
	run("Invert LUT");
	run("Clear Results");
	run("Duplicate...",  "title=nucleus_mask duplicate");
	selectWindow("nucleus_mask");
	saveAs("Tiff", path1+i+"_nucleus_mask.tif");
	selectWindow("Mask of "+title);
	run("Duplicate...", "duplicate");
	selectWindow("Mask of "+title+"-1");
	run("Duplicate...", "duplicate");
	selectWindow("Mask of "+title+"-2");
	run("Options...", "iterations=5 count=1 black do=Erode stack");
	imageCalculator("Subtract create stack", "Mask of "+title+"-1","Mask of "+title+"-2");
	selectWindow("Result of Mask of "+title+"-1");
	saveAs("Tiff", path1+i+"_nucleus_rim_mask.tif");
	selectWindow("Mask of "+title+"-1");
	run("Analyze Particles...", "size=2000-Infinity pixel show=Outlines display exclude clear include add stack");
	roiManager("Save", path1+i+"_nuclear_periphery_roiset.zip");
	selectWindow("Drawing of Mask of "+title+"-1");
	close();
	selectWindow("Mask of "+title+"-2");
	close();
	selectWindow("Mask of "+title+"-1");
	close();
	selectWindow("Mask of "+title);
	close();
	selectWindow(i+"_nucleus_mask.tif");
	close();
	selectWindow(i+"_nucleus_rim_mask.tif");
	close();
	selectWindow(title);
	close();
	run("Bio-Formats Importer", "open="+path1+fname+" autoscale color_mode=Default rois_import=[ROI manager] split_timepoints split_channels view=Hyperstack stack_order=XYCZT");
	selectWindow(title);
	run("Gaussian Blur...", "sigma=2 stack");
	run("Clear Results");
	run("Auto Threshold", "method=MaxEntropy white stack");
	run("Open", "stack");
	run("Analyze Particles...", "size=50-Infinity pixel show=Masks display exclude clear include add stack");
	saveAs("Results", path1+i+"_DAPI_dense_regions.csv");
	selectWindow("Mask of "+title);
	run("Invert LUT");
	saveAs("Tiff", path1+i+"_DAPI_dense_regions.tif");
	close();
	selectWindow(title);
	close();
	run("Clear Results");
}


path1="/Users/masahiro/Desktop/Experiments/IF/dapi/210822/image_analysis/GSCLC/";
for(i=1; i<25;i++){
	fname=""+i+"_SR_high.czi";
	run("Bio-Formats Importer", "open="+path1+fname+" autoscale color_mode=Default rois_import=[ROI manager] split_timepoints split_channels view=Hyperstack stack_order=XYCZT");
	title=getTitle();
	run("Duplicate...",  "title=image duplicate");
	saveAs("Tiff", path1+i+".tif");
	close();
	selectWindow(title);
	run("Gaussian Blur...", "sigma=2 stack");
	setAutoThreshold("Default dark");
	//run("Threshold...");
	setAutoThreshold("Otsu dark");
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Analyze Particles...", "size=2000-Infinity pixel show=Masks display exclude clear include stack");
	run("Invert LUT");
	run("Clear Results");
	run("Duplicate...",  "title=nucleus_mask duplicate");
	selectWindow("nucleus_mask");
	saveAs("Tiff", path1+i+"_nucleus_mask.tif");
	selectWindow("Mask of "+title);
	run("Duplicate...", "duplicate");
	selectWindow("Mask of "+title+"-1");
	run("Duplicate...", "duplicate");
	selectWindow("Mask of "+title+"-2");
	run("Options...", "iterations=5 count=1 black do=Erode stack");
	imageCalculator("Subtract create stack", "Mask of "+title+"-1","Mask of "+title+"-2");
	selectWindow("Result of Mask of "+title+"-1");
	saveAs("Tiff", path1+i+"_nucleus_rim_mask.tif");
	selectWindow("Mask of "+title+"-1");
	run("Analyze Particles...", "size=2000-Infinity pixel show=Outlines display exclude clear include add stack");
	roiManager("Save", path1+i+"_nuclear_periphery_roiset.zip");
	selectWindow("Drawing of Mask of "+title+"-1");
	close();
	selectWindow("Mask of "+title+"-2");
	close();
	selectWindow("Mask of "+title+"-1");
	close();
	selectWindow("Mask of "+title);
	close();
	selectWindow(i+"_nucleus_mask.tif");
	close();
	selectWindow(i+"_nucleus_rim_mask.tif");
	close();
	selectWindow(title);
	close();
	run("Bio-Formats Importer", "open="+path1+fname+" autoscale color_mode=Default rois_import=[ROI manager] split_timepoints split_channels view=Hyperstack stack_order=XYCZT");
	selectWindow(title);
	run("Gaussian Blur...", "sigma=2 stack");
	run("Clear Results");
	run("Auto Threshold", "method=MaxEntropy white stack");
	run("Open", "stack");
	run("Analyze Particles...", "size=50-Infinity pixel show=Masks display exclude clear include add stack");
	saveAs("Results", path1+i+"_DAPI_dense_regions.csv");
	selectWindow("Mask of "+title);
	run("Invert LUT");
	saveAs("Tiff", path1+i+"_DAPI_dense_regions.tif");
	close();
	selectWindow(title);
	close();
	run("Clear Results");
}
