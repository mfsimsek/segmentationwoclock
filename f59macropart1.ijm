input = getDirectory("Input directory"); //select the folder including nd2 images of the data to be analyzed
//~~~ opens all files in the folder ~~~
//for nd2 files FIJI uses "BioFormats Importer" Plugin, which opens a window dialogue for each file on how to open the file (split channels, composite, grayscale, etc.) 
//once you are set with how to open a file, you can eliminate that pop-up window by defaulting options (Plugins>Bio-Formats>Bio-Formats Plugins Configuration>select nd2 and check Windowless)
list = getFileList(input);
ll=list.length;
for (i = 0; i < ll; i++) {
	open(list[i]);
}
input=getInfo("image.Directory");
folder=input+"pFakMasks/"; //creates a folder for pFak mask files
File.makeDirectory(folder);

//~~~ if code gets interrupted in the middle re-run from here ~~~
ll=nImages;
//~~~ goes through each image open ~~~
for (im=1; im<ll+1; im++){
input=getInfo("image.Directory");
folder=input+"pFakMasks/"; 
title=getTitle;
subtitle=substring(title,1,lastIndexOf(title,"."));
input=getInfo("image.Directory");

run("Duplicate...", "title=pFak.tif duplicate stack channels=2");
run("Subtract Background...", "rolling=15 stack");
run("Remove Outliers...", "radius=2 threshold=1500 which=Bright stack");
run("Z Project...", "projection=[Standard Deviation]");
run("8-bit");
run("Enhance Local Contrast (CLAHE)", "blocksize=63 histogram=1024 maximum=3 mask=*None*");

setTool("line");
	beep();  
	
	waitForUser("Draw a horizontal line between the start and end positions, then click OK.");
//	roiManager("Add");
	getLine(x1, y1, x2, y2, lineWidth);
	angle=floor(-180*atan((y2-y1)/(x2-x1))/PI);
run("Rotate... ", "angle="+toString(angle)+" grid=1 interpolation=Bicubic fill enlarge");	
setTool("freehand");
	beep();  
	waitForUser("Select a bounding line for the region of boundaries to be analyzed, then click OK.");

run("Duplicate...", "duplicate");
run("8-bit");
run("Select All");
run("Auto Threshold", "method=MaxEntropy white");
run("Remove Outliers...", "radius=2 threshold=50 which=Bright");
setOption("BlackBackground", true);
//run("Convert to Mask");
run("Restore Selection");
run("Make Inverse");
run("Colors...", "foreground=black background=black selection=yellow");
run("Fill", "slice");
run("Select All");
run("Convert to Mask");
//run("Invert");

saveAs("Tiff", folder+"Mask_"+subtitle+".tif");
selectWindow("STD_pFak.tif");
run("Duplicate...", " ");
imageCalculator("Min create", "STD_pFak-1.tif","Mask_"+subtitle+".tif");
selectWindow("Result of STD_pFak-1.tif");
saveAs("Tiff", folder+"Masked_"+subtitle+".tif");

selectWindow(title);
	close();
selectWindow("pFak.tif");
	close();
selectWindow("STD_pFak.tif");
	close();
selectWindow("STD_pFak-1.tif");
	close();
selectWindow("Mask_"+subtitle+".tif");
	close();
selectWindow("Masked_"+subtitle+".tif");
	close();
}
//run("Simple Neurite Tracer", " ");
