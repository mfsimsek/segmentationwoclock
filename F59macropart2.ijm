
//2nd Part after saving 5 ROIs
//run("Invert");
//run("Erode");

title=getTitle;
subtitle=substring(title,5,lastIndexOf(title,"."));
selectWindow("Mask_"+subtitle+".tif");
run("Open Next");
titlenext=getTitle;
input=getInfo("image.Directory");
open(input+title);
run("Dilate");

roiManager("Deselect");
roiManager("Save", input+"roiset"+subtitle+".zip");
selectWindow(title);
run("Set Measurements...", "area mean standard redirect=None decimal=3");
//roiManager("Select", newArray(0,1,2,3,4,5));
roiManager("Measure");
saveAs("Results", input+"sample"+subtitle+".csv");
run("Close All");
wlist = getList("window.titles"); 
for (i=0; i<wlist.length; i++){      
	window = wlist[i]; 
	selectWindow(window); 
	run("Close"); 
	} 

open(input+titlenext);
//run("Change Tracing Image", "choice="+titlenext+" validatecalibration=false");
//run("Simple Neurite Tracer", " ");
//runMacro("C:/Users/sim4yx/OneDrive - cchmc/FIJI_Macros/xirpmacropart1.ijm");