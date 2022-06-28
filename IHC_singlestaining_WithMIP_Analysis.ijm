input = getDirectory("Input directory"); //select the folder including nd2 images of the data to be analyzed
list = getFileList(input);
ll=list.length;

  Dialog.create("Staining Details");
  Dialog.addCheckbox("1", 0);  
  Dialog.addCheckbox("2", 1); 
  Dialog.addCheckbox("3", 0); 
  Dialog.addCheckbox("4", 0);  
  Dialog.addString("Name", "protein"); 
  beep();  
  Dialog.show();  
  for (i = 0; i < 4; i++) {
  if(Dialog.getCheckbox())
  	channel=i+1;
  }
  abody1=Dialog.getString();
  
for (i = 0; i < ll; i++) {
	open(list[i]);
}
//channel=2;
//abody1="ppERK";
//ll=46;
for (im=1; im<ll+1; im++){
	
input=getInfo("image.Directory");
title=getTitle;
extract=".nd2";
subtitle=substring(title,0,indexOf(title,extract));
folder=input+subtitle+"/";
File.makeDirectory(folder);
folder1=folder+"/"+abody1+"/";
File.makeDirectory(folder);
File.makeDirectory(folder1);

run("Duplicate...", "title=NewStack.tif duplicate"); 
Stack.getDimensions(width, height, channels, slices, frames); 
getPixelSize(unit, pixelWidth, pixelHeight);
run("Despeckle");

run("Duplicate...", "duplicate channels=1");
run("32-bit");
run("ROF Denoise", "theta=25");
run("16-bit");
for (k=1;k<=slices;k++){
	setSlice(k);
run("Subtract Background...", "sliding=20");
run("Enhance Local Contrast (CLAHE)", "blocksize=31 histogram=512 maximum=3 mask=*None*");
}
setOption("BlackBackground", true);
run("Make Binary", "method=IJ_IsoData background=Dark calculate black");
run("16-bit");
a=4095/255;
run("Multiply...", "value=a stack");
saveAs("Tiff", folder1+"Mask.tif");

run("Duplicate...", "duplicate");
run("Invert", "stack");
saveAs("Tiff", folder1+"Cyto_Mask.tif");

selectWindow("NewStack.tif"); 
selectWindow("NewStack.tif"); 
selectWindow("NewStack.tif"); 
run("Duplicate...", "duplicate channels="+toString(channel));
run("Remove Outliers...", "radius=2 threshold=50 which=Bright stack");
saveAs("Tiff", folder1+abody1+".tif");
imageCalculator("Min create stack", abody1+".tif","Mask.tif");
selectWindow("Result of "+abody1+".tif");
selectWindow("Result of "+abody1+".tif");
selectWindow("Result of "+abody1+".tif");
saveAs("Tiff", folder1+abody1+"_Masked.tif");

imageCalculator("Min create stack", abody1+".tif","Cyto_Mask.tif");
selectWindow("Result of "+abody1+".tif");
selectWindow("Result of "+abody1+".tif");
selectWindow("Result of "+abody1+".tif");
saveAs("Tiff", folder1+abody1+"_Cyto_Masked.tif");
			
selectWindow("NewStack.tif"); 
selectWindow("NewStack.tif"); 
selectWindow("NewStack.tif"); 
close();

//a0=0;
side=0;
selectWindow(abody1+".tif");
selectWindow(abody1+".tif");
selectWindow(abody1+".tif");
resetMinAndMax;
  Dialog.create("How is that sample mounted?");
  Dialog.addChoice("Mounting:", newArray("Flat","Lateral","Skip"));
  beep();  
  Dialog.show();
  mount = Dialog.getChoice();
  if (mount=="Flat") {
  sides=2;
  }
  if (mount=="Lateral") {
  sides=1;
  }
  if (mount=="Skip") {
  sides=0;
  }
  while (side<sides) {
for (j=2; j<slices; j++) { 
	//a=floor((j-3)/3)+1;	
	selectWindow(abody1+"_Masked.tif");
	selectWindow(abody1+"_Masked.tif");
	selectWindow(abody1+"_Masked.tif");
	resetMinAndMax;
	if (j==2){
		setTool("polyline"); 
		beep();  
				waitForUser("Draw Line To Analyze, then click OK.");
				run("Properties... ", "  width="+toString(Math.floor(61/pixelWidth)));
		roiManager("Add");
		roiManager("save", folder1+"LOI"+toString(side)+".zip");
		}
		
	//if (a>a0) {
		//str=a*3-1;
		//stp=a*3+1;
		str=j-1;
		stp=j+1;
		
		//if (stp<slices+1) {
			selectWindow("Mask.tif"); 
			selectWindow("Mask.tif"); 
			selectWindow("Mask.tif"); 
			run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			selectWindow("MAX_Mask.tif"); 
			selectWindow("MAX_Mask.tif"); 
			selectWindow("MAX_Mask.tif"); 
			roiManager("Select", 0);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_Mask_"+toString(j)+toString(side)+".csv");
			//saveAs("Measurements", folder1+"Values_"+abody1+"_"+toString(j-1)+toString(side)+".csv");
			//saveAs("Measurements", folder1+"Values_"+abody1+"_"+toString(j-2)+toString(side)+".csv");
			selectWindow("MAX_Mask.tif"); 
			close();

			selectWindow("Cyto_Mask.tif"); 
			selectWindow("Cyto_Mask.tif"); 
			selectWindow("Cyto_Mask.tif"); 
			run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			selectWindow("MAX_Cyto_Mask.tif"); 
			selectWindow("MAX_Cyto_Mask.tif"); 
			selectWindow("MAX_Cyto_Mask.tif"); 
			roiManager("Select", 0);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_Cyto_Mask_"+toString(j)+toString(side)+".csv");			
			//saveAs("Measurements", folder1+"Values_Cyto_Mask_"+abody1+"_"+toString(j-1)+toString(side)+".csv");
			//saveAs("Measurements", folder1+"Values_Cyto_Mask_"+abody1+"_"+toString(j-2)+toString(side)+".csv");
			selectWindow("MAX_Cyto_Mask.tif"); 
			close();
						
			selectWindow(abody1+"_Masked.tif");	
			selectWindow(abody1+"_Masked.tif");	
			selectWindow(abody1+"_Masked.tif");	
			run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			selectWindow("MAX_"+abody1+"_Masked.tif"); 
			selectWindow("MAX_"+abody1+"_Masked.tif"); 
			selectWindow("MAX_"+abody1+"_Masked.tif"); 
			roiManager("Select", 0);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_"+abody1+"_"+toString(j)+toString(side)+".csv");
			//saveAs("Measurements", folder1+"Values_"+abody1+"_"+toString(j-1)+toString(side)+".csv");
			//saveAs("Measurements", folder1+"Values_"+abody1+"_"+toString(j-2)+toString(side)+".csv");
			selectWindow("MAX_"+abody1+"_Masked.tif"); 
			close();

			selectWindow(abody1+"_Cyto_Masked.tif");	
			selectWindow(abody1+"_Cyto_Masked.tif");
			selectWindow(abody1+"_Cyto_Masked.tif");
			run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			selectWindow("MAX_"+abody1+"_Cyto_Masked.tif"); 
			selectWindow("MAX_"+abody1+"_Cyto_Masked.tif"); 
			selectWindow("MAX_"+abody1+"_Cyto_Masked.tif"); 			
			roiManager("Select", 0);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_"+abody1+"_Cyto_"+toString(j)+toString(side)+".csv");
			//saveAs("Measurements", folder1+"Values_"+abody1+"_Cyto_"+toString(j-1)+toString(side)+".csv");
			//saveAs("Measurements", folder1+"Values_"+abody1+"_Cyto_"+toString(j-2)+toString(side)+".csv");			
			selectWindow("MAX_"+abody1+"_Cyto_Masked.tif"); 
			close();
			

			selectWindow(abody1+".tif");
			selectWindow(abody1+".tif");
			selectWindow(abody1+".tif");	
			run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			selectWindow("MAX_"+abody1+".tif");  
			selectWindow("MAX_"+abody1+".tif"); 
			selectWindow("MAX_"+abody1+".tif"); 
			roiManager("Select", 0);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_"+abody1+"unmasked_"+toString(j)+toString(side)+".csv");
			//saveAs("Measurements", folder1+"Values_"+abody1+"unmasked_"+toString(j-1)+toString(side)+".csv");
			//saveAs("Measurements", folder1+"Values_"+abody1+"unmasked_"+toString(j-2)+toString(side)+".csv");
			selectWindow("MAX_"+abody1+".tif"); 
			close();
				
			//}
		//}
	//a0=floor((j-3)/3)+1;	
	}
	roiManager("Delete");	
	side=side+1;
}
	selectWindow(title);
	close();
	selectWindow("Mask.tif");
	selectWindow("Mask.tif");
	selectWindow("Mask.tif");
	close();
	selectWindow("Cyto_Mask.tif");
	selectWindow("Cyto_Mask.tif");
	selectWindow("Cyto_Mask.tif");
	close();
	selectWindow(abody1+".tif");
	selectWindow(abody1+".tif");
	selectWindow(abody1+".tif");
	close();
	selectWindow(abody1+"_Masked.tif");
	selectWindow(abody1+"_Masked.tif");
	selectWindow(abody1+"_Masked.tif");
	close();
	selectWindow(abody1+"_Cyto_Masked.tif");
	selectWindow(abody1+"_Cyto_Masked.tif");
	selectWindow(abody1+"_Cyto_Masked.tif");
	close();

	
wlist = getList("window.titles"); 
for (i=0; i<wlist.length; i++){      
	window = wlist[i]; 
	selectWindow(window); 
	selectWindow(window); 
	selectWindow(window); 
	run("Close"); 
	} 
	} 
