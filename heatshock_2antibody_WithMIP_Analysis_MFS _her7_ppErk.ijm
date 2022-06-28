input = getDirectory("Input directory"); //select the folder including nd2 images of the data to be analyzed
list = getFileList(input);
ll=list.length;
for (i = 0; i < ll; i++) {
	open(list[i]);
}
for (im=1; im<ll+1; im++){
input=getInfo("image.Directory");
title=getTitle;
extract=".nd2";
subtitle=substring(title,0,indexOf(title,extract));
i1=indexOf(title,"_")+1;
iend=lastIndexOf(title,"_");
abodyall=substring(title,i1,iend);
ibody=indexOf(abodyall,"_");
abody1=substring(abodyall,0,ibody);
abody2=substring(abodyall,ibody+1);
folder=input+subtitle+"/";
File.makeDirectory(folder);
folder1=folder+"/"+abody1+"/";
folder2=folder+"/"+abody2+"/";
File.makeDirectory(folder);
File.makeDirectory(folder1);
File.makeDirectory(folder2);

run("Duplicate...", "title=NewStack.tif duplicate"); 
Stack.getDimensions(width, height, channels, slices, frames); 
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
saveAs("Tiff", folder2+"Mask.tif");

Mask=Mask*4095/255;
run("Duplicate...", "duplicate");
//run("32-bit");
//run("Subtract...", "value=255 stack");
//run("Abs", "stack");
//run("8-bit");
run("Invert", "stack");
saveAs("Tiff", folder1+"Cyto_Mask.tif");
saveAs("Tiff", folder2+"Cyto_Mask.tif");

selectWindow("NewStack.tif"); 
run("Duplicate...", "duplicate channels=2");
run("Remove Outliers...", "radius=2 threshold=50 which=Bright stack");
saveAs("Tiff", folder1+abody1+".tif");
imageCalculator("Min create stack", abody1+".tif","Mask.tif");
selectWindow("Result of "+abody1+".tif");
saveAs("Tiff", folder1+abody1+"_Masked.tif");

imageCalculator("Min create stack", abody1+".tif","Cyto_Mask.tif");
selectWindow("Result of "+abody1+".tif");
saveAs("Tiff", folder1+abody1+"_Cyto_Masked.tif");

selectWindow("NewStack.tif"); 
run("Duplicate...", "duplicate channels=3");
saveAs("Tiff", folder2+abody2+".tif");
imageCalculator("Min create stack", abody2+".tif","Mask.tif");
selectWindow("Result of "+abody2+".tif");
saveAs("Tiff", folder2+abody2+"_Masked.tif");

imageCalculator("Min create stack", abody2+".tif","Cyto_Mask.tif");
selectWindow("Result of "+abody2+".tif");
saveAs("Tiff", folder2+abody2+"_Cyto_Masked.tif");


//selectWindow("NewStack.tif"); 
//run("Duplicate...", "duplicate channels=2");
//saveAs("Tiff", folder1+"GFP.tif");
//saveAs("Tiff", folder2+"GFP.tif");
//imageCalculator("Min create stack", "GFP.tif","Mask.tif");
//selectWindow("Result of GFP.tif");
//saveAs("Tiff", folder1+"GFP_Masked.tif");
//saveAs("Tiff", folder2+"GFP_Masked.tif");


			
selectWindow("NewStack.tif"); 
close();

a0=0;
side=0;

  Dialog.create("How is that sample mounted?");
  Dialog.addChoice("Mounting:", newArray("Lateral","Flat","Skip"));
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
for (j=1; j<slices+1; j++) { 
	a=floor((j-3)/3)+1;	
	selectWindow(abody1+"_Cyto_Masked.tif");
	if (j==1){
		setTool("polyline"); 
		beep();  
				waitForUser("Select Line To Analyze, then click OK.");
		roiManager("Add");
		}

			selectWindow(abody2+"_Masked.tif");	
			//run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			//selectWindow("MAX_"+abody+"_Masked.tif"); 
			setTool("polyline");  
			roiManager("Select", 0);
			setSlice(j);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder2+"Values_"+abody2+"_"+toString(j)+toString(side)+".csv");
			//selectWindow("MAX_"+abody+"_Masked.tif"); 
			//close();

			selectWindow(abody2+"_Cyto_Masked.tif");	
			//run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			//selectWindow("MAX_"+abody+"_Cyto_Masked.tif"); 
			setTool("polyline");  
			roiManager("Select", 0);
			setSlice(j);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder2+"Values_"+abody2+"_Cyto_"+toString(j)+toString(side)+".csv");
			//selectWindow("MAX_"+abody+"_Cyto_Masked.tif"); 
			//close();
			

			selectWindow(abody2+".tif");	
			//run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			//selectWindow("MAX_"+abody+".tif");  
			roiManager("Select", 0);
			setSlice(j);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder2+"Values_"+abody2+"unmasked_"+toString(j)+toString(side)+".csv");
			//selectWindow("MAX_"+abody+".tif"); 
			//close();

			
			selectWindow("Mask.tif"); 
			//run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			//selectWindow("MAX_Mask.tif"); 
			roiManager("Select", 0);
			setSlice(j);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_Mask_"+toString(j)+toString(side)+".csv");
			saveAs("Measurements", folder2+"Values_Mask_"+toString(j)+toString(side)+".csv");
			//selectWindow("MAX_Mask.tif"); 
			//close();

			selectWindow("Cyto_Mask.tif"); 
			//run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			//selectWindow("MAX_Cyto_Mask.tif"); 
			roiManager("Select", 0);
			setSlice(j);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_Cyto_Mask_"+toString(j)+toString(side)+".csv");
			saveAs("Measurements", folder2+"Values_Cyto_Mask_"+toString(j)+toString(side)+".csv");
			//selectWindow("MAX_Cyto_Mask.tif"); 
			//close();

		
	if (a>a0) {
		str=a*3-1;
		stp=a*3+1;
		if (stp<slices+1) {
			selectWindow(abody1+"_Masked.tif");	
			run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			selectWindow("MAX_"+abody1+"_Masked.tif"); 
			setTool("polyline");  
			roiManager("Select", 0);
			//setSlice(j);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_"+abody1+"_"+toString(j)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_"+abody1+"_"+toString(j-1)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_"+abody1+"_"+toString(j-2)+toString(side)+".csv");
			selectWindow("MAX_"+abody1+"_Masked.tif"); 
			close();

			selectWindow(abody1+"_Cyto_Masked.tif");	
			run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			selectWindow("MAX_"+abody1+"_Cyto_Masked.tif"); 
			setTool("polyline");  
			roiManager("Select", 0);
			//setSlice(j);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_"+abody1+"_Cyto_"+toString(j)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_"+abody1+"_Cyto_"+toString(j-1)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_"+abody1+"_Cyto_"+toString(j-2)+toString(side)+".csv");			
			selectWindow("MAX_"+abody1+"_Cyto_Masked.tif"); 
			close();
			

			selectWindow(abody1+".tif");	
			run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			selectWindow("MAX_"+abody1+".tif");  
			roiManager("Select", 0);
			//setSlice(j);
			run("Clear Results");
			profile = getProfile();
			for (i=0; i<profile.length; i++)
				setResult("Value", i, profile[i]);
			updateResults();
			saveAs("Measurements", folder1+"Values_"+abody1+"unmasked_"+toString(j)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_"+abody1+"unmasked_"+toString(j-1)+toString(side)+".csv");
			saveAs("Measurements", folder1+"Values_"+abody1+"unmasked_"+toString(j-2)+toString(side)+".csv");
			selectWindow("MAX_"+abody1+".tif"); 
			close();

			
			
			
			//selectWindow("GFP_Masked.tif"); 
			//run("Z Project...", "start=&str stop=&stp projection=[Max Intensity]");
			//selectWindow("MAX_GFP_Masked.tif"); 
			//roiManager("Select", 0);
			//setSlice(j);
			//run("Clear Results");
			//profile = getProfile();
			//for (i=0; i<profile.length; i++)
			//	setResult("Value", i, profile[i]);
			//updateResults();
			//saveAs("Measurements", folder1+"Values_GFP_"+toString(j)+toString(side)+".csv");
			//saveAs("Measurements", folder2+"Values_GFP_"+toString(j)+toString(side)+".csv");
			//selectWindow("MAX_GFP_Masked.tif"); 
			//close();			
			}
		}
	a0=floor((j-3)/3)+1;
	
	if(j==slices){
		roiManager("Delete");
		
	}
	}
	side=side+1;
}
	selectWindow(title);
	close();
	selectWindow("Mask.tif");
	close();
	selectWindow("Cyto_Mask.tif");
	close();
	selectWindow(abody1+".tif");
	close();
	selectWindow(abody1+"_Masked.tif");
	close();
	selectWindow(abody1+"_Cyto_Masked.tif");
	close();
	selectWindow(abody2+".tif");
	close();
	selectWindow(abody2+"_Masked.tif");
	close();
	selectWindow(abody2+"_Cyto_Masked.tif");


			saveAs("Measurements", folder1+"Values_"+abody1+"_"+toString(100)+".csv");
			saveAs("Measurements", folder1+"Values_"+abody1+"_"+toString(110)+".csv");
			saveAs("Measurements", folder1+"Values_"+abody1+"_Cyto_"+toString(100)+".csv");
			saveAs("Measurements", folder1+"Values_"+abody1+"_Cyto_"+toString(110)+".csv");
			saveAs("Measurements", folder1+"Values_"+abody1+"_unmasked_"+toString(100)+".csv");
			saveAs("Measurements", folder1+"Values_"+abody1+"_unmasked_"+toString(110)+".csv");
	
	close();
	//selectWindow("GFP.tif");
	//close();
	//selectWindow("GFP_Masked.tif");
	//close();
	
wlist = getList("window.titles"); 
for (i=0; i<wlist.length; i++){      
	window = wlist[i]; 
	selectWindow(window); 
	run("Close"); 
	} 
	} 
