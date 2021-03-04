
//Crop yo Files
fs=File.separator;
//dirStor=getDirectory("Select Source Directory");
dirStor="/home/sdenecke/Applications/Custom_Applications/Wiggle_Index_chrysalida/Wiggle_Analyze/"
dirSav=dirStor+"zCroppedTiffStackLL";
File.makeDirectory(dirSav);
setBatchMode(true);
lis=getFileList(dirStor);
ds=lengthOf(dirStor);
hi=substring(dirStor,(ds-2),(ds-1));
parent=substring(dirStor,0,(ds-3));

tableTitle="Wiggle Index_LowLeft";
tableTitle2="["+tableTitle+"]";
if (isOpen("Wiggle Index_LowLeft")){
	print("Table already open");}
	
	else{
		
		run("Table...", "name="+tableTitle2+" width=600 height=250");
		print(tableTitle2, "\\Headings:Image Name\tWiggle Index\tImage Scale\tImage Blur\tFrames Averaged\tProjection Type\tThreshold Cutoff");}

for (b=0; b<lis.length; b++){
  folderName2=dirStor+lis[b];
  //folderName2=replace(folderName,"/",fs);
  list2=getFileList(folderName2);
  fileName=folderName2+list2[0];
  if (endsWith(fileName, ".jpg")) {
    run("Image Sequence...", "open=["+ fileName + "] number=250 sort use");
 	makeRectangle(0, 950, 525, 600);
	run("Crop");
	saveName=getTitle();
    	saveAs("Tiff", dirSav+fs+saveName+"LowLeft");
    	close();
    	call("java.lang.System.gc");
	wait(1000);	
  }



//Get yo Wiggle on
fs=File.separator;
dirStore=dirSav+fs;
dirSave=dirStore+fs+"Output";
File.makeDirectory(dirSave);



scaleNumber=1;
gaussianNumber=2;
averageNumber=150;
projectionType="Standard Deviation";
thresholdCutoff=15;
setBatchMode(true);

list=getFileList(dirStore);

for (u=0; u<list.length; u++){
	fileName=dirStore+list[u];
	
		if (endsWith(fileName, ".tif")) {
			open(fileName);
			nameStore=getTitle();
			selectWindow(nameStore);
				if(gaussianNumber>1){
					selectWindow(nameStore);
					run("Gaussian Blur...", "sigma="+gaussianNumber+" stack");
}

			selectWindow(nameStore);
			if (scaleNumber==1){
				run("Duplicate...", "title=For_Processing duplicate range=1-10000");
			}

			if (scaleNumber<1){
			run("Scale...", "x="+scaleNumber+" y="+scaleNumber+" z=1.0 width=74 height=74 depth=450 interpolation=Bicubic process create title=For_Processing");
			}
			selectWindow("For_Processing");
			getDimensions(tempWidth, tempHeight, c, tempSlices, f);

			newImage("Wiggle Stack", "32-bit Black", tempWidth, tempHeight, 1);
			wait(100);
			//setBatchMode(true);
			selectWindow("For_Processing");
			setSlice(1);


			for (i=0; i<tempSlices-1; i++){
				selectWindow("For_Processing");
				currentSlice=getSliceNumber();
				run("Duplicate...", "title=[temp stack] duplicate range="+currentSlice+"-"+(currentSlice+(averageNumber-1))+"");
				selectWindow("temp stack");
				run("Z Project...", "start=1 stop="+averageNumber+" projection=[Standard Deviation]");
				selectWindow("STD_temp stack");
				run("Copy");
				selectWindow("Wiggle Stack");
				run("Add Slice");
				run("Paste");
				selectWindow("temp stack");
				close();
				selectWindow("STD_temp stack");
				close();
				selectWindow("For_Processing");
				setSlice(currentSlice+1);
	
}
			selectWindow("For_Processing");
			close();
			selectWindow(nameStore);
			close();



			selectWindow("Wiggle Stack");
			setSlice(1);
			run("Delete Slice");


			if(projectionType=="Maximum"){
				selectWindow("Wiggle Stack");
				run("Z Project...", "start=1 stop="+tempSlices+" projection=[Max Intensity]");
				selectWindow("MAX_Wiggle Stack");
				rename("Wiggle Project");
				selectWindow("Wiggle Stack");
				close();
}

			if(projectionType=="Average"){
				selectWindow("Wiggle Stack");
				run("Z Project...", "start=1 stop="+tempSlices+" projection=[Average Intensity]");
				selectWindow("AVG_Wiggle Stack");
				rename("Wiggle Project");
				selectWindow("Wiggle Stack");
				close();
}

			if(projectionType=="Standard Deviation"){
				selectWindow("Wiggle Stack");
				run("Z Project...", "start=1 stop="+tempSlices+" projection=[Standard Deviation]");
				selectWindow("STD_Wiggle Stack");
				rename("Wiggle Project");
				selectWindow("Wiggle Stack");
				close();
}



			selectWindow("Wiggle Project");
			setMinAndMax(0,200);
			run("8-bit");
			run("16_colors");
			run("Enhance Contrast", "saturated=0.00");
			setMinAndMax(thresholdCutoff,255);
			run("Apply LUT");

			run("Set Measurements...", "  mean limit display redirect=None decimal=3");

			run("Measure");
			resetThreshold();
			run("16_colors");
			run("Enhance Contrast", "saturated=0.00");
			wiggleIndex=getResult("Mean");
			rename(nameStore);
			print(tableTitle2, nameStore + "\t" + wiggleIndex + "\t" + scaleNumber + "\t" + gaussianNumber + "\t" + averageNumber + "\t" + projectionType + "\t" + thresholdCutoff);
			saveName=getTitle();
			saveAs("Tiff", dirSave+fs+saveName);
			
			if (isOpen("Results")) { 
       			selectWindow("Results"); 
       			run("Close"); 
   }
			if (isOpen("Log")) { 
       			selectWindow("Log"); 
       			run("Close"); 
   }
			selectWindow("Wiggle Index_LowLeft");
			saveAs("Text", dirSave+fs+"Wiggle Index LowLeft.tsv");
			File.delete(fileName);
   			call("java.lang.System.gc");
			wait(1000);
   			
		}}
			
}
//selectWindow("Wiggle Index_LowLeft");
//saveAs("Text", parent+fs+hi+"Wiggle Index LowLeft.tsv");
run("Quit");
