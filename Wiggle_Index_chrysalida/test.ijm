
//Crop yo Files
fs=File.separator;
//dirStor=getDirectory("Select Source Directory");
dirStor="/home/sdenecke/Wiggle_Analyze/"
dirSav=dirStor+"zCroppedTiffStackLL";
File.makeDirectory(dirSav);
setBatchMode(true);
lis=getFileList(dirStor);
//ds=lengthOf(dirStor);
//hi=substring(dirStor,(ds-2),(ds-1));
//parent=substring(dirStor,0,(ds-3));



tableTitle="Wiggle Index_LowLeft";
tableTitle2="["+tableTitle+"]";
if (isOpen("Wiggle Index_LowLeft")){
	print("Table already open");}
	
	else{
		
		run("Table...", "name="+tableTitle2+" width=600 height=250");
		print(tableTitle2, "\\Headings:Image Name\tWiggle Index\tImage Scale\tImage Blur\tFrames Averaged\tProjection Type\tThreshold Cutoff");}





folderNameTest=dirStor+lis[1];
listTest=getFileList(folderNameTest);
fileNameTest=folderNameTest+listTest[1];
run("Image Sequence...", "open=["+ fileNameTest + "] number=250 sort use");
makeRectangle(0, 100, 0, 100);
run("Crop");
saveName=getTitle();
saveAs("Tiff", dirSav+fs+saveName+"LowLeft");
close();
call("java.lang.System.gc");
wait(1000);	

//setTool("rectangle");
//waitForUser("Waiting for user to draw a rectangle around Low Left...");
//rect2=Roi.getBounds();
//makeRectangle(rect);
//close();
//call("java.lang.System.gc");
//wait(1000);	
