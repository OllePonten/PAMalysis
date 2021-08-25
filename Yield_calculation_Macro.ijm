// Calculate Yield-original script from Kolbowski
//
// Example for calculation from PamWin multitiff export.
// Command on all the images in a stack.


  imagestack = getTitle();
  open("01-pim.lut"); // set LUT

  setBatchMode(true);
  input = getImageID();
  run("Duplicate...", "title=Ft.tif duplicate");
  Ftid = getImageID();
  run("Duplicate...", "title=Fs.tif duplicate");
  Fsid = getImageID();

// Remove Fm and absorbtivity 
  selectImage(Ftid);
  for (i=1; i<=4; i++) 
     run("Delete Slice");
  n = nSlices();
  for (i=n; i>=2; i = i-2) {
     setSlice(i);
     run("Delete Slice");
  }
// Remove Fo and absorbtivity
  selectImage(Fsid);
  for (i=1; i<=4; i++) 
     run("Delete Slice");
  n = nSlices();
  for (i=n-1; i>=1; i = i-2) {
     setSlice(i);
     run("Delete Slice");
  }

// Create Ft mask and eliminate outliers
  n = nSlices();
  selectImage(Ftid);
  run("Duplicate...", "title=Ftmask.tif duplicate");
  setThreshold(10,255);

  newImage("Maskstack", "8-Bit", 640, 480, n);
  stack = getImageID;

  for (i=1; i<=n; i++) {
    selectImage("Ftmask.tif");
    setSlice(i);  
    run("Create Mask");
    run("Remove Outliers...","radius=1 threshold=1 which=Dark");
    run("Copy");

    selectImage(stack);
    setSlice(i);  
    run("Paste"); 
  }

  imageCalculator("Substract create stack", "Fs.tif", "Ft.tif"); 
  rename("dF.tif");
  imageCalculator("AND stack", "dF.tif", "Maskstack"); 

  imageCalculator("Divide create 32-bit stack", "dF.tif", "Fs.tif"); 
  rename("yield_" + imagestack);
  setMinAndMax(0,1);
  setBatchMode(false);


// Make mask from the original yield images

run("Duplicate...", "duplicate");
setAutoThreshold("Huang dark");
run("Convert to Mask", "method=Huang background=Dark calculate black");
rename("Mask")

//This creates the stack of the mask images + NIR reference image, make sure that image1 refers to the mask image and NOT the yield image

selectWindow("Mask"); // change name here
run("Make Substack...", "  slices=3");
run("8-bit");
run("Concatenate...", "  title=[Concatenated Stacks] image1=[Mask] image2=[Substack (3)] "); // change name here image3=[--None--]

// Select ROIs in ROI manager on file= concantenated stacks (=mask), add one selection of ROIs by hand ("add(t)") at a time. I usually select on the NIR image and then select the desired masek image (see below_)

run("ROI Manager...");
//roiManager("Delete");
waitForUser("Draw ROI, then hit OK");   
setTool("rectangle");

// Select desired mask image (e.g. max quant yield =1, effective quantum yield =8)

//run("Analyze Particles...", "add in_situ slice");
maxsize = getNumber("Max size in pixels for particles", 150);
run("Analyze Particles...", "size=3-"+maxsize+" pixel add in_situ slice"); // size has to be defined for each cell type to omit spurious particles but also capture all cells
selectWindow("yield_" + imagestack); // change name here
waitForUser("Select image to measure, then hit OK");   
run("Set Measurements...", "area mean limit display redirect=None decimal=3");
roiManager("multi-measure measure_all");
//saveAs("Results", "C:\\Users\\Olle\\Desktop\\Strain421_pH8_2.csv"); // change file save name here

// The resulting CSV file can be imported into Origin/Excel and hereafter filtered by area to omit spurious -small- measurements. Alternatively this parameter can also be adjusted to liking
// in the Analyze particle settings
