i = getTitle();
setBatchMode(true);
setBatchMode("hide");
corr_cutoff = 0.5; //The cutoff ratio of "max_corr - min_corr" that will be defined as the boundary of the interline stimulation
stim_padding = 1; //Number of pixels to add to the boundary of the interline to ensure all interline is excluded
stim_bin_half_width = 5; //± # of slices before and after the stim period to add to autocorrelation analysis
//Make sure slices are on the Z-axis
selectWindow(i);
Stack.getDimensions(stack_width, stack_height, stack_channels, stack_slices, stack_frames);
if(stack_slices == 1 && stack_channels > 1) run("Re-order Hyperstack ...", "channels=[Slices (z)] slices=[Channels (c)] frames=[Frames (t)]");
else if (stack_slices == 1 && stack_frames > 1) run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
if(isOpen("")){
	close(i);
	selectWindow("");
	rename(i);
}

interline_boundary = findInterlineBoundary(i);
stim_period = findStimPeriod(i, interline_boundary);
registerStack(i, interline_boundary);
autocorrelation("Reg_Stack", interline_boundary, stim_period, stack_width);
dfStack("Reg_Stack", interline_boundary, stack_width);

//Restore original boundary
selectWindow("Reg_Stack");
Stack.getDimensions(reg_width, reg_height, reg_channels, reg_slices, reg_frames);
run("Canvas Size...", "width=" + reg_width + interline_boundary[0] + " height=" + reg_height + " position=Center-Right zero"); 
run("Canvas Size...", "width=" + stack_width + " height=" + reg_height + " position=Center-Left zero");
makeRectangle(interline_boundary[0], interline_boundary[1], interline_boundary[2], interline_boundary[3]);
run("Make Inverse");
run("Set...", "value=NaN stack");
run("Select None");
setBatchMode("exit and display");

function autocorrelation(i, interline_boudary, mean_stim_period, original_width){
	//Create variance stack
	selectWindow(i);
	run("32-bit");
	Stack.getDimensions(stack_width, stack_height, stack_channels, stack_slices, stack_frames);
	newImage("Autocorrelation", "32-bit black", stack_width, stack_height, stack_slices-1);
	selectWindow(i);
	run("Z Project...", "projection=[Standard Deviation]");
	selectWindow("STD_" + i);
	run("Square");
	rename("Variance");
	
	//Create mean stack
	selectWindow(i);
	run("Z Project...", "projection=[Average Intensity]");
	selectWindow("AVG_" + i);
	rename("Mean");
	
	//Create reference and phase shift stack
	selectWindow(i);
	run("Duplicate...", "title=Ref duplicate");
	run("Duplicate...", "title=Shift duplicate");
	imageCalculator("Subtract 32-bit stack", "Ref","Mean");
	imageCalculator("Subtract 32-bit stack", "Shift","Mean");

	//Generate autocorrelation stack
	for(k=0; k<stack_slices-1; k++){
		showProgress(k/(stack_slices-1));
	
		//Increment phase shift by 1
		if(k != 0){
			selectWindow("Ref");
			Stack.setSlice(nSlices);
			run("Delete Slice");
			selectWindow("Shift");
			Stack.setSlice(1);
			run("Delete Slice");
		}
		
		stim_phase_shift = (k%round(mean_stim_period));
		if((stim_phase_shift <= stim_bin_half_width || (mean_stim_period - stim_phase_shift) <= stim_bin_half_width) && round(k/mean_stim_period) > 0){
			//Generate normalization image
			selectWindow("Variance");
			run("Duplicate...", "title=Normalization");
			run("Multiply...", "value=" + stack_slices-k);
		
			//Calculate autocorrelation
			imageCalculator("Multiply create 32-bit stack", "Ref","Shift");
			selectWindow("Result of Ref");
			run("Z Project...", "projection=[Sum Slices]");
			imageCalculator("Divide create 32-bit", "SUM_Result of Ref","Normalization");
		
			//Copy result to autocorrelation stack
			run("Select All");
			run("Copy");
			selectWindow("Autocorrelation");
			Stack.setSlice(k+1);
			run("Paste");
			run("Select None");
		
			//Add phase position relative to stim period to metadata
			stim_string = "Stimulation period " + round(k/mean_stim_period) + " ";
			if(stim_phase_shift > (mean_stim_period / 2)){
				stim_phase_shift = round(mean_stim_period-stim_phase_shift);
				stim_string += "- " + stim_phase_shift + " slices.";
			}
			else{
				stim_string += "+ " + stim_phase_shift + " slices.";
			}
			setMetadata("label",  stim_string);
		
			//Close images
			close("Result of SUM_Result of Ref");
			close("SUM_Result of Ref");
			close("Result of Ref");
			close("Normalization");
		}
	}
	close("Ref");
	close("Shift");
	close("Mean");
	close("Variance");
	
	//Remove blank slices from autocorrelation stack
	selectWindow("Autocorrelation");
	slice = 1;
	while(slice <= nSlices){
		setSlice(slice++);
		getStatistics(dummy, dummy, dummy, max);
		if(max == 0){
			run("Delete Slice");
			slice--;
		}
	}
	
	//Create autocorrelation image
	selectWindow("Autocorrelation");
	run("Canvas Size...", "width=" + stack_width + interline_boundary[0] + " height=" + stack_height + " position=Center-Right zero"); //Restore original boundary
	run("Canvas Size...", "width=" + original_width + " height=" + stack_height + " position=Center-Left zero");
	makeRectangle(interline_boundary[0], interline_boundary[1], interline_boundary[2], interline_boundary[3]);
	run("Make Inverse");
	run("Set...", "value=NaN stack");
	run("Select None");
	run("Z Project...", "projection=[Average Intensity]");
	selectWindow("AVG_Autocorrelation");
	reds = newArray(255, 253, 251, 249, 247, 245, 243, 241, 239, 237, 235, 233, 231, 229, 227, 225, 223, 221, 219, 217, 215, 213, 211, 209, 207, 205, 203, 201, 199, 197, 195, 193, 191, 189, 187, 185, 183, 181, 179, 177, 175, 173, 171, 169, 167, 165, 163, 161, 159, 157, 155, 153, 151, 149, 147, 145, 143, 141, 139, 137, 135, 133, 131, 129, 127, 125, 123, 121, 119, 117, 115, 113, 111, 109, 107, 105, 103, 101, 99, 97, 95, 93, 91, 89, 87, 85, 83, 81, 79, 77, 75, 73, 71, 69, 67, 65, 63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	greens = newArray(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119, 121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143, 145, 147, 149, 151, 153, 155, 157, 159, 161, 163, 165, 167, 169, 171, 173, 175, 177, 179, 181, 183, 185, 187, 189, 191, 193, 195, 197, 199, 201, 203, 205, 207, 209, 211, 213, 215, 217, 219, 221, 223, 225, 227, 229, 231, 233, 235, 237, 239, 241, 243, 245, 247, 249, 251, 253);
	blues = newArray(255, 253, 251, 249, 247, 245, 243, 241, 239, 237, 235, 233, 231, 229, 227, 225, 223, 221, 219, 217, 215, 213, 211, 209, 207, 205, 203, 201, 199, 197, 195, 193, 191, 189, 187, 185, 183, 181, 179, 177, 175, 173, 171, 169, 167, 165, 163, 161, 159, 157, 155, 153, 151, 149, 147, 145, 143, 141, 139, 137, 135, 133, 131, 129, 127, 125, 123, 121, 119, 117, 115, 113, 111, 109, 107, 105, 103, 101, 99, 97, 95, 93, 91, 89, 87, 85, 83, 81, 79, 77, 75, 73, 71, 69, 67, 65, 63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11, 9, 7, 5, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	setLut(reds, greens, blues);
	getStatistics(dummy, dummy, min, max);
	if(abs(min) > max) setMinAndMax(min, abs(min));
	else setMinAndMax(-1*max, max);
	close("Autocorrelation");
}

function registerStack(i, interline_boundary){
	//Crop off interline stim
	selectWindow(i);
	run("Duplicate...", "title=crop_stack duplicate");
	i = "crop_stack";
	selectWindow(i);
	makeRectangle(interline_boundary[0], interline_boundary[1], interline_boundary[2], interline_boundary[3]);
	run("Crop");
	slices = nSlices;
	setSlice(1);
	run("Duplicate...", "title=[Ref_Slice]");
	selectWindow("Ref_Slice");
	run("Duplicate...", "title=[Reg_Stack]");
	for(a=2; a<=slices; a++){
		showProgress(a, slices);
		selectWindow(i);
		setSlice(a);
		run("Duplicate...", "title=[Reg_Slice]");
		run("Concatenate...", "  title=Reg_Pair open image1=Ref_Slice image2=Reg_Slice image3=[-- None --]");
		selectWindow("Reg_Pair");
		setSlice(1);
		run("StackReg", "transformation=Affine");
		selectWindow("Reg_Pair");

		run("Re-order Hyperstack ...", "channels=[Frames (t)] slices=[Slices (z)] frames=[Channels (c)]");
		run("Split Channels");
		selectWindow("C1-Reg_Pair");
		rename("Ref_Slice");
		run("Concatenate...", "  title=Reg_Stack open image1=Reg_Stack image2=C2-Reg_Pair image3=[-- None --]");
	
		//run("Stack to Images");
		//run("Concatenate...", "  title=Reg_Stack open image1=Reg_Stack image2=Reg_Slice image3=[-- None --]");
	}
	close(i);
	close("Ref_Slice");

	//Make sure slices are on the Z-axis
	selectWindow("Reg_Stack");
	Stack.getDimensions(stack_width, stack_height, stack_channels, stack_slices, stack_frames);
	if(stack_slices == 1 && stack_channels > 1) run("Re-order Hyperstack ...", "channels=[Slices (z)] slices=[Channels (c)] frames=[Frames (t)]");
	else if (stack_slices == 1 && stack_frames > 1) run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
	if(isOpen("")){
		close("Reg_Stack");
		selectWindow("");
		rename("Reg_Stack");
	}
}

function findStimPeriod(i, interline_boudary){	
	//Find stim events
	selectWindow(i);
	getDimensions(width, height, dummy, stack_slices, frames);
	makeRectangle(1, 0, width-2, height);
	run("Make Inverse"); //Select the interline region
	Stack.getStatistics(dummy, stack_mean, dummy, dummy, dummy);
	z_profile = newArray(stack_slices);
	max_stim_slices = -1;
	stim_start = false;
	index = 0;
	for(a=1; a<=stack_slices; a++){
		Stack.setSlice(a);
		getStatistics(dummy, mean, dummy, dummy, dummy);
		if(mean > stack_mean && !stim_start){
			stim_start = true;
			z_profile[index++] = a;
		}
		else if(mean < stack_mean && stim_start){
			stim_start = false;
			z_profile[index++] = a;
			cur_stim_slices = a - z_profile[index-2];
			if(cur_stim_slices > max_stim_slices){
				max_stim_slices = cur_stim_slices;
			}
		}
	}
	run("Select None");
	z_profile = Array.slice(z_profile,0,index);
	
	//Calculate stim period
	mean_stim_period = 0;
	for(a=0; a<z_profile.length-2; a++){
		mean_stim_period += z_profile[a+2] - z_profile[a];
	}
	mean_stim_period /= z_profile.length-2;
	run("Select None");
	return mean_stim_period;
}

function findInterlineBoundary(i){
	//Get correlation of image with interline pattern
	selectWindow(i);
	getDimensions(width, height, dummy, dummy, dummy);
	ref_array = newArray(height);
	corr_array = newArray(width);
	for(a=0; a<ref_array.length; a++) ref_array[a] = a%2; //Populate ref array with interline pattern
	Array.getStatistics(ref_array, dummy, dummy, ref_mean, dummy);
	run("Z Project...", "projection=[Max Intensity]");
	selectWindow("MAX_" + i);
	for(x=0; x<width; x++){
		makeRectangle(x, 0, 1, height);
		getStatistics(dummy, px_mean, dummy, dummy, dummy);
		num_corr = 0;
		px_corr = 0;
		ref_corr = 0;
		for(y=0; y<height; y++){
			px = getPixel(x,y);
			num_corr += (px - px_mean)*(ref_array[y] - ref_mean);
			px_corr += (px - px_mean)*(px - px_mean);
			ref_corr += (ref_array[y] - ref_mean)*(ref_array[y] - ref_mean);
		}
		corr_array[x] = abs(num_corr/(sqrt(px_corr)*sqrt(ref_corr)));
	}
	close("MAX_" + i);
	Array.getStatistics(corr_array, corr_min, corr_max, corr_mean, corr_std);
	x=-1;
	w = -1;
	cutoff = corr_cutoff*(corr_max - corr_min);
	interline_boundary = newArray(4);
	interline_boundary[0] = -1;
	interline_boundary[1] = 0;
	interline_boundary[2] = -1;
	interline_boundary[3] = height;
	for(a=0; a<corr_array.length; a++){
		if(interline_boundary[0]<0 && corr_array[a] < corr_cutoff){
			interline_boundary[0] = a + stim_padding;
		}
		else if(interline_boundary[0]>0 && interline_boundary[2]<0 && corr_array[a] > corr_cutoff){
			interline_boundary[2] = a - stim_padding - 1 - interline_boundary[0];
		}
	}
	return interline_boundary;
}


function dfStack(i, interline_boudary, original_width){
	selectWindow(i);
	
	run("Select None");
	run("Duplicate...", "title=dF-F duplicate");
	selectWindow("dF-F");
	run("32-bit");
	run("Median...", "radius=0.5 stack"); //Denoise
	run("Z Project...", "projection=Median"); //Get baseline reference
	imageCalculator("Divide stack", "dF-F","MED_dF-F");
	close("MED_dF-F");
	selectWindow("dF-F");
	run("Subtract...", "value=1 stack"); //Normalize to 0
	run("Duplicate...", "title=dF_Mask duplicate");
	selectWindow("dF_Mask");
	setThreshold(0.000000000, 1000000000000000000000000000000.000000000);
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=Default background=Dark");
	run("Divide...", "value=255 stack");
	run("32-bit");
	imageCalculator("Multiply stack", "dF-F","dF_Mask");
	close("dF_Mask");

	//Restore original boundary
	selectWindow("dF-F");
	Stack.getDimensions(stack_width, stack_height, stack_channels, stack_slices, stack_frames);
	run("Canvas Size...", "width=" + stack_width + interline_boundary[0] + " height=" + stack_height + " position=Center-Right zero"); 
	run("Canvas Size...", "width=" + original_width + " height=" + stack_height + " position=Center-Left zero");
	makeRectangle(interline_boundary[0], interline_boundary[1], interline_boundary[2], interline_boundary[3]);
	run("Make Inverse");
	run("Set...", "value=NaN stack");
	run("Select None");	
}


//selectWindow("AVG_Autocorrelation");
//run("Duplicate...", "title=pixelcount");
//selectWindow("pixelcount");
//run("Threshold...");
//setThreshold(0.1131, 1000000000000000000000000000000.0000);
//setOption("BlackBackground", false);
//run("Make Binary");
//run("Set Measurements...", "area limit redirect=None decimal=3");
//run("Measure");
//close("Threshold");

selectWindow("Reg_Stack");
run("Z Project...", "projection=[Average Intensity]");

close("pixelcount");
//close("Reg_Stack");
//close("dF-F");
