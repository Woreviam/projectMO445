#include "ift.h"

/* Author: Alexandre Xavier Falcão (September 10th, 2023)

   Description: Detects objects (e.g., parasite eggs) from the decoded
   saliency maps.

*/

//It was a intended changed of the original function in order to compute components that no necesarily the pixels are equal value but different than zero.
iftImage *iftFastLabelComp2(const iftImage *bin, const iftAdjRel *Ain) {
    iftAdjRel *A = NULL;
    if (Ain == NULL) {
        if (iftIs3DImage(bin))
            A = iftSpheric(1.74);
        else A = iftCircular(1.45);
    }
    else A = iftCopyAdjacency(Ain);

  iftImage *label=NULL;
  int i,p,q,l=1, *cost;
  iftVoxel u,v;
  iftGQueue *Q;

  label  = iftCreateImageFromImage(bin);
  cost   = iftAllocIntArray(bin->n);
  Q      = iftCreateGQueue(2,bin->n,cost);

  for (p=0; p < bin->n; p++){
    if ((bin->val[p]!=0)&&(label->val[p]==0)){
      cost[p] = 1;
      iftInsertGQueue(&Q,p);
    }
  }

  while(!iftEmptyGQueue(Q)){
    p = iftRemoveGQueue(Q);
    if (cost[p]==1){ /* p is the first in a component */
      cost[p]=0;
      label->val[p] = l; l++;
    }
    u = iftGetVoxelCoord(bin,p);
    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(bin,v)){
	q = iftGetVoxelIndex(bin,v);
	if ((bin->val[p] > 0 && bin->val[q] > 0)&&(label->val[q] == 0)){
	      label->val[q] = label->val[p];
	      iftRemoveGQueueElem(Q,q);
	      cost[q]=0;
	      iftInsertGQueue(&Q,q);
	}
      }
    }
  }

  iftDestroyAdjRel(&A);

  iftDestroyGQueue(&Q);
  iftFree(cost);
  return(label);
}

//We added a parameter to play with the scale bounding box factor.
iftImage *DrawBoxes(iftImage *orig, iftImage *comp, iftColor YCbCr, float scale, int id)
{
  iftImage  *box     = iftCopyImage(orig);

  if (iftMaximumValue(comp)==0)
    return(box);

  iftAdjRel *A       = iftCircular(1.5);
  iftImage  *label   = iftFastLabelComp(comp, A);
  int        nlabels = iftMaximumValue(label);
  
  //printf("number of components in image %d is: %d\n", id, nlabels);
  
  for (int l=1; l <= nlabels; l++) {
    iftImage *bin     = iftExtractObject(label,l);
    iftVoxel center;
    iftBoundingBox bb = iftMinBoundingBox(bin, &center);
    iftDestroyImage(&bin);
    int xsize = (bb.end.x - bb.begin.x) * scale;
    int ysize = (bb.end.y - bb.begin.y) * scale;
    
    if (xsize*ysize > 0){
      iftAdjRel *B = iftRectangular(xsize,ysize);
      iftAdjRel *C = iftAdjacencyBoundaries(B,A);    
      iftDestroyAdjRel(&B);
      for (int i=1; i < C->n; i++) {
	iftVoxel v = iftGetAdjacentVoxel(C,center,i);
	if (iftValidVoxel(box,v)){
	  int q       = iftGetVoxelIndex(box,v);
	  box->val[q] = YCbCr.val[0];
	  box->Cb[q]  = YCbCr.val[1];
	  box->Cr[q]  = YCbCr.val[2];
	}
      }
      iftDestroyAdjRel(&C);
    }
  }
  iftDestroyAdjRel(&A);
  iftDestroyImage(&label);
  
  return(box);
}

//IoU avaliation also handles false positive cases, we work with two different scale bounding box
double IoU_avaliation(iftImage *orig, iftImage *comp, iftImage *gt, float scale1, float scale2){

	if(iftMaximumValue(comp) == 0 && iftMaximumValue(gt) == 0)return 1.0;
	if(iftMaximumValue(comp) == 0 || iftMaximumValue(gt) == 0)return 0;
		
	/* Initialization of frequency table */
	iftImage *frequency = iftCopyImage(orig);
	for(int i = 0; i < frequency -> n; i++)frequency->val[i] = 0;
	
	
	iftAdjRel *A       = iftCircular(1.5);
  	iftImage  *label   = iftFastLabelComp(comp, A);
  	int        nlabels = iftMaximumValue(label);
  	
  	for (int l=1; l <= nlabels; l++) {
    	iftImage *bin     = iftExtractObject(label,l);
    	iftVoxel center;
    	iftBoundingBox bb = iftMinBoundingBox(bin, &center);
    	iftDestroyImage(&bin);
    	int xsize = (bb.end.x - bb.begin.x)*scale1;
    	int ysize=(bb.end.y - bb.begin.y)*scale1;
    	if (xsize*ysize > 0){
    		
      		iftAdjRel *B = iftRectangular(xsize,ysize);
      		for (int i=1; i < B->n; i++) {
				iftVoxel v = iftGetAdjacentVoxel(B,center,i);
				if (iftValidVoxel(frequency,v)){
	  				int q       = iftGetVoxelIndex(frequency,v);
	  				frequency->val[q]++;
				}
      		}	
      		iftDestroyAdjRel(&B);
    	}
  	}
  	
  	label   = iftFastLabelComp(gt, A);
  	nlabels = iftMaximumValue(label);
  	
  	for (int l=1; l <= nlabels; l++) {
    	iftImage *bin     = iftExtractObject(label,l);
    	iftVoxel center;
    	iftBoundingBox bb = iftMinBoundingBox(bin, &center);
    	iftDestroyImage(&bin);
    	int xsize = (bb.end.x - bb.begin.x) * scale2;
    	int ysize=(bb.end.y - bb.begin.y) * scale2;
    	if (xsize*ysize > 0){
    		
      		iftAdjRel *B = iftRectangular(xsize,ysize);
      		for (int i=1; i < B->n; i++) {
				iftVoxel v = iftGetAdjacentVoxel(B,center,i);
				if (iftValidVoxel(frequency,v)){
	  				int q       = iftGetVoxelIndex(frequency,v);
	  				frequency->val[q]++;
				}
      		}	
      		iftDestroyAdjRel(&B);
    	}
  	}
  	
  	iftDestroyAdjRel(&A);
  	iftDestroyImage(&label);
	
	int area_intersection = 0, area_union = 0;
	for(int i = 0; i < frequency -> n; i++){
		
		if((frequency->val[i]) > 0)area_union++;
		if((frequency->val[i]) > 1)area_intersection++;
	}
	
	return (area_intersection * 1.0) / area_union;
}

//DICE avaliation also handles false positive cases, we work with two different scale bounding box
double DICE_avaliation(iftImage *orig, iftImage *comp, iftImage *gt, float scale){

	if(iftMaximumValue(comp) == 0 && iftMaximumValue(gt) == 0)return 1.0;
	if(iftMaximumValue(comp) == 0 || iftMaximumValue(gt) == 0)return 0;
		
	/* Initialization of frequency table */
	iftImage *frequency = iftCopyImage(orig);
	for(int i = 0; i < frequency -> n; i++)frequency->val[i] = 0;
	
	
	iftAdjRel *A       = iftCircular(1.5);
  	iftImage  *label   = iftFastLabelComp(comp, A);
  	int        nlabels = iftMaximumValue(label);
  	
  	for (int l=1; l <= nlabels; l++) {
    	iftImage *bin     = iftExtractObject(label,l);
    	iftVoxel center;
    	iftBoundingBox bb = iftMinBoundingBox(bin, &center);
    	iftDestroyImage(&bin);
    	int xsize = (bb.end.x - bb.begin.x)*scale;
    	int ysize=(bb.end.y - bb.begin.y)*scale;
    	if (xsize*ysize > 0){
    		
      		iftAdjRel *B = iftRectangular(xsize,ysize);
      		for (int i=1; i < B->n; i++) {
				iftVoxel v = iftGetAdjacentVoxel(B,center,i);
				if (iftValidVoxel(frequency,v)){
	  				int q       = iftGetVoxelIndex(frequency,v);
	  				frequency->val[q]++;
				}
      		}	
      		iftDestroyAdjRel(&B);
    	}
  	}
  	
  	label   = iftFastLabelComp(gt, A);
  	nlabels = iftMaximumValue(label);
  	
  	for (int l=1; l <= nlabels; l++) {
    	iftImage *bin     = iftExtractObject(label,l);
    	iftVoxel center;
    	iftBoundingBox bb = iftMinBoundingBox(bin, &center);
    	iftDestroyImage(&bin);
    	int xsize = (bb.end.x - bb.begin.x);
    	int ysize=(bb.end.y - bb.begin.y);
    	if (xsize*ysize > 0){
    		
      		iftAdjRel *B = iftRectangular(xsize,ysize);
      		for (int i=1; i < B->n; i++) {
				iftVoxel v = iftGetAdjacentVoxel(B,center,i);
				if (iftValidVoxel(frequency,v)){
	  				int q       = iftGetVoxelIndex(frequency,v);
	  				frequency->val[q]++;
				}
      		}	
      		iftDestroyAdjRel(&B);
    	}
  	}
  	
  	iftDestroyAdjRel(&A);
  	iftDestroyImage(&label);
	
	int area_intersection = 0, area_union = 0;
	for(int i = 0; i < frequency -> n; i++){
		
		area_union += frequency->val[i];
		if((frequency->val[i]) > 1)area_intersection++;
	}
	
	return (2 * area_intersection * 1.0) / area_union;
}

//Computes an iftImage with only one single region, the one with hightest density pixel values.
iftImage *iftComponentDensityArea(iftImage *bin, iftImage *salie, iftAdjRel *A)
{
  iftImage *label;
  int      *size, ncomps, p;
  float *density;
  
  label  = iftFastLabelComp(bin,A);
  ncomps = iftMaximumValue(label);
  
  size   = iftAllocIntArray(ncomps+1);
  density   = iftAllocFloatArray(ncomps+1);
  //printf("testing density[0]: %.4lf\n", density[0]);
  
  for (p=0; p < label->n; p++)
    if (label->val[p] != 0){
      size[label->val[p]]++;
      density[label->val[p]] += salie->val[p];
    }
  
  //for(int i = 1; i <= ncomps; i++)printf("pair testing: %d %.4lf\n", size[i], density[i]);
  
  for (p=0; p < label->n; p++)
    label->val[p] = (density[label->val[p]] * 1.0)/size[label->val[p]];

  iftFree(size);
  iftCopyVoxelSize(bin,label);

  return(label);
}

/*
1) Initially we filter the components with area >= thre_min and <= thres_max
2) We compute the density of the components by dividing the sum of values of pixesl in a component by the area of it
3) We chose the component with greatest density for the final chosen region
4) At the end there will be at most one component to draw the bounding box
*/
iftImage *iftSelectMaxCompDensity(iftImage *bin, iftImage *salie, iftAdjRel *A, int thres_min, int thres_max, int id)
{
  iftImage *area  = iftComponentArea(bin, A);
  iftImage *filtered_areas = iftCopyImage(area);
  //bool flag = true;
  for (int p=0; p < area->n; p++){
    
    /*if(flag && bin->val[p] > 0 && bin->val[p] < 255.0){
    	
    	flag = false;
    	printf("pixel found: %d\n", bin->val[p]);
    }*/
    if ((filtered_areas->val[p] >= thres_min) && (filtered_areas->val[p] <= thres_max)){
    	
    	//if(area->val[p] < min_area)min_area = area->val[p];
    	filtered_areas->val[p] = bin->val[p];
    }
    else filtered_areas->val[p] = 0;
  }
  
  iftCopyVoxelSize(bin, area);
  iftCopyVoxelSize(bin, filtered_areas);
  //printf("max area in image %d: %d\n", id, min_area);
  
  iftImage *density  = iftComponentDensityArea(filtered_areas, salie, A);
  float max_density = 0, max_area = 0;
  
  for (int p=0; p < density->n; p++){
    if (density->val[p] > max_density){
    	
    	max_density = density->val[p];
    	max_area = area->val[p];
    }
    else if(max_density == density->val[p] && max_area < area->val[p])max_area = area->val[p];
  }
  
  //printf("density in image %d: %.4lf\n", id, max_density);
  
  for (int p=0; p < density->n; p++){
    if ((density->val[p] == max_density && max_area == area->val[p]) )
      density->val[p] = bin->val[p];
    else
      density->val[p] = 0;
  }
  
  //printf("%.4lf\n", max_density);
  //iftImage *label  = iftFastLabelComp2(density,A);
  //printf("detection-number of components in image %d: %d\n", id, iftMaximumValue(label));
  
  iftCopyVoxelSize(bin,density);
  return(density);
}



int main(int argc, char *argv[])
{

  /* Example: detection salie 1 boxes */
  
  if (argc!=4){ 
    iftError("Usage: detection <P1> <P2> <P3>\n"
	     "[1] folder with the salience maps\n"
	     "[2] layer (1,2,...) to create the results\n"
	     "[3] output folder with the resulting bounding boxes\n",	 
	     "main");
  }
  
  timer *tstart = iftTic();

  char *filename     = iftAllocCharArray(512);
  int layer          = atoi(argv[2]);
  char suffix[12];
  sprintf(suffix,"_layer%d.png",layer);
  iftFileSet *fs     = iftLoadFileSetFromDirBySuffix(argv[1], suffix, true);
  char *output_dir   = argv[3];
  iftMakeDir(output_dir);
  iftColor RGB, YCbCr;
	
  float detection_error = 0;
  int testCases = 0;
  
  for(int i = 0; i < fs->n; i++) {

    //printf("Processing image %d of %ld\r", i + 1, fs->n);
    char *basename1   = iftFilename(fs->files[i]->path,suffix);      
    char *basename2   = iftFilename(fs->files[i]->path,".png");      
    iftImage *salie   = iftReadImageByExt(fs->files[i]->path);
    sprintf(filename,"./truelabels/%s.png",basename1);
    iftImage *gt      = iftReadImageByExt(filename);
    sprintf(filename,"./images/%s.png",basename1);
    iftImage *orig    = iftReadImageByExt(filename);
    int Imax          = iftNormalizationValue(iftMaximumValue(orig));
 
    /* Detect parasite */

    iftAdjRel *A   = iftCircular(1.5);
    /*
    bool flag = true;
    for(int j = 0; j < salie->n; j++)
    	if(salie->val[j] > 0 && salie->val[j] < 255.0 && flag)printf("pixel tested: %d\n", salie->val[j]), flag = true;
    */
    
    
    
    iftImage *bin  = iftThreshold(salie,iftOtsu(salie),IFT_INFINITY_INT,255);
    
    
    
    iftImage *comp = iftSelectMaxCompDensity(bin, salie, A, 2000, 9000, i);
    
    //iftImage *label  = iftFastLabelComp(comp,A);
    //int ncomps = iftMaximumValue(label);
    //printf("* number of components: %d\n", ncomps);
    
    //iftImage *comp = iftSelectMaxCompDensity(salie, A, 400, 4000, i);
    
    //printf("* availiation %d: %.3lf\n", i, avaliation(orig, comp, gt, 1.5));
    detection_error += IoU_avaliation(orig, comp, gt, 1.6, 1.5);
    //printf("detection_error in image %d: %.4lf\n", i + 1, detection_error);
    //detection_error += DICE_avaliation(orig, comp, gt, 1.5);
    
    RGB.val[0]     = Imax;
    RGB.val[1]     = 0;
    RGB.val[2]     = Imax;
    YCbCr          = iftRGBtoYCbCr(RGB,Imax);
    iftImage *img1 = DrawBoxes(orig, comp, YCbCr, 1.6, i);
    RGB.val[0]     = 0;
    RGB.val[1]     = Imax;
    RGB.val[2]     = Imax; 
    YCbCr          = iftRGBtoYCbCr(RGB,Imax);
    iftImage *img2 = DrawBoxes(img1, gt, YCbCr, 1.5, i);

    iftDestroyImage(&bin);
    iftDestroyImage(&comp);
    iftDestroyAdjRel(&A);
    iftDestroyImage(&salie);
    iftDestroyImage(&gt);
    iftDestroyImage(&orig);
    iftDestroyImage(&img1);

    /* save resulting image */

    sprintf(filename,"%s/%s.png",output_dir,basename2);
    iftWriteImageByExt(img2,filename);

    iftDestroyImage(&img2);
    iftFree(basename1);
    iftFree(basename2);
    testCases++;
  }
  
  //printf("average detection_error: %.4lf\n",  detection_error*1.0/testCases);
  
  FILE *fp;
  fp = fopen("detection_result.txt", "w");
  fprintf(fp, "%.4lf\n", detection_error/(1.0 * testCases));  
  fclose(fp);
  
  //printf("%.4lf\n", detection_error);
  iftFree(filename);
  iftDestroyFileSet(&fs);
  
  printf("\nDone ... %s\n", iftFormattedTime(iftCompTime(tstart, iftToc())));
  
  return (0);
}
