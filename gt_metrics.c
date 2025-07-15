#include "ift.h"

/* Author: Alexandre Xavier Falcão (September 10th, 2023)

   Description: Detects objects (e.g., parasite eggs) from the decoded
   saliency maps.

*/


iftImage *DrawBoxes(iftImage *orig, iftImage *comp, iftColor YCbCr, float scale, int type, int id)
{
  iftImage  *box     = iftCopyImage(orig);

  if (iftMaximumValue(comp)==0)
    return(box);

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
    
    
    //if(type == 0)printf("  result: %d xsize = %d, ysize = %d\n", id, xsize, ysize);
    //else printf("  ground-truth: %d xsize = %d, ysize = %d\n", id, xsize, ysize);
    
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

float min_area = 1e20, max_area = 0;
void update_metrics(iftImage *gt, int imageID){

	if(iftMaximumValue(gt) == 0)return;
	
	printf("size: %d %d\n", gt->xsize, gt->ysize);
	/* computing area */
	
  	
  	
  	iftAdjRel *A       = iftCircular(1.5);
  	
  	iftImage *label   = iftFastLabelComp(gt, A);
  	int nlabels = iftMaximumValue(label);
  	
  	for (int l=1; l <= nlabels; l++) {
  		
	    	iftImage *bin     = iftExtractObject(label,l);
	    	int ones = 0, zeros = 0;
	    	for(int i = 0; i < bin->n; i++){
	    		if((bin->val[i]) == 1)ones++;
	    		else zeros++;
	    	}
	    	
	    	int area = ones;
  		if(area > zeros)area = zeros;
  		if(area < 5)continue;
  		printf("area of image %d: %d\n", imageID, area);
  		min_area = fmin(min_area, 1.0*area);
  		max_area = fmax(max_area, 1.0*area);
  	}
  	
  	return;
}



int main(int argc, char *argv[])
{

  /* Example: gt_metrics */
  
 
  timer *tstart = iftTic();
  char *filename     = iftAllocCharArray(512);

  
  for(int i = 0; i < 96; i++) {
  
    sprintf(filename,"./truelabels/%06d.png",i + 1);
    iftImage *gt      = iftReadImageByExt(filename);
    update_metrics(gt, i + 1);
  }
  
  printf("min_area: %.4lf and max_area: %.4lf\n", min_area, max_area);
  iftFree(filename);
  printf("\nDone ... %s\n", iftFormattedTime(iftCompTime(tstart, iftToc())));
  
  return (0);
}
