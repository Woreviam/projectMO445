#include "ift.h"

/* Author: Alexandre Xavier Falcão (September 10th, 2023)

   Description: Delineates objects (e.g., parasite eggs) from the
   decoded saliency maps.

*/

/* Complete the code below to delineate objects using Dynamic
   Trees. Use Equation 16 from articles/2019_IFT_DGCI.pdf */

iftImage *DynamicTrees(iftImage *orig, iftImage *seeds_in, iftImage *seeds_out)
{
  iftMImage  *mimg     = iftImageToMImage(orig,LAB_CSPACE);
  iftImage   *pathval  = NULL, *label = NULL, *root = NULL;
  float *tree_L = NULL;
  float *tree_A = NULL;
  float *tree_B = NULL;
  int *number_of_nodes = NULL;
  int intMax = iftRound(sqrtf(3.0)*
  		iftMax(iftMax(iftMMaximumValue(mimg, 0),
	  		iftMMaximumValue(mimg, 1)),
	  		iftMMaximumValue(mimg, 2)
  		));
  iftGQueue *Q = NULL;
  iftAdjRel *adjacency = iftCircular(1.0);
  iftVoxel u,v;
  int i,p,q,r,tmp;
  //Initialization
  pathval = iftCreateImage(orig->xsize,orig->ysize,orig->zsize);
  label = iftCreateImage(orig->xsize,orig->ysize,orig->zsize);
  root = iftCreateImage(orig->xsize,orig->ysize,orig->zsize);

  tree_L = iftAllocFloatArray(orig->n);
  tree_A = iftAllocFloatArray(orig->n);
  tree_B = iftAllocFloatArray(orig->n);
  number_of_nodes = iftAllocIntArray(orig->n);

  Q = iftCreateGQueue(intMax + 1, orig->n, pathval->val);

//Setting things up
  for(p = 0; p < orig->n; p++){
    pathval->val[p] = IFT_INFINITY_INT;
    if(seeds_in->val[p] != 0){
      root->val[p] = p;
      label->val[p] = seeds_in->val[p];
      pathval->val[p] = 0;
    }
    else{
      if(seeds_out->val[p] != 0){
        root->val[p] = p;
        label->val[p] = 0;
        pathval->val[p] = 0;
      }
    }
    iftInsertGQueue(&Q, p);
  }
//The IFT itself
  while(!iftEmptyGQueue(Q)){
    p = iftRemoveGQueue(Q);
    r = root->val[p];
    tree_L[r] += mimg->val[p][0];
    tree_A[r] += mimg->val[p][1];
    tree_B[r] += mimg->val[p][2];
    number_of_nodes[r] += 1;
    
    u = iftGetVoxelCoord(orig,p);
    for(i = 1; i < adjacency->n; i++){
      v = iftGetAdjacentVoxel(adjacency, u, i);
      if(iftValidVoxel(orig,v)){
        q = iftGetVoxelIndex(orig,v);

        if(Q->L.elem[q].color != IFT_BLACK){
          int WL = iftRound(sqrtf(
            powf((mimg->val[q][0]-tree_L[r]/number_of_nodes[r]), 2.0) +
            powf((mimg->val[q][1]-tree_A[r]/number_of_nodes[r]), 2.0) +
            powf((mimg->val[q][2]-tree_B[r]/number_of_nodes[r]), 2.0)
          ));

          tmp = iftMax(pathval->val[p], WL);

          if(tmp < pathval->val[q]){
            if(Q->L.elem[q].color == IFT_GRAY){
              iftRemoveGQueueElem(Q,q);
            }
            label->val[q] = label->val[p];
            root->val[q] = root->val[q];
            pathval->val[q] = tmp;
            iftInsertGQueue(&Q, q);

          }
        }
      }
    }

  }

  iftDestroyAdjRel(&adjacency);
  iftDestroyGQueue(&Q);
  iftDestroyImage(&pathval);
  return (label);
}

iftImage *ImageGradient(iftImage *img, iftAdjRel *A)
{
  iftImage *gradI = iftCreateImage(img->xsize,img->ysize,img->zsize);
  float    *mag   = iftAllocFloatArray(A->n); 
  float    *gx    = iftAllocFloatArray(3);
  float    *gy    = iftAllocFloatArray(3);
  float    *gz    = iftAllocFloatArray(3);
  
  for (int i=0; i < A->n; i++)
    mag[i]=sqrt(A->dx[i]*A->dx[i]+A->dy[i]*A->dy[i]+A->dz[i]*A->dz[i]);
  
  for (ulong p=0; p < img->n; p++){
    iftVoxel u  = iftGetVoxelCoord(img,p);

    for (int b=0; b < 3; b++){
      gx[b] = 0; gy[b] = 0; gz[b] = 0;
    }

    for (int i=1; i < A->n; i++) {
      iftVoxel v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(img,v)){
	int q = iftGetVoxelIndex(img,v);	    
	for (int b=0; b < 3; b++){
	  gx[b] += ((float)img->val[q]-(float)img->val[p])*A->dx[i]/mag[i];
	  gy[b] += ((float)img->Cb[q]-(float)img->Cb[p])*A->dy[i]/mag[i];
	  gz[b] += ((float)img->Cr[q]-(float)img->Cr[p])*A->dz[i]/mag[i];
	}
      }
    }
    float Gx=0.0, Gy=0.0, Gz=0.0;
    for (int b=0; b < 3; b++){
      gx[b] = gx[b] / (A->n-1);
      gy[b] = gy[b] / (A->n-1);
      gz[b] = gz[b] / (A->n-1);
      Gx += gx[b]; Gy += gy[b]; Gz += gz[b]; 
    }
    Gx /= 3; Gy /= 3; Gz /= 3;
    gradI->val[p] = iftRound(sqrtf(Gx*Gx + Gy*Gy + Gz*Gz));
  }

  iftFree(mag);
  iftFree(gx);
  iftFree(gy);
  iftFree(gz);
  
  return(gradI);
}

iftImage *Watershed(iftImage *gradI, iftImage *seeds_in, iftImage *seeds_out)
{
  iftImage   *pathval = NULL, *label = NULL;
  iftGQueue  *Q = NULL;
  int         i, p, q, tmp;
  iftVoxel    u, v;
  iftAdjRel     *A = iftCircular(1.0);

  // Initialization
  
  pathval    = iftCreateImage(gradI->xsize, gradI->ysize, gradI->zsize);
  label      = iftCreateImage(gradI->xsize, gradI->ysize, gradI->zsize);
  Q          = iftCreateGQueue(iftMaximumValue(gradI)+1, gradI->n, pathval->val);

  /* Initialize costs */
    
  for (p = 0; p < gradI->n; p++)
  {
    pathval->val[p] = IFT_INFINITY_INT;
    if (seeds_in->val[p] != 0){
      label->val[p]   = seeds_in->val[p];
      pathval->val[p] = 0;
    }else{
      if (seeds_out->val[p] != 0){
	label->val[p] = 0;
	pathval->val[p] = 0;
      }
    }    
    iftInsertGQueue(&Q, p);
  }      
  
  /* Propagate Optimum Paths by the Image Foresting Transform */

  while (!iftEmptyGQueue(Q))
  {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(gradI, p);

    for (i = 1; i < A->n; i++)
    {
      v = iftGetAdjacentVoxel(A, u, i);

      if (iftValidVoxel(gradI, v))
      {
        q = iftGetVoxelIndex(gradI, v);

	if (Q->L.elem[q].color != IFT_BLACK) {
	  
          tmp = iftMax(pathval->val[p], gradI->val[q]);

          if (tmp < pathval->val[q])  {
	    iftRemoveGQueueElem(Q,q);
            label->val[q]    = label->val[p];

            pathval->val[q]  = tmp;
	    iftInsertGQueue(&Q,q); 	    
	  }
        }
      }
    }
  }
  
  iftDestroyAdjRel(&A);
  iftDestroyGQueue(&Q);
  iftDestroyImage(&pathval);

  return (label);
}

//IoU avaliation also handles false positive cases, we work with two different scale bounding box
double IoU_avaliation(iftImage *orig, iftImage *comp, iftImage *gt, float scale){

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
		
		if((frequency->val[i]) > 0)area_union++;
		if((frequency->val[i]) > 1)area_intersection++;
	}
	
	return (area_intersection * 1.0) / area_union;
}

//DICE avaliation also handles false positive cases, we work with two different scale bounding box
double DICE_avaliation(iftImage *orig, iftImage *comp, iftImage *gt, float scale1, float scale2){

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
		
		if((frequency->val[i]) > 0)area_union += frequency->val[i];
		if((frequency->val[i]) > 1)area_intersection++;
	}
	
	return (2 * area_intersection * 1.0) / area_union;
}

int main(int argc, char *argv[])
{

  /* Example: delineation salie 1 objs */
  
  if (argc!=4){ 
    iftError("Usage: delineation <P1> <P2> <P3>\n"
	     "[1] folder with the salience maps\n"
	     "[2] layer (1,2,...) to create the results\n"
	     "[3] output folder with the resulting images\n",	 
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
  iftColorTable *ctb = iftCreateRandomColorTable(10);
  iftAdjRel *A       = iftCircular(3.0); 
  iftAdjRel *B       = iftCircular(1.5);      
  iftAdjRel *C       = iftCircular(1.0);

  float delineation_error = 0;
  
  for(int i = 0; i < fs->n; i++) {
    
    printf("Processing image %d of %ld\r", i + 1, fs->n);
    
    char *basename1   = iftFilename(fs->files[i]->path,suffix); // 000001_layer1.png - _layer1.png = 000001  
    char *basename2   = iftFilename(fs->files[i]->path,".png"); // 000001_layer1.png - png = 000001_layer1      
    iftImage *salie   = iftReadImageByExt(fs->files[i]->path);
    
    sprintf(filename,"./truelabels/%s.png",basename1);
    iftImage *gt      = iftReadImageByExt(filename);
    
    sprintf(filename,"./images/%s.png",basename1); //preparing the ith image file name
    iftImage *orig    = iftReadImageByExt(filename); // reading the ith image 
    
    /* Delineate parasite */

    iftImage *gradI = ImageGradient(orig,A);
    iftImage  *bin  = iftThreshold(salie,iftOtsu(salie),IFT_INFINITY_INT,255);
    iftImage  *seeds_in = iftSelectCompInAreaInterval(bin, NULL, 1000, 4472);
    iftDestroyImage(&bin);
    iftImage  *img      = NULL;

    if (iftMaximumValue(seeds_in)==0)
      img = iftCopyImage(orig);
    else{
      iftSet *S       = NULL;
      bin             = iftDilateBin(seeds_in,&S,15.0);
      iftDestroySet(&S);
      iftImage *seeds_out = iftComplement(bin);
      iftDestroyImage(&bin);
      bin             = iftFastLabelComp(seeds_in,NULL);
      iftDestroyImage(&seeds_in);
      seeds_in        = bin;

       iftImage *label = Watershed(gradI,seeds_in,seeds_out);
      //iftImage  *label  = DynamicTrees(orig,seeds_in,seeds_out);
      iftFImage *weight = iftSmoothWeightImage(gradI, 0.5);
      iftImage *smooth_label = iftFastSmoothObjects(label, weight, 5);
      iftDestroyImage(&label);
      label = smooth_label;
      iftDestroyFImage(&weight);
      img              = iftCopyImage(orig);
      
      //delineation_error += IoU_avaliation(img, label, gt, 1.6);
      delineation_error += DICE_avaliation(img, label, gt, 1.8, 1.5);
      
      iftDrawBorders(img, label, C, ctb->color[1], B);
      iftDestroyImage(&label);
      iftDestroyImage(&seeds_out);
    }

    iftDestroyImage(&seeds_in);
    iftDestroyImage(&salie);
    iftDestroyImage(&orig);
    iftDestroyImage(&gradI);
      
    /* save resulting image */

    sprintf(filename,"%s/%s.png",output_dir,basename2);
    iftWriteImageByExt(img,filename);

    iftDestroyImage(&img);
    iftFree(basename1);
    iftFree(basename2);
  }
  
  FILE *fp;
  fp = fopen("delineation_result.txt", "w");
  fprintf(fp, "%.4lf\n", delineation_error/(1.0 * (fs->n)));  
  fclose(fp);
  
  iftDestroyColorTable(&ctb);
  iftFree(filename);
  iftDestroyFileSet(&fs);
  iftDestroyAdjRel(&A);
  iftDestroyAdjRel(&B);
  iftDestroyAdjRel(&C);
  
  printf("\nDone ... %s\n", iftFormattedTime(iftCompTime(tstart, iftToc())));
  
  return (0);
}
