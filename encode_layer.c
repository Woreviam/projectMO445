#include "ift.h"

/* Author: Alexandre Xavier Falcão (September, 10th 2023) 

   Description: Executes a convolutional block to encode the current
   layer using the model of each training image. 

*/

float *LoadBias(char *basepath)
{
  int number_of_kernels;
  char filename[200];
  FILE *fp;
  float *bias;
  
  sprintf(filename, "%s-bias.txt", basepath);
  fp = fopen(filename, "r");
  fscanf(fp, "%d", &number_of_kernels);
  bias = iftAllocFloatArray(number_of_kernels);
  for (int k = 0; k < number_of_kernels; k++) {
    fscanf(fp, "%f ", &bias[k]);
  }
  fclose(fp);

  return(bias);
}

iftAdjRel *GetPatchAdjacency(iftMImage *mimg, iftFLIMLayer layer)
{
  iftAdjRel *A;

  if (iftIs3DMImage(mimg)){
    A = iftCuboidWithDilationForConv(layer.kernel_size[0],
				     layer.kernel_size[1],
				     layer.kernel_size[2],
				     layer.dilation_rate[0],
				     layer.dilation_rate[1],
				     layer.dilation_rate[2]);
  }else{
    A = iftRectangularWithDilationForConv(layer.kernel_size[0],
					  layer.kernel_size[1],
					  layer.dilation_rate[0],
					  layer.dilation_rate[1]);    
  }

  return(A);
}

int *LoadTruelabels(char *filename)
{
  int    number_of_kernels;
  FILE  *fp;
  int   *truelabel;
  
  fp = fopen(filename, "r");
  fscanf(fp, "%d", &number_of_kernels);
  truelabel = iftAllocIntArray(number_of_kernels);
  for (int k = 0; k < number_of_kernels; k++) {
    fscanf(fp, "%d ", &truelabel[k]);
  }
  fclose(fp);

  return(truelabel);
}

int main(int argc, char *argv[]) {
    timer *tstart;

    /* Example: encode_layer arch.json 1 flim_models */
    
    if (argc!=4)
      iftError("Usage: encode_layer <P1> <P2> <P3>\n"
	       "[1] architecture of the network (.json) \n"
	       "[2] layer number (1, 2, 3) \n"
	       "[3] folder with the models \n",
	       "main");

    tstart = iftTic();

    iftFLIMArch *arch   = iftReadFLIMArch(argv[1]);
    int          layer  = atoi(argv[2]);
    char    *model_dir  = argv[3];
    char    *filename   = iftAllocCharArray(512);
    char     input_dir[20], output_dir[20]; 
    
    sprintf(input_dir,"layer%d",layer-1);
    sprintf(output_dir,"layer%d",layer);
    iftMakeDir(output_dir);

    iftFileSet *fs = iftLoadFileSetFromDirBySuffix(input_dir, ".mimg", true);

    for (int i=0; i < fs->n; i++) {
      iftMImage *mimg  = iftReadMImage(fs->files[i]->path);
      char *basename   = iftFilename(fs->files[i]->path, ".mimg");
      sprintf(filename,"%s/%s-conv%d-kernels.npy",model_dir,basename,layer);

      if (iftFileExists(filename)){ /* encode layer using its kernels */
	iftMatrix *K     = iftReadMatrix(filename);	       
	sprintf(filename,"%s/%s-conv%d",model_dir,basename,layer);
	float     *bias  = LoadBias(filename);
	iftAdjRel *A     = GetPatchAdjacency(mimg, arch->layer[layer-1]);
	
	/* Complete the code below to compute convolution with a kernel
	   bank followed by bias */
	iftMatrix *XI = iftMImageToFeatureMatrix(mimg, A, NULL);
	iftMatrix *XJ = iftMultMatrices(XI, K);
	
	
	iftMImage *activ = iftMatrixToMImage(XJ, mimg->xsize, mimg->ysize, mimg->zsize, K->ncols, 'c');
	
	iftDestroyMatrix(&XI);
	iftDestroyMatrix(&XJ);
	iftDestroyAdjRel(&A);
	
	iftSetMImage(activ, 0);
	
	/* ReLU in place */
	if(arch->layer[layer-1].relu){
		for (int p = 0; p < activ->n; p++){
			for (int b = 0; b < activ->m; b++){
				activ->val[p][b] += bias[b];
				if(activ->val[p][b] < 0)
					activ->val[p][b] = 0;
			}
		}
	}
	else{
		for (int p = 0; p < activ->n; p++){
			for (int b = 0; b < activ->m; b++){
				activ->val[p][b] += bias[b];
			}
		}
	
	}
	
	/* Pooling */
	
	if (strcmp(arch->layer[layer-1].pool_type, "no_pool") != 0){
	  iftMImage *pool = NULL;
	  if (strcmp(arch->layer[layer-1].pool_type, "avg_pool") == 0) {
	    pool = iftFLIMAtrousAveragePooling(activ,
					       arch->layer[layer-1].pool_size[0],
					       arch->layer[layer-1].pool_size[1],
					       arch->layer[layer-1].pool_size[2],
					       1,
					       arch->layer[layer-1].pool_stride);
	    iftDestroyMImage(&activ);
	    activ = pool;
	  } else {
	    if (strcmp(arch->layer[layer-1].pool_type, "max_pool") == 0) { 
	      pool = iftFLIMAtrousMaxPooling(activ,
					     arch->layer[layer-1].pool_size[0],
					     arch->layer[layer-1].pool_size[1],
					     arch->layer[layer-1].pool_size[2],
					     1,
					     arch->layer[layer-1].pool_stride);
	      iftDestroyMImage(&activ);
	      activ = pool;
	    } else {
	      iftError("Invalid pooling in layer %d","main",layer);
	    }
	  }
	}
	
	sprintf(filename,"%s/%s.mimg",output_dir,basename);
	iftWriteMImage(activ,filename);
	iftDestroyMatrix(&K);
	iftDestroyMImage(&activ);
	iftFree(bias);
      }
      
      iftFree(basename);
      iftDestroyMImage(&mimg);
    }
    
    
    iftFree(filename);
    iftDestroyFileSet(&fs);
    iftDestroyFLIMArch(&arch);
    
    printf("Done ... %s\n", iftFormattedTime(iftCompTime(tstart, iftToc())));
    return (0);
}
