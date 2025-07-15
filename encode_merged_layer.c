#include "ift.h"

/* Author: Alexandre Xavier Falcão (September, 10th 2023) 

   Description: Executes a convolutional block to encode the current
   layer using the consolidated model of all training images.

*/


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

void LoadMergedModel(char *basepath, iftMatrix **K, float **bias)
{
  char filename[200];
  FILE *fp;
  int nkernels;
  
  sprintf(filename, "%s-kernels.npy", basepath);
  *K = iftReadMatrix(filename);
  sprintf(filename, "%s-bias.txt", basepath);
  fp = fopen(filename, "r");
  fscanf(fp, "%d", &nkernels);
  *bias = iftAllocFloatArray(nkernels);
  for (int k = 0; k < nkernels; k++) {
    fscanf(fp, "%f ", &((*bias)[k]));
  }  
  fclose(fp);
}

int main(int argc, char *argv[]) {
    timer *tstart;

    /* Example: encode_merged_layer arch.json 1 flim_models */
    
    if (argc!=4)
      iftError("Usage: encode_merged_layer <P1> <P2> <P3> \n"
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
    
    iftMatrix *Kmerged=NULL;
    float     *bias_merged=NULL;

    sprintf(filename,"%s/conv%d",model_dir,layer);
    LoadMergedModel(filename,&Kmerged,&bias_merged);
    
    /* Encode layer with merged model */      
	
    for (int i=0; i < fs->n; i++) {
      iftMImage *mimg  = iftReadMImage(fs->files[i]->path);
      char *basename   = iftFilename(fs->files[i]->path, ".mimg");
      iftAdjRel *A     = GetPatchAdjacency(mimg, arch->layer[layer-1]);
      
      /* Complete the code below to compute convolution with a kernel
	 bank followed by bias */
	
      iftMatrix *XI = iftMImageToFeatureMatrix(mimg, A, NULL);
      iftMatrix *XJ = iftMultMatrices(XI, Kmerged);
      iftMImage *activ = iftMatrixToMImage(XJ, mimg->xsize, mimg->ysize, mimg->zsize, Kmerged->ncols, 'c');
      iftDestroyMatrix(&XI);
      iftDestroyMatrix(&XJ);
      iftDestroyAdjRel(&A);
      
      /* ReLU in place */
	if(arch->layer[layer-1].relu){
		for (int p = 0; p < activ->n; p++){
			for (int b = 0; b < activ->m; b++){
				activ->val[p][b] += bias_merged[b];
				if(activ->val[p][b] < 0)
					activ->val[p][b] = 0; 
			}
		}
	}
	else{
		for (int p = 0; p < activ->n; p++){
			for (int b = 0; b < activ->m; b++){
				activ->val[p][b] += bias_merged[b];
			}
		}
	
	}
      
      /* pooling */
      
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
      iftDestroyMImage(&activ);
      iftFree(basename);
      iftDestroyMImage(&mimg);
    }

    iftDestroyMatrix(&Kmerged);
    iftFree(bias_merged);
    iftFree(filename);
    iftDestroyFileSet(&fs);
    iftDestroyFLIMArch(&arch);
    
    printf("Done ... %s\n", iftFormattedTime(iftCompTime(tstart, iftToc())));
    return (0);
}
