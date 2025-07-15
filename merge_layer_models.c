#include "ift.h"

/* Author: Alexandre Xavier Falcão (September 10th, 2023)

   After you have created one model per training image, this code
   merges the kernels, biases, and weights from all models. The
   consolidated model is saved in the same folder of all models.

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

float *LoadKernelWeights(char *filename)
{
  int    number_of_kernels;
  FILE  *fp;
  float *weight;
  
  fp = fopen(filename, "r");
  fscanf(fp, "%d", &number_of_kernels);
  weight = iftAllocFloatArray(number_of_kernels);
  for (int k = 0; k < number_of_kernels; k++) {
    fscanf(fp, "%f ", &weight[k]);
  }
  fclose(fp);

  return(weight);
}

void SaveMergedData(char *basepath, iftMatrix *K, float *bias, float *weight)
{
  char filename[200];
  FILE *fp;

  sprintf(filename, "%s-kernels.npy", basepath);
  iftWriteMatrix(K,filename);
  sprintf(filename, "%s-bias.txt", basepath);
  fp = fopen(filename, "w");
  fprintf(fp, "%d\n", K->ncols);
  for (int k = 0; k < K->ncols; k++) {
    fprintf(fp, "%f ", bias[k]);
  }  
  fclose(fp);
  sprintf(filename, "%s-weights.txt", basepath);
  fp = fopen(filename, "w");
  fprintf(fp, "%d\n", K->ncols);
  for (int k = 0; k < K->ncols; k++) {
    fprintf(fp, "%f ", weight[k]);
  }  
  fclose(fp);
}

int main(int argc, char *argv[]) {
    timer *tstart;

    /* Example: merge_layer_models arch.json 1 flim_models */

    if (argc!=4)
      iftError("Usage: merge_layer_models <P1> <P2> <P3> \n"
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
    float     *weight_merged=NULL;

    /* Compute the total number of kernels and the feature vector
       size */      

    int ncols=0, nrows=0;
    for (int i=0; i < fs->n; i++) {
      char *basename   = iftFilename(fs->files[i]->path, ".mimg");
      sprintf(filename,"%s/%s-conv%d-kernels.npy",model_dir,basename,layer);
      if (iftFileExists(filename)){ /* encode layer using its kernels */
	iftMatrix *K   = iftReadMatrix(filename);	       
	ncols         += K->ncols;
	if (nrows == 0){
	  nrows        = K->nrows;
	}else{
	  if (nrows != K->nrows){	      
	    iftWarning("Cannot merge kernels with different weight vector sizes. Adopting one kernel per marked image.","main");
	  }
	}
      }
	iftFree(basename);
      }

    /* Merge data */
      
    Kmerged          = iftCreateMatrix(ncols,nrows);
    bias_merged      = iftAllocFloatArray(ncols);
    weight_merged    = iftAllocFloatArray(ncols);
    int  col         = 0;
    for (int i=0; i < fs->n; i++) {
      char *basename   = iftFilename(fs->files[i]->path, ".mimg");
      sprintf(filename,"%s/%s-conv%d-kernels.npy",model_dir,basename,layer);
      if (iftFileExists(filename)){ /* encode layer using its kernels */
	iftMatrix *K   = iftReadMatrix(filename);	       
	sprintf(filename,"%s/%s-conv%d",model_dir,basename,layer);
	float     *bias  = LoadBias(filename);
	sprintf(filename,"%s/%s-conv%d-weights.txt",
		model_dir,basename,layer);
	float    *weight = LoadKernelWeights(filename);
	
	for (int c=0; c < K->ncols; c++, col++) { /* copy kernel, bias and
						     weights */
	  bias_merged[col]   = bias[c];
	  weight_merged[col] = weight[c];
	  for (int r=0; r < K->nrows; r++) {
	    iftMatrixElem(Kmerged, col, r) = iftMatrixElem(K, c, r);
	  }
	}
	iftDestroyMatrix(&K);
	iftFree(bias);
	iftFree(weight);
      }
      iftFree(basename);
    }

    sprintf(filename,"%s/conv%d",model_dir,layer);
    SaveMergedData(filename,Kmerged,bias_merged,weight_merged);
      
    iftDestroyMatrix(&Kmerged);
    iftFree(bias_merged);
    iftFree(weight_merged);    
    iftFree(filename);
    iftDestroyFileSet(&fs);
    iftDestroyFLIMArch(&arch);
    
    printf("Done ... %s\n", iftFormattedTime(iftCompTime(tstart, iftToc())));
    return (0);
}
