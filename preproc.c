#include "ift.h"

/* Author: Alexandre Xavier Falcão (September 10th 2023) 

   Description: Preprocess the original images by a median filter to
   remove noise

*/

int main(int argc, char *argv[])
{

  /* Example: preproc images 1.5 filtered_images */


  if (argc!=4){ 
    iftError("Usage: preproc <P1> <P2> <P3>\n"
	     "[1] folder with the original images\n"
	     "[2] adjacency radius of the filter\n"
	     "[3] folder with the filtered images\n",	 
	     "main");
  }
  
  timer *tstart = iftTic();

  iftFileSet *fs     = iftLoadFileSetFromDirBySuffix(argv[1],".png", true);
  char *filename     = iftAllocCharArray(512);
  float adj_radius   = atof(argv[2]);
  iftAdjRel  *A      = iftCircular(adj_radius);
  char *output_dir   = argv[3];
  iftMakeDir(output_dir);

  for(int i = 0; i < fs->n; i++) {
    printf("Processing image %d of %ld\r", i + 1, fs->n);
    char *basename    = iftFilename(fs->files[i]->path,".png");      
    sprintf(filename,"%s/%s.png",output_dir,basename);
    iftImage *orig    = iftReadImageByExt(fs->files[i]->path);
    iftImage *median  = iftMedianFilter(orig,A);
    iftDestroyImage(&orig);
    iftWriteImageByExt(median,filename);
    iftDestroyImage(&median);
    iftFree(basename);
  }

  iftFree(filename);
  iftDestroyFileSet(&fs);
  iftDestroyAdjRel(&A);
  
  printf("\nDone ... %s\n", iftFormattedTime(iftCompTime(tstart, iftToc())));
  
  return (0);
}
