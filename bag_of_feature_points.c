#include "ift.h"

/* 
   Author: Alexandre Xavier Falcão (September 10th, 2023)

   Description: This code assumes that labeled markers have been drawn on
   a few images and saved in a separated folder. You should draw
   markers with label 1 for object and label 0 for background. You may
   use FLIMBuilder in libmo445/demo/Qt/FLIMBuilder/FLIMBuilder to do
   that. The code extracts 3x3 patches from marker pixels, computes N
   clusters and selects N feature points (cluster centers) per
   marker. The number N is an input parameter. The code also adds a
   marker id to each marker pixel to facilitate the process. 

*/


void AddMarkerLabels(iftImage *img, iftLabeledSet *S) {
    iftLabeledSet *seed    = S;
    iftImage      *markers = iftCreateImage(img->xsize, img->ysize, img->zsize);
    iftAdjRel     *A = NULL;

    if (iftIs3DImage(img))
      A = iftSpheric(sqrtf(3.0));
    else
      A = iftCircular(sqrtf(2.0));

    while (seed != NULL) {
      int p = seed->elem;
      markers->val[p] = 255;
      seed = seed->next;
    }

    iftImage *lbmarkers = iftLabelComp(markers, A);
    iftDestroyAdjRel(&A);
    iftDestroyImage(&markers);

    seed = S;
    while (seed != NULL) {
        int p = seed->elem;
	if (lbmarkers->val[p] != 0)
	  seed->marker = lbmarkers->val[p];
        seed = seed->next;
    }

    iftDestroyImage(&lbmarkers);
}

/* It creates a dataset in which voxels of a same marker are assigned
   to the same truelabel -- the marker id. Such samples are
   represented by 3x3 patches centered at the marker voxels. */

iftDataSet *ExtractMarkerDataSet(iftMImage *img, iftLabeledSet *S, iftAdjRel *A) {
  int nsamples    = iftLabeledSetSize(S);
  int tensor_size = img->m * A->n;
  iftDataSet *Z = iftCreateDataSet(nsamples, tensor_size);
  int ninput_channels = img->m;

  Z->nclasses = 0;
  int s       = 0;
  iftLabeledSet *seed = S;
  while (seed != NULL) {
    int p = seed->elem;
    Z->sample[s].id = seed->elem;
    Z->sample[s].truelabel = seed->marker;
    if (Z->sample[s].truelabel > Z->nclasses)
      Z->nclasses = Z->sample[s].truelabel;
    iftVoxel u = iftMGetVoxelCoord(img, p);
    int j = 0;
    for (int k = 0; k < A->n; k++) {
      iftVoxel v = iftGetAdjacentVoxel(A, u, k);
      if (iftMValidVoxel(img, v)) {
	int q = iftMGetVoxelIndex(img, v);
	for (int b = 0; b < ninput_channels; b++) {
	  Z->sample[s].feat[j] = img->val[q][b];
	  j++;
	}
      } else {
	for (int b = 0; b < img->m; b++) {
	  Z->sample[s].feat[j] = 0;
	  j++;
	}
      }
    }
    s++;
    seed = seed->next;
  }
  
  iftSetStatus(Z, IFT_TRAIN);
  iftAddStatus(Z, IFT_SUPERVISED);
  
  return (Z);
}

iftLabeledSet *SelectClusterCenters(iftImage *img, iftLabeledSet *S, float adj_radius, int nclusters_per_marker)
{
  iftImage      *markers = iftCreateImage(img->xsize, img->ysize, img->zsize);
  iftImage      *labels  = iftCreateImage(img->xsize, img->ysize, img->zsize);
  iftLabeledSet *seed    = S, *T = NULL;

  /* Store marker id and label of each marker voxel */
  while (seed != NULL) {
    markers->val[seed->elem]=seed->marker;
    labels->val[seed->elem]=seed->label;
    seed = seed->next;
  }

  /* Extract the marker dataset: the marker id is the truelabel of the
     sample */
  
  iftMImage *mimg = iftImageToMImage(img, LABNorm2_CSPACE);
  iftAdjRel *A    = NULL;
  if (iftIs3DImage(img)){
    A = iftSpheric(adj_radius);
  }else{
    A = iftCircular(adj_radius);
  }  
  iftDataSet *Z = ExtractMarkerDataSet(mimg,S,A);

  /* Estimate the given number of cluster centers per marker */
  
  for (int c = 1; c <= Z->nclasses; c++) {
    iftDataSet *Z1  = iftExtractSamplesFromClass(Z, c);
    int ngroups = iftMin(Z1->nsamples, nclusters_per_marker);
    iftKMeans(Z1, ngroups, 100, 0.001, NULL, NULL, iftEuclideanDistance);
    for (int s = 0; s < Z1->nsamples; s++) {
      if (iftHasSampleStatus(Z1->sample[s], IFT_PROTOTYPE)) {
	int p = Z1->sample[s].id;
	iftInsertLabeledSetMarkerAndHandicap(&T,p,
					     labels->val[p],
					     markers->val[p],
					     0);
      }
    }
    iftDestroyDataSet(&Z1);
  }
  
  iftDestroyImage(&markers);
  iftDestroyImage(&labels);
  iftDestroyAdjRel(&A);
  iftDestroyMImage(&mimg);
  iftDestroyDataSet(&Z);
  
  return(T);
}

int main(int argc, char *argv[]) {
    timer *tstart;

    /* Example: 

       bag_of_feature_points filtered_images markers 1 bag 


       If you use disks as markers, 1 feature point per
       marker is fine. However, if you draw strokes, try from 3-5
       feature points per marker. Recall that each feature point will
       create one kernel. 
    */ 
    
    if (argc!=5)
      iftError("Usage: bag_of_feature_points <P1> <P2> <P3> <P4>\n"
	       "P1: image folder \n"
	       "P2: marker folder \n"
	       "P3: number of feature points per marker \n"
	       "P4: new marker folder, representing a bag of feature points\n",
	       "main");

    tstart = iftTic();

    /* Read image filenames and get the files' extension */
    
    iftFileSet *fs = iftLoadFileSetFromDir(argv[1], 0);
    char ext[10];
    sprintf(ext,"%s",iftFileExt(fs->files[0]->path));
    iftMakeDir("layer0");
    iftMakeDir(argv[4]);
    char *filename = iftAllocCharArray(512);
    int nfpts_per_marker = atoi(argv[3]);
    
    /* Convert images to multichannel images and save them as layer0 */
    
    for (int i = 0; i < fs->n; i++) {
      char *basename  = iftFilename(fs->files[i]->path, ext);
      //printf("dir: %s\n", fs->files[i]->path);
      iftImage *img   = iftReadImageByExt(fs->files[i]->path);
      iftMImage *mimg = NULL;	
      mimg = iftImageToMImage(img, LABNorm2_CSPACE); /* even grayscale
							images are
							converted to
							three-channel
							images */

      sprintf(filename,"layer0/%s.mimg",basename);
      iftWriteMImage(mimg,filename);
      iftFree(basename);
      iftDestroyImage(&img);
      iftDestroyMImage(&mimg);
    }
    iftDestroyFileSet(&fs);

    /* Select points from markers and save them as a bag of feature
       points */
    
    fs = iftLoadFileSetFromDirBySuffix(argv[2],"-seeds.txt", 1);  
    
    for (int i=0; i < fs->n; i++) {
      char *basename   = iftFilename(fs->files[i]->path,"-seeds.txt");
      sprintf(filename,"%s/%s%s",argv[1],basename,ext);
      printf("filename >>>>>> %s\n", filename);
      iftImage *img    = iftReadImageByExt(filename);
      iftLabeledSet *S = iftReadSeeds(img,fs->files[i]->path);
      AddMarkerLabels(img, S);
      float radius;
      if (iftIs3DImage(img))
	radius=1.75;
      else
	radius=1.42; 
      
      iftLabeledSet *T = SelectClusterCenters(img, S, radius, nfpts_per_marker);
      iftDestroyLabeledSet(&S);
      S = T;
      sprintf(filename,"%s/%s-fpts.txt",argv[4],basename);
      iftWriteLabeledSet(S,img->xsize,img->ysize,img->zsize,filename);
      iftFree(basename);
      iftDestroyImage(&img);
      iftDestroyLabeledSet(&S);
    }

    iftFree(filename);    
    iftDestroyFileSet(&fs);
    
    printf("Done ... %s\n", iftFormattedTime(iftCompTime(tstart, iftToc())));
    return (0);
}
