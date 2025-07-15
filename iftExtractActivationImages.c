#include <ift.h>

int main(int argc, char *argv[]) {

    if (argc != 2) {
        iftError("usage: iftExtractActivationImages [layer_path]", "iftExtractActivationImages");
    }

    char *path = argv[1];
    iftFileSet *fset = iftLoadFileSetFromDir(path, 1);

    for (int i = 0; i < fset->n; i++) {
        char *file = fset->files[i]->path;
        iftMImage *mimg = iftReadMImage(file);
        for (int b = 0; b < mimg->m; b++) {
            iftImage *img = iftMImageToImage(mimg, 255, b);
            iftWriteImageByExt(img, "%s/%03d.png", iftBasename(file), b);
            iftDestroyImage(&img);
        }
        iftDestroyMImage(&mimg);
    }

    iftDestroyFileSet(&fset);
}
