#include <err.h>
#include "image_transformation/image/image_processing.h"
#include "detection/processing.h"
#include <sys/stat.h>
#include <sys/types.h>

int main(int argc, char** argv)
{
  if (argc != 2)
        errx(EXIT_FAILURE, "Usage: image-file");

  SDL_Surface* surface = binarize_and_rotate(argv);
  IMG_SavePNG(surface, "binarized.png");
  // - Create a surface from the colored image.
  /*
  SDL_Surface* surface = load_image(argv[1]);
  //SDL_Surface* new_one = load_image(argv[1]);
  
  // - Create an array to implement the histogram of the grayscale image
  int histo[256];
  for (int i = 0; i < 256; i++)
    histo[i] = 0;
  
  // - Convert the surface into simple binarization
  //surface_to_grayscale(surface);


  IMG_SavePNG(surface, "Gauss.png");
  */

  surface = sobel(surface);

  IMG_SavePNG(surface, "sobel.png");

  surface = dilate(surface);
  surface = dilate(surface);
  IMG_SavePNG(surface, "dilated.png");

  surface = erode(surface);
  IMG_SavePNG(surface, "eroded.png");



  /*
  SDL_Surface* new_one = SDL_ConvertSurface(surface, surface->format, surface->flags);
  
  // - Fill the histogram with values
  histogram(new_one, histo);
  
  // - Equalize the histogram to enhance contrasts
  equalized(new_one, histo);

  // - Calculate the threshold of the image used to binarize
  Uint8 threshold = Otsusmethod(new_one, histo);
  //printf("%u\n", threshold);
  
  simple_binarize(surface, threshold); //128 == threshold
  */
  
  //auto_rotate(surface);
  //IMG_SavePNG(surface, "rotated.png");

  //componentLabeling(surface);
  //IMG_SavePNG(surface, "component.png");

  houghLines(surface);
  //IMG_SavePNG(surface, "hough2.png");
  //auto_rotate(surface);
  //SDL_FreeSurface(new_one);

  int res = mkdir("cases", 0777);
  if (res == 1)
    errx(1, "could not create : 'cases/'");

  cut(surface);
  IMG_SavePNG(surface, "last_process.png");
  SDL_FreeSurface(surface);
  SDL_Quit();
  return 1;
}
