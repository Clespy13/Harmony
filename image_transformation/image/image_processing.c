#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <math.h>
#include "image_processing.h"
#include "../../detection/processing.h"


SDL_Surface* load_image(const char* path)
{
    if (path == NULL)
        errx(EXIT_FAILURE, "%s", SDL_GetError());
    SDL_Surface* tmp = IMG_Load(path);
    SDL_Surface* surface = SDL_ConvertSurfaceFormat(tmp, SDL_PIXELFORMAT_RGB888, 0);
    SDL_FreeSurface(tmp);
    return surface;
}



SDL_Surface* binarize_and_rotate(char** argv)
{
  // - Create a surface from the colored image.
  SDL_Surface* surface = load_image(argv[1]);
  SDL_Surface* new_one = load_image(argv[1]);
  
  // - Create an array to implement the histogram of the grayscale image
  int histo[256];
  for (int i = 0; i < 256; i++)
    histo[i] = 0;
  
  // - Fill the histogram with values
  histogram(new_one, histo);
  
  // - Equalize the histogram to enhance contrasts
  equalized(new_one, histo);
  
  // - Calculate the threshold of the image used to binarize
  Uint8 threshold = Otsusmethod(new_one, histo);
  //printf("%u\n", threshold);
  
  // - Convert the surface into simple binarization
  surface_to_grayscale(surface);
  //surface = gaussianBlur(surface);
  simple_binarize(surface, threshold); //128 == threshold
  auto_rotate(surface);
  SDL_FreeSurface(new_one);
  return surface;
}

/*

int main(int argc, char** argv)
{
  if (argc != 2)
        errx(EXIT_FAILURE, "Usage: image-file");

  SDL_Surface* surface = binarize_and_rotate(argv);
  IMG_SavePNG(surface, "binarized_and_rotated.png");
  SDL_FreeSurface(surface);
  SDL_Quit();
  return 1;
}

*/
