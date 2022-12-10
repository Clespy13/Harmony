#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

#ifndef IMAGE_PROCESSING_H
#define IMAGE_PROCESSING_H

Uint32 pixel_to_grayscale(Uint32 pixel_color, SDL_PixelFormat* format);
void surface_to_grayscale(SDL_Surface* surface);
Uint32 simple(Uint32 pixel_color, SDL_PixelFormat* format, Uint32 threshold);
void simple_binarize(SDL_Surface* surface, Uint32 threshold);
void histogram(SDL_Surface* surface, int* histo);
void equalized(SDL_Surface* surface, int* histo);
Uint8 Otsusmethod(SDL_Surface* surface, int* histo);
int hough_transform(SDL_Surface* src_surface);
int auto_rotate(SDL_Surface* src_surface);
void pixel_rotation(double x,double y,double angle,double* rotated_x,double* rotated_y,double center_image_x,double center_image_y);
void image_rotation(SDL_Surface* image,double angle);
double degrees_to_rad(double degrees);
void put_pixel(SDL_Surface *surface, unsigned x_position, unsigned y_position, Uint32 pixel);
Uint32 get_pixel(SDL_Surface *surface, unsigned x_position, unsigned y_position);
void scale(SDL_Surface* src_surface);
SDL_Surface* load_image(const char* path);
SDL_Surface* binarize_and_rotate(char** argv);

#endif
