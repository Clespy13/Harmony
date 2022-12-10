#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>

#ifndef IMAGE_ROTATION_H
#define IMAGE_ROTATION_H

void pixel_rotation(double x,double y,double angle,double* rotated_x,double* rotated_y,double center_image_x,double center_image_y);

void image_rotation(SDL_Surface* image,double angle);

double degrees_to_rad(double degrees);

void put_pixel(SDL_Surface *surface, unsigned x_position, unsigned y_position, Uint32 pixel);

Uint32 get_pixel(SDL_Surface *surface, unsigned x_position, unsigned y_position);

#endif
