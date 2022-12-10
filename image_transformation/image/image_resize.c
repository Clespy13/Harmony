#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
//#include "image_rotation.h"
//#include "image_resize.h"
#include "image_processing.h"


#define getByte(value, n) (value >> (n*8) & 0xFF)

float max(float a, float b) { return (a < b) ? a : b; };

float lerp(float s, float e, float t){return s+(e-s)*t;}

float blerp(float c00, float c10, float c01, float c11, float tx, float ty)
{
    return lerp(lerp(c00, c10, tx), lerp(c01, c11, tx), ty);
}

void scale(SDL_Surface*src_surface)
{
    SDL_Surface* dest_surface;
    int newWidth =28;
    int newHeight=28;

	dest_surface = SDL_CreateRGBSurface(SDL_SWSURFACE, newWidth, newHeight,
		src_surface->format->BitsPerPixel, src_surface->format->Rmask,
		src_surface->format->Gmask, src_surface->format->Bmask,
		src_surface->format->Amask);

	SDL_FillRect(dest_surface,
		NULL, SDL_MapRGB(dest_surface->format, 255, 255, 255));

	if (dest_surface == NULL)
	{
		errx(EXIT_FAILURE,"%s", SDL_GetError());
	}

	SDL_LockSurface(src_surface);
	SDL_LockSurface(dest_surface);

    int x, y;
    for(x= 0, y=0; y < newHeight; x++)
    {
        if(x > newWidth)
        {
            x = 0; 
            y++;
        }
        // Image should be clamped at the edges and not scaled.
        float gx = max(x / (float)(newWidth) * (src_surface->w) - 0.5f, src_surface->w - 1);
        float gy = max(y / (float)(newHeight) * (src_surface->h) - 0.5, src_surface->h - 1);
        int gxi = (int)gx;
        int gyi = (int)gy;
        uint32_t result=0;
        uint32_t c00 = get_pixel(src_surface, gxi, gyi);
        uint32_t c10 = get_pixel(src_surface, gxi+1, gyi);
        uint32_t c01 = get_pixel(src_surface, gxi, gyi+1);
        uint32_t c11 = get_pixel(src_surface, gxi+1, gyi+1);
        uint8_t i;
        for(i = 0; i < 3; i++)
        {
            //((uint8_t*)&result)[i] = blerp( ((uint8_t*)&c00)[i], ((uint8_t*)&c10)[i], ((uint8_t*)&c01)[i], ((uint8_t*)&c11)[i], gxi - gx, gyi - gy); // this is shady
            result |= (uint8_t)blerp(getByte(c00, i), getByte(c10, i), getByte(c01, i), getByte(c11, i), gx - gxi, gy -gyi) << (8*i);
        }
        put_pixel(dest_surface,x, y, result);
    }
    SDL_UnlockSurface(src_surface);
	SDL_UnlockSurface(dest_surface);
    *src_surface = *dest_surface;
    free(dest_surface);
}
