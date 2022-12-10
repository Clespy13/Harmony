#include <stdio.h>
#include <err.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL_ttf.h>

#define WINDOW_SIZE 1152
#define BOX_SIZE (WINDOW_SIZE / 9)

SDL_Surface* draw(SDL_Surface* surface)
{
  SDL_Rect rect = {0, 0, WINDOW_SIZE, WINDOW_SIZE};
  SDL_LockSurface(surface);
  SDL_FillRect(surface, &rect, SDL_MapRGB(surface->format, 255, 255, 255));
  for (int i = 1; i < 9; i++)
    {
      if (i % 3 == 0)
	{
	  rect.x = BOX_SIZE * i - 2;
	  rect.w = 4;
	}
      else
	{
	  rect.x = BOX_SIZE * i - 1;
	  rect.w = 2;
	}
      SDL_FillRect(surface, &rect, SDL_MapRGB(surface->format, 0, 0, 0));
    }
  rect.w = WINDOW_SIZE;
  rect.x = 0;
  for (int i = 1; i < 9; i++)
    {
      if (i % 3 == 0)
	{
	  rect.y = BOX_SIZE * i - 2;
	  rect.h = 4;
	}
      else
	{
	  rect.y = BOX_SIZE * i - 1;
	  rect.h = 2;
	}
      SDL_FillRect(surface, &rect, SDL_MapRGB(surface->format, 0, 0, 0));
    }
  
  SDL_UnlockSurface(surface);
  return surface;
}


int main(int argc, char** argv)
{
  if (argc > 2)
    errx(1, "Arguments invalides");

  //Accessing the wanted file
  FILE* unsolved = NULL;
  unsolved = fopen(argv[1], "r");

  //If the file is not defined, error
  if (unsolved == NULL)
    errx(1, "La grille donnée n'existe pas");

  //TTF_Init();
  
  if (TTF_Init() == -1)
    {
      errx(1, "Erreur d'initialisation de TTF_Init %s", TTF_GetError());
    }

  //printf("%s\n", SDL_GetBasePath());
  TTF_Font* font = NULL;
  font = TTF_OpenFont("arial.ttf", 25);
  //if (font == 0)
  //  errx(1, "%s", TTF_GetError());
  SDL_Color color = {0, 0, 0, 0};
  SDL_Surface* text = TTF_RenderText_Solid(font, "Hello World!", color);
  SDL_Surface* surface = SDL_CreateRGBSurface(0, WINDOW_SIZE, WINDOW_SIZE, 32, 0x00ff0000, 0x0000ff00, 0x000000ff, 0xff000000);
  surface = draw(surface);
  SDL_BlitSurface(text, NULL, surface, NULL);
  IMG_SavePNG(surface, "white_screen.png");
  //IMG_SavePNG(text, "text.png");
  TTF_CloseFont(font);
  SDL_FreeSurface(text);
  SDL_FreeSurface(surface);
  TTF_Quit();
  SDL_Quit();
  
}
