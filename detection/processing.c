#include "processing.h"
#include "../image_transformation/image/image_processing.h"

/// @brief Puts the pixel at pos x,y in the surface
/// @param surface the surface to put the pixel in
/// @param x the x coordinate to put the pixel in
/// @param y the y coordinate to put the pixel in
/// @param pixel the pixel to put in the surface
void putPixel(SDL_Surface* surface, int x, int y, Uint32 pixel) {
	Uint32* pixels = surface->pixels;
	if (y*surface->pitch/4+x > 0 && y*surface->pitch/4+x < surface->w*surface->h)
		pixels[ ( y * surface->pitch/4 ) + x ] = pixel;
}

/// @brief Retrieve the pixel at position x,y in surface
/// @param surface the surface we want to retrieve the pixel from
/// @param x the x coordinate where we want to get the pixel
/// @param y the y coordinate where we want to get the pixel
/// @return the pixel at the x,y coordinate
Uint32 getPixel(SDL_Surface* surface, int x, int y) {
	Uint32* pixels = surface->pixels;
	if (y*surface->pitch/4+x > 0 && y*surface->pitch/4+x < surface->w*surface->h)
		return pixels[ ( y * surface->pitch/4 ) + x ];
	return 0;
}

/*

/// @brief Transforms a pixel into it's grayscaled value
/// @param pixelColor the pixel to modify
/// @param format the format of the pixel
/// @return the pixel grayscaled
Uint32 pixelToGrayscale(Uint32 pixelColor, SDL_PixelFormat* format) {
	Uint8 r, g, b;
	SDL_GetRGB(pixelColor, format, &r, &g, &b);
	Uint8 average = 0.3*r + 0.59*g + 0.11*b;
	return SDL_MapRGB(format, average, average, average);
}

/// @brief Converts an entire surface to grayscale
/// @param surface the surface to transform
void surfaceToGrayscale(SDL_Surface* surface) {
	Uint32* pixels = surface->pixels;
	
	int len = surface->w * surface->h;
	
	SDL_PixelFormat* format = surface->format;
	SDL_LockSurface(surface);
	
	for (int i = 0; i < len; i++)
		pixels[i] = pixelToGrayscale(pixels[i], format);
	
	SDL_UnlockSurface(surface);
}

*/

/// @brief Gaussian Blur on Image using a kernel of 5x5 to smooth out the image for
/// better binarization later on
/// @param surface the surface we want to edit
/// @return The modified surface blurred out 
SDL_Surface* gaussianBlur(SDL_Surface* surface) {
	SDL_PixelFormat* format = surface->format;
	SDL_Surface* tempSurface = SDL_ConvertSurface(surface, format, surface->flags);

	int w = surface->w;
	int h = surface->h;

	int kernel[] = {1,  4,  7,  4, 1, 
					4, 16, 26, 16, 4,
					7, 26, 41, 26, 7,
					4, 16, 26, 16, 4,
					1,  4,  7,  4, 1,};

	int gaussWidth = 5;
	int halfGaussWidth = gaussWidth / 2;
	int gaussSum = 0, c = 0;
	for (int i = 0; i < gaussWidth * gaussWidth; i++)
		gaussSum += kernel[i];

	Uint8 sumr, sumg, sumb;
	Uint8 r, g, b;
	int flag = 0;
	
	SDL_LockSurface(tempSurface);
	SDL_LockSurface(surface);

	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			sumr = 0, sumg = 0, sumb = 0, flag = 0, c = 0;
			for (int i = -halfGaussWidth; i < halfGaussWidth+1; i++)
			{
				for (int j = -halfGaussWidth; j < halfGaussWidth+1; j++)
				{
					if (getPixel(surface, x+j, y+i)) {
						flag = 1;
						SDL_GetRGB(getPixel(surface, x+j, y+i), format, &r, &g, &b);
						sumr += r * kernel[c] / gaussSum;
						sumg += g * kernel[c] / gaussSum;
						sumb += b * kernel[c] / gaussSum;
					}
					c++;
				}
			}

			if (flag) putPixel(tempSurface, x, y, SDL_MapRGB(format, sumr, sumg, sumb));
		}		
		
	}

	SDL_UnlockSurface(tempSurface);
	SDL_UnlockSurface(surface);

	return tempSurface;
	
}

void DrawLine(SDL_Surface *img, int x0, int y0, int x1, int y1, float wd, Uint32 pixel_color) {
	int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
	int dy = abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
	int err = dx - dy, e2, x2, y2;
	float ed = dx + dy == 0 ? 1 : sqrt((float) dx * dx + (float) dy * dy);


	Uint32 pixel = pixel_color;

	for (wd = (wd + 1) / 2 ; ; )
	{
		if (x0 >= 0 && y0 >= 0 && x0 < img -> h && y0 < img -> w)
		{
			putPixel(img, y0, x0, pixel);
		}

		e2 = err;
		x2 = x0;

		if (2 * e2 >= -dx)
		{
			for (e2 += dy, y2 = y0; e2 < ed * wd && (y1 != y2 || dx > dy); e2 += dx)
			{
				if (x0 >= 0 && x0 < img -> h && (y2 + sy) >= 0 && (y2 + sy) < img -> w)
				{
					putPixel(img, (y2 += sy), x0, pixel);
				}
			}

			if (x0 == x1)
			{
				break;
			}

			e2 = err;
			err -= dy;
			x0 += sx;
		}

		if (2 * e2 <= dy)
		{
			for (e2 = dx - e2; e2 < ed * wd && (x1 != x2 || dx < dy); e2 += dy)
			{
				if ((x2 + sx >= 0 && x2 + sx < img -> h) && (y0 >= 0 && y0 < img -> w))
				{
					putPixel(img, y0, x2 += sx, pixel);
				}
			}

			if (y0 == y1)
			{
				break;
			}

			err += dx;
			y0 += sy;
		}


	}
}

typedef struct point {
	double x, y;
} point;

typedef struct line {
	point p1;
	point p2;
	struct line* next;
} line;

void vectors(line* L) {
	while (L) {
		L->p1.x = L->p1.x-L->p2.x;
		L->p1.y = L->p1.y-L->p2.y;
		L->p2.x = 0;
		L->p2.y = 0;
		L = L->next;
	}
}

double findAngle(int ux, int uy, int vx, int vy)
{
    int mult, prod;
    
    mult = sqrt(ux*ux+uy*uy)*sqrt(vx*vx+vy*vy);
    prod = ux*vx + uy*vy;
    if (prod < 0)
        prod = -prod;

    if (mult == 0)      //division per 0 exception
        return 1.5708;    //value of radian right angle
    //compute the angle between u and v
    return acos(((double)prod)/((double)mult));
}

point findIntersection(line l1, line l2)
{
	double a1 = l1.p2.y - l1.p1.y;
	double b1 = l1.p1.x - l1.p2.x;
	double c1 = a1 * (l1.p1.x) + b1 * (l1.p1.y);


	double a2 = l2.p2.y - l2.p1.y;
	double b2 = l2.p1.x - l2.p2.x;
	double c2 = a2 * (l2.p1.x) + b2 * (l2.p1.y);

	double determinant = a1*b2 - a2*b1; // 3*-1 - -4*-3 = -3 -12 = -15
	//printf("a1: %f b1: %f c1: %f\na2: %f b2: %f c2: %f\ndet: %f\n", a1, b1, c1, a2, b2, c2, determinant);
	point p;
	if (determinant == 0) {
		p.x = 0;
		p.y = 0;
	} else {
		p.x = (b2 * c1 - b1 * c2) / determinant; // -1 * 0 - -3*-5 = -15/-15
		//printf("x: %f y: %f\n", (b2*c1 - b1*c2) / determinant, (a1*c2 - a2*c1) / determinant);
		p.y = (a1 * c2 - a2 * c1) / determinant; // 3*-5 - -4*0 = -15/-15
	}
	return p;
}

void freeMat(double **mat, int n)
{
    for (int i = 0; i < n; i++)
    {
        free(mat[i]);
    }
    free(mat);
}

double **allocMat(int size1, int size2)
{
    double **mat = calloc(size1, sizeof(double *));
    if (mat == NULL)
    {
        printf("Error allocating memory for matrix\n");
        exit(1);
    }
    for (int i = 0; i < size2; i++)
    {
        mat[i] = calloc(size2, sizeof(double));
        if (mat[i] == NULL)
        {
            printf("Error allocating memory for matrix\n");
            exit(1);
        }
    }
    return mat;
}

void houghLines(SDL_Surface* surface) {
	/*line l1, l2;
	l1.p1.x = 0;
	l1.p1.y = 1;
	l1.p2.x = 0;
	l1.p2.y = 4;
	l2.p1.x = 1;
	l2.p1.y = 8;
	l2.p2.x = 1;
	l2.p2.y = 4;
	point p = findIntersection(l1, l2);
	printf("p(%f, %f)", p.x, p.y);*/
	int h = surface->h;
	int w = surface->w;
	SDL_PixelFormat* format = surface->format;


	SDL_Surface* tempSurface = SDL_ConvertSurface(surface, format, surface->flags);
	int diag = sqrt(w*w+h*h);
    int (*accumulator)[181] = malloc(sizeof(int[2*diag][181]));

	double maxTheta = 180.0, minTheta = 0.0;
	double maxRho = diag, minRho = -diag;
	double nbRho = 2 * diag, nbTheta = nbRho;

	double rhoStep = (maxRho - minRho) / nbRho;
	double *rhos = calloc(nbRho+1, sizeof(double));
	double step = (maxTheta - minTheta) / nbTheta;
    double *thetas = calloc(nbTheta + 1, sizeof(double));

	int index = 0;
    for (double val = minRho; val <= maxRho && index < nbTheta;
         val += rhoStep, index++)
    {
        rhos[index] = val;
    }

	index = 0;

    for (double val = minTheta; val <= maxTheta && index < nbTheta;
         val += step, index++)
    {
        thetas[index] = val;
    }

	double *cos_thetas = calloc(nbTheta + 1, sizeof(double));
    double *sin_thetas = calloc(nbTheta + 1, sizeof(double));

	for (int theta = 0; theta < nbTheta; theta++)
	{
		thetas[theta] = theta * (M_PI / 180);

		double t = thetas[theta];
		cos_thetas[theta] = cos(t);
		sin_thetas[theta] = sin(t);
		step++;
	}

	int max = 0;
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			Uint8 r, g, b;
            Uint32 color = getPixel(surface, x, y);
            SDL_GetRGB(color, surface->format, &r, &g, &b);
            float average = 0.3*r+0.59*g+0.11*b; 
			if (average >= 200) {
				for (int theta = 0; theta <= 180; theta++)
                {
                    int rho = (int)(x * cos_thetas[theta] + y * sin_thetas[theta]);
					accumulator[theta][rho]++;

					if (accumulator[theta][rho] > max)
                    {
                        max = accumulator[theta][rho];
                    }
                }
			}
		}
	}

	float threshold = 0.45 * max;
	line* L = calloc(200, sizeof(line));

	for (int rho = 0; rho < 2*diag; rho++)
	{
		for (int theta = 0; theta <= 180; theta++)
		{
			if (abs(theta - 45) < 10 || abs(theta - 135) < 20) // skip diagonal lines
					continue;
			if (abs(theta - 90) > 10 && abs(theta) > 10) // skewed line => skip 
				continue;

			if (accumulator[rho][theta] >= threshold) {
				for (int y = 0; y < h; y++) {
					for (int x = 0; x < w; x++) {
						if (rho == (int)(x*cos_thetas[theta]+y*sin_thetas[theta]))
							putPixel(surface, x, y, SDL_MapRGB(format, 255, 0, 0));
					}
				}
			}

			/*int val = accumulator[rho][theta];

			if (val >= threshold) {
				double r = rhos[prev_rho], t = thetas[prev_theta];

				if (theta > tempMaxTheta) {
					tempMaxTheta = t;
                    rounded_angle = (unsigned int) (t * (180 / M_PI));
					if (rounded_angle < 181)
                    	histogram[rounded_angle]++;
				}

				double c = cos(t), s = sin(t);
				point p1, p2, p3;
				p1.x = c*r;
				p1.y = s*r;

				p2.x = p1.x + diag * -s;
				p2.y = p1.y + diag * c;

				p3.x = p1.x - diag * -s;
				p3.y = p1.y - diag * c;

				line* l = calloc(1, sizeof(line));
				l->p1.x = p2.x;
				l->p1.y = p2.y;
				l->p2.x = p3.x;
				l->p2.y = p3.y;
				l->next = L->next;
				L->next = l;
				DrawLine(tempSurface, l->p1.y, l->p1.x, l->p2.y, l->p2.x, 0.000003 * w, SDL_MapRGB(format, 255, 0, 0));
			}*/
		}
	}

	/*double maxTheta = 180.0, minTheta = 0.0;
	double maxRho = diag, minRho = -diag;
	double nbRho = 2 * diag, nbTheta = nbRho;

	double rhoStep = (maxRho - minRho) / nbRho;
	double *rhos = calloc(nbRho+1, sizeof(double));
	double step = (maxTheta - minTheta) / nbTheta;
    double *thetas = calloc(nbTheta + 1, sizeof(double));

	int index = 0;
    for (double val = minRho; val <= maxRho && index < nbTheta;
         val += rhoStep, index++)
    {
        rhos[index] = val;
    }

	index = 0;

    for (double val = minTheta; val <= maxTheta && index < nbTheta;
         val += step, index++)
    {
        thetas[index] = val;
    }

	double *cos_thetas = calloc(nbTheta + 1, sizeof(double));
    double *sin_thetas = calloc(nbTheta + 1, sizeof(double));


	Uint8 r, g, b;

	for (int theta = 0; theta < nbTheta; theta++)
	{
		thetas[theta] = theta * (M_PI / 180);

		double t = thetas[theta];
		cos_thetas[theta] = cos(t);
		sin_thetas[theta] = sin(t);
		step++;
	}

	int **accumulator = allocMat(nbTheta + 1, nbRho + 1);

	int theta, rho;
	int croppedRho, max = 0;
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
            Uint32 color = getPixel(surface, x, y);
            SDL_GetRGB(color, format, &r, &g, &b);
            float average = 0.3*r+0.59*g+0.11*b;
			if (average >= 200) {
				for (int theta = 0; theta < nbTheta; theta++)
                {
					rho = (x * cos_thetas[theta]) + (y * sin_thetas[theta]);
					croppedRho = rho + diag;
					accumulator[croppedRho][theta]++;

					if (accumulator[croppedRho][theta] > max)
                    {
                        max = accumulator[croppedRho][theta];
                    }
                }
			}
		}
	}

	//printf("max: %i", maxi);
	// define threshold 
	float threshold = 0.4 * max;

	// every value >= threshold in accu is saved


	int first = 1;

	line* L = calloc(200, sizeof(line));

	double tempMaxTheta = 0.0;
    unsigned int histogram[181] = { 0 };
    unsigned int rounded_angle;

    int prev = accumulator[0][0];
    int prev_theta = 0, prev_rho = 0;
    int boolIsIncreasing = 1;

	for (int theta = 0; theta <= nbTheta; theta++)
	{
		for (int rho = 0; rho <= nbRho; rho++)
		{
			if (abs(theta - 45) < 10 || abs(theta - 135) < 20) // skip diagonal lines
					continue;
			if (abs(theta - 90) > 10 && abs(theta) > 10) // skewed line => skip 
				continue;

			int val = accumulator[rho][theta];

			if (val >= prev) {
				prev = val;
                prev_rho = rho;
                prev_theta = theta;
                boolIsIncreasing = 1;
                continue;
			} else if (val < prev && boolIsIncreasing)
				boolIsIncreasing = 0;
			else if (val < prev) {
				prev = val;
                prev_rho = rho;
                prev_theta = theta;
                continue;
			}

			if (val >= threshold) {
				double r = rhos[prev_rho], t = thetas[prev_theta];

				if (t > tempMaxTheta) {
					tempMaxTheta = t;
                    rounded_angle = (unsigned int) (t * (180 / M_PI));
					if (rounded_angle < 181)
                    	histogram[rounded_angle]++;
				}

				double c = cos(t), s = sin(t);
				point p1, p2, p3;
				p1.x = c*r;
				p1.y = s*r;

				p2.x = p1.x + diag * -s;
				p2.y = p1.y + diag * c;

				p3.x = p1.x - diag * -s;
				p3.y = p1.y - diag * c;

				line* l = calloc(1, sizeof(line));
				l->p1.x = p2.x;
				l->p1.y = p2.y;
				l->p2.x = p3.x;
				l->p2.y = p3.y;
				l->next = L->next;
				L->next = l;
				DrawLine(tempSurface, l->p1.y, l->p1.x, l->p2.y, l->p2.x, 0.000003 * w, SDL_MapRGB(format, 255, 0, 0));
			}

			if (accumulator[r][t] > threshold) {
				theta = thetas[t];
				first = 1;
				
				if (abs(theta - 45) < 10 || abs(theta - 135) < 20) // skip diagonal lines
					continue;
				if (abs(theta - 90) > 10 && abs(theta) > 10) // skewed line => skip 
					continue;

				line* l = calloc(1, sizeof(line));
				for (int y = 0; y < h; y++)
				{
					for (int x = 0; x < w; x++)
					{
						if (r == ceil(x * cos_thetas[t] + y * sin_thetas[t])) {
							//putPixel(surface, x, y, SDL_MapRGB(format, 255, 0, 0));
							//printf("%i %i ", i, j);
							
							if (first)
							{
								l->p1.x = y;
								l->p1.y = x;
								l->p2.x = y;
								l->p2.y = x;
								first = 0;
							}
							l->p2.x = y;
							l->p2.y = x;
						}
					}
				}
				l->next = L->next;
				L->next = l;
			}
		}	
	}
	free(thetas);
	free(rhos);
	freeMat(accumulator, nbTheta+1);*/

  	IMG_SavePNG(tempSurface, "lines.png");

	//vectors(L);

	/*line* tmp = L->next;
	while (L)
	{
		//float angle = findAngle(L->p1.x, L->p1.y, tmp->p1.x, tmp->p1.y);
		//printf("angle: %f\n", angle);
		//if (angle >= 1) {
		if (tmp) {
			point p = findIntersection(*L, *tmp);
			if (p.x != 0 && p.y != 0)
				putPixel(surface, p.x, p.y, SDL_MapRGB(format, 255, 0, 0));
			tmp = tmp->next;
		}
		L = L->next;
	}*/

}


/*void hough_transform2(SDL_Surface* surface)
{
	int height = surface->h;
	int width = surface->w;
	SDL_PixelFormat* format = surface->format;
    // Biggest possible accumulator.
    int diag = sqrt(width*width + height*height);
	int nbThetas = 181;

    int (*accumulator)[nbThetas] = malloc(sizeof(int[2*diag+1][nbThetas]));
	memset(accumulator, 0, sizeof(*accumulator));

	Uint8 r, g, b;
	double thetas[nbThetas];
	double sin_thetas[nbThetas];
	double cos_thetas[nbThetas];
	
	int step = 0;
	for (int i = 0; i < nbThetas; i++)
	{
		thetas[i] = step;

		sin_thetas[i] = sin(thetas[i] * (M_PI / 180));
		cos_thetas[i] = cos(thetas[i] * (M_PI / 180));
		step++;
	}

	for (int i = 0; i < 2 * diag + 1; i++)
	{
		for (int j = 0; j < nbThetas; j++)
		{
			accumulator[i][j] = 0;
		}
	}

	int theta, rho;

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
            Uint32 color = getPixel(surface, x, y);
            SDL_GetRGB(color, format, &r, &g, &b);
            float average = 0.3*r+0.59*g+0.11*b;
			if (average >= 200) {
				for (int i = 0; i < nbThetas; i++)
                {
					rho = ceil((x * cos_thetas[i]) + (y * sin_thetas[i]));
					theta = thetas[i];
					accumulator[rho][theta]++;
                }
			}
		}
	}
	// define threshold genre 0.5*max accu
	// parcours accu tt >= threshold = good
	// parcours image si dans save alors line = good
	
	// get max of accu
	int maxi = accumulator[0][0];
	for (int y = 0; y < 2*diag+1; y++)
		for (int x = 0; x < nbThetas; x++)
			if (accumulator[y][x] >= maxi) 
				maxi = accumulator[y][x];
	//printf("max: %i", maxi);
	// define threshold 
	float threshold = 0.45 * maxi;

	// every value >= threshold in accu is saved


	int coords[height][width][4];

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			coords[i][j][0] = 0;
			coords[i][j][1] = 0;
			coords[i][j][2] = 0;
			coords[i][j][3] = 0;
		}
	}

	int first = 1;

	for (int r = 0; r < 2*diag+1; r++)
	{
		for (int t = 0; t < nbThetas; t++)
		{
			if (accumulator[r][t] > threshold) {
				theta = thetas[t];
				first = 1;
				
				if (abs(theta - 45) < 10 || abs(theta - 135) < 20) // skip diagonal lines
					continue;
				if (abs(theta - 90) > 10 && abs(theta) > 10) // skewed line => skip 
					continue;

				for (int i = 0, y = 0; y < height; y++, i++)
				{
					for (int j = 0, x = 0; x < width; x++, j++)
					{
						if (r == ceil(x * cos_thetas[t] + y * sin_thetas[t])) {
							//putPixel(surface, x, y, SDL_MapRGB(format, 255, 0, 0));
							//printf("%i %i ", i, j);
							if (first)
							{
								coords[i][j][0] = x;
								coords[i][j][1] = y;
								coords[i][j][2] = x;
								coords[i][j][3] = y;
								first = 0;
							}
							coords[i][j][2] = x;
							coords[i][j][3] = y;
						}
					}
				}
			}
		}	
	}

	free(accumulator);

	// for each coordinates saved in coords get the vectors of the lines
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int x1 = coords[y][x][0];
			int y1 = coords[y][x][1];
			int x2 = coords[y][x][2];
			int y2 = coords[y][x][3];
			//putPixel(surface, x, y, SDL_MapRGB(format, 255, 0, 0));
			coords[y][x][0] = abs(x1-x2);
			coords[y][x][1] = abs(y1-y2);
			coords[y][x][2] = 0;
			coords[y][x][3] = 0;
		}
//		printf("\n");
	}

	//printf("Here\n");
	float alpha = 0;
	int v1x = 0, v1y = 0;
	printf("angle: %f ", alpha);
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (!v1x) {
				v1x = coords[y][x][0];
				v1y = coords[y][x][1];
			}
			else {
				int v2x = coords[y][x][0], v2y = coords[y][x][1];
				//float alpha = 0;
				//float alpha = findAngle(v1x, v1y, v2x, v2y);
				//v1x = 0; v1y = 0;
				//if (alpha) {
				//	putPixel(surface, v1x, v1y, SDL_MapRGB(format, 255, 0, 0));
				//}
				//v1x = 0; v1y = 0;
			}
			//if (coords[y][x][0] || coords[y][x][1]) {
			//	putPixel(surface, coords[y][x][0], coords[y][x][1], SDL_MapRGB(format, 0, 0, 255));
			//}
		}
	}

	while (vL) {
		hLine* temp = hL;
		while (temp) {
			if (temp->x == vL->x && temp->y == vL->y) {
				//coords[temp->y][temp->x] = 1;
				putPixel(surface, temp->x, temp->y, SDL_MapRGB(format, 0, 255, 0));
			}
			temp = temp->next;
		}
		vL = vL->next;
	}
	

	//free(accumulator);
	// Processing every saved value there was in the accu to get the lines
	vLine* vL = calloc(200, sizeof(vLine));
	hLine* hL = calloc(200, sizeof(hLine));
	int countVL = 0, countHL = 0;

	int i = 0;

	for (int r = 0; r < width; r++) {
		for (int t = 0; t < height; t++) {
		  i = 0;
			if (savedAccu[t][r] != 0)
				for (int y = 0; y < height; y++, i++) {
					for (int x = 0; x < width; x++) {
						if (r == ceil(x * cos(degrees_to_rad2(t)) 
									+ y * sin(degrees_to_rad2(t)))) {
							// to see the lines
							if (sin(degrees_to_rad2(t)) > 0.85) { // horizontal
								hLine* l = calloc(1, sizeof(hLine));
								l->x = x;
								l->y = y;
								l->r = r;
								l->t = t;
								l->next = hL->next;
								hL->next = l;
								countHL++;
								//putPixel(surface, x, y, SDL_MapRGB(format, 255, 0, 0));
							} else if (sin(degrees_to_rad2(t)) < 0.2) { // vertical
								vLine* l = calloc(1, sizeof(vLine));
								l->x = x;
								l->y = y;
								l->r = r;
								l->t = t;
								l->next = vL->next;
								vL->next = l;
								countVL++;
								//putPixel(surface, x, y, SDL_MapRGB(format, 0, 255, 0));
							}
						}
					}
			}
		}
	}
	free(savedAccu);

	vLine* filtered = calloc(countVL, sizeof(vLine));
	vLine* temp = vL;
	int lasty = 0, lastx = 0;
	int value = 1500;
	printf("vl: %i, hl: %i\n", countVL, countHL);
	int c = 0;
	
	while (vL) {
		while (temp)
		{
			if (temp->y != lasty)
				lastx = temp->x;

			if (temp->x >= vL->x-value && temp->x <= vL->x+value) {
				vLine* l = calloc(1, sizeof(vLine));
				l->x = lastx;
				l->y = temp->y;
				l->r = temp->r;
				l->t = temp->t;
				l->next = filtered->next;
				filtered->next = l;
				lasty = temp->y;
				lastx = temp->x;
				c++;
				//putPixel(surface, temp->x, temp->y, SDL_MapRGB(format, 0, 255, 0));
			}
			temp = temp->next;
		}
		vL = vL->next;
	}

	printf("c: %i\n", c);

	while (filtered)
	{
		putPixel(surface, filtered->x, filtered->y, SDL_MapRGB(format, 0, 255, 0));
		filtered = filtered->next;
	}

	//int (*coords)[diag] = calloc(diag, sizeof(*coords));
	

	// Quick Method to get intersection points between horizontal and vertical lines
	while (vL) {
		hLine* temp = hL;
		while (temp) {
			if (temp->x == vL->x && temp->y == vL->y) {
				coords[temp->y][temp->x] = 1;
				//putPixel(surface, temp->x, temp->y, SDL_MapRGB(format, 0, 255, 0));
			}
			temp = temp->next;
		}
		vL = vL->next;
	}
	
	int x1 = 0, x2 = 0, x3 = 0, x4 = 0, 
		y1 = 0, y2 = 0, y3 = 0, y4 = 0;

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (coords[y][x]) {
				x1 = x;
				y1 = y;
				x = width;
				y = height;
			}
		}
	}

	for (int y = 0; y < height; y++) {
		for (int x = width; x > 0; x--) {
			if (coords[y][x] && abs(x1-x) > 1000) {
				x2 = x;
				y2 = y;
				x = 0;
				y = height;
			}
		}
	}

	for (int y = height; y > 0; y--) {
		for (int x = 0; x < width; x++) {
			if (coords[y][x] && abs(y2-y) > 100) {
				x3 = x;
				y3 = y;
				x = width;
				y = 0;
			}
		}
	}


	for (int y = height; y > 0; y--) {
		for (int x = width; x > 0; x--) {
			if (coords[y][x] && abs(x3-x) > 100 && abs(y3-y) > 100) {
				x4 = x;
				y4 = y;
				x = 0;
				y = 0;
			}
		}
	}

	for (int i = 0; i < 50; i++)
	{
		putPixel(surface, x1+i, y1+i, SDL_MapRGB(format, 255, 0, 0));
		putPixel(surface, x2+i, y2+i, SDL_MapRGB(format, 255, 0, 0));
		putPixel(surface, x3+i, y3-i, SDL_MapRGB(format, 0, 255, 0));
		putPixel(surface, x4+i, y4-i, SDL_MapRGB(format, 0, 255, 0));
	}


	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (coords[y][x])
			{
				int lastx = width;
				int lasty = y;
				while (lastx > 0 && lastx != x && !coords[lasty][lastx])
				{
					lastx--;
				}
					//if (coords[lasty+i][lastx]) break;
					//if (lastx == 0)
						//lastx = width;
				lasty += i;
				//printf("x(%i, %i) y(%i, %i)\n", x, lastx, y, lasty);
				for (int i = y; i <= lasty; i++)
				{
					for (int j = x; j <= lastx; j++)
					{
						coords[i][j] = 1;
						putPixel(surface, j, i, SDL_MapRGB(format, 255, 0, 0));
					}
					
				}
				
				//printf("same height? y: %i, lasty: %i\n", y, lasty);
			}
			
		}
	}


	// ========= temp ==========
	//for (int k = 0; k < width; k++)
	//	printf("x1(%i,%i) x2(%i,%i) ", coords[0][k][0], coords[0][k][1], coords[0][k][2], coords[0][k][3]);
	// ========= temp ==========

	// for each coordinates saved in coords get the vectors of the lines
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int x1 = coords[y][x][0];
			int y1 = coords[y][x][1];
			int x2 = coords[y][x][2];
			int y2 = coords[y][x][3];
			//putPixel(surface, x, y, SDL_MapRGB(format, 255, 0, 0));
			coords[y][x][0] = abs(x1-x2);
			coords[y][x][1] = abs(y1-y2);
			coords[y][x][2] = 0;
			coords[y][x][3] = 0;
		}
//		printf("\n");
	}

	// goal is to get the intersection of two vectors and it should give every intersection of every line
	// then we can save only the 4 corners and crop the image to the grid
	int v1x = 0, v1y = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (!v1x) {
				v1x = coords[y][x][0];
				v1y = coords[y][x][1];
			}
			else {
				printf("Here\n");
				int v2x = coords[y][x][0], v2y = coords[y][x][1];
				float alpha = (v1x*v2x+v1y*v2y);
				printf("angle1: %f ", alpha);
				alpha /= (sqrt(v1x*v1x+v1y*v1y)*sqrt(v2x*v2x+v2y*v2y));
				printf("angle: %f\n", alpha);
				//v1x = 0; v1y = 0;
				if (alpha) {
					putPixel(surface, v1x, v1y, SDL_MapRGB(format, 0, 0, 255));
				}
				v1x = 0; v1y = 0;
			}
			//if (coords[y][x][0] || coords[y][x][1]) {
			//	putPixel(surface, coords[y][x][0], coords[y][x][1], SDL_MapRGB(format, 0, 0, 255));
			//}
		}
	}
	
	//free(coords);
}*/

/// @brief cut the image in 9x9
/// @param surface  the surface to cut
void cut(SDL_Surface* surface) {
	int height = surface->h;
	int width = surface->w;
	int currWidth = 0, currHeight = 0, count = 0;
	for (int i = 9; i > 0; i--)
	{
		currWidth = 0;
		for (int j = 9; j > 0; j--)
		{
			SDL_Surface* tempSurface = SDL_CreateRGBSurface(surface->flags, width/9, height/9, surface->format->BitsPerPixel, surface->format->Rmask, surface->format->Gmask, surface->format->Bmask, surface->format->Amask);
			SDL_Rect rect = {currWidth, currHeight, width/9, height/9};
			currWidth += width/9;
			SDL_BlitSurface(surface, &rect, tempSurface, 0);
			char* name;
			if (count < 9)
				name = "cases/0%i.png";
			else
				name = "cases/%i.png";
			char res[100];
			count++;
			sprintf(res, name, count);
			SDL_SaveBMP(tempSurface, res);
		}
		currHeight += height/9;
	}
}

/// @brief simple truncate function
/// @param x the value to truncate
/// @return value maxed at 255
float truncate(float x) {
	return x > 255 ? 255 : x;
}


/// @brief apply sobel algorithm on image
/// @param surface the surface to apply algorithm on
/// @return the new surface where sobel was applied
SDL_Surface* sobel(SDL_Surface* surface) {
	SDL_PixelFormat* format = surface->format;
	SDL_Surface* tempSurface = SDL_ConvertSurface(surface, format, surface->flags);
	int w = surface->w;
	int h = surface->h;

	Uint8 r, g, b;

	float sumr, sumg, sumb;
	float resr, resg, resb;
	float tempr, tempg, tempb;
	int c = 0;

	int xKernel[] = {-1, 0, 1,
					 -2, 0, 2,
					 -1, 0, 1,};

	int yKernel[] = {-1,  -2,  -1,
					  0,   0, 	0,
					  1,   2, 	1,};

	SDL_LockSurface(tempSurface);
	SDL_LockSurface(surface);
	
	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			sumr = 0, sumg = 0, sumb = 0, 
			tempr = 0, tempg = 0, tempb = 0, c = 0;
			for (int i = -1; i < 2; i++)
			{
				for (int j = -1; j < 2; j++)
				{
					if (getPixel(surface, x+j, y+i)) {
						SDL_GetRGB(getPixel(surface, x+j, y+i), format, &r, &g, &b);
						tempr += r * xKernel[c];
						tempg += g * xKernel[c];
						tempb += b * xKernel[c]; // horizontal convolution
						sumr += r * yKernel[c];
						sumg += g * yKernel[c];
						sumb += b * yKernel[c]; // vertical convolution
					}
					c++;
				}
			}
			resr = truncate(sqrt(tempr * tempr + sumr * sumr));
			resg = truncate(sqrt(tempg * tempg + sumg * sumg));
			resb = truncate(sqrt(tempb * tempb + sumb * sumb)); // magnitude
			putPixel(tempSurface, x, y, SDL_MapRGB(format, resr, resg, resb));
		}
	}
	
	SDL_UnlockSurface(tempSurface);
	SDL_UnlockSurface(surface);

	return tempSurface;
}

SDL_Surface* dilate(SDL_Surface* surface) {
	SDL_PixelFormat* format = surface->format;
	SDL_Surface* tempSurface = SDL_ConvertSurface(surface, format, surface->flags);
	int w = surface->w;
	int h = surface->h;

	Uint8 r, g, b;

	SDL_LockSurface(tempSurface);
	SDL_LockSurface(surface);
	
	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			SDL_GetRGB(getPixel(surface, x, y), format, &r, &g, &b);
			float average = 0.3*r+0.59*g+0.11*b; 
			if (average < 200) { // if it is a black pixel
				SDL_GetRGB(getPixel(surface, x-1, y), format, &r, &g, &b); //back pixel
				average = 0.3*r+0.59*g+0.11*b; 
				if (average >= 200)
					putPixel(tempSurface, x, y, SDL_MapRGB(format, 255, 255, 255));
				else {
					SDL_GetRGB(getPixel(surface, x, y-1), format, &r, &g, &b); //top pixel
					average = 0.3*r+0.59*g+0.11*b;
					if (average >= 200)
						putPixel(tempSurface, x, y, SDL_MapRGB(format, 255, 255, 255));
					else {
						SDL_GetRGB(getPixel(surface, x+1, y), format, &r, &g, &b); //right pixel
						average = 0.3*r+0.59*g+0.11*b;
						if (average >= 200)
							putPixel(tempSurface, x, y, SDL_MapRGB(format, 255, 255, 255));
						else {
							SDL_GetRGB(getPixel(surface, x, y+1), format, &r, &g, &b); //top pixel
							average = 0.3*r+0.59*g+0.11*b;
							if (average >= 200)
								putPixel(tempSurface, x, y, SDL_MapRGB(format, 255, 255, 255));
						}
					}
				}
			}
		}
	}

	SDL_UnlockSurface(tempSurface);
	SDL_UnlockSurface(surface);

	return tempSurface;
}

SDL_Surface* erode(SDL_Surface* surface) {
	SDL_PixelFormat* format = surface->format;
	SDL_Surface* tempSurface = SDL_ConvertSurface(surface, format, surface->flags);
	int w = surface->w;
	int h = surface->h;

	Uint8 r, g, b;
	int c = 0;
	float sumr, sumg, sumb;

	int kernel[] = { 0, 1, 0,
					 1, 1, 1,
					 0, 1, 0,};

	SDL_LockSurface(tempSurface);
	SDL_LockSurface(surface);
	
	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			sumr = 0, sumg = 0, sumb = 0, c = 0;
			SDL_GetRGB(getPixel(surface, x, y), format, &r, &g, &b);
			float average = 0.3*r+0.59*g+0.11*b; 
			if (average > 200) { // if it is a white pixel
				for (int i = -1; i < 2; i++)
				{
					for (int j = -1; j < 2; j++)
					{
						if (getPixel(surface, x+j, y+i)) {
							SDL_GetRGB(getPixel(surface, x+j, y+i), format, &r, &g, &b);
							sumr += r * kernel[c];
							sumg += g * kernel[c];
							sumb += b * kernel[c];
						}
						c++;
					}
				}

				average = 0.3*sumr/5+0.59*sumg/5+0.11*sumb/5;
				if (average <= 120)//if (sumr >= 2295)
					putPixel(tempSurface, x, y, SDL_MapRGB(format, 0, 0, 0));
			}
		}
	}

	SDL_UnlockSurface(tempSurface);
	SDL_UnlockSurface(surface);

	return tempSurface;
}

/*void secondPass(SDL_Surface* surface, int* map[], int* equivalenceList) {
	SDL_PixelFormat* format = surface->format;
	int w = surface->w;
	int h = surface->h;

	Uint8 r, g, b;
	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			SDL_GetRGB(getPixel(surface, x, y), format, &r, &g, &b);
			float average = 0.3*r+0.59*g+0.11*b;
			if (average >= 200) {
				if (map[y][x] != equivalenceList[map[y][x]]) {
					map[y][x] = equivalenceList[map[y][x]];	
					//secondPass(surface, map, equivalenceList);
				}
			}
		}
	}
}*/

void fill(SDL_Surface *surface_to_fill, int step, int **M, int rank, 
      int maxi, int width, int height)
{
    int pos;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            pos = step < 0 ? width-1-j : j;
            if (M[i][pos] != maxi)
            {
                putPixel(surface_to_fill, pos, i, SDL_MapRGB(surface_to_fill->format, 200, 0, 0));
            }
            else
            {
                j = width;
            }
        }
    }
}

void next_nbh(int *src)
{
    int x = src[0];
    int y = src[1];

    if (x && !y) src[1] = -x;
    else if (x && x == -y) src[0] = 0;
    else if (!x && y) src[0] = y;
    else src[1] = 0;
}

int moore(SDL_Surface *surface, int ind, int width, int height, 
        int **M, int *s, int prev_x, int prev_y)
{
    //printf("moore \n");
    int maxi = 0, x = s[0], y = s[1], c0, c1;
    int x_min = x, x_max = x, y_min = y, y_max = y, cmp = 0;


    Uint32 pixel, a;
    Uint8 r, g, bl;

    int *p, *b, *c;
    p = malloc(2*sizeof(int));
    c = malloc(2*sizeof(int));
    b = malloc(2*sizeof(int));

    p[0] = x; p[1] = y;
    b[0] = prev_x; b[1] = prev_y;
    c[0] = prev_x - x; c[1] = prev_y - y; 

    next_nbh(c);

    char test = 1;
    while (c[0] != s[0] || c[1] != s[1])
    {
        while ((p[0] + c[0] < 0 || p[1] + c[1] < 0 
                || p[0] + c[0] >= height || p[1] + c[1] >= width) && cmp < 10)
        {
            next_nbh(c);
            cmp++;
        }

        if (cmp > 9)
            break;

        pixel = getPixel(surface, c[1]+p[1],c[0]+p[0]);
        SDL_GetRGB(pixel, surface->format, &r, &g, &bl);
        a = 0.3*r+0.59*g+0.11*bl;
        if (a> 225)
        {
            putPixel(surface, p[1],p[0], SDL_MapRGB(surface->format, 0, 0, 0));

            c0 = c[0]; c1 = c[1];
            M[p[0]+c0][p[1]+c1] = ind;
            b[0] = p[0] ; b[1] = p[1];
            p[0] += c0 ; p[1] += c1;
            c[0] = b[0] - p[0]; c[1] = b[1] - p[1];
            x_min = p[0] < x_min? p[0] : x_min;
            x_max = p[0] > x_max? p[0] : x_max;
            y_min = p[1] < y_min? p[1] : y_min;
            y_max = p[1] > y_max? p[1] : y_max;
            next_nbh(c);
            cmp = 0;
        }
        else
        {
            b[0] = p[0] + c[0]; b[1] = p[1] + c[1];
            next_nbh(c);
            cmp++;
        }
    }
    maxi = sqrt((x_max-x_min)*(x_max-x_min) + (y_max-y_min)*(y_max-y_min));

    return maxi;
}

void componentLabeling(SDL_Surface *surface)
{
    //definir M[h][w] ou comme pixels
    //mise à 0
    //pointeur ?
    //
    //Choisir le sens de rotation
    //
    //--Pendant le process de moore :--
    //Si on tombe sur un px de couleur,
    //Si c'est le pixel de départ, on arrête 
    //Si il est marqué, on annule et on va tout droit
    //Sinon on y va, on le marque et on repart

    
    int width = surface->w;
    int height = surface->h;

    //printf("width = %i ; height = %i\n",width, height);
    //TODO : fill all image contours in black
    /*
    for (int i = 0; i<width;i++)
    {
        put_pixel(surface, 0, i, 200);
        put_pixel(surface, height-1, i, 200);
        //fill black pixel[0][i]
        //fill black pixel[height-1][width]
    }
    for (int i = 0; i<height;i++)
    {
        put_pixel(surface, i, 0, 200);
        put_pixel(surface, i, width-1, 200);
        //fill black pixel[i][0]
        //fill_black pixel[i][width-1]
    }
    */

    int *s = calloc(2, sizeof(int));

    int **M;
    M = calloc(height, sizeof(int *));
    for (int i = 0; i < height; i++)
    {
        M[i] = calloc(width, sizeof(int));
    }

    int j;

    Uint32 pixel;
    Uint8 r, g, b;
    Uint32 a;

    int moore_res;
    int maxi = 0;
    int ind_max = 0; //index of the maximum form
    int ind = 0; //current filling index for each contour


    for (int i = 0; i < height; i++)
    {
        j = 0;
        while(j < width)
        {
	        pixel = getPixel(surface,i,j);
	    
	        SDL_GetRGB(pixel, surface->format, &r, &g, &b);
	        a = 0.3*r+0.59*g+0.11*b;

	        if (a>225)
	        {   
	            if (M[i][j] == 0)
                {
                    //printf("\nind = %i\n",ind);
                    ind++;
                    M[i][j] = ind;
                    //printf("calling moore\n");
                    s[0] = i;
                    s[1] = j;
                    moore_res = moore(surface, ind, width, height, M, s, i-1, j);
                    
                    //printf("moore res = %i ; round %i\n", moore_res, i);
                    if (moore_res > maxi)
                    {
                        maxi = moore_res;
                        ind_max = ind;
                    }


                }
                j = width;
            }
            else
                j++;
        }
    }

    fill(surface, 1, M, 1, ind_max, width, height);
    fill(surface, -1, M, 1, ind_max, width, height);

    /*for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (M[i][j] == 0)
                printf("_");
            else
                printf("%i", M[i][j]);
        }
        printf("\n");
    }
    printf("%i\n", ind_max);*/
    //maybe call fill in the main function rather than here

}



/*void componentLabeling(SDL_Surface* surface) {
	SDL_PixelFormat* format = surface->format;
	int w = surface->w;
	int h = surface->h;

	Uint8 r, g, b;
	int label = 0;
	int *equivalenceList = malloc(sizeof(int[h*w]));
	int (*map)[w] = malloc(sizeof(int[h][w]));
	//memset(equivalenceList, 0, sizeof(*equivalenceList));
	memset(map, 0, sizeof(*map));

	for (int y = 0; y < h*w; y++)
	{
		equivalenceList[y] = 0;
	}
	
	int kernel[2][3] = {{1, 1, 1,},
						{1, 0, 0,}};

	int flag = 1;
	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			SDL_GetRGB(getPixel(surface, x, y), format, &r, &g, &b);
			float average = 0.3*r+0.59*g+0.11*b;
			flag = 1;
			if (average > 200) {
				for (int i = -1; i < 1; i++)
				{
					for (int j = -1; j < 2; j++)
					{
						if (kernel[i+1][j+1] && getPixel(surface, x+j, y+i)) {
							SDL_GetRGB(getPixel(surface, x+j, y+i), format, &r, &g, &b);
							average = 0.3*r+0.59*g+0.11*b;
							if (average > 200) {
								if (map[y+i][x+j] != label) {
									if (label > map[y+i][x+j]) {
										map[y][x] = map[y+i][x+j];
										equivalenceList[label] = map[y][x];
									} else {
										map[y][x] = map[y+i][x+j];
										equivalenceList[map[y+i][x+j]] = map[y+i][x+j];
									}
								}
								else
									map[y][x] = map[y+i][x+j];
								flag = 0;
							} else {
								flag &= 1;
							}
						}
					}
				}

				if (flag) {				
					label++;
					map[y][x] = label;
					equivalenceList[label] = label;
				}
			}
		}
	}

	//secondPass(surface, map, equivalenceList);
	printf("l: %i", label);
	for (int y = 0; y < h; y++)
	{	
		for (int x = 0; x < w; x++)
		{
			//printf("%i ", map[y][x]);
			//if (map[y][x])
			SDL_GetRGB(getPixel(surface, x, y), format, &r, &g, &b);
			float average = 0.3*r+0.59*g+0.11*b;
			if (average >= 200) {
				//if (map[y][x] != equivalenceList[map[y][x]])
				//	map[y][x] = equivalenceList[map[y][x]];
				//int smallestLabel = label;
				//for (int k = 0; k < 50; k++) {
				//		printf("(%i %i) \n", map[0][x], equivalenceList[map[0][x]][i]);
						//if (smallestLabel > equivalenceList[map[y][x]][k])
						//	smallestLabel = equivalenceList[map[y][x]][k];
				//}
				//printf("small: %i\n", smallestLabel);

				SDL_GetRGB(getPixel(surface, x, y), format, &r, &g, &b);
				float average = 0.3*r+0.59*g+0.11*b;
				if (average >= 200) {
					while (map[y][x] != equivalenceList[map[y][x]] && equivalenceList[map[y][x]]) {
						map[y][x] = equivalenceList[map[y][x]];	
						//secondPass(surface, map, equivalenceList);
					}
				}

			}
		}
	}
	
	int labels[label+1];
	for	(int i = 0; i < label+1; i++) {
		labels[i] = 0;
	}
	
	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			if (map[y][x])
				labels[map[y][x]]++;
		}
	}

	int maxLabel = -1;
	for (int i = 0; i < label+1; i++) {
		printf("i: %i\n", labels[i]);
		if (labels[i] >= maxLabel)
			maxLabel = i;
	}

	printf("\nmax: %i", maxLabel);
	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			SDL_GetRGB(getPixel(surface, x, y), format, &r, &g, &b);
			float average = 0.3*r+0.59*g+0.11*b;
			if (average >= 200)
				putPixel(surface, x, y, SDL_MapRGB(format, 10 * map[y][x], 20 * map[y][x], 15 * map[y][x]));
		}
	}


	//printf("\n");

	free(*map);
	free(equivalenceList);


}*/

/*

// =========== OTSU GREGOIRE =========

Uint32 simple(Uint32 pixel_color, SDL_PixelFormat* format, Uint32 threshold)
{
  Uint8 r, g, b;
  SDL_GetRGB(pixel_color, format, &r, &g, &b);
  Uint32 color;
  if (r > threshold)
    color = SDL_MapRGB(format, 255, 255, 255);
  else
    color = SDL_MapRGB(format, 0, 0, 0);
  return color;
}


void simple_binarize(SDL_Surface* surface, Uint32 threshold)
{
  if (surface == NULL)
    errx(EXIT_FAILURE, "%s", SDL_GetError());
  Uint32* pixels = surface->pixels;
  int len = surface->w * surface->h;
  SDL_LockSurface(surface);
  SDL_PixelFormat* format = surface->format;
  for (int i = 0; i < len; i++)
    {
      pixels[i] = simple(pixels[i], format, threshold);
    }
  SDL_UnlockSurface(surface);
}


void histogram(SDL_Surface* surface, int* histo)
{
  Uint8 r, g, b;
  int len = surface->w * surface->h;
  Uint32* pixels = surface->pixels;
  SDL_PixelFormat* format = surface->format;
  SDL_LockSurface(surface);
  
  //creates a histogram of the grayscale of the image
  for (int i = 0; i < len; i++)
    {
      //pixels[i] = pixelToGrayscale(pixels[i],format);
      SDL_GetRGB(pixels[i],format,&r,&g,&b);
      histo[b]++;
    }
  SDL_UnlockSurface(surface);
}


void equalized(SDL_Surface* surface, int* histo)
{
  int sum = 0;
  int len = surface->w * surface->h;
  for (int i = 0; i < 256; i++)
    {
      sum += histo[i];
      histo[i] = (sum * 255) / len;
    }
  Uint32* pixels = surface->pixels;
  SDL_PixelFormat* format = surface->format;
  SDL_LockSurface(surface);
  int out;
  Uint8 r, g, b;
  for (int i = 0; i < len; i++)
    {
      SDL_GetRGB(pixels[i],format,&r,&g,&b);
      out = histo[b];
      pixels[i] = SDL_MapRGB(format, out, out, out);
    }
  SDL_UnlockSurface(surface);
  for (int i = 0; i < 256; i++)
    histo[i] = 0;
  histogram(surface, histo);
}


Uint8 Otsusmethod(SDL_Surface* surface, int* histo)
{
  int len = surface->w * surface->h;
  
  //calculating the total sum of the histogram
  unsigned long sum0 = 0, sum1 = 0;
  for (int i = 0; i < 256; i++)
    {
      sum0 += histo[i] * i;
    }

  unsigned long p0 = 0, p1 = 0;
  unsigned long mean0 = 0, mean1 = 0;
  float between = 0, maxi = 0;
  Uint8 threshold0 = 0, threshold1 = 0;
  
  for (int i = 0; i < 256; i++)
    {
      p0 += histo[i];
      p1 = len - p0;
      if (p0 > 0 && p1 > 0)
	{
	  sum1 += histo[i] * i;
	  mean0 = sum1 / p0;
	  mean1 = (sum0 - sum1) / p1;
	  between = p0 * p1 * (mean0 - mean1) * (mean0 - mean1);
	  if (between >= maxi)
	    {
	      threshold0 = i;
	      if (between > maxi)
	      threshold1 = i;
	      maxi = between;
	    }
	}
    }
  return (threshold0 + threshold1) / 2;
}

*/

// ============ TEMP =================
/*
void largestConnectedComponent(SDL_Surface* surface) {
	SDL_PixelFormat* format = surface->format;
	int w = surface->w;
	int h = surface->h;

	Uint8 r, g, b;

	int zone = 0;putPixel(tempSurface, x, y, SDL_MapRGB(format, 255, 255, 255));

	int equivalenceList[w];
	float (*map)[w] = malloc(sizeof(float[h][w]));
	memset(equivalenceList, 0, sizeof(equivalenceList));
	memset(map, 0, sizeof(*map));

	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
				if (x-1 > 0 && y-1 > 0) { // if top and back pixel exist
				SDL_GetRGB(getPixel(surface, x-1, y), format, &r, &g, &b); // get values of back pixel
				if (r == 255) { // check if that pixel is white
					SDL_GetRGB(getPixel(surface, x, y-1), format, &r, &g, &b); // get values of top pixel
					if (r == 255) { // check if top pixel is white
						int zoneToTake = map[y-1][x] < map[y][x-1] ? map[y-1][x] : map[y][x-1];
						equivalenceList[zone] = zoneToTake;
						map[y][x] = zoneToTake;
					} else { // if only back pixel is white
						map[y][x] = map[y][x-1];
					}
				} else {
					SDL_GetRGB(getPixel(surface, x, y-1), format, &r, &g, &b); // get values of top pixel
					if (r == 255) { // check if top pixel is white
						map[y][x] = map[y-1][x];
					} else { // if only back pixel is white
						zone++;
						map[y][x] = zone;
						equivalenceList[zone] = zone;
			 		}
				}
			}
			else if (y-1 > 0) {
				SDL_GetRGB(getPixel(surface, x, y-1), format, &r, &g, &b); // get values of top pixel
				if (r == 255) { // check if top pixel is white
					map[y][x] = map[y-1][x];
			 	}
			}
			else if (x-1 > 0) {
				SDL_GetRGB(getPixel(surface, x-1, y), format, &r, &g, &b); // get values of top pixel
				if (r == 255) { // check if top pixel is white
					map[y][x] = map[y][x-1];
				}
			}
		}}
	}

	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			map[y][x] = equivalenceList[(int)map[y][x]];
		}
	}
	

	for (int j = 0; j < w; j++)
 	{
 		printf("%i ", equivalenceList[j]);
 	}

	int freq[zone];
	memset(freq, 0, sizeof(freq));
	printf(" zones: %i ", zone);
	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			int value = map[y][x];
			if (value != 0) freq[value]++;
		}
	}

	

    int max = freq[0];
	int savedIndex = 0;
	for (int i = 1; i < zone; i++) {
		if (freq[i] > max) {
            max = freq[i];
			savedIndex = i;
		}
	}

	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			if(map[y][x] !=0 && map[y][x] != savedIndex) {// 
				putPixel(surface, x, y, SDL_MapRGB(format, 0+10*map[y][x], 0, 0));
			}
		}
	}


}


void largestConnectedComponent(SDL_Surface* surface) {
	SDL_PixelFormat* format = surface->format;
	int w = surface->w;
	int h = surface->h;

	Uint8 r, g, b;

	int zone = 0;
	int equivalenceList[w];
	float (*map)[w] = malloc(sizeof(float[h][w]));
	memset(map, 0, sizeof(*map));
	memset(equivalenceList, 0, sizeof(equivalenceList));

	// map every pack of pixel in the map to get the largest
	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			SDL_GetRGB(getPixel(surface, x, y), format, &r, &g, &b);
			if (r == 255) {
				


				if (x-1 > 0 && y-1 > 0) { // if top and back pixel exist
					SDL_GetRGB(getPixel(surface, x-1, y), format, &r, &g, &b); // get values of back pixel
					if (r == 255) { // check if that pixel is white
						SDL_GetRGB(getPixel(surface, x, y-1), format, &r, &g, &b); // get values of top pixel
						if (r == 255) { // check if top pixel is white
							int zoneToTake = map[y-1][x] < map[y][x-1] ? map[y-1][x] : map[y][x-1];
							equivalenceList[zone] = zoneToTake;
							map[y][x] = zoneToTake;
						} else { // if only back pixel is white
							map[y][x] = map[y][x-1];
						}
					} else {
						SDL_GetRGB(getPixel(surface, x, y-1), format, &r, &g, &b); // get values of top pixel
						if (r == 255) { // check if top pixel is white
							map[y][x] = map[y-1][x];
						} else { // if only back pixel is white
							zone++;
							map[y][x] = zone;
							equivalenceList[zone] = zone;
						}
					}
				}
				else if (y-1 > 0) {
					SDL_GetRGB(getPixel(surface, x, y-1), format, &r, &g, &b); // get values of top pixel
					if (r == 255) { // check if top pixel is white
						map[y][x] = map[y-1][x];
					}
				}
				else if (x-1 > 0) {
					SDL_GetRGB(getPixel(surface, x-1, y), format, &r, &g, &b); // get values of top pixel
					if (r == 255) { // check if top pixel is white
						map[y][x] = map[y][x-1];
					}
				}
			}
		}
	}

	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			map[y][x] = equivalenceList[(int)map[y][x]];
		}
	}
	

	for (int j = 0; j < w; j++)
 	{
 		printf("%i ", equivalenceList[j]);
 	}

	int freq[zone];
	memset(freq, 0, sizeof(freq));
	printf(" zones: %i ", zone);
	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			int value = map[y][x];
			if (value != 0) freq[value]++;
		}
	}

	

    int max = freq[0];
	int savedIndex = 0;
	for (int i = 1; i < zone; i++) {
		if (freq[i] > max) {
            max = freq[i];
			savedIndex = i;
		}
	}

	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			if(map[y][x] !=0 && map[y][x] != savedIndex) {// 
				putPixel(surface, x, y, SDL_MapRGB(format, 0+10*map[y][x], 0, 0));
			}
		}
	}


	// ====== TO REMOVE ========
	// for (int i = 0; i < h; i++)
	// {
	// 	for (int j = 0; j < w; j++)
	// 	{
	// 		printf("%i ", map[i][j]);
	// 	}
	//  	printf("\n");
	// }
}
*/
