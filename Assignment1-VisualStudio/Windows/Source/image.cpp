#include "image.h"
#include "bmp.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iostream>


/**
 * Image
 **/
Image::Image (int width_, int height_)
{
    assert(width_ > 0);
    assert(height_ > 0);

    width           = width_;
    height          = height_;
    num_pixels      = width * height;
    pixels          = new Pixel[num_pixels];
    sampling_method = IMAGE_SAMPLING_POINT;

    assert(pixels != NULL);
}


Image::Image (const Image& src)
{
    width           = src.width;
    height          = src.height;
    num_pixels      = width * height;
    pixels          = new Pixel[num_pixels];
    sampling_method = IMAGE_SAMPLING_POINT;

    assert(pixels != NULL);
    memcpy(pixels, src.pixels, src.width * src.height * sizeof(Pixel));
}


Image::~Image ()
{
    delete [] pixels;
    pixels = NULL;
}

/*
void Image::AddNoise (double factor)
{

}
*/

void Image::Brighten (double factor)
{
  /* Your Work Here  (section 3.2.1 of assignment)*/
	for (int i = 0; i < num_pixels; i++)
	{
		pixels[i].SetClamp(pixels[i].r * factor,
						   pixels[i].g * factor,
						   pixels[i].b * factor);
	}
}


void Image::ChangeContrast (double factor)
{
  /* Your Work Here (section 3.2.2) */
	/* get the greyscale */
	int lum = 0;
	for (int i = 0; i < num_pixels; i++)
	{
		lum += pixels[i].Luminance();
	}
	lum /= num_pixels;

	for (int i = 0; i < num_pixels; i++)
	{
		pixels[i].SetClamp((1.0 - factor) * lum + pixels[i].r * factor,
						   (1.0 - factor) * lum + pixels[i].g * factor,
						   (1.0 - factor) * lum + pixels[i].b * factor);
	}
}


void Image::ChangeSaturation(double factor)
{
  /* Your Work Here (section 3.2.3) */
	/* get the greyscale */
	for (int i = 0; i < num_pixels; i++)
	{
		int tmp = pixels[i].Luminance();
		pixels[i].SetClamp((1.0 - factor) * tmp + pixels[i].r * factor,
						   (1.0 - factor) * tmp + pixels[i].g * factor,
						   (1.0 - factor) * tmp + pixels[i].b * factor);
	}
}

void Image::ChangeGamma(double factor)
{
  /* Your Work Here (section 3.2.4) */
	for (int i = 0; i < num_pixels; i++)
	{
		pixels[i].SetClamp(pow(pixels[i].r / 255.0, 1.0 / factor) * 255.0,
						   pow(pixels[i].g / 255.0, 1.0 / factor) * 255.0,
						   pow(pixels[i].b / 255.0, 1.0 / factor) * 255.0);
	}

}

Image* Image::Crop(int x, int y, int w, int h)
{
  /* Your Work Here (section 3.2.5) */
	Image * ret = new Image(w, h);
	Pixel * pixs = new Pixel[w*h];
	for (int i = x; i < x + w; i++)
	{
		for (int j = y; j < y + h; j++)
		{
			pixs[(j-y)*w + (i-x)] = pixels[j*width + i];
		}
	}
	ret->pixels = pixs;
	return ret;
}

/*
void Image::ExtractChannel(int channel)
{
  // For extracting a channel (R,G,B) of image.  
  // Not required for the assignment
}
*/

void Image::Quantize (int nbits)
{
  /* Your Work Here (Section 3.3.1) */
	for (int i = 0; i < num_pixels; i++) 
	{
		double num = pow(2.0, nbits);
		float pr = pixels[i].r / 256.0;
		float pg = pixels[i].g / 256.0;
		float pb = pixels[i].b / 256.0;
		pixels[i].SetClamp(floor((255.0) * floor(pr * num) / (num - 1)),
						   floor((255.0) * floor(pg * num) / (num - 1)),
						   floor((255.0) * floor(pb * num) / (num - 1)));
	}
}


void Image::RandomDither (int nbits)
{
  /* Your Work Here (Section 3.3.2) */
	for (int i = 0; i < num_pixels; i++)
	{
		double num = pow(2.0, nbits);
		float pr = pixels[i].r / 256.0;
		float pg = pixels[i].g / 256.0;
		float pb = pixels[i].b / 256.0;
		float ran = (-0.5f + rand() / (float)(RAND_MAX));
		pixels[i].SetClamp(floor((255.0) * floor(pr * num + ran) / (num - 1)),
			floor((255.0) * floor(pg * num + ran) / (num - 1)),
			floor((255.0) * floor(pb * num + ran) / (num - 1)));
	}

}


/* Matrix for Bayer's 4x4 pattern dither. */
/* uncomment its definition if you need it */

/*
static int Bayer4[4][4] =
{
    {15, 7, 13, 5},
    {3, 11, 1, 9},
    {12, 4, 14, 6},
    {0, 8, 2, 10}
};


void Image::OrderedDither(int nbits)
{
  // For ordered dithering
  // Not required for the assignment
}

*/

/* Error-diffusion parameters for Floyd-Steinberg*/
const double
    ALPHA = 7.0 / 16.0, //right
    BETA  = 3.0 / 16.0, //bottomleft
    GAMMA = 5.0 / 16.0, //bottom
    DELTA = 1.0 / 16.0; //botomright

void Image::FloydSteinbergDither(int nbits)
{
  /* Your Work Here (Section 3.3.3) */
	double num = pow(2.0, nbits);
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			float pr = GetPixel(j,i).r / 256.0;
			float pg = GetPixel(j, i).g / 256.0;
			float pb = GetPixel(j, i).b / 256.0;
			pr = floor((255.0) * floor(pr * num) / (num - 1));
			pg = floor((255.0) * floor(pg * num) / (num - 1));
			pb = floor((255.0) * floor(pb * num) / (num - 1));
			double r_err = GetPixel(j, i).r - pr;
			double g_err = GetPixel(j, i).g - pg;
			double b_err = GetPixel(j, i).b - pb;
			GetPixel(j, i).SetClamp(ComponentClamp(pr), ComponentClamp(pg), ComponentClamp(pb));

			if (ValidCoord(j+1, i))
			{
				GetPixel(j + 1, i).SetClamp(ComponentClamp(GetPixel(j + 1, i).r + ALPHA * r_err),
							ComponentClamp(GetPixel(j + 1, i).g + ALPHA * g_err),
							ComponentClamp(GetPixel(j + 1, i).b + ALPHA * b_err));
			}

			if (ValidCoord(j - 1, i+1))
			{
				GetPixel(j - 1, i + 1).SetClamp(ComponentClamp(GetPixel(j - 1, i + 1).r + BETA * r_err),
								ComponentClamp(GetPixel(j - 1, i + 1).g + BETA * g_err),
								ComponentClamp(GetPixel(j - 1, i + 1).b + BETA * b_err));
			}

			if (ValidCoord(j, i + 1))
			{
				GetPixel(j, i + 1).SetClamp(ComponentClamp(GetPixel(j, i + 1).r + GAMMA * r_err),
								ComponentClamp(GetPixel(j, i+1).g + GAMMA * g_err),
								ComponentClamp(GetPixel(j, i+1).b + GAMMA * b_err));
			}

			if (ValidCoord(j + 1, i + 1))
			{
				GetPixel(j + 1, i + 1).SetClamp(ComponentClamp(GetPixel(j + 1, i + 1).r + DELTA * r_err),
								ComponentClamp(GetPixel(j + 1, i+1).g + DELTA * g_err),
								ComponentClamp(GetPixel(j + 1, i+1).b + DELTA * b_err));
			}

			//if (j != width - 1) // not rightest
			//{
			//	if (j == 0) // left most
			//	{
			//		if (i != height - 1) // not bottom
			//		{
			//			//left most not bottom
			//			//right
			//			GetPixel(j + 1, i).SetClamp(ComponentClamp(GetPixel(j + 1, i).r + ALPHA * r_err),
			//				ComponentClamp(GetPixel(j + 1, i).g + ALPHA * g_err),
			//				ComponentClamp(GetPixel(j + 1, i).b + ALPHA * b_err));
			//			//bottomleft
			//			GetPixel(width - 1, i+1).SetClamp(ComponentClamp(GetPixel(width - 1, i + 1).r + BETA * r_err),
			//				ComponentClamp(GetPixel(width - 1, i + 1).g + BETA * g_err),
			//				ComponentClamp(GetPixel(width - 1, i + 1).b + BETA * b_err));
			//			//bottom
			//			GetPixel(j, i+1).SetClamp(ComponentClamp(GetPixel(j, i + 1).r + GAMMA * r_err),
			//				ComponentClamp(GetPixel(j, i + 1).g + GAMMA * g_err),
			//				ComponentClamp(GetPixel(j, i + 1).b + GAMMA * b_err));
			//			//bottomright
			//			GetPixel(j+1, i + 1).SetClamp(ComponentClamp(GetPixel(j+1, i + 1).r + DELTA * r_err),
			//				ComponentClamp(GetPixel(j+1, i + 1).g + DELTA * g_err),
			//				ComponentClamp(GetPixel(j+1, i + 1).b + DELTA * b_err));

			//		}
			//		else // bottom leftmost
			//		{
			//			//right
			//			GetPixel(j + 1, i).SetClamp(ComponentClamp(GetPixel(j + 1, i).r + ALPHA * r_err),
			//				ComponentClamp(GetPixel(j + 1, i).g + ALPHA * g_err),
			//				ComponentClamp(GetPixel(j + 1, i).b + ALPHA * b_err));
			//			//bottomleft
			//			GetPixel(width - 1, 0).SetClamp(ComponentClamp(GetPixel(width - 1, 0).r + BETA * r_err),
			//				ComponentClamp(GetPixel(width - 1, 0).g + BETA * g_err),
			//				ComponentClamp(GetPixel(width - 1, 0).b + BETA * b_err));
			//			//bottom
			//			GetPixel(j, 0).SetClamp(ComponentClamp(GetPixel(j, 0).r + GAMMA * r_err),
			//				ComponentClamp(GetPixel(j, 0).g + GAMMA * g_err),
			//				ComponentClamp(GetPixel(j, 0).b + GAMMA * b_err));
			//			//bottomright
			//			GetPixel(j + 1, 0).SetClamp(ComponentClamp(GetPixel(j + 1, 0).r + DELTA * r_err),
			//				ComponentClamp(GetPixel(j + 1, 0).g + DELTA * g_err),
			//				ComponentClamp(GetPixel(j + 1, 0).b + DELTA * b_err));
			//		}
			//	}
			//	else // middle
			//	{
			//		if (i != height - 1) // not bottom
			//		{
			//			// middle not bottom
			//			//right
			//			GetPixel(j + 1, i).SetClamp(ComponentClamp(GetPixel(j + 1, i).r + ALPHA * r_err),
			//				ComponentClamp(GetPixel(j + 1, i).g + ALPHA * g_err),
			//				ComponentClamp(GetPixel(j + 1, i).b + ALPHA * b_err));
			//			//bottomleft
			//			GetPixel(j-1, i+1).SetClamp(ComponentClamp(GetPixel(j - 1, i + 1).r + BETA * r_err),
			//				ComponentClamp(GetPixel(j - 1, i + 1).g + BETA * g_err),
			//				ComponentClamp(GetPixel(j - 1, i + 1).b + BETA * b_err));
			//			//bottom
			//			GetPixel(j, i+1).SetClamp(ComponentClamp(GetPixel(j, i+1).r + GAMMA * r_err),
			//				ComponentClamp(GetPixel(j, i+1).g + GAMMA * g_err),
			//				ComponentClamp(GetPixel(j, i+1).b + GAMMA * b_err));
			//			//bottomright
			//			GetPixel(j + 1, i+1).SetClamp(ComponentClamp(GetPixel(j + 1, i+1).r + DELTA * r_err),
			//				ComponentClamp(GetPixel(j + 1, i+1).g + DELTA * g_err),
			//				ComponentClamp(GetPixel(j + 1, i+1).b + DELTA * b_err));
			//		}
			//		else // bottom middle
			//		{
			//			//right
			//			GetPixel(j + 1, i).SetClamp(ComponentClamp(GetPixel(j + 1, i).r + ALPHA * r_err),
			//				ComponentClamp(GetPixel(j + 1, i).g + ALPHA * g_err),
			//				ComponentClamp(GetPixel(j + 1, i).b + ALPHA * b_err));
			//			//bottomleft
			//			GetPixel(j - 1, 0).SetClamp(ComponentClamp(GetPixel(j - 1, 0).r + BETA * r_err),
			//				ComponentClamp(GetPixel(j - 1, 0).g + BETA * g_err),
			//				ComponentClamp(GetPixel(j - 1, 0).b + BETA * b_err));
			//			//bottom
			//			GetPixel(j, 0).SetClamp(ComponentClamp(GetPixel(j, 0).r + GAMMA * r_err),
			//				ComponentClamp(GetPixel(j, 0).g + GAMMA * g_err),
			//				ComponentClamp(GetPixel(j, 0).b + GAMMA * b_err));
			//			//bottomright
			//			GetPixel(j + 1, 0).SetClamp(ComponentClamp(GetPixel(j + 1, 0).r + DELTA * r_err),
			//				ComponentClamp(GetPixel(j + 1, 0).g + DELTA * g_err),
			//				ComponentClamp(GetPixel(j + 1, 0).b + DELTA * b_err));
			//		}
			//	}
			//}
			//else // rightest
			//{
			//	if (i != height - 1) // not bottom
			//	{
			//		// rightest not bottom
			//		//right
			//		GetPixel(0, i).SetClamp(ComponentClamp(GetPixel(0, i).r + ALPHA * r_err),
			//			ComponentClamp(GetPixel(0, i).g + ALPHA * g_err),
			//			ComponentClamp(GetPixel(0, i).b + ALPHA * b_err));
			//		//bottomleft
			//		GetPixel(j - 1, i+1).SetClamp(ComponentClamp(GetPixel(j - 1, i+1).r + BETA * r_err),
			//			ComponentClamp(GetPixel(j - 1, i+1).g + BETA * g_err),
			//			ComponentClamp(GetPixel(j - 1, i+1).b + BETA * b_err));
			//		//bottom
			//		GetPixel(j, i+1).SetClamp(ComponentClamp(GetPixel(j, i+1).r + GAMMA * r_err),
			//			ComponentClamp(GetPixel(j, i+1).g + GAMMA * g_err),
			//			ComponentClamp(GetPixel(j, i+1).b + GAMMA * b_err));
			//		//bottomright
			//		GetPixel(0, i+1).SetClamp(ComponentClamp(GetPixel(0, i+1).r + DELTA * r_err),
			//			ComponentClamp(GetPixel(0, i+1).g + DELTA * g_err),
			//			ComponentClamp(GetPixel(0, i+1).b + DELTA * b_err));
			//	}
			//	else // bottom rightest
			//	{
			//		//right
			//		GetPixel(0, i).SetClamp(ComponentClamp(GetPixel(0, i).r + ALPHA * r_err),
			//			ComponentClamp(GetPixel(0, i).g + ALPHA * g_err),
			//			ComponentClamp(GetPixel(0, i).b + ALPHA * b_err));
			//		//bottomleft
			//		GetPixel(j - 1, 0).SetClamp(ComponentClamp(GetPixel(j - 1, 0).r + BETA * r_err),
			//			ComponentClamp(GetPixel(j - 1, 0).g + BETA * g_err),
			//			ComponentClamp(GetPixel(j - 1, 0).b + BETA * b_err));
			//		//bottom
			//		GetPixel(j, 0).SetClamp(ComponentClamp(GetPixel(j, 0).r + GAMMA * r_err),
			//			ComponentClamp(GetPixel(j, 0).g + GAMMA * g_err),
			//			ComponentClamp(GetPixel(j, 0).b + GAMMA * b_err));
			//		//bottomright
			//		GetPixel(0, 0).SetClamp(ComponentClamp(GetPixel(0, 0).r + DELTA * r_err),
			//			ComponentClamp(GetPixel(0, 0).g + DELTA * g_err),
			//			ComponentClamp(GetPixel(0, 0).b + DELTA * b_err));
			//	}
			//}
		}
	}
}

void ImageComposite(Image *bottom, Image *top, Image *result)
{
  // Extra Credit (Section 3.7).
  // This hook just takes the top image and bottom image, producing a result
  // You might want to define a series of compositing modes as OpenGL does
  // You will have to use the alpha channel here to create Mattes
  // One idea is to composite your face into a famous picture
}

void Image::Convolve(int *filter, int n, int normalization, int absval) {
  // This is my definition of an auxiliary function for image convolution 
  // with an integer filter of width n and certain normalization.
  // The absval param is to consider absolute values for edge detection.
  
  // It is helpful if you write an auxiliary convolve function.
  // But this form is just for guidance and is completely optional.
  // Your solution NEED NOT fill in this function at all
  // Or it can use an alternate form or definition
}

void Image::Blur(int n)
{
  /* Your Work Here (Section 3.4.1) */
	double * coef = new double[n*n];
	float sigma = floor(n / 2.0) / 2.0;
	int mean = n / 2;
	double sum = 0;
	double PI = 3.1415926;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			coef[j * n + i] = exp(-0.5f * (pow((i - mean) / sigma, 2) + pow((j - mean) / sigma, 2))) / (2.0f * PI * pow(sigma, 2));
			sum += coef[j * n + i];
		}
	}

	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
			//int ctr = 0;
			double r = 0;
			double g = 0;
			double b = 0;
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
					//ctr++;
					int x = w - (j - mean);
					int y = h - (k - mean);
					x = (x >= 0) ? ((x >= width)? 2 * width - x - 1: x): -x - 1;
					y = (y >= 0) ? ((y >= height) ? 2 * height - y - 1 : y) : -y - 1;
					
					if (ValidCoord(x, y))
					{
						r = r + GetPixel(x, y).r * coef[k * n + j];
						g = g + GetPixel(x, y).g * coef[k * n + j];
						b = b + GetPixel(x, y).b * coef[k * n + j];
					}
				}
			}
			GetPixel(w, h).SetClamp(r/sum, g/sum, b/sum);
		}
	}
}

void Image::Sharpen() 
{
  /* Your Work Here (Section 3.4.2) */
	int * coef = new int[9];
	coef[0] = -1;
	coef[1] = -2;
	coef[2] = -1;
	coef[3] = -2;
	coef[4] = 19;
	coef[5] = -2;
	coef[6] = -1;
	coef[7] = -2;
	coef[8] = -1;
	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
			//int ctr = 0;
			int r = 0;
			int g = 0;
			int b = 0;
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					//ctr++;
					int x = w - (j - 1);
					int y = h - (k - 1);
					x = (x >= 0) ? ((x >= width) ? 2 * width - x - 1 : x) : -x - 1;
					y = (y >= 0) ? ((y >= height) ? 2 * height - y - 1 : y) : -y - 1;

					if (ValidCoord(x, y))
					{
						r = r + GetPixel(x, y).r * coef[k * 3 + j];
						g = g + GetPixel(x, y).g * coef[k * 3 + j];
						b = b + GetPixel(x, y).b * coef[k * 3 + j];
					}
				}
			}
			GetPixel(w, h).SetClamp(r / 7,g / 7,b / 7);
		}
	}
}

void Image::EdgeDetect(int threshold)
{
  /* Your Work Here (Section 3.4.3) */
	int * coefH = new int[9];
	coefH[0] = -1;
	coefH[1] = 0;
	coefH[2] = 1;
	coefH[3] = -2;
	coefH[4] = 0;
	coefH[5] = 2;
	coefH[6] = -1;
	coefH[7] = 0;
	coefH[8] = 1;

	Pixel * hor_pix = new Pixel[num_pixels];
	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
			//int ctr = 0;
			int r = 0;
			int g = 0;
			int b = 0;
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					//ctr++;
					int x = w - (j - 2);
					int y = h - (k - 2);
					x = (x >= 0) ? ((x >= width)? width - 1: x): 0;
					y = (y >= 0) ? ((y >= height) ? height - 1 : y) : 0;

					if (ValidCoord(x, y))
					{
						r = r + GetPixel(x, y).r * coefH[k * 3 + j];
						g = g + GetPixel(x, y).g * coefH[k * 3 + j];
						b = b + GetPixel(x, y).b * coefH[k * 3 + j];
					}
				}
			}
			hor_pix[h*width + w].SetClamp(abs(r), abs(g), abs(b));
		}
	}

	int * coefV = new int[9];
	coefV[0] = 1;
	coefV[1] = 2;
	coefV[2] = 1;
	coefV[3] = 0;
	coefV[4] = 0;
	coefV[5] = 0;
	coefV[6] = -1;
	coefV[7] = -2;
	coefV[8] = -1;

	Pixel * ver_pix = new Pixel[num_pixels];
	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
			//int ctr = 0;
			int r = 0;
			int g = 0;
			int b = 0;
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					//ctr++;
					int x = w - (j - 2);
					int y = h - (k - 2);
					x = (x >= 0) ? ((x >= width) ? width - 1 : x) : 0;
					y = (y >= 0) ? ((y >= height) ? height - 1 : y) : 0;
					
					if (ValidCoord(x, y))
					{
						r = r + GetPixel(x, y).r * coefV[k * 3 + j];
						g = g + GetPixel(x, y).g * coefV[k * 3 + j];
						b = b + GetPixel(x, y).b * coefV[k * 3 + j];
					}
				}
			}
			ver_pix[h*width + w].SetClamp(abs(r), abs(g), abs(b));
		}
	}


	double gradient;
	for (int i = 0; i < num_pixels; i++)
	{
		gradient = sqrt(pow(hor_pix[i].Luminance(),2) + pow(ver_pix[i].Luminance(),2));
		(gradient > threshold) ? pixels[i].SetClamp(255, 255, 255) : pixels[i].SetClamp(0, 0, 0);
	}
}


Image* Image::Scale(int sizex, int sizey)
{
  /* Your Work Here (Section 3.5.1) */
	Image * new_img = new Image(sizex,sizey);
	double xScale = (double)sizex / (double)width;
	double yScale = (double)sizey / (double)height;
	for (int i = 0; i < sizey; i++)
	{
		for (int j = 0; j < sizex; j++)
		{
			double u = j / xScale;
			double v = i / yScale;
			new_img->GetPixel(j,i)=Sample(u,v,xScale,yScale,1);
		}
	}
  return new_img;
}

void Image::Shift(double sx, double sy)
{
  /* Your Work Here (Section 3.5.2) */
	Pixel * pixs = new Pixel[num_pixels];
	bool used = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			int newX, newY;
			newX = (sx > 0) ? width - 1 - j : j;
			newY = (sy > 0) ? height - 1 - i : i;
			if (sx == floor(sx) && sy == floor(sy))
			{
				GetPixel(newX, newY) = (ValidCoord(newX - sx, newY - sy)) ? GetPixel(newX - sx, newY - sy) : Pixel(0, 0, 0);
			}
			else
			{
				//GetPixel(j,i) = Sample(j,i,sx,sy,0);
				used = 1;
				pixs[j + i * width] = Sample(j, i, sx, sy, 0);
			}
		}
	}
	if (used)
	{
		memcpy(pixels, pixs, width * height * sizeof(Pixel));
	}
	delete[] pixs;
}


/*
Image* Image::Rotate(double angle)
{
  // For rotation of the image
  // Not required in the assignment
  // But you can earn limited extra credit if you fill it in
  // (It isn't really that hard) 

    return NULL;
}
*/
int getMedian(int * arr) {
	int n = 9;
	bool swapped = true;
	while (swapped)
	{
		swapped = false;
		for (int i = 1; i <= n-1; i++)
		{
			if (arr[i-1] > arr[i])
			{
				int tmp = arr[i - 1];
				arr[i - 1] = arr[i];
				arr[i] = tmp;
				swapped = true;
			}
		}
		n = n - 1;
	}
	return arr[(int)(9 / 2)];
}


void Image::Fun()
{
    /* Your Work Here (Section 3.6) */
	int median = 1;
	for (int i = median; i < height - median; i++)
	{
		for (int j = median; j < width - median; j++)
		{
			int * window_r = new int[9];
			int * window_g = new int[9];
			int * window_b = new int[9];
			for (int x = 0; x < 3; x++)
			{
				for (int y = 0; y < 3; y++)
				{
					window_r[x + y * 3] = GetPixel(j + x - 1, i + y - 1).r;
					window_g[x + y * 3] = GetPixel(j + x - 1, i + y - 1).g;
					window_b[x + y * 3] = GetPixel(j + x - 1, i + y - 1).b;
				}
			}
			GetPixel(j, i).SetClamp(getMedian(window_r), getMedian(window_g), getMedian(window_b));
		}
	}
}


Image* ImageMorph (Image* I0, Image* I1, int numLines, Line* L0, Line* L1, double t)
{
  /* Your Work Here (Section 3.7) */
  // This is extra credit.
  // You can modify the function definition. 
  // This definition takes two images I0 and I1, the number of lines for 
  // morphing, and a definition of corresponding line segments L0 and L1
  // t is a parameter ranging from 0 to 1.
  // For full credit, you must write a user interface to join corresponding 
  // lines.
  // As well as prepare movies 
  // An interactive slider to look at various morph positions would be good.
  // From Beier-Neely's SIGGRAPH 92 paper

    return NULL;
}


/**
 * Image Sample
 **/
void Image::SetSamplingMethod(int method)
{
  // Sets the filter to use for Scale and Shift
  // You need to implement point sampling, hat filter and mitchell

    assert((method >= 0) && (method < IMAGE_N_SAMPLING_METHODS));
    sampling_method = method;
}

Pixel Image::Sample (double u, double v, double sx, double sy, int type)
{
  // To sample the image in scale and shift
  // This is an auxiliary function that it is not essential you fill in or 
  // you may define it differently.
  // u and v are the floating point coords of the points to be sampled.
  // sx and sy correspond to the scale values. 
  // In the assignment, it says implement MinifyX MinifyY MagnifyX MagnifyY
  // separately.  That may be a better way to do it.
  // This hook is primarily to get you thinking about that you have to have 
  // some equivalent of this function.
  if (sampling_method == IMAGE_SAMPLING_POINT) {
    // Your work here
	  float x, y;
	  x = (type) ? u : u - sx;
	  y = (type) ? v : v - sy;
	  //x = (x - floor(x) >= 0.5) ? ceil(x) : floor(x);
	  //y = (y- floor(y) >= 0.5) ? ceil(y) : floor(y);
	  x = floor(x);
	  y = floor(y);
	  if (x < 0 || x >= width || y < 0 || y >= height)
	  {
		  return Pixel(0, 0, 0);
	  }

	  return GetPixel((int)x,(int)y);
  }

  else if (sampling_method == IMAGE_SAMPLING_HAT) {
    // Your work here
	  float x, y;
	  float xw, yw;
	  x = (type) ? u: u-sx;
	  y = (type) ? v : v - sy;
	  xw = (type) ? ((sx < 1) ? 1.0f / sx : 1.0) : 1.0;
	  yw = (type) ? ((sy < 1) ? 1.0f / sy : 1.0) : 1.0;
	  if (x < 0 || x >= width || y < 0 || y >=height)
	  {
		  return Pixel(0, 0, 0);
	  }
	  float r, g, b;
	  float nom;
	  r = 0;
	  g = 0;
	  b = 0;
	  nom = 0;
	  for (int i = ceil(x - xw); i < ceil(x + xw); i++)
	  {
		  for (int j = ceil(y - yw); j < ceil(y + yw); j++)
		  {
			  if (ValidCoord(i, j))
			  {
				  float hatX = abs((i - x) / xw);
				  float hatY = abs((j - y) / yw);
				  hatX = (hatX <= 1) ? 1.0 - hatX : 0.0;
				  hatY = (hatY <= 1) ? 1.0 - hatY : 0.0;
				  float ret = hatX * hatY;
				  nom += ret;
				  r += ret * GetPixel(i, j).r;
				  g += ret * GetPixel(i, j).g;
				  b += ret * GetPixel(i, j).b;
			  }
		  }
	  }
	  r /= nom;
	  g /= nom;
	  b /= nom;
	  r = (r > 255) ? 255 : ((r < 0) ? 0 : r);
	  g = (g > 255) ? 255 : ((g < 0) ? 0 : g);
	  b = (b > 255) ? 255 : ((b < 0) ? 0 : b);
	  return Pixel(r, g, b);
  }

  else if (sampling_method == IMAGE_SAMPLING_MITCHELL) {
    // Your work here
	  float x, y;
	  float xw, yw;
	  x = (type) ? u : u - sx;
	  y = (type) ? v : v - sy;
	  xw = (type) ? ((sx < 1) ? 1.0f / sx : 1.0) : 1.0;
	  yw = (type) ? ((sy < 1) ? 1.0f / sy : 1.0) : 1.0;
	  if (x < 0 || x >= width || y < 0 || y >= height)
	  {
		  return Pixel(0, 0, 0);
	  }
	  float r, g, b;
	  float nom;
	  r = 0;
	  g = 0;
	  b = 0;
	  nom = 0;
	  for (double i = ceil(x - xw); i < ceil(x + xw); i++)
	  {
		  for (double j = ceil(y - yw); j < ceil(y + yw); j++)
		  {
			  if (ValidCoord(i, j)) {
				  float mitX = fabs((i - x) / xw);
				  float mitY = fabs((j - y) / yw);
				  if (mitX < 1)
				  {
					  mitX = ((7.0f * pow(mitX, 3)) - (12.0f * pow(mitX, 2)) + (16.0f / 3.0f)) / 6.0f;
				  }
				  else if (mitX < 2 && mitX >= 1)
				  {
					  mitX = (((-7.0f / 3.0f) * pow(mitX, 3)) + (12.0f * pow(mitX, 2)) - (20.0f * mitX) + (32.0f / 3.0f)) / 6.0f;
				  }
				  if (mitY < 1)
				  {
					  mitY = ((7.0f * pow(mitY, 3)) - (12.0f * pow(mitY, 2)) + (16.0f / 3.0f)) / 6.0f;
				  }
				  else if (mitY < 2 && mitY >= 1)
				  {
					  mitY = (((-7.0f / 3.0f) * pow(mitY, 3)) + (12.0f * pow(mitY, 2)) - (20.0f * mitY) + (32.0f / 3.0f)) / 6.0f;
				  }
				  float ret = mitX * mitY;
				  nom += ret;
				  r += ret * (float)GetPixel(i, j).r;
				  g += ret * (float)GetPixel(i, j).g;
				  b += ret * (float)GetPixel(i, j).b;
			  }
		  }
	  }
	  r /= nom;
	  g /= nom;
	  b /= nom;
	  r = (r > 255) ? 255 : ((r < 0) ? 0 : r);
	  g = (g > 255) ? 255 : ((g < 0) ? 0 : g);
	  b = (b > 255) ? 255 : ((b < 0) ? 0 : b);
	  return Pixel(r, g, b);

  }

  else {
    fprintf(stderr,"I don't understand what sampling method is used\n") ;
    exit(1) ;
  }

  return Pixel() ;
}

