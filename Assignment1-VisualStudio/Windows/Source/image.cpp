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
	for (int i = 0; i < num_pixels; i++)
	{
		double num = pow(2.0, nbits);
		double r_err = pixels[i].r - (255.0 / num) * floor(pixels[i].r * num / 256.0);
		double g_err = pixels[i].g - (255.0 / num) * floor(pixels[i].g * num / 256.0);
		double b_err = pixels[i].b - (255.0 / num) * floor(pixels[i].b * num / 256.0);

		pixels[i].SetClamp(pixels[i].r - (255.0 / num) * floor(pixels[i].r * num / 256.0),
			(255.0 / num) * floor(pixels[i].g * num / 256.0),
			(255.0 / num) * floor(pixels[i].b * num / 256.0));

		// left most
		if (i % width == 0)
		{
			if (i >= (height - 1) * width)
			{

			}
			else
			{

			}
		}
		else if (i % (width - 1) == 0) // right most
		{
			if (i >= (height - 1) * width)
			{

			}
			else
			{

			}
		}
		else
		{
			// bottom
			if (i >= (height - 1) * width)
			{

			}
			else
			{
				//right
				pixels[i + 1].SetClamp(pixels[i + 1].r + r_err * ALPHA, pixels[i + 1].g + g_err * ALPHA, pixels[i + 1].b + b_err * ALPHA);
				//bottom
				pixels[i + width].SetClamp(pixels[i + width].r + r_err * GAMMA, pixels[i + width].g + g_err * GAMMA, pixels[i + width].b + b_err * GAMMA);
				//bottom left
				pixels[i + width + 1].SetClamp(pixels[i + width + 1].r + r_err * BETA, pixels[i + width + 1].g + g_err * BETA, pixels[i + width + 1].b + b_err * BETA);
				//bottom right
				pixels[i + width - 1].SetClamp(pixels[i + width - 1].r + r_err * DELTA, pixels[i + width - 1].g + g_err * DELTA, pixels[i + width - 1].b + b_err * DELTA);
			}
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
/*
void Image::Convolve(int *filter, int n, int normalization, int absval) {
  // This is my definition of an auxiliary function for image convolution 
  // with an integer filter of width n and certain normalization.
  // The absval param is to consider absolute values for edge detection.
  
  // It is helpful if you write an auxiliary convolve function.
  // But this form is just for guidance and is completely optional.
  // Your solution NEED NOT fill in this function at all
  // Or it can use an alternate form or definition
}*/

void Image::Blur(int n)
{
  /* Your Work Here (Section 3.4.1) */
	float sigma = floor(n / 2.0) / 2.0;
	float constan = 1.0 / (2.0 * 3.14 * pow(sigma, 2));
	for (int i = 0; i < num_pixels; i++)
	{
		
		Pixel tmp_pix = pixels[i];
		tmp_pix.SetClamp(0,0,0);
		int col = i % width;
		int row = floor(i / width);
		//std::vector<std::vector<int>> arr;
		for (int j = -1; j < (n-1); j++)
		{
			for (int k = -1; k < (n-1); k++)
			{
				
				if (ValidCoord(row + j, col + k))
				{
					int mul = constan * exp(pow((j+1), 2) + pow((k+1), 2) / (-2.0 * pow(sigma, 2)));

					tmp_pix.SetClamp(tmp_pix.r + GetPixel(row + j, col + k).r * (mul),
						tmp_pix.g + GetPixel(row + j, col + k).g * (mul),
							tmp_pix.b + GetPixel(row + j, col + k).b * (mul));
				}
				
			}
		}
		pixels[i] = tmp_pix;
	}
}
int LefRem( int x ,int lim) {
	if (x < 0)
	{
		return -x - 1;
	}
	else if (x >= lim)
	{
		return 2 * lim - x - 1;
	}
	else
	{
		return x;
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
					int x = LefRem(w - (j - 1), width);
					int y = LefRem(h - (k - 1), width);
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


int Lim(int x, int lim) {
	if (x < 0)
		return 0;
	else if (x >= lim)
		return lim - 1;
	else
		return x;
}
void Image::Convolve(int *filter, int n, int normalization, int absval) {
	// This is my definition of an auxiliary function for image convolution 
	// with an integer filter of width n and certain normalization.
	// The absval param is to consider absolute values for edge detection.

	// It is helpful if you write an auxiliary convolve function.
	// But this form is just for guidance and is completely optional.
	// Your solution NEED NOT fill in this function at all
	// Or it can use an alternate form or definition
	int sumr, sumg, sumb;
	sumr = sumg = sumb = 0;

	for (int h = 0; h < height; h++) {
		for (int w = 0; w < width; w++) {
			for (int p = 0; p < n; p++) {
				for (int q = 0; q < n; q++) {
					int x = 0;
					int y = 0;
					int mp = 0;
					if (absval) {
						mp = 2;
						x = Lim(w - (p - mp), width);
						y = Lim(h - (q - mp), height);
					}
					else {
						mp = n / 2;
						x = LefRem(w - (p - mp), width);
						y = LefRem(h - (q - mp), height);
					}
					if (!ValidCoord(x, y)) continue;
					Pixel curr = GetPixel(x, y);
					int filt = filter[q * n + p];
					sumr += (int)curr.r * filt;
					sumg += (int)curr.g * filt;
					sumb += (int)curr.b * filt;
				}
			}
			if (absval) {
				sumr = abs(sumr);
				sumg = abs(sumg);
				sumb = abs(sumb);
			}
			GetPixel(w, h).SetClamp(sumr / normalization, sumg / normalization, sumb / normalization);
			sumr = sumg = sumb = 0;
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
					int x = Lim(w - (j - 2), width);
					int y = Lim(h - (k - 2), height);
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
					int x = Lim(w - (j - 2), width);
					int y = Lim(h - (k - 2), height);
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

  return NULL ;
}

void Image::Shift(double sx, double sy)
{
  /* Your Work Here (Section 3.5.2) */
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


void Image::Fun()
{
    /* Your Work Here (Section 3.6) */
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

Pixel Image::Sample (double u, double v, double sx, double sy)
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
  }

  else if (sampling_method == IMAGE_SAMPLING_HAT) {
    // Your work here
  }

  else if (sampling_method == IMAGE_SAMPLING_MITCHELL) {
    // Your work here
  }

  else {
    fprintf(stderr,"I don't understand what sampling method is used\n") ;
    exit(1) ;
  }

  return Pixel() ;
}

