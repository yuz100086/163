#include "image.h"
#include "bmp.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>


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
	if (factor == 1)
	{
		return;
	}

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
						   (1.0 - factor) * lum + pixels[i].r * factor,
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
		pixels[i].SetClamp(pow(pixels[i].r, 1.0/factor),
						   pow(pixels[i].g, 1.0 / factor),
						   pow(pixels[i].b, 1.0 / factor));
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
		pixels[i].SetClamp((255.0 / num) * floor(pixels[i].r * num / 256.0),
						   (255.0 / num) * floor(pixels[i].g * num / 256.0),
						   (255.0 / num) * floor(pixels[i].b * num / 256.0));
	}
}


void Image::RandomDither (int nbits)
{
  /* Your Work Here (Section 3.3.2) */
	srand(time_t(NULL));
	for (int i = 0; i < num_pixels; i++)
	{
		double num = pow(2.0, nbits);
		pixels[i].SetClamp((255.0 / num) * floor(pixels[i].r * num / 256.0 + 0.5 - static_cast <float> (rand()) / static_cast <float> (RAND_MAX)),
			(255.0 / num) * floor(pixels[i].g * num / 256.0 + 0.5 - static_cast <float> (rand()) / static_cast <float> (RAND_MAX)),
			(255.0 / num) * floor(pixels[i].b * num / 256.0 + 0.5 - static_cast <float> (rand()) / static_cast <float> (RAND_MAX)));
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
}

void Image::Sharpen() 
{
  /* Your Work Here (Section 3.4.2) */
}

void Image::EdgeDetect(int threshold)
{
  /* Your Work Here (Section 3.4.3) */
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

