/*
 * PNGIO.hpp
 *
 *  Created on: 31.01.2011
 *      Author: jagiella
 */

#ifndef PNGIO_HPP_
#define PNGIO_HPP_

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#define PNG_DEBUG 3
#include <png.h>

/** \addtogroup IO
 *  @{
 */

class PNG{
	int x, y;

	int width, height;
	png_byte color_type;
	png_byte bit_depth;

	png_structp png_ptr;
	png_infop info_ptr;
	int number_of_passes;
	png_bytep * row_pointers;
	png_byte *image_data;

public:
	bool read_png_file(char* file_name);
	void write_png_file(char* file_name);

	void readpng_cleanup(int free_image_data);
	void writepng_cleanup(int free_image_data);

	void process_file(void);
	png_byte getPixel( int x, int y, int channel);
	float getPixelIntensity( int x, int y, int channel);

	int getHeight();
	int getWidth();
};

/** @}*/


#endif /* PNGIO_HPP_ */
