/*
 * PNGIO.cpp
 *
 *  Created on: 31.01.2011
 *      Author: jagiella
 */

#include "PNG.hpp"


void abort_(const char * s, ...)
{
        va_list args;
        va_start(args, s);
        vfprintf(stderr, s, args);
        fprintf(stderr, "\n");
        va_end(args);
        abort();
}


/*void PNG::read_png_file(char* file_name)
{
	png_byte header[8];    // 8 is the maximum size that can be checked

        / * open file and test for it being a png
        FILE *fp = fopen(file_name, "rb");
        if (!fp)
                abort_("[read_png_file] File %s could not be opened for reading", file_name);
        fread(header, 1, 8, fp);
        if( png_sig_cmp(header, 0, 8))
                abort_("[read_png_file] File %s is not recognized as a PNG file", file_name);


        / * initialize stuff
        png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

        if (!png_ptr)
                abort_("[read_png_file] png_create_read_struct failed");

        info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr)
                abort_("[read_png_file] png_create_info_struct failed");

        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[read_png_file] Error during init_io");

        png_init_io(png_ptr, fp);
        png_set_sig_bytes(png_ptr, 8);

        png_read_info(png_ptr, info_ptr);

        width = png_get_image_width(png_ptr, info_ptr);
        height = png_get_image_height(png_ptr, info_ptr);
        color_type = png_get_color_type(png_ptr, info_ptr);
        bit_depth = png_get_bit_depth(png_ptr, info_ptr);

        number_of_passes = png_set_interlace_handling(png_ptr);
        png_read_update_info(png_ptr, info_ptr);


        // read file
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[read_png_file] Error during read_image");

        // ROW-WISE READING
       / * row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
        for (y=0; y<height; y++)
                row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr,info_ptr));
* /
        // IMAGE-WISE READING
        unsigned long rowbytes = png_get_rowbytes(png_ptr, info_ptr);
        //png_byte *image_data; // vector: png_byte * rowbytes * height
    	if ((image_data = (png_byte *) malloc(rowbytes * height)) == NULL) {
    		png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    		return;
    	}
    	// matrix: png_bytep * height x rowbytes
    	row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    	for (int i = 0; i < height; i++)
    		row_pointers[i] = image_data + i * rowbytes;
    		//row_pointers[i] = &image_data[i * rowbytes];


        png_read_image(png_ptr, row_pointers);
        png_read_end(png_ptr, NULL);

        fclose(fp);
}*/

bool PNG::read_png_file(char* file_name)
{
	png_byte header[8];    // 8 is the maximum size that can be checked

        /* open file and test for it being a png */
        FILE *fp = fopen(file_name, "rb");
        if (!fp){
                //abort_("[read_png_file] File %s could not be opened for reading", file_name);
        	fprintf(stderr, "[read_png_file] File %s could not be opened for reading\n", file_name);
        	return false;
        }
        fread(header, 1, 8, fp);
        if( png_sig_cmp(header, 0, 8))
                abort_("[read_png_file] File %s is not recognized as a PNG file", file_name);


        /* initialize stuff */
        png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

        if (!png_ptr)
                abort_("[read_png_file] png_create_read_struct failed");

        info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr)
                abort_("[read_png_file] png_create_info_struct failed");

        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[read_png_file] Error during init_io");

        png_init_io(png_ptr, fp);
        png_set_sig_bytes(png_ptr, 8);

        png_read_info(png_ptr, info_ptr);

        width = png_get_image_width(png_ptr, info_ptr);
        //fprintf(stderr, "width=%li\n", width);
        height = png_get_image_height(png_ptr, info_ptr);
        //fprintf(stderr, "height=%li\n", height);
        color_type = png_get_color_type(png_ptr, info_ptr);
        if(color_type==PNG_COLOR_TYPE_GRAY)
        	fprintf(stderr, "color_type=PNG_COLOR_TYPE_GRAY\n");
        else if(color_type==PNG_COLOR_TYPE_RGB)
        	fprintf(stderr, "color_type=PNG_COLOR_TYPE_RGB\n");
        else if(color_type==PNG_COLOR_TYPE_PALETTE)
        	fprintf(stderr, "color_type=PNG_COLOR_TYPE_PALETTE\n");
        else
           	fprintf(stderr, "color_type=%li\n", color_type);
        bit_depth = png_get_bit_depth(png_ptr, info_ptr);
        fprintf(stderr, "bit_depth=%li\n", bit_depth);

        number_of_passes = png_set_interlace_handling(png_ptr);
        png_read_update_info(png_ptr, info_ptr);


        /* read file */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[read_png_file] Error during read_image");

        // ROW-WISE READING
       /* row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
        for (y=0; y<height; y++)
                row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr,info_ptr));
*/
        // IMAGE-WISE READING
        unsigned long rowbytes = png_get_rowbytes(png_ptr, info_ptr);
        //fprintf(stderr, "row bytes=%li\n", rowbytes);
        //png_byte *image_data; // vector: png_byte * rowbytes * height
    	if ((image_data = (png_byte *) malloc(rowbytes * height)) == NULL) {
    		png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    		return false;
    	}
    	// matrix: png_bytep * height x rowbytes
    	row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    	for (int i = 0; i < height; i++)
    		row_pointers[i] = image_data + i * rowbytes;
    		//row_pointers[i] = &image_data[i * rowbytes];


        png_read_image(png_ptr, row_pointers);
        png_read_end(png_ptr, NULL);

        fclose(fp);

        return true;
}

void PNG::write_png_file(char* file_name)
{
        /* create file */
        FILE *fp = fopen(file_name, "wb");
        if (!fp)
                abort_("[write_png_file] File %s could not be opened for writing", file_name);


        /* initialize stuff */
        png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

        if (!png_ptr)
                abort_("[write_png_file] png_create_write_struct failed");

        info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr)
                abort_("[write_png_file] png_create_info_struct failed");

        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during init_io");

        png_init_io(png_ptr, fp);


        /* write header */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during writing header");

        png_set_IHDR(png_ptr, info_ptr, width, height,
                     bit_depth, color_type, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

        png_write_info(png_ptr, info_ptr);


        /* write bytes */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during writing bytes");

        png_write_image(png_ptr, row_pointers);


        /* end write */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during end of write");

        png_write_end(png_ptr, NULL);

        /* cleanup heap allocation */
        /*for (y=0; y<height; y++)
                free(row_pointers[y]);
        free(row_pointers);*/
        //free(image_data);

        fclose(fp);
}


void PNG::process_file(void)
{
        if (png_get_color_type(png_ptr, info_ptr) == PNG_COLOR_TYPE_RGB)
                abort_("[process_file] input file is PNG_COLOR_TYPE_RGB but must be PNG_COLOR_TYPE_RGBA "
                       "(lacks the alpha channel)");

        if (png_get_color_type(png_ptr, info_ptr) != PNG_COLOR_TYPE_RGBA)
                abort_("[process_file] color_type of input file must be PNG_COLOR_TYPE_RGBA (%d) (is %d)",
                       PNG_COLOR_TYPE_RGBA, png_get_color_type(png_ptr, info_ptr));

        for (y=0; y<height; y++) {
                png_byte* row = row_pointers[y];
                for (x=0; x<width; x++) {
                        png_byte* ptr = &(row[x*4]);
                        printf("Pixel at position [ %d - %d ] has RGBA values: %d - %d - %d - %d\n",
                               x, y, ptr[0], ptr[1], ptr[2], ptr[3]);

                        /* set red value to 0 and green value to the blue one */
                        ptr[0] = 0;
                        ptr[1] = ptr[2];
                }
        }
}

png_byte PNG::getPixel( int x, int y, int channel)
{
	png_byte* row = row_pointers[y];

	if( this->color_type==PNG_COLOR_TYPE_GRAY){
		return row[x];
	}else if( this->color_type==PNG_COLOR_TYPE_RGB){
		png_byte* ptr = &(row[x*3]);
		return ptr[channel];
	}

	fprintf( stderr, "ERROR: getPixel() not defined for color_type %i\n", color_type);
	return 0;
}

void binaer(unsigned short i, int size) {

	if (size > 1) {
		binaer(i >> 1, size - 1);
	}

	fprintf( stderr, "%d", (i & 1));
}

float PNG::getPixelIntensity( int x, int y, int channel)
{
	void* row = row_pointers[y];

	if( this->color_type==PNG_COLOR_TYPE_GRAY){
		if( this->bit_depth == 8){
			unsigned char value = ((unsigned char*)row)[x];
			return ((float)value / 256.);
		}else if( this->bit_depth == 16){
			//fprintf( stderr, "[%i]\n", ((unsigned short*)row)[x]);
			//binaer( ((unsigned short*)row)[x], 16);
			//fprintf( stderr, "\n");
			//int test = ((unsigned char*)row)[x*3];
			printf( "%i %i %i %i %lf\n", x,y,((unsigned char*)row)[x*2], ((unsigned char*)row)[x*2+1], ((unsigned short*)row)[x]/65536.);
			unsigned short value = ((unsigned short*)row)[x];
			return ((float)value / 65536.);
		}
	}else if( this->color_type==PNG_COLOR_TYPE_RGB){
		if( this->bit_depth == 8){
			unsigned char value = ((unsigned char*)row)[x*3 + channel];
			return ((float)value / 256.);
		}else if( this->bit_depth == 16){
			unsigned short value = ((unsigned short*)row)[x*3 + channel];
			return ((float)value / 65536.);
		}
	}

	fprintf( stderr, "ERROR: getPixelIntensity(%i,%i) not defined for color_type %i and/or bit_depth %i\n", x,y,color_type, bit_depth);
	return 0;
}

int PNG::getHeight()
{
	return this->height;
}

int PNG::getWidth()
{
	return this->width;
}


void PNG::readpng_cleanup(int free_image_data)
{
	if (png_ptr && info_ptr) {
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        png_ptr = NULL;
        info_ptr = NULL;
    }
    if (free_image_data && image_data) {
        free(image_data);
        image_data = NULL;
        free( row_pointers);
    }
}


void PNG::writepng_cleanup(int free_image_data)
{
	if (png_ptr && info_ptr) {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        png_ptr = NULL;
        info_ptr = NULL;
    }
    if (free_image_data && image_data) {
        free(image_data);
        image_data = NULL;
        free( row_pointers);
   }
}
