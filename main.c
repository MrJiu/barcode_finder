/***************************************************************************
 *   v4l2grab Version 0.2                                                  *
 *   Copyright (C) 2009 by Tobias Müller                                   *
 *   Tobias_Mueller@twam.info                                              *
 *                                                                         *
 *   based on V4L2 Specification, Appendix B: Video Capture Example        *
 *   (http://v4l2spec.bytesex.org/spec/capture-example.html)               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


// use ITU-R float conversion for YUV422toRGB888 by default
#if !defined(ITU_R_FLOAT) && !defined(ITU_R_INT) && !defined(NTSC)
#define ITU_R_FLOAT
#endif

#if ((defined(ITU_R_FLOAT)) && (defined(ITU_R_INT)) && (defined(NTSC))) || ((defined(ITU_R_FLOAT)) && (defined(ITU_R_INT))) || ((defined(ITU_R_FLOAT)) && (defined(NTSC))) ||  ((defined(ITU_R_INT)) && (defined(NTSC)))
#error Only one conversion for YUV422toRGB888 is allowed!
#endif

// compile with all three access methods
#if !defined(IO_READ) && !defined(IO_MMAP) && !defined(IO_USERPTR)
#define IO_READ
#define IO_MMAP
#define IO_USERPTR
#endif

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <assert.h>
#include <getopt.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <malloc.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/ioctl.h>
#include <asm/types.h>
#include <linux/videodev2.h>
#include <jpeglib.h>

//#define NUM_FFT 64
#include "fft/fft_brin.h"

#define CLEAR(x) memset (&(x), 0, sizeof (x))

typedef enum {
#ifdef IO_READ
        IO_METHOD_READ,
#endif
#ifdef IO_MMAP
        IO_METHOD_MMAP,
#endif
#ifdef IO_USERPTR
        IO_METHOD_USERPTR,
#endif
} io_method;

struct buffer {
        void *                  start;
        size_t                  length;
};

struct surface {
	unsigned char *buf;
	unsigned int width;
	unsigned int height;
	unsigned int bytes_per_pixel;
};

struct point {
	unsigned int x;
	unsigned int y;
};

static io_method        io              = IO_METHOD_MMAP;
static int              fd              = -1;
struct buffer *         buffers         = NULL;
static unsigned int     n_buffers       = 0;

// global settings
static unsigned int width = 640;
static unsigned int height = 480;
static unsigned char jpegQuality = 70;
static char* jpegFilename = NULL;
static int swap_pixels = 0;
static int do_fft = 0;
static int find_barcode = 0;
static int barcode_algorithm;
static char* deviceName = "/dev/video0";
static char use_mjpeg = 1; // try to use MJPEG data, if that fails fall back to YUV

static void captureStop(void);
static void deviceUninit(void);
static void deviceClose(void);

void sighandler(int signum)
{
	printf("Caught signal %d, clean up and exit. TODO: set event flag and clean up in appropriate places...\n", signum);
	captureStop();
	deviceUninit();
	deviceClose();
	exit(EXIT_SUCCESS);
}

unsigned char jpg_header[] = {
	0xff, 0xd8,                   // SOI
	0xff, 0xe0,                   // APP0
	0x00, 0x10,                   // APP0 Hdr size
	0x4a, 0x46, 0x49, 0x46, 0x00, // ID string: 'JFIF' + null terminator
	0x01, 0x01,                   // Version
	0x00,                         // Bits per type
	0x00, 0x00,                   // X density
	0x00, 0x00,                   // Y density
	0x00,                         // X Thumbnail size
	0x00                          // Y Thumbnail size
};

/* JPEG DHT Segment for YCrCb omitted from MJPEG data */
static unsigned char jpeg_dht_seg[] = {
	0xff, 0xc4, 0x01, 0xa2,

	0x00, 0x00, 0x01, 0x05, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
	0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b,

	0x01, 0x00, 0x03, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00,
	0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b,

	0x10, 0x00, 0x02, 0x01, 0x03, 0x03, 0x02, 0x04, 0x03, 0x05, 0x05, 0x04, 0x04, 0x00, 0x00, 0x01, 0x7d,
	0x01, 0x02, 0x03, 0x00, 0x04, 0x11, 0x05, 0x12, 0x21, 0x31, 0x41, 0x06, 0x13, 0x51, 0x61, 0x07,
	0x22, 0x71, 0x14, 0x32, 0x81, 0x91, 0xa1, 0x08, 0x23, 0x42, 0xb1, 0xc1, 0x15, 0x52, 0xd1, 0xf0,
	0x24, 0x33, 0x62, 0x72, 0x82, 0x09, 0x0a, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x25, 0x26, 0x27, 0x28,
	0x29, 0x2a, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3a, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48, 0x49,
	0x4a, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x5a, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69,
	0x6a, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79, 0x7a, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89,
	0x8a, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9a, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7,
	0xa8, 0xa9, 0xaa, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xc2, 0xc3, 0xc4, 0xc5,
	0xc6, 0xc7, 0xc8, 0xc9, 0xca, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7, 0xd8, 0xd9, 0xda, 0xe1, 0xe2,
	0xe3, 0xe4, 0xe5, 0xe6, 0xe7, 0xe8, 0xe9, 0xea, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8,
	0xf9, 0xfa,

	0x11, 0x00, 0x02, 0x01, 0x02, 0x04, 0x04, 0x03, 0x04, 0x07, 0x05, 0x04, 0x04, 0x00, 0x01, 0x02, 0x77,
	0x00, 0x01, 0x02, 0x03, 0x11, 0x04, 0x05, 0x21, 0x31, 0x06, 0x12, 0x41, 0x51, 0x07, 0x61, 0x71,
	0x13, 0x22, 0x32, 0x81, 0x08, 0x14, 0x42, 0x91, 0xa1, 0xb1, 0xc1, 0x09, 0x23, 0x33, 0x52, 0xf0,
	0x15, 0x62, 0x72, 0xd1, 0x0a, 0x16, 0x24, 0x34, 0xe1, 0x25, 0xf1, 0x17, 0x18, 0x19, 0x1a, 0x26,
	0x27, 0x28, 0x29, 0x2a, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3a, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48,
	0x49, 0x4a, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x5a, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68,
	0x69, 0x6a, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79, 0x7a, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
	0x88, 0x89, 0x8a, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9a, 0xa2, 0xa3, 0xa4, 0xa5,
	0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xc2, 0xc3,
	0xc4, 0xc5, 0xc6, 0xc7, 0xc8, 0xc9, 0xca, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7, 0xd8, 0xd9, 0xda,
	0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7, 0xe8, 0xe9, 0xea, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8,
	0xf9, 0xfa
};


/*
 * Convert mjpeg image buffer (from webcam) to raw pixel buffer (RGB888)
 *
 * To make a valid .jpeg structure from a webcam MPEG frame (so that libjpeg
 * can decompress it), do the following.
 * 1. Start with (fixed) jpg header
 * 2. Append (fixed) huffman table
 * 3. Append MJPEG frame data, but skip the two byte start-of-image (SOI) marker
 *
 * @param out_buf points to malloc'ed data (if call succeeds), user is responsible for free()
 * @return 0 success
 * @return -1 error
 */
int MJPEGtoRGB888(unsigned char *mjpeg_frame, size_t mjpeg_frame_len, int width, int height, unsigned char **out_buf)
{
	int x;
	unsigned char r, g, b;
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	JSAMPARRAY pJpegBuffer; /* Output row buffer */
	int row_stride;         /* physical row width in output buffer */
	unsigned char jpeg_file_buffer[1024*1024*1];
	int offset = 0;

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);

	// To make a valid .jpeg structure from a webcam MPEG frame, do the following.

	// 1. Start with (fixed) jpg header
	memcpy(jpeg_file_buffer, jpg_header, sizeof jpg_header);

	// 2. Append (fixed) huffman table
	memcpy(jpeg_file_buffer + sizeof jpg_header, jpeg_dht_seg, sizeof jpeg_dht_seg);

	// 3. Append MJPEG frame data, but skip the two byte start-of-image (SOI) marker
	memcpy(jpeg_file_buffer + sizeof jpg_header + sizeof jpeg_dht_seg,
			mjpeg_frame + 2,
			mjpeg_frame_len - 2);

	jpeg_mem_src(&cinfo, jpeg_file_buffer, (unsigned long)mjpeg_frame_len + sizeof jpg_header + sizeof jpeg_dht_seg);

	if (JPEG_HEADER_OK != jpeg_read_header(&cinfo, 1)) {
		printf("jpeg_read_header: error\n");
		return -1;
	}

	jpeg_start_decompress(&cinfo);
	*out_buf = malloc(sizeof (unsigned char) * width * height * 3); // 3 bytes per pixel for RGB888

	if (!*out_buf){
		printf("%s: out of memory\n", __func__);
		return -1;
	}
	row_stride = width * cinfo.output_components ;
	pJpegBuffer = (*cinfo.mem->alloc_sarray)((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

	while (cinfo.output_scanline < cinfo.output_height) {
		jpeg_read_scanlines(&cinfo, pJpegBuffer, 1);
		for (x = 0; x < width; x++) {
			r = pJpegBuffer[0][cinfo.output_components * x];
			if (cinfo.output_components > 2) {
				g = pJpegBuffer[0][cinfo.output_components * x + 1];
				b = pJpegBuffer[0][cinfo.output_components * x + 2];
			} else {
				g = r;
				b = r;
			}
			*(*out_buf + offset++) = r;
			*(*out_buf + offset++) = g;
			*(*out_buf + offset++) = b;
		}
	}

	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	return 0;
}


/**
	Convert from YUV422 format to RGB888. Formulae are described on http://en.wikipedia.org/wiki/YUV

	\param width width of image
	\param height height of image
	\param src source
	\param dst destination
*/
static void YUV422toRGB888(int width, int height, unsigned char *src, unsigned char *dst)
{
	int line, column;
	unsigned char *py, *pu, *pv;
	unsigned char *tmp = dst;

	/* In this format each four bytes is two pixels. Each four bytes is two Y's, a Cb and a Cr.
	   Each Y goes to one of the pixels, and the Cb and Cr belong to both pixels. */
	py = src;
	pu = src + 1;
	pv = src + 3;

	#define CLIP(x) ( (x)>=0xFF ? 0xFF : ( (x) <= 0x00 ? 0x00 : (x) ) )

	for (line = 0; line < height; ++line) {
		for (column = 0; column < width; ++column) {
#ifdef ITU_R_FLOAT
			// ITU-R float
			*tmp++ = CLIP((double)*py + 1.402*((double)*pv-128.0));
			*tmp++ = CLIP((double)*py - 0.344*((double)*pu-128.0) - 0.714*((double)*pv-128.0));
			*tmp++ = CLIP((double)*py + 1.772*((double)*pu-128.0));
#endif

#ifdef ITU_R_INT
			// ITU-R integer
			*tmp++ = CLIP( *py + (*pv-128) + ((*pv-128) >> 2) + ((*pv-128) >> 3) + ((*pv-128) >> 5) );
			*tmp++ = CLIP( *py - (((*pu-128) >> 2) + ((*pu-128) >> 4) + ((*pu-128) >> 5)) - (((*pv-128) >> 1) + ((*pv-128) >> 3) + ((*pv-128) >> 4) + ((*pv-128) >> 5)) );  // 52 58 
			*tmp++ = CLIP( *py + (*pu-128) + ((*pu-128) >> 1) + ((*pu-128) >> 2) + ((*pu-128) >> 6) );
#endif

#ifdef NTSC
			// NTSC integer
			*tmp++ = CLIP( (298*(*py-16) + 409*(*pv-128) + 128) >> 8 );
			*tmp++ = CLIP( (298*(*py-16) - 100*(*pu-128) - 208*(*pv-128) + 128) >> 8 );
			*tmp++ = CLIP( (298*(*py-16) + 516*(*pu-128) + 128) >> 8 );
#endif
			// increase py every time
			py += 2;

			// increase pu,pv every second time
			if ((column & 1)==1) {
				pu += 4;
				pv += 4;
			}
		}
	}
}

/**
	Print error message and terminate programm with EXIT_FAILURE return code.

	\param s error message to print
*/
static void errno_exit(const char* s)
{
	fprintf(stderr, "%s error %d, %s\n", s, errno, strerror (errno));
	exit(EXIT_FAILURE);
}

/**
	Do ioctl and retry if error was EINTR ("A signal was caught during the ioctl() operation."). Parameters are the same as on ioctl.

	\param fd file descriptor
	\param request request
	\param argp argument
	\returns result from ioctl
*/
static int xioctl(int fd, int request, void* argp)
{
	int r;

	do r = ioctl(fd, request, argp);
	while (-1 == r && EINTR == errno);

	return r;
}

/**
	Write image to jpeg file.

	\param img image to write
	\param width image width
	\param height image height
	\param filename filename to write
*/
static void jpegWrite(unsigned char *img, unsigned int width, unsigned int height,
		const char *filename)
{
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	
	JSAMPROW row_pointer[1];
	FILE *outfile = fopen(filename, "wb" );

	// try to open file for saving
	if (!outfile) {
		errno_exit("jpeg");
	}

	// create jpeg data
	cinfo.err = jpeg_std_error( &jerr );
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, outfile);

	// set image parameters
	cinfo.image_width = width;	
	cinfo.image_height = height;
	cinfo.input_components = 3;
	cinfo.in_color_space = JCS_RGB;

	// set jpeg compression parameters to default
	jpeg_set_defaults(&cinfo);
	// and then adjust quality setting
	jpeg_set_quality(&cinfo, jpegQuality, TRUE);

	// start compress 
	jpeg_start_compress(&cinfo, TRUE);

	// feed data
	while (cinfo.next_scanline < cinfo.image_height) {
		row_pointer[0] = &img[cinfo.next_scanline * cinfo.image_width *  cinfo.input_components];
		jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}

	// finish compression
	jpeg_finish_compress(&cinfo);

	// destroy jpeg data
	jpeg_destroy_compress(&cinfo);

	// close output file
	fclose(outfile);
}


int clip256(int val)
{
	if (val < 0) {
		return 0;
	} else if (val > 256) {
		return 256;
	} else {
		return val;
	}
}


int fetch_pixel(struct surface *surf, unsigned int x, unsigned int y,
		unsigned char *r, unsigned char *g, unsigned char *b)
{
	unsigned int offset;

	if (surf->bytes_per_pixel != 3) {
		return -1;
	}

	offset = surf->width * surf->bytes_per_pixel * y + surf->bytes_per_pixel * x;
	*r = surf->buf[offset + 0];
	*g = surf->buf[offset + 1];
	*b = surf->buf[offset + 2];

	return 0;
}


int get_luminance(struct surface *surf, struct point *p0, struct point *p1,
		int *max, int *min, int *avg)
{
	unsigned int i, j;
	unsigned int luminance;
	unsigned long total_luminance = 0;
	unsigned char r, g, b;

	*max = 0;
	*min = 255;

	if (p0->x > surf->width || p0->y > surf->height
			|| p1->x > surf->width || p1->y > surf->height) {
		return -1;
	}

	for (i = p0->y; i < p1->y; i++) {
		for (j = p0->x; j < p1->x; j++) {
			fetch_pixel(surf, j, i, &r, &g, &b);
			luminance = (r + g + b) / 3;
			total_luminance += luminance;

			if (luminance > *max) {
				*max = luminance;
			}

			if (luminance < *min) {
				*min = luminance;
			}
		}
	}

	unsigned int pixels = (p1->x - p0->x) * (p1->y - p0->y);
	*avg = total_luminance / pixels;
	return 0;
}


int draw_pixel(struct surface *surf, unsigned int x, unsigned int y,
		unsigned char r, unsigned char g, unsigned char b)
{
	unsigned int offset;

	if (surf->bytes_per_pixel != 3) {
		return -1;
	}

	offset = surf->width * surf->bytes_per_pixel * y + surf->bytes_per_pixel * x;
	surf->buf[offset + 0] = r;
	surf->buf[offset + 1] = g;
	surf->buf[offset + 2] = b;
	return 0;
}

// http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
int draw_line(struct surface *surf,
		unsigned int x0, unsigned int y0,
		unsigned int x1, unsigned int y1,
		unsigned char r, unsigned char g, unsigned char b)
{
	int dx, dy, sx, sy, err, e2;

	if (surf->bytes_per_pixel != 3) {
		return -1;
	}

	dx = abs(x1 - x0);
	dy = abs(y1 - y0);
	sx = x0 < x1 ? 1 : -1;
	sy = y0 < y1 ? 1 : -1;
	err = dx - dy;

	while (1) {
		draw_pixel(surf, x0, y0, r, g, b);
		if (x0 == x1 && y0 == y1) {
			break;
		}
		e2 = 2 * err;
		if (e2 > -dy) {
			err = err - dy;
			x0 = x0 + sx;
		}
		if (e2 < dx) {
			err = err + dx;
			y0 = y0 + sy;
		}
	}

	return 0;
}


// Scan a line in surface *surf and count the number of contiguous sections
// that has a luminance below the given threshold (+ hysteresis).
int contiguous_dark_sections(struct surface *surf, unsigned int line,
		unsigned char luminance_threshold, unsigned int hysteresis,
		unsigned int *num_sections)
{
	enum state {NO_MATCH, MATCH} state = NO_MATCH;
	unsigned char r, g, b;
	unsigned char luminance;
	int i;
	int bar_width = 0;

	*num_sections = 0;

	if (line >= surf->height) {
		return -1;
	}

	for (i = 0; i < surf->width; i++) {
		fetch_pixel(surf, i, line, &r, &g, &b);
		luminance = (r + g + b) / 3;

		if (luminance < luminance_threshold - hysteresis) {
			if (state != MATCH) {
				state = MATCH;
				bar_width = 1;
				(*num_sections)++;
				draw_pixel(surf, i, line, 255, 0, 0);
			} else {
				bar_width++;
			}
		} else if (luminance >= luminance_threshold + hysteresis) {
			if (state == MATCH) {
				state = NO_MATCH;
				bar_width = 0;
				draw_pixel(surf, i, line, 0, 255, 0);
			}
		} else {
			if (state == MATCH) {
				bar_width++;
			}
		}
	}

	return 0;
}


// locate global maximum of the values in buf, move to both sides until
// value<(max-threshold), where threshold=(max-min)*ratio, and return the
// front and back edges of the peak area.
int find_peak_area(int *buf, unsigned int len, int *front_edge, int *back_edge)
{
	int min, max, max_pos;
	int front_edge_index, back_edge_index;
	int i;
	int threshold;

	min = max = buf[0];
	max_pos = 0;

	for (i = 0; i < len; i++) {
		if (buf[i] < min) {
			min = buf[i];
		}
		if (buf[i] > max) {
			max = buf[i];
			max_pos = i;
		}
	}

	threshold = (max - min) * 0.4;

	for (int i=max_pos; i<len; i++) {
		if (buf[i] < (max - threshold)) {
			break;
		} else {
			back_edge_index = i;
		}
	}

	for (int i=max_pos; i>=0; i--) {
		if (buf[i] < (max - threshold)) {
			break;
		} else {
			front_edge_index = i;
		}
	}

	if (front_edge && back_edge) {
		*front_edge = front_edge_index;
		*back_edge = back_edge_index;
	}

	return 0;
}

// from https://github.com/xuphys/peakdetect
int detect_peak(
        const int *data, /* the data */
        int data_count, /* row count of data */
        int* emi_peaks, /* emission peaks will be put here */
        int* num_emi_peaks, /* number of emission peaks found */
        int max_emi_peaks, /* maximum number of emission peaks */
        int* absop_peaks, /* absorption peaks will be put here */
        int* num_absop_peaks, /* number of absorption peaks found */
        int max_absop_peaks, /* maximum number of absorption peaks
*/
        double delta, /* delta used for distinguishing peaks */
        int emi_first /* should we search emission peak first of
absorption peak first? */
        )
{
    int i;
    double mx;
    double mn;
    int mx_pos = 0;
    int mn_pos = 0;
    int is_detecting_emi = emi_first;


    mx = data[0];
    mn = data[0];

    *num_emi_peaks = 0;
    *num_absop_peaks = 0;

    for(i = 1; i < data_count; ++i)
    {
        if(data[i] > mx)
        {
            mx_pos = i;
            mx = data[i];
        }
        if(data[i] < mn)
        {
            mn_pos = i;
            mn = data[i];
        }

        if(is_detecting_emi &&
                data[i] < mx - delta)
        {
            if(*num_emi_peaks >= max_emi_peaks) /* not enough spaces */
                return 1;

            emi_peaks[*num_emi_peaks] = mx_pos;
            ++ (*num_emi_peaks);

            is_detecting_emi = 0;

            i = mx_pos - 1;

            mn = data[mx_pos];
            mn_pos = mx_pos;
        }
        else if((!is_detecting_emi) &&
                data[i] > mn + delta)
        {
            if(*num_absop_peaks >= max_absop_peaks)
                return 2;

            absop_peaks[*num_absop_peaks] = mn_pos;
            ++ (*num_absop_peaks);

            is_detecting_emi = 1;

            i = mn_pos - 1;

            mx = data[mn_pos];
            mx_pos = data[mn_pos];
        }
    }

    return 0;
}


int do_find_barcode3(struct surface *surf)
{
	int y, x;
	int data[surf->width];
        int emi_peaks[surf->width];        /* emission peaks will be put here */
        int num_emi_peaks;                 /* number of emission peaks found */
        int max_emi_peaks = surf->width;   /* maximum number of emission peaks */
        int absop_peaks[surf->width];      /* absorption peaks will be put here */
        int num_absop_peaks;               /* number of absorption peaks found */
        int max_absop_peaks = surf->width; /* maximum number of absorption peaks */
        double delta = 40;                 /* delta used for distinguishing peaks */
        int emi_first = 1;                 /* should we search emission peak first of absorption peak first? */
	unsigned char r, g, b, luminance;

	int bars[surf->height];
	int start_line, end_line;
	int barcode_height;
	int max_bars = 0;
	int center;

	for (y=0; y<surf->height; y+=1) {
		// build luminance array for this scan line
		for (x=0; x<surf->width; x+=1) {
			fetch_pixel(surf, x, y, &r, &g, &b);
			luminance = (r + g + b) / 3;
			data[x] = luminance;
		}

		// locate peaks
		detect_peak(data,
				sizeof data / sizeof data[0],
				emi_peaks, /* emission peaks will be put here */
				&num_emi_peaks, /* number of emission peaks found */
				max_emi_peaks, /* maximum number of emission peaks */
				absop_peaks, /* absorption peaks will be put here */
				&num_absop_peaks, /* number of absorption peaks found */
				max_absop_peaks, /* maximum number of absorption peaks */
				delta, /* delta used for distinguishing peaks */
				emi_first /* should we search emission peak first of absorption peak first? */
			   );

		//printf("line peaks: %d %d\n", y, num_emi_peaks);
		bars[y] = num_emi_peaks;
		if (max_bars < num_emi_peaks) {
			max_bars = num_emi_peaks;
		}

		// draw
		for (x=0; x<num_emi_peaks; x++) {
			draw_pixel(surf, emi_peaks[x], y, 255, 0, 0);
		}
		for (x=0; x<num_absop_peaks; x++) {
			draw_pixel(surf, absop_peaks[x], y, 0, 255, 0);
		}
	}

	find_peak_area(bars, sizeof bars / sizeof (bars[0]), &start_line, &end_line);
	barcode_height = end_line - start_line;
	if (barcode_height <= 0 || max_bars < 40) {
		barcode_height = 0;
		start_line = -1;
		end_line = -1;
	}

	center = start_line + barcode_height / 2;
	//printf("max_bars %d\n", max_bars);
	printf("%d pixel high barcode found centered around line %d (start %d, end %d)\n",
				barcode_height, center, start_line, end_line);
	draw_line(surf, 1, start_line, 5, start_line, 255, 255, 255);
	draw_line(surf, 1, end_line, 5, end_line, 255, 255, 255);
	return 0;
}


int do_find_barcode2(struct surface *surf)
{
	struct point p0, p1;
	int max, min, avg;
	unsigned int hysteresis = 5;
	int barcode_height;
	int start_line, end_line;
	int center;
	unsigned int dark_bars;
	int max_bars = 0;
	int bars[surf->height];
	int i;

	memset(bars, 0, sizeof bars);

	p0.x = 0;
	p1.x = width;

	for (i=0; i<surf->height-1; i+=1) {
		p0.y = i;
		p1.y = i+1;
		get_luminance(surf, &p0, &p1, &max, &min, &avg);
		contiguous_dark_sections(surf, i, avg, hysteresis, &dark_bars);
		bars[i] = dark_bars;
		if (dark_bars > max_bars) {
			max_bars = dark_bars;
		}
	}

	//for (i=0; i < sizeof bars / sizeof (bars[0]); i++) {
		//printf("bars line %d %d\n", i, bars[i]);
	//}

	//printf("max_bars %d\n", max_bars);
	find_peak_area(bars, sizeof bars / sizeof (bars[0]), &start_line, &end_line);
	barcode_height = end_line - start_line;
	if (barcode_height <= 0 || max_bars < 30) {
		barcode_height = 0;
		start_line = -1;
		end_line = -1;
	}

	center = start_line + barcode_height / 2;
	printf("%d pixel high barcode found centered around line %d (start %d, end %d)\n",
				barcode_height, center, start_line, end_line);
	draw_line(surf, 1, start_line, 5, start_line, 255, 255, 255);
	draw_line(surf, 1, end_line, 5, end_line, 255, 255, 255);

	return 0;
}


int do_find_barcode1(struct surface *surf)
{
	signed int line_where_barcode_begins, line_where_barcode_ends;
	unsigned int barcount[surf->height];
	unsigned int bar_begin, bar_end;
	struct point p0, p1;
	int max, min, avg;
	unsigned int luminance_threshold;
	int i, j;
	unsigned char r, g, b;
	unsigned int luminance, bar_height, height_tmp;

	line_where_barcode_begins = height-1;
	line_where_barcode_ends = 0;

	p0.x = 0;
	p0.y = 0;
	p1.x = width;
	p1.y = height;

	get_luminance(surf, &p0, &p1, &max, &min, &avg);
	//printf("luminance (max/min/avg): %d %d %d\n", max, min, avg);
	luminance_threshold = 3 * (avg - (avg - min) * 0.5);

	for (i=0; i<surf->height; i+=20) {
		memset(barcount, 0, sizeof barcount);
		for (j=0; j<surf->width; j++) {
			fetch_pixel(surf, j, i, &r, &g, &b);
			luminance = r + g + b;
			//printf("luminance: %d at (%d, %d)\n", luminance, j, i);
			if (luminance < luminance_threshold) {
				// look for contiguous vertical bar
				bar_height = 0;
				bar_begin = bar_end = i;
				/*printf("luminance below luminance_threshold at (%d, %d)\n", j, i);*/

				// scan downwards
				for (height_tmp=i; height_tmp<height; height_tmp++) {
					bar_height++;
					bar_end = height_tmp;
					fetch_pixel(surf, j, height_tmp, &r, &g, &b);
					luminance = r + g + b;
					if (luminance > luminance_threshold) {
						break;
					} else {
						draw_pixel(surf, j, height_tmp, 0, 0, 255);
					}
				}

				// scan upwards
				for (height_tmp=i; height_tmp>0; height_tmp--) {
					bar_height++;
					bar_begin = height_tmp;
					fetch_pixel(surf, j, height_tmp-1, &r, &g, &b);
					luminance = r + g + b;
					if (luminance > luminance_threshold) {
						break;
					} else {
						draw_pixel(surf, j, height_tmp, 0, 255, 0);
					}
				}

				if (bar_height > 50) {
					//draw_pixel(surf, j, height_tmp, 255, 0, 0);
					for (int i=bar_begin; i<=bar_end; i++) {
						barcount[i]++;
					}
					draw_line(surf, j, bar_begin, j, bar_end, 255, 0, 0);
					//printf("bar_height: %d at (%d, %d)\n", bar_height, j, i);
				}
			}
		}

		if (barcount[i] > 20) {
			for (int k=i; k>0; k--) {
				if (barcount[k] > 20) {
					if (line_where_barcode_begins > k) {
						line_where_barcode_begins = k;
					}
				}
			}

			for (int k=i; k<height; k++) {
				if (barcount[k] > 20) {
					if (line_where_barcode_ends < k) {
						line_where_barcode_ends = k;
					}
				}
			}
		}
	}

	int barcode_height = line_where_barcode_ends - line_where_barcode_begins;
	if (barcode_height <= 0) {
		barcode_height = 0;
		line_where_barcode_begins = -1;
		line_where_barcode_ends = -1;
	} else {
		draw_line(surf, 0, line_where_barcode_begins, 0, line_where_barcode_ends, 255, 255, 0);
	}

	int center = line_where_barcode_begins + barcode_height / 2;
	printf("%d pixel high barcode found centered around line %d (start %d, end %d)\n",
			barcode_height, center, line_where_barcode_begins, line_where_barcode_ends);
	if (barcode_height > 90) {
		//printf("WARN: very high barcode: %d\n", barcode_height);
	}

	return 0;
}

/**
	process image read
*/
static int imageProcess(const void *p, size_t buflen)
{
	static int skip_frames = 1; // the first frame(s) is usually bad, so skip it
	unsigned char* src = (unsigned char*)p;
	unsigned char *dst, *dst_copy;
	unsigned int i, line, offset, offset_in_line;
	char tmp;
	struct surface surf;
	int ret;

	if (skip_frames-- > 0) {
		return 0;
	}

	if (use_mjpeg) {
		ret = MJPEGtoRGB888(src, buflen, width, height, &dst);
		if (ret != 0) {
			return -1;
		}
	} else {
		dst = malloc(width*height*3*sizeof(char));
		if (dst == NULL) {
			perror("malloc");
			return -1;
		}

		YUV422toRGB888(width, height, src, dst);
	}

	// swap pixels
	if (swap_pixels) {
		for (i=0; i<width*height*3; i+=6) {
			tmp = dst[i+0];
			dst[i+0] = dst[i+3+0];
			dst[i+3+0] = tmp;

			tmp = dst[i+1];
			dst[i+1] = dst[i+3+1];
			dst[i+3+1] = tmp;

			tmp = dst[i+2];
			dst[i+2] = dst[i+3+2];
			dst[i+3+2] = tmp;
		}
	}

	// make a copy we can mess up
	dst_copy = malloc(width*height*3*sizeof(char));
	if (dst_copy == NULL) {
		perror("malloc");
		goto free_dst;
	}
	memcpy(dst_copy, dst, width*height*3*sizeof(char));

	surf.bytes_per_pixel = 3;
	surf.width = width;
	surf.height = height;
	surf.buf = dst_copy;


	if (find_barcode) {
		switch (barcode_algorithm) {
			case 1: do_find_barcode1(&surf); break;
			case 2: do_find_barcode2(&surf); break;
			case 3: do_find_barcode3(&surf); break;
			default:
				printf("No such barcode algorithm: %d (valid [1..3])\n",
						barcode_algorithm);
				break;
		}
	}

	if (jpegFilename) {
		//printf("  jpegWrite()\n");
		// write jpeg
		jpegWrite(dst_copy, width, height, jpegFilename);
	}

	if (!do_fft) {
		free(dst);
		free(dst_copy);
		return 1;
	}

	offset_in_line = (3 * sizeof (char)) * 40;

	//printf("calculating FFT\n");
	fflush(NULL);
	for (line=0; line<height; line++) {
		offset = (width * 3 * sizeof (char)) * line + offset_in_line;
		for (i=0; i<NUM_FFT; i++) {
			Real[i] = dst[i*3 + offset + 0] + dst[i*3 + offset + 1] + dst[i*3 + offset + 2];
			Real[i] *= 42;
			//Real[i] = 128 * src[i*2 + offset];   multiply by N to get better FFT resolution
		}

		//// print the Real array
		////printf("%03d:", line);
		//for (i=0; i<NUM_FFT; i++) {
			//printf("%d\t", Real[i]);
		//}
		//printf("\n");

		CLEAR(Imag);
		// fft_brin is hardcoded for NUM_FFT samples defined in project
		// and it must be 2^n (We'll run a fft on a longer window
		// simply by not using the leftover piece)
		fft_brin(Real);

		// Convert result to magnitudes (in-place)
		fft_mag(NUM_FFT/2, Real);

		// print the FFT array
		//printf("%03d:", line);
		for (i=0; i<NUM_FFT/2; i++) {
			printf("%d\t", Real[i]);
		}
		printf("\n");
	}

	free(dst_copy);
free_dst:
	free(dst);
	//printf("imageProcess() end\n");
	return 1;
}

/**
	read single frame
*/
static int frameRead(void)
{
	int ret = 0;
	struct v4l2_buffer buf;
#ifdef IO_USERPTR
	unsigned int i;
#endif

	//printf("frameRead()\n");
	fflush(NULL);

	switch (io) {
#ifdef IO_READ
		case IO_METHOD_READ:
			if (-1 == read (fd, buffers[0].start, buffers[0].length)) {
				switch (errno) {
					case EAGAIN:
						return 0;

					case EIO:
						// Could ignore EIO, see spec.
						// fall through

					default:
						errno_exit("read");
				}
			}

			ret = imageProcess(buffers[0].start, buffers[0].length);
			break;
#endif

#ifdef IO_MMAP
		case IO_METHOD_MMAP:
			CLEAR(buf);

			buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
			buf.memory = V4L2_MEMORY_MMAP;

			if (-1 == xioctl(fd, VIDIOC_DQBUF, &buf)) {
				switch (errno) {
					case EAGAIN:
						return 0;

					case EIO:
						// Could ignore EIO, see spec
						// fall through

					default:
						errno_exit("VIDIOC_DQBUF");
				}
			}

			assert(buf.index < n_buffers);

			ret = imageProcess(buffers[buf.index].start, buffers[buf.index].length);

			if (-1 == xioctl(fd, VIDIOC_QBUF, &buf))
				errno_exit("VIDIOC_QBUF");

			break;
#endif

#ifdef IO_USERPTR
			case IO_METHOD_USERPTR:
				CLEAR (buf);

				buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
				buf.memory = V4L2_MEMORY_USERPTR;

				if (-1 == xioctl(fd, VIDIOC_DQBUF, &buf)) {
					switch (errno) {
						case EAGAIN:
							return 0;

						case EIO:
							// Could ignore EIO, see spec.
							// fall through

						default:
							errno_exit("VIDIOC_DQBUF");
					}
				}

				for (i = 0; i < n_buffers; ++i)
					if (buf.m.userptr == (unsigned long) buffers[i].start && buf.length == buffers[i].length)
						break;

				assert (i < n_buffers);

				ret = imageProcess((void *)buf.m.userptr, buffers[i].length);

				if (-1 == xioctl(fd, VIDIOC_QBUF, &buf))
					errno_exit("VIDIOC_QBUF");
				break;
#endif
	}

	return ret;
}

/** 
	mainloop: read frames and process them
*/
static void mainLoop(void)
{
	unsigned int count;
	unsigned int numberOfTimeouts;

	numberOfTimeouts = 0;
	count = 1;

	while (count-- > 0) {
		for (;;) {
			fd_set fds;
			struct timeval tv;
			int r;

			FD_ZERO(&fds);
			FD_SET(fd, &fds);

			/* Timeout. */
			tv.tv_sec = 1;
			tv.tv_usec = 0;

			r = select(fd + 1, &fds, NULL, NULL, &tv);

			if (-1 == r) {
				if (EINTR == errno)
					continue;

				errno_exit("select");
			}

			if (0 == r) {
				if (numberOfTimeouts <= 0) {
					count++;
				} else {
					fprintf(stderr, "select timeout\n");
					exit(EXIT_FAILURE);
				}
			}

			if (frameRead())
				break;

			/* EAGAIN - continue select loop. */
		}
	}
}

/**
	stop capturing
*/
static void captureStop(void)
{
	enum v4l2_buf_type type;

	switch (io) {
#ifdef IO_READ
		case IO_METHOD_READ:
			/* Nothing to do. */
			break;
#endif

#ifdef IO_MMAP
		case IO_METHOD_MMAP:
#endif
#ifdef IO_USERPTR
		case IO_METHOD_USERPTR:
#endif
#if defined(IO_MMAP) || defined(IO_USERPTR)
			type = V4L2_BUF_TYPE_VIDEO_CAPTURE;

			if (-1 == xioctl(fd, VIDIOC_STREAMOFF, &type))
			errno_exit("VIDIOC_STREAMOFF");

			break;
#endif 
	}
}

/**
  start capturing
*/
static void captureStart(void)
{
	unsigned int i;
	enum v4l2_buf_type type;

	switch (io) {
#ifdef IO_READ    
		case IO_METHOD_READ:
			/* Nothing to do. */
			break;
#endif

#ifdef IO_MMAP
		case IO_METHOD_MMAP:
			for (i = 0; i < n_buffers; ++i) {
				struct v4l2_buffer buf;

				CLEAR(buf);

				buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
				buf.memory = V4L2_MEMORY_MMAP;
				buf.index = i;

				if (-1 == xioctl(fd, VIDIOC_QBUF, &buf))
					errno_exit("VIDIOC_QBUF");
				}

			type = V4L2_BUF_TYPE_VIDEO_CAPTURE;

			if (-1 == xioctl(fd, VIDIOC_STREAMON, &type))
				errno_exit("VIDIOC_STREAMON");

			break;
#endif

#ifdef IO_USERPTR
		case IO_METHOD_USERPTR:
			for (i = 0; i < n_buffers; ++i) {
				struct v4l2_buffer buf;

			CLEAR (buf);

			buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
			buf.memory = V4L2_MEMORY_USERPTR;
			buf.index = i;
			buf.m.userptr = (unsigned long) buffers[i].start;
			buf.length = buffers[i].length;

			if (-1 == xioctl(fd, VIDIOC_QBUF, &buf))
				errno_exit("VIDIOC_QBUF");
			}

			type = V4L2_BUF_TYPE_VIDEO_CAPTURE;

			if (-1 == xioctl(fd, VIDIOC_STREAMON, &type))
				errno_exit("VIDIOC_STREAMON");

			break;
#endif
	}
}

static void deviceUninit(void)
{
	unsigned int i;

	switch (io) {
#ifdef IO_READ
		case IO_METHOD_READ:
			free(buffers[0].start);
			break;
#endif

#ifdef IO_MMAP
		case IO_METHOD_MMAP:
			for (i = 0; i < n_buffers; ++i)
				if (-1 == munmap (buffers[i].start, buffers[i].length))
					errno_exit("munmap");
			break;
#endif

#ifdef IO_USERPTR
		case IO_METHOD_USERPTR:
			for (i = 0; i < n_buffers; ++i)
				free (buffers[i].start);
			break;
#endif
	}

	free(buffers);
}

#ifdef IO_READ
static void readInit(unsigned int buffer_size)
{
	buffers = calloc(1, sizeof(*buffers));

	if (!buffers) {
		fprintf(stderr, "Out of memory\n");
		exit(EXIT_FAILURE);
	}

	buffers[0].length = buffer_size;
	buffers[0].start = malloc(buffer_size);

	if (!buffers[0].start) {
		fprintf (stderr, "Out of memory\n");
		exit(EXIT_FAILURE);
	}
}
#endif

#ifdef IO_MMAP
static void mmapInit(void)
{
	struct v4l2_requestbuffers req;

	CLEAR(req);

	req.count = 4;
	req.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
	req.memory = V4L2_MEMORY_MMAP;

	if (-1 == xioctl(fd, VIDIOC_REQBUFS, &req)) {
		if (EINVAL == errno) {
			fprintf(stderr, "%s does not support memory mapping\n", deviceName);
			exit(EXIT_FAILURE);
		} else {
			errno_exit("VIDIOC_REQBUFS");
		}
	}

	if (req.count < 2) {
		fprintf(stderr, "Insufficient buffer memory on %s\n", deviceName);
		exit(EXIT_FAILURE);
	}

	buffers = calloc(req.count, sizeof(*buffers));

	if (!buffers) {
		fprintf(stderr, "Out of memory\n");
		exit(EXIT_FAILURE);
	}

	for (n_buffers = 0; n_buffers < req.count; ++n_buffers) {
		struct v4l2_buffer buf;

		CLEAR(buf);

		buf.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
		buf.memory = V4L2_MEMORY_MMAP;
		buf.index = n_buffers;

		if (-1 == xioctl(fd, VIDIOC_QUERYBUF, &buf))
			errno_exit("VIDIOC_QUERYBUF");

		buffers[n_buffers].length = buf.length;
		buffers[n_buffers].start = mmap (NULL /* start anywhere */, buf.length, PROT_READ | PROT_WRITE /* required */, MAP_SHARED /* recommended */, fd, buf.m.offset);

		if (MAP_FAILED == buffers[n_buffers].start)
			errno_exit("mmap");
	}
}
#endif

#ifdef IO_USERPTR
static void userptrInit(unsigned int buffer_size)
{
	struct v4l2_requestbuffers req;
	unsigned int page_size;

	page_size = getpagesize();
	buffer_size = (buffer_size + page_size - 1) & ~(page_size - 1);

	CLEAR(req);

	req.count = 4;
	req.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
	req.memory = V4L2_MEMORY_USERPTR;

	if (-1 == xioctl(fd, VIDIOC_REQBUFS, &req)) {
		if (EINVAL == errno) {
			fprintf(stderr, "%s does not support user pointer i/o\n", deviceName);
			exit(EXIT_FAILURE);
		} else {
			errno_exit("VIDIOC_REQBUFS");
		}
	}

	buffers = calloc(4, sizeof (*buffers));

	if (!buffers) {
		fprintf(stderr, "Out of memory\n");
		exit(EXIT_FAILURE);
	}

	for (n_buffers = 0; n_buffers < 4; ++n_buffers) {
		buffers[n_buffers].length = buffer_size;
		buffers[n_buffers].start = memalign (/* boundary */ page_size, buffer_size);

		if (!buffers[n_buffers].start) {
			fprintf(stderr, "Out of memory\n");
			exit(EXIT_FAILURE);
		}
	}
}
#endif

/**
	initialize device
*/
static void deviceInit(void)
{
	struct v4l2_capability cap;
	struct v4l2_cropcap cropcap;
	struct v4l2_crop crop;
	struct v4l2_format fmt;
	unsigned int min;

	if (-1 == xioctl(fd, VIDIOC_QUERYCAP, &cap)) {
		if (EINVAL == errno) {
			fprintf(stderr, "%s is no V4L2 device\n",deviceName);
			exit(EXIT_FAILURE);
		} else {
			errno_exit("VIDIOC_QUERYCAP");
		}
	}

	if (!(cap.capabilities & V4L2_CAP_VIDEO_CAPTURE)) {
		fprintf(stderr, "%s is no video capture device\n",deviceName);
		exit(EXIT_FAILURE);
	}

	switch (io) {
#ifdef IO_READ
		case IO_METHOD_READ:
			if (!(cap.capabilities & V4L2_CAP_READWRITE)) {
				fprintf(stderr, "%s does not support read i/o\n",deviceName);
				exit(EXIT_FAILURE);
			}
			break;
#endif

#ifdef IO_MMAP
		case IO_METHOD_MMAP:
#endif
#ifdef IO_USERPTR
		case IO_METHOD_USERPTR:
#endif
#if defined(IO_MMAP) || defined(IO_USERPTR)
      			if (!(cap.capabilities & V4L2_CAP_STREAMING)) {
				fprintf(stderr, "%s does not support streaming i/o\n",deviceName);
				exit(EXIT_FAILURE);
			}
			break;
#endif
	}


	/* Select video input, video standard and tune here. */
	CLEAR(cropcap);

	cropcap.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;

	if (0 == xioctl(fd, VIDIOC_CROPCAP, &cropcap)) {
		crop.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
		crop.c = cropcap.defrect; /* reset to default */

		if (-1 == xioctl(fd, VIDIOC_S_CROP, &crop)) {
			switch (errno) {
				case EINVAL:
					/* Cropping not supported. */
					break;
				default:
					/* Errors ignored. */
					break;
			}
		}
	} else {        
		/* Errors ignored. */
	}

	CLEAR(fmt);

	// v4l2_format
	fmt.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
	fmt.fmt.pix.width = width; 
	fmt.fmt.pix.height = height;
	fmt.fmt.pix.field = V4L2_FIELD_SEQ_TB;

	if (use_mjpeg) {
		fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_MJPEG;
		if (-1 == xioctl(fd, VIDIOC_S_FMT, &fmt)) {
			printf("info: MJPEG not supported, falling back to YUV\n");
			fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_YUYV;
			use_mjpeg = 0;
		}
	} else {
		fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_YUYV;
		use_mjpeg = 0;
	}

	if (-1 == xioctl(fd, VIDIOC_S_FMT, &fmt))
		errno_exit("VIDIOC_S_FMT");

	/* Note VIDIOC_S_FMT may change width and height. */
	if (width != fmt.fmt.pix.width) {
		width = fmt.fmt.pix.width;
		fprintf(stderr,"Image width set to %i by device %s.\n",width,deviceName);
	}

	if (height != fmt.fmt.pix.height) {
		height = fmt.fmt.pix.height;
		fprintf(stderr,"Image height set to %i by device %s.\n",height,deviceName);
	}

	/* Buggy driver paranoia. */
	min = fmt.fmt.pix.width * 2;
	if (fmt.fmt.pix.bytesperline < min)
		fmt.fmt.pix.bytesperline = min;
	min = fmt.fmt.pix.bytesperline * fmt.fmt.pix.height;
	if (fmt.fmt.pix.sizeimage < min)
		fmt.fmt.pix.sizeimage = min;

	switch (io) {
#ifdef IO_READ
		case IO_METHOD_READ:
			readInit(fmt.fmt.pix.sizeimage);
			break;
#endif

#ifdef IO_MMAP
		case IO_METHOD_MMAP:
			mmapInit();
			break;
#endif

#ifdef IO_USERPTR
		case IO_METHOD_USERPTR:
			userptrInit(fmt.fmt.pix.sizeimage);
			break;
#endif
	}
}

/**
	close device
*/
static void deviceClose(void)
{
	if (-1 == close (fd))
		errno_exit("close");

	fd = -1;
}

/**
	open device
*/
static void deviceOpen(void)
{
	struct stat st;

	// stat file
	if (-1 == stat(deviceName, &st)) {
		fprintf(stderr, "Cannot identify '%s': %d, %s\n", deviceName, errno, strerror (errno));
		exit(EXIT_FAILURE);
	}

	// check if its device
	if (!S_ISCHR (st.st_mode)) {
		fprintf(stderr, "%s is no device\n", deviceName);
		exit(EXIT_FAILURE);
	}

	// open device
	fd = open(deviceName, O_RDWR /* required */ | O_NONBLOCK, 0);

	// check if opening was successfull
	if (-1 == fd) {
		fprintf(stderr, "Cannot open '%s': %d, %s\n", deviceName, errno, strerror (errno));
		exit(EXIT_FAILURE);
	}
}

/**
	print usage information
*/
static void usage(FILE* fp, int argc, char** argv)
{
	fprintf (fp,
		"Usage: %s [options]\n\n"
		"Options:\n"
		"-d | --device name   Video device name [/dev/video0]\n"
		"-h | --help          Print this message\n"
		"-o | --output        JPEG output filename\n"
		"-s | --swap-pixels   swap pixels (needed on Atmel ISI)\n"
		"-f | --fft           FFT data to stdout\n"
		"-b | --barcode N     use barcode detection algorithm N\n"
		"-q | --quality       JPEG quality (0-100) [%d]\n"
		"-m | --mmap          Use memory mapped buffers\n"
		"-r | --read          Use read() calls\n"
		"-u | --userptr       Use application allocated buffers\n"
		"-W | --width         width [%d]\n"
		"-H | --height        height [%d]\n"
		"",
		argv[0],
		jpegQuality,
		width,
		height);
	}

static const char short_options [] = "d:ho:sfb:q:mruW:H:";

static const struct option
long_options [] = {
	{ "device",     required_argument,      NULL,           'd' },
	{ "help",       no_argument,            NULL,           'h' },
	{ "output",     required_argument,      NULL,           'o' },
	{ "fft",        no_argument,            NULL,           'f' },
	{ "barcode",    required_argument,      NULL,           'b' },
	{ "quality",    required_argument,      NULL,           'q' },
	{ "mmap",       no_argument,            NULL,           'm' },
	{ "read",       no_argument,            NULL,           'r' },
	{ "userptr",    no_argument,            NULL,           'u' },
	{ "width",      required_argument,      NULL,           'W' },
	{ "height",     required_argument,      NULL,           'H' },
	{ 0, 0, 0, 0 }
};

int main(int argc, char **argv)
{
	char *endptr;

	for (;;) {
		int index, c = 0;

		c = getopt_long(argc, argv, short_options, long_options, &index);

		if (-1 == c)
			break;

		switch (c) {
			case 0: /* getopt_long() flag */
				break;

			case 'd':
				deviceName = optarg;
				break;

			case 'h':
				// print help
				usage(stdout, argc, argv);
				exit(EXIT_SUCCESS);

			case 'o':
				// set jpeg filename
				jpegFilename = optarg;
				break;

			case 's':
				swap_pixels = 1;
				break;

			case 'f':
				do_fft = 1;
				break;

			case 'b':
				find_barcode = 1;
				barcode_algorithm = strtol(optarg, NULL, 10);
				if (barcode_algorithm < 1 || barcode_algorithm > 3) {
					fprintf(stderr, "algorithm must be either 1, 2 or 3\n");
					return 1;
				}
				break;

			case 'q':
				// set jpeg quality
				jpegQuality = strtol(optarg, &endptr, 10);
				if ((jpegQuality == 0 && endptr == optarg) || jpegQuality > 100) {
					fprintf(stderr, "quality must be an integer between 0 and 100\n");
					return 1;
				}
				break;

			case 'm':
#ifdef IO_MMAP
				io = IO_METHOD_MMAP;
#else
				fprintf(stderr, "You didn't compile for mmap support.\n");
				exit(EXIT_FAILURE);         
#endif
				break;

			case 'r':
#ifdef IO_READ
				io = IO_METHOD_READ;
#else
				fprintf(stderr, "You didn't compile for read support.\n");
				exit(EXIT_FAILURE);         
#endif
				break;

			case 'u':
#ifdef IO_USERPTR
				io = IO_METHOD_USERPTR;
#else
				fprintf(stderr, "You didn't compile for userptr support.\n");
				exit(EXIT_FAILURE);         
#endif
				break;

			case 'W':
				// set width
				width = strtol(optarg, &endptr, 10);
				if (endptr == optarg) {
					fprintf(stderr, "width must be an integer\n");
					return 1;
				}
				break;

			case 'H':
				// set height
				height = strtol(optarg, &endptr, 10);
				if (endptr == optarg) {
					fprintf(stderr, "width must be an integer\n");
					return 1;
				}
				break;

			default:
				usage(stderr, argc, argv);
				exit(EXIT_FAILURE);
		}
	}

	signal(SIGINT, sighandler);

	// open and initialize device
	deviceOpen();
	deviceInit();

	// start capturing
	captureStart();

	// process frames
	mainLoop();

	// stop capturing
	captureStop();

	// close device
	deviceUninit();
	deviceClose();

	exit(EXIT_SUCCESS);

	return 0;
}
