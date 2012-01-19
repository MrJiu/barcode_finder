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
static char* deviceName = "/dev/video0";

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
*/
static void jpegWrite(unsigned char* img)
{
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	
	JSAMPROW row_pointer[1];
	FILE *outfile = fopen( jpegFilename, "wb" );

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


int get_luminance(struct surface *surf, int *max, int *min, int *avg)
{
	unsigned int i, j;
	unsigned int luminance;
	unsigned long total_luminance = 0;
	unsigned char r, g, b;

	*max = 0;
	*min = 255;

	for (i=0; i<surf->height; i++) {
		for (j=0; j<surf->width; j++) {
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

	*avg = total_luminance / (width * height);
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


/**
	process image read
*/
static void imageProcess(const void* p)
{
	unsigned char* src = (unsigned char*)p;
	unsigned char *dst, *dst_copy;
	unsigned int i, j, line, offset, offset_in_line;
	unsigned int luminance, bar_height, height_tmp;
	char tmp;
	unsigned int luminance_threshold;
	signed int line_where_barcode_begins, line_where_barcode_ends;
	struct surface surf;
	unsigned char r, g, b;
	unsigned int bar_begin, bar_end;
	unsigned int barcount[height];
	int max, min, avg;

	line_where_barcode_begins = height-1;
	line_where_barcode_ends = 0;

	dst = malloc(width*height*3*sizeof(char));
	if (dst == NULL) {
		perror("malloc");
		exit(1);
	}

	//printf("imageProcess() start\n");
	//fflush(NULL);
	//printf("  convert from YUV422 to RGB888\n");
	// convert from YUV422 to RGB888
	YUV422toRGB888(width,height,src,dst);

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
		exit(1);
	}
	memcpy(dst_copy, dst, width*height*3*sizeof(char));

	surf.bytes_per_pixel = 3;
	surf.width = width;
	surf.height = height;
	surf.buf = dst_copy;

	get_luminance(&surf, &max, &min, &avg);
	//printf("luminance (max/min/avg): %d %d %d\n", max, min, avg);
	luminance_threshold = avg * 3 * 0.7;

	// find vertical bar(s)
	if (find_barcode) {
		for (i=0; i<height; i+=20) {
			memset(barcount, 0, sizeof barcount);
			for (j=0; j<width; j++) {
				fetch_pixel(&surf, j, i, &r, &g, &b);
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
						fetch_pixel(&surf, j, height_tmp, &r, &g, &b);
						luminance = r + g + b;
						if (luminance > luminance_threshold) {
							break;
						} else {
							draw_pixel(&surf, j, height_tmp, 0, 0, 255);
						}
					}

					// scan upwards
					for (height_tmp=i; height_tmp>0; height_tmp--) {
						bar_height++;
						bar_begin = height_tmp;
						fetch_pixel(&surf, j, height_tmp-1, &r, &g, &b);
						luminance = r + g + b;
						if (luminance > luminance_threshold) {
							break;
						} else {
							draw_pixel(&surf, j, height_tmp, 0, 255, 0);
						}
					}

					if (bar_height > 50) {
						//draw_pixel(&surf, j, height_tmp, 255, 0, 0);
						for (int i=bar_begin; i<=bar_end; i++) {
							barcount[i]++;
						}
						draw_line(&surf, j, bar_begin, j, bar_end, 255, 0, 0);
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
			draw_line(&surf, 0, line_where_barcode_begins, 0, line_where_barcode_ends, 255, 255, 0);
		}

		int center = line_where_barcode_begins + barcode_height / 2;
		printf("%d pixel high barcode found centered around line %d (start %d, end %d)\n",
				barcode_height, center, line_where_barcode_begins, line_where_barcode_ends);
		if (barcode_height > 90) {
			//printf("WARN: very high barcode: %d\n", barcode_height);
		}
	}

	if (jpegFilename) {
		//printf("  jpegWrite()\n");
		// write jpeg
		jpegWrite(dst_copy);
	}

	if (!do_fft) {
		free(dst);
		free(dst_copy);
		return;
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

	//printf("imageProcess() end\n");
	free(dst);
	free(dst_copy);
}

/**
	read single frame
*/
static int frameRead(void)
{
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

			imageProcess(buffers[0].start);
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

			imageProcess(buffers[buf.index].start);

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

				imageProcess((void *)buf.m.userptr);

				if (-1 == xioctl(fd, VIDIOC_QBUF, &buf))
					errno_exit("VIDIOC_QBUF");
				break;
#endif
	}

	return 1;
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
	fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_YUYV; //V4L2_PIX_FMT_RGB555X; //V4L2_PIX_FMT_RGB24;//  V4L2_PIX_FMT_YUYV
	fmt.fmt.pix.field = V4L2_FIELD_SEQ_TB;

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
		"-b | --barcode       detect barcode\n"
		"-q | --quality       JPEG quality (0-100)\n"
		"-m | --mmap          Use memory mapped buffers\n"
		"-r | --read          Use read() calls\n"
		"-u | --userptr       Use application allocated buffers\n"
		"-W | --width         width\n"
		"-H | --height        height\n"
		"",
		argv[0]);
	}

static const char short_options [] = "d:ho:sfbq:mruW:H:";

static const struct option
long_options [] = {
	{ "device",     required_argument,      NULL,           'd' },
	{ "help",       no_argument,            NULL,           'h' },
	{ "output",     required_argument,      NULL,           'o' },
	{ "fft",        no_argument,            NULL,           'f' },
	{ "barcode",    no_argument,            NULL,           'b' },
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
				break;

			case 'q':
				// set jpeg quality
				jpegQuality = atoi(optarg);
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
				width = atoi(optarg);
				break;

			case 'H':
				// set height
				height = atoi(optarg);
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
