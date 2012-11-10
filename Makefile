CC = $(CROSS_COMPILE)gcc
CFLAGS += -g -std=c99 -Wall -I ./fft/ -DITU_R_INT -DNUM_FFT=512
LDFLAGS += -lm -ljpeg

SRC = main.c fft/fft_brin.c fft/FFT_Code_Tables.c

all: build

build: barcode_finder

barcode_finder: $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $@ $(LDFLAGS)

clean:
	rm -f barcode_finder

rebuild: clean build

.PHONY: clean rebuild
