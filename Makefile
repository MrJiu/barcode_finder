CC = $(CROSS_COMPILE)gcc
CFLAGS += -std=c99 -Wall -I ./fft/ -DITU_R_INT -DNUM_FFT=512
LDFLAGS += -lm

build: barcode_finder

barcode_finder: main.c fft/fft_brin.c fft/FFT_Code_Tables.c
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f barcode_finder

rebuild: clean build

.PHONY: clean rebuild
