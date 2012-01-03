#ifndef FFT_BRIN  // FFT_BRIN
#define FFT_BRIN

#ifndef NUM_FFT
#error You have to define NUM_FFT before including the FFT lib
#endif

#include "FFT_Code_Tables.h"
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
//#define NO_IMAG

//-----------------------------------------------------------------------------
// Global CONSTANTS and Variable Type Definitions
//-----------------------------------------------------------------------------
#define NUM_BITS 16 // Number of Bits in Data
// cause the program to stop after one
// data set.
typedef union IBALONG { // Integer or Byte-addressable LONG
  long l; // long: Var.l
  uint16_t i[2]; // u int16_t: Var.i[0]:Var.i[1]
  uint8_t b[4]; // u char: Var.b[0]:Var.b[1]:
  // Var.b[2]:Var.b[3]
} IBALONG;

typedef union BAINT { // Byte-addressable INT
  int16_t i; // int16_t: Var.i
  uint8_t b[2]; // u char: Var.b[0]:Var.b[1]
} BAINT;

//-----------------------------------------------------------------------------
// Function PROTOTYPES
//-----------------------------------------------------------------------------
void fft_ (int16_t Win_Array[], uint8_t SE_data);
void fft_Int_FFT(int16_t ReArray[], int16_t ImArray[]);
void fft_Bit_Reverse(int16_t BR_Array[]);
void fft_WindowCalc(int16_t Win_Array[], uint8_t SE_data);
void fft_mag(uint16_t L, int16_t * result);
void fft_phase(uint16_t L, int16_t * result);
void fft_magphase(uint16_t L, int16_t * mag_out, int16_t * phase_out);

//uint16_t BinNum;

extern int16_t Real[NUM_FFT];
extern int16_t Imag[NUM_FFT];

// bit Conversion_Set_Complete; // This indicates when the data has been
// stored, and is ready to be processed
// using the FFT routines

void fft_mag(uint16_t L, int16_t * result);
void fft_Int_FFT(int16_t ReArray[], int16_t ImArray[]);
void fft_Bit_Reverse(int16_t BR_Array[]);
void fft_WindowCalc(int16_t Win_Array[], uint8_t SE_data);
void fft_brin(int16_t samples[]);


#endif  // End define FFT_BRIN
