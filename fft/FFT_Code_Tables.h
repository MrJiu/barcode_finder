//-----------------------------------------------------------------------------
// FFT_Code_Tables.h
//-----------------------------------------------------------------------------
// Copyright 2003 Cygnal Integrated Products, Inc.
//
// AUTH: BD
// DATE: 30 JAN 03
//
// This header file is used to provide Sine, Bit Reversal, and Window tables
// for calculating an FFT. The tables are stored in FLASH memory (code space).
// All of the tables are conditionally used at compile time, so only one table
// of each type is used in the software.
//
// Target: C8051F12x
// Tool chain: KEIL C51 6.03
//
//#define NUM_FFT 128 // Length of FFT to process Define this outside!
// Must be 2^N, where N is an integer >= 2
//#define WINDOW_TYPE 4 // WINDOW_TYPE specifies the window to use on the data
// The available window functions are:
// 0 = No Window
// 1 = Triangle Window
// 2 = Hanning Window
// 3 = Hamming Window
// 4 = Blackman Window
// SinTable[] - SIN Tables are first 1/4 of a SIN wave - used to perform
// complex math functions. These are encoded such that a value of 1.0
// corresponds to 32768, and a value of -1.0 corresponds to -32768.
// BRTable[] - Bit Reversal tables are used to bit-reverse sort the data and
// perform other indexing functions. The Bit Reversal tables are stored
// as 1/2 of their actual value, and the real value is computed at
// runtime.
// WindowFunc[] - Tables used to window data. These are encoded such that
// 1.0 corresponds to 65536, and 0.0 corresponds to 0.
//-----------------------------------------------------------------------------
// SIN and BR Tables for NUM_FFT = 1024
//-----------------------------------------------------------------------------
#ifndef FFT_SINTAB_H
#define FFT_SINTAB_H
#if (NUM_FFT == 1024)
extern int SinTable[256];
extern unsigned int BRTable[512];
#endif

//-----------------------------------------------------------------------------
// SIN and BR Tables for NUM_FFT = 512
//-----------------------------------------------------------------------------
#if (NUM_FFT == 512)
extern int SinTable[128];
extern unsigned char BRTable[256];
#endif

//-----------------------------------------------------------------------------
// SIN and BR Tables for NUM_FFT = 256
//-----------------------------------------------------------------------------
#if (NUM_FFT == 256)
extern int SinTable[64];
extern unsigned char BRTable[128];
#endif

//-----------------------------------------------------------------------------
// SIN and BR Tables for NUM_FFT = 128
//-----------------------------------------------------------------------------
#if (NUM_FFT == 128)
extern int SinTable[32];
extern unsigned char BRTable[64];
#endif

//-----------------------------------------------------------------------------
// SIN and BR Tables for NUM_FFT = 64
//-----------------------------------------------------------------------------
#if (NUM_FFT == 64)
extern int SinTable[16];
extern unsigned char BRTable[32];
#endif

//-----------------------------------------------------------------------------
// SIN and BR Tables for NUM_FFT = 32
//-----------------------------------------------------------------------------
#if (NUM_FFT == 32)
extern int SinTable[8];
extern unsigned char BRTable[16];
#endif

//-----------------------------------------------------------------------------
// SIN and BR Tables for NUM_FFT = 16
//-----------------------------------------------------------------------------
#if (NUM_FFT == 16)
extern int SinTable[4];
extern unsigned char BRTable[8];
#endif

//-----------------------------------------------------------------------------
// SIN and BR Tables for NUM_FFT = 8
//-----------------------------------------------------------------------------
#if (NUM_FFT == 8)
extern int SinTable[2];
extern unsigned char BRTable[4];
#endif

//-----------------------------------------------------------------------------
// SIN and BR Table for NUM_FFT = 4
//-----------------------------------------------------------------------------
#if (NUM_FFT == 4)
extern int SinTable[1];
extern unsigned char BRTable[2];
#endif

//-----------------------------------------------------------------------------
// Window Functions for NUM_FFT = 1024
//-----------------------------------------------------------------------------
#if (NUM_FFT == 1024)
#if (WINDOW_TYPE == 1)
// Triangle Window
extern unsigned int WindowFunc[512];
#endif
#if (WINDOW_TYPE == 2)
// Hanning Window
extern unsigned int WindowFunc[512];
#endif
#if (WINDOW_TYPE == 3)
// Hamming Window
extern unsigned int WindowFunc[512];
#endif
#if (WINDOW_TYPE == 4)
// Blackman Window
extern unsigned int WindowFunc[512];
#endif
#endif

//-----------------------------------------------------------------------------
// Window Functions for NUM_FFT = 512
//-----------------------------------------------------------------------------
#if (NUM_FFT == 512)
#if (WINDOW_TYPE == 1)
// Triangle Window
extern unsigned int WindowFunc[256];
#endif
#if (WINDOW_TYPE == 2)
// Hanning Window
extern unsigned int WindowFunc[256];
#endif
#if (WINDOW_TYPE == 3)
// Hamming Window
extern unsigned int WindowFunc[256];
#endif
#if (WINDOW_TYPE == 4)
// Blackman Window
extern unsigned int WindowFunc[256];
#endif
#endif

//-----------------------------------------------------------------------------
// Window Functions for NUM_FFT = 256
//-----------------------------------------------------------------------------
#if (NUM_FFT == 256)
#if (WINDOW_TYPE == 1)
// Triangle Window
extern unsigned int WindowFunc[128];
#endif
#if (WINDOW_TYPE == 2)
// Hanning Window
extern unsigned int WindowFunc[128];
#endif
#if (WINDOW_TYPE == 3)
// Hamming Window
extern unsigned int WindowFunc[128];
#endif
#if (WINDOW_TYPE == 4)
// Blackman Window
extern unsigned int WindowFunc[128];
#endif
#endif

//-----------------------------------------------------------------------------
// Window Functions for NUM_FFT = 128
//-----------------------------------------------------------------------------
#if (NUM_FFT == 128)
#if (WINDOW_TYPE == 1)
// Triangle Window
extern unsigned int WindowFunc[64];
#endif
#if (WINDOW_TYPE == 2)
// Hanning Window
extern unsigned int WindowFunc[64];
#endif
#if (WINDOW_TYPE == 3)
// Hamming Window
extern unsigned int WindowFunc[64];
#endif
#if (WINDOW_TYPE == 4)
// Blackman Window
extern unsigned int WindowFunc[64];
#endif
#endif

//-----------------------------------------------------------------------------
// Window Functions for NUM_FFT = 64
//-----------------------------------------------------------------------------
#if (NUM_FFT == 64)
#if (WINDOW_TYPE == 1)
// Triangle Window
extern unsigned int WindowFunc[32];
#endif
#if (WINDOW_TYPE == 2)
// Hanning Window
extern unsigned int WindowFunc[32];
#endif
#if (WINDOW_TYPE == 3)
// Hamming Window
extern unsigned int WindowFunc[32];
#endif
#if (WINDOW_TYPE == 4)
// Blackman Window
extern unsigned int WindowFunc[32];
#endif
#endif

//-----------------------------------------------------------------------------
// Window Functions for NUM_FFT = 32
//-----------------------------------------------------------------------------
#if (NUM_FFT == 32)
#if (WINDOW_TYPE == 1)
// Triangle Window
extern unsigned int WindowFunc[16];
#endif
#if (WINDOW_TYPE == 2)
// Hanning Window
extern unsigned int WindowFunc[16];
#endif
#if (WINDOW_TYPE == 3)
// Hamming Window
extern unsigned int WindowFunc[16];
#endif
#if (WINDOW_TYPE == 4)
// Blackman Window
extern unsigned int WindowFunc[16];
#endif
#endif

//-----------------------------------------------------------------------------
// Window Functions for NUM_FFT = 16
//-----------------------------------------------------------------------------
#if (NUM_FFT == 16)
#if (WINDOW_TYPE == 1)
// Triangle Window
extern unsigned int WindowFunc[8];
#endif
#if (WINDOW_TYPE == 2)
// Hanning Window
extern unsigned int WindowFunc[8];
#endif
#if (WINDOW_TYPE == 3)
// Hamming Window
extern unsigned int WindowFunc[8];
#endif
#if (WINDOW_TYPE == 4)
// Blackman Window
extern unsigned int WindowFunc[8];
#endif
#endif

//-----------------------------------------------------------------------------
// Window Functions for NUM_FFT = 8
//-----------------------------------------------------------------------------
#if (NUM_FFT == 8)
#if (WINDOW_TYPE == 1)
// Triangle Window
extern unsigned int WindowFunc[4];
#endif
#if (WINDOW_TYPE == 2)
// Hanning Window
extern unsigned int WindowFunc[4];
#endif
#if (WINDOW_TYPE == 3)
// Hamming Window
extern unsigned int WindowFunc[4];
#endif
#if (WINDOW_TYPE == 4)
// Blackman Window
extern unsigned int WindowFunc[4];
#endif
#endif

//-----------------------------------------------------------------------------
// Window Functions for NUM_FFT = 4
//-----------------------------------------------------------------------------
#if (NUM_FFT == 4)
#if (WINDOW_TYPE == 1)
// Triangle Window
extern unsigned int WindowFunc[2];
#endif
#if (WINDOW_TYPE == 2)
// Hanning Window
extern unsigned int WindowFunc[2];
#endif
#if (WINDOW_TYPE == 3)
// Hamming Window
extern unsigned int WindowFunc[2];
#endif
#if (WINDOW_TYPE == 4)
// Blackman Window
extern unsigned int WindowFunc[2];
#endif
#endif

#endif  // End define FFT_SINTAB_H
