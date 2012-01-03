/*
 * fft_brin.c
 *
 * Created: 06.12.2011 17:09:35
 *  Author: kya
 */ 

#include "fft_brin.h"

//-----------------------------------------------------------------------------
// Global Variables
//-----------------------------------------------------------------------------
// storage of FFT: requires NUM_FFT*4 Bytes after DATA_BEGIN address
//__no_init 
int16_t Real[NUM_FFT]; // = {0,2187,4248,6066,7546,8621,9258,9463,9276,8770,8046,7221,6420,5764,5359,5289,5607,6330,7438,8876,10559,12376,14204,15917,17394,18533,19258,19526,19331,18703,17706,16433,15000,13531,12150,10973,10092,9575,9450,9714,10323,11203,12253,13356,14385,15220,15753,15899,15607,14858,13673,12107,10247,8202,6095,4056,2205,644,-550,-1331,-1693,-1664,-1308,-715,0,715,1308,1664,1693,1331,550,-644,-2205,-4056,-6095,-8202,-10247,-12107,-13673,-14858,-15607,-15899,-15753,-15220,-14385,-13356,-12253,-11203,-10323,-9714,-9450,-9575,-10092,-10973,-12150,-13531,-15000,-16433,-17706,-18703,-19331,-19526,-19258,-18533,-17394,-15917,-14204,-12376,-10559,-8876,-7438,-6330,-5607,-5289,-5359,-5764,-6420,-7221,-8046,-8770,-9276,-9463,-9258,-8621,-7546,-6066,-4248,-2187}; // Testdata already in here!
//__no_init 
int16_t Imag[NUM_FFT]; //
// NUM_FFT is defined in the “FFT_Code_Tables.h” header file

#if (NUM_FFT >= 256)
uint16_t fft_index, index, ADC_Index;
#endif

#if (NUM_FFT < 256)
uint8_t fft_index, index, ADC_Index;
#endif


void fft_brin(int16_t samples[])
{
  fft_WindowCalc(samples, 0);    // Window Real Data, do not convert from single ended to signed
  fft_Bit_Reverse(samples);      // Sort Real (Input) Data in bit-reverse order
  fft_Int_FFT(samples, Imag);    // Perform FFT on data
}


////-----------------------------------------------------------------------------
//// WindowCalc
////-----------------------------------------------------------------------------
////
//// Uses the values in WindowFunc[] to window the stored data.
////
//// The WindowFunc[] Array contains window coefficients between samples
//// 0 and (NUM_FFT/2)-1, and samples from NUM_FFT/2 to NUM_FFT-1 are the mirror
//// image of the other side.
//// Window values are interpreted as a fraction of 1 (WindowFunc[x]/65536).
//// The value at NUM_FFT/2 is assumed to be 1.0 (65536).
////
//// If SE_data = 1, the input data is assumed to be single-ended, and is
//// converted to a 2’s complement, differential representation, to cancel the DC
//// offset.
////
//void fft_WindowCalc(int16_t Win_Array[], uint8_t SE_data)
//{
//#if (WINDOW_TYPE != 0) // Use this section if a window has been specified
  //IBALONG NewVal;
  //if (SE_data) // If data is single-ended,
    //Win_Array[0] ^= 0x8000; // convert it to differential
  //NewVal.l = (long)Win_Array[0] * WindowFunc[0];
  //if ((NewVal.l < 0)&&(NewVal.i[0]))      // Swapped Kyrre
    //Win_Array[0] = NewVal.i[1] + 1;       // Swapped Kyrre
  //else Win_Array[0] = NewVal.i[1];        // Swapped Kyrre
  //if (SE_data) // If data is single-ended,
    //Win_Array[NUM_FFT/2] ^= 0x8000; // convert it to differential
  //for (fft_index = 1; fft_index < NUM_FFT/2; fft_index++)
  //{
    //// Array positions 1 to (NUM_FFT/2 - 1)
    //if (SE_data) // If data is single-ended,
      //Win_Array[fft_index] ^= 0x8000; // convert it to differential
    //NewVal.l = (long)Win_Array[fft_index] * WindowFunc[fft_index];
    //if ((NewVal.l < 0)&&(NewVal.i[0]))      // Swapped Kyrre
      //Win_Array[fft_index] = NewVal.i[1] + 1;   // Swapped Kyrre
    //else Win_Array[fft_index] = NewVal.i[1];    // Swapped Kyrre
    //// Array positions (NUM_FFT/2 + 1) to (NUM_FFT - 1)
    //if (SE_data) // If data is single-ended,
      //Win_Array[NUM_FFT-fft_index] ^= 0x8000; // convert it to differential
    //NewVal.l = (long)Win_Array[NUM_FFT-fft_index] * WindowFunc[fft_index];
    //if ((NewVal.l < 0)&&(NewVal.i[0]))              // Swapped Kyrre
      //Win_Array[NUM_FFT-fft_index] = NewVal.i[1] + 1;   // Swapped Kyrre
    //else Win_Array[NUM_FFT-fft_index] = NewVal.i[1];    // Swapped Kyrre
  //}
//#endif
//#if (WINDOW_TYPE == 0) // Compile this if no window has been specified
  //if (SE_data) // If data is single-ended,
  //{ // convert it to differential
    //for (fft_index = 0; fft_index < NUM_FFT; fft_index++)
    //{
      //Win_Array[fft_index] ^= 0x8000; // XOR MSB with ‘1’ to invert
    //}
  //}
//#endif
//} // END WindowCalc

//-----------------------------------------------------------------------------    
// WindowCalc    
//-----------------------------------------------------------------------------    
//    
// Uses the values in WindowFunc[] to window the stored data.    
//    
// The WindowFunc[] Array contains window coefficients between samples    
// 0 and (NUM_FFT/2)-1, and samples from NUM_FFT/2 to NUM_FFT-1 are the mirror    
// image of the other side.    
// Window values are interpreted as a fraction of 1 (WindowFunc[x]/65536).    
// The value at NUM_FFT/2 is assumed to be 1.0 (65536).    
//    
// If SE_data = 1, the input data is assumed to be single-ended, and is    
// converted to a 2's complement, differential representation, to cancel the DC    
// offset.    
//    
void fft_WindowCalc(int16_t Win_Array[], uint8_t SE_data)
{   
   
#if (WINDOW_TYPE != 0)  // Use this section if a window has been specified    
     
    IBALONG NewVal;   
   
   if (SE_data)                              // If data is single-ended,    
      Win_Array[0] ^= 0x8000;                // convert it to differential    
   NewVal.l = (long)Win_Array[0] * WindowFunc[0];   
     
   if ((NewVal.l < 0)&&(NewVal.i[1]))   
     Win_Array[0] = NewVal.i[0] + 1;   
   else Win_Array[0] = NewVal.i[0];   
   
   if (SE_data)                              // If data is single-ended,    
      Win_Array[NUM_FFT/2] ^= 0x8000;        // convert it to differential    
   
  for (index = 1; index < NUM_FFT/2; index++)   
  {   
      // Array positions 1 to (NUM_FFT/2 - 1)    
      if (SE_data)                           // If data is single-ended,    
         Win_Array[index] ^= 0x8000;         // convert it to differential    
      NewVal.l = (long)Win_Array[index] * WindowFunc[index];   
      if ((NewVal.l < 0)&&(NewVal.i[1]))   
        Win_Array[index] = NewVal.i[0] + 1;   
   
 else Win_Array[index] = NewVal.i[0];   
   
      // Array positions (NUM_FFT/2 + 1) to (NUM_FFT - 1)    
      if (SE_data)                           // If data is single-ended,    
         Win_Array[NUM_FFT-index] ^= 0x8000; // convert it to differential    
      NewVal.l = (long)Win_Array[NUM_FFT-index] * WindowFunc[index];   
       
      if ((NewVal.l < 0)&&(NewVal.i[1]))   
         Win_Array[NUM_FFT-index] = NewVal.i[0] + 1;   
      else Win_Array[NUM_FFT-index] = NewVal.i[0];   
   
   }   
   
#endif    
   
#if (WINDOW_TYPE == 0)  // Compile this if no window has been specified    
   
   if (SE_data)                              // If data is single-ended,    
   {                                         // convert it to differential    
   
      for (index = 0; index < NUM_FFT; index++)   
      {   
         Win_Array[index] ^= 0x8000;         // XOR MSB with '1' to invert    
      }   
   }   
   
#endif    
   
}  // END WindowCalc    
   



//-----------------------------------------------------------------------------
// Bit_Reverse
//-----------------------------------------------------------------------------
//
// Sorts data in Bit Reversed Address order
//
// The BRTable[] array is used to find which values must be swapped. Only
// half of this array is stored, to save code space. The second half is
// assumed to be a mirror image of the first half.
//
void fft_Bit_Reverse(int16_t BR_Array[])
{
#if (NUM_FFT >= 512)
  uint16_t swapA, swapB, sw_cnt; // Swap Indices
#endif
#if (NUM_FFT <= 256)
  uint8_t swapA, swapB, sw_cnt; // Swap Indices
#endif
  int16_t TempStore;
  // Loop through locations to swap
  for (sw_cnt = 1; sw_cnt < NUM_FFT/2; sw_cnt++)
  {
    swapA = sw_cnt; // Store current location
    swapB = BRTable[sw_cnt] * 2; // Retrieve bit-reversed index
    if (swapB > swapA) // If the bit-reversed index is
    { // larger than the current index,
      TempStore = BR_Array[swapA]; // the two data locations are
      BR_Array[swapA] = BR_Array[swapB]; // swapped. Using this comparison
      BR_Array[swapB] = TempStore; // ensures that locations are only
    } // swapped once, and never with
    // themselves
    swapA += NUM_FFT/2; // Now perform the same operations
    swapB++; // on the second half of the data
    if (swapB > swapA)
    {
      TempStore = BR_Array[swapA];
      BR_Array[swapA] = BR_Array[swapB];
      BR_Array[swapB] = TempStore;
    }
  }
} // END Bit Reverse Order Sort


////-----------------------------------------------------------------------------
//// Int_FFT
////-----------------------------------------------------------------------------
////
//// Performs a Radix-2 Decimation-In-Time FFT on the input array ReArray[]
////
//// During each stage of the FFT, the values are calculated using a set of
//// “Butterfly” equations, as listed below:
////
//// Re1 = Re1 + (Cos(x)*Re2 + Sin(x)*Im2)
//// Re2 = Re1 - (Cos(x)*Re2 + Sin(x)*Im2)
//// Im1 = Im1 + (Cos(x)*Im2 - Sin(x)*Re2)
//// Im2 = Im1 - (Cos(x)*Im2 - Sin(x)*Re2)
////
//// The routine implements this calculation using the following values:
////
//// Re1 = ReArray[indexA], Re2 = ReArray[indexB]
//// Im1 = ImArray[indexA], Im2 = ImArray[indexB]
//// x = the angle: 2*pi*(sin_index/NUM_FFT), in radians. The necessary values
//// are stored in code space in the SinTable[] array.
////
////
//// Key Points for using this FFT routine:
////
//// 1) It expects REAL data (in ReArray[]), in 2’s complement, 16-bit binary
//// format and assumes a value of 0 for all imaginary locations
//// (in ImArray[]).
////
//// 2) It expects the REAL input data to be sorted in bit-reversed index order.
////
//// 3) SIN and COS values are retrieved and calculated from a table consisting
//// of 1/4 of a period of a SIN function.
////
//// 4) It is optimized to use integer math only (no floating-point operations),
//// and for storage space. The input, all intermediate stages, and the
//// output of the FFT are stored as 16-bit INTEGER values. This limits the
//// precision of the routine. When using input data of less than 16-bits,
//// the best results are produced by left-justifying the data prior to
//// windowing and performing the FFT.
////
//// 5) The algorithm is a Radix-2 type, meaning that the number of samples must
//// be 2^N, where N is an integer. The minimum number of samples to process
//// is 4. The constant NUM_FFT contains the number of samples to process.
////
////
//void fft_Int_FFT(int16_t ReArray[], int16_t ImArray[])
//{
//#if (NUM_FFT >= 512)
//  uint16_t sin_index, g_cnt, s_cnt; // Keeps track of the proper index
//  uint16_t indexA, indexB; // locations for each calculation
//#endif
//#if (NUM_FFT <= 256)
//  uint8_t sin_index, g_cnt, s_cnt; // Keeps track of the proper index
//  uint8_t indexA, indexB; // locations for each calculation
//#endif
//  uint16_t group = NUM_FFT/4, stage = 2;
//  long CosVal, SinVal;
//  long TempImA, TempImB, TempReA, TempReB, TempReA2, TempReB2;
//  IBALONG ReTwid, ImTwid, TempL;
//  // FIRST STAGE - optimized for REAL input data only. This will set all
//  // Imaginary locations to zero.
//  //
//  // Shortcuts have been taken to remove unnecessary multiplications during this
//  // stage. The angle “x” is 0 radians for all calculations at this point, so
//  // the SIN value is equal to 0.0 and the COS value is equal to 1.0.
//  // Additionally, all Imaginary locations are assumed to be ‘0’ in this stage of
//  // the algorithm, and are set to ‘0’.
//  indexA = 0;
//  for (g_cnt = 0; g_cnt < NUM_FFT/2; g_cnt++)
//  {
//    indexB = indexA + 1;
//    TempReA = ReArray[indexA];
//    TempReB = ReArray[indexB];
//    // Calculate new value for ReArray[indexA]
//    TempL.l = (long)TempReA + TempReB;
//    if ((TempL.l < 0)&&(0x01 & TempL.b[0]))   // Swapped Kyrre
//      TempReA2 = (TempL.l >> 1) + 1;
//    else TempReA2 = TempL.l >> 1;
//    // Calculate new value for ReArray[indexB]
//    TempL.l = (long)TempReA - TempReB;
//    if ((TempL.l < 0)&&(0x01 & TempL.b[0]))   // Swapped Kyrre
//      ReArray[indexB] = (TempL.l >> 1) + 1;
//    else ReArray[indexB] = TempL.l >> 1;
//    ReArray[indexA] = TempReA2;
//    ImArray[indexA] = 0; // set Imaginary locations to ‘0’
//    ImArray[indexB] = 0;
//    indexA = indexB + 1;
//  }
//  // END OF FIRST STAGE
//  while (stage <= NUM_FFT/2)
//  {
//    indexA = 0;
//    sin_index = 0;
//    for (g_cnt = 0; g_cnt < group; g_cnt++)
//    {
//      for (s_cnt = 0; s_cnt < stage; s_cnt++)
//      {
//        indexB = indexA + stage;
//        TempReA = ReArray[indexA];
//        TempReB = ReArray[indexB];
//#ifndef NO_IMAG
//        TempImA = ImArray[indexA];
//        TempImB = ImArray[indexB];
//#endif
//        // The following first checks for the special cases when the angle “x” is
//        // equal to either 0 or pi/2 radians. In these cases, unnecessary
//        // multiplications have been removed to improve the processing speed.
//        if (sin_index == 0) // corresponds to “x” = 0 radians
//        {
//          // Calculate new value for ReArray[indexA]
//          TempL.l = (long)TempReA + TempReB;
//          if ((TempL.l < 0)&&(0x01 & TempL.b[0]))   // Swapped Kyrre
//            TempReA2 = (TempL.l >> 1) + 1;
//          else TempReA2 = TempL.l >> 1;
//          // Calculate new value for ReArray[indexB]
//          TempL.l = (long)TempReA - TempReB;
//          if ((TempL.l < 0)&&(0x01 & TempL.b[0]))   // Swapped Kyrre
//            TempReB2 = (TempL.l >> 1) + 1;
//          else TempReB2 = TempL.l >> 1;
//#ifndef NO_IMAG
//          // Calculate new value for ImArray[indexB]
//          TempL.l = (long)TempImA - TempImB;
//          if ((TempL.l < 0)&&(0x01 & TempL.b[0]))   // Swapped Kyrre
//            TempImB = (TempL.l >> 1) + 1;
//          else TempImB = TempL.l >> 1;
//          // Calculate new value for ImArray[indexA]
//          TempL.l = (long)TempImA + TempImB;
//          if ((TempL.l < 0)&&(0x01 & TempL.b[0]))   // Swapped Kyrre
//            TempImA = (TempL.l >> 1) + 1;
//          else TempImA = TempL.l >> 1;
//#endif
//        }
//        else if (sin_index == NUM_FFT/4) // corresponds to “x” = pi/2 radians
//        {
//          // Calculate new value for ReArray[indexB]
//          TempL.l = (long)TempReA - TempImB;
//          if ((TempL.l < 0)&&(0x01 & TempL.b[0]))   // Swapped Kyrre
//            TempReB2 = (TempL.l >> 1) + 1;
//          else TempReB2 = TempL.l >> 1;
//          // Calculate new value for ReArray[indexA]
//          TempL.l = (long)TempReA + TempImB;
//          if ((TempL.l < 0)&&(0x01 & TempL.b[0]))   // Swapped Kyrre
//            TempReA2 = (TempL.l >> 1) + 1;
//          else TempReA2 = TempL.l >> 1;
//#ifndef NO_IMAG
//          // Calculate new value for ImArray[indexB]
//          TempL.l = (long)TempImA + TempReB;
//          if ((TempL.l < 0)&&(0x01 & TempL.b[0]))   // Swapped Kyrre
//            TempImB = (TempL.l >> 1) + 1;
//          else TempImB = TempL.l >> 1;
//          // Calculate new value for ImArray[indexA]
//          TempL.l = (long)TempImA - TempReB;
//          if ((TempL.l < 0)&&(0x01 & TempL.b[0]))   // Swapped Kyrre
//            TempImA = (TempL.l >> 1) + 1;
//          else TempImA = TempL.l >> 1;
//#endif
//        }
//        else
//        {
//          // If no multiplication shortcuts can be taken, the SIN and COS
//          // values for the Butterfly calculation are fetched from the
//          // SinTable[] array.
//          if (sin_index > NUM_FFT/4)
//          {
//            SinVal = SinTable[(NUM_FFT/2) - sin_index];
//            CosVal = -SinTable[sin_index - (NUM_FFT/4)];
//          }
//          else
//          {
//            SinVal = SinTable[sin_index];
//            CosVal = SinTable[(NUM_FFT/4) - sin_index];
//          }
//          // The SIN and COS values are used here to calculate part of the
//          // Butterfly equation
//          ReTwid.l = ((long)TempReB * CosVal) +
//            ((long)TempImB * SinVal);
//          ImTwid.l = ((long)TempImB * CosVal) -
//            ((long)TempReB * SinVal);
//          // Using the values calculated above, the new variables
//          // are computed
//          // Calculate new value for ReArray[indexA]
////          TempL.i[1] = 0;
////          TempL.i[0] = TempReA;
//          TempL.i[0] = 0;                     // Swapped Kyrre
//          TempL.i[1] = TempReA;               // Swapped Kyrre
//          TempL.l = TempL.l >> 1;
//          ReTwid.l += TempL.l;
////          if ((ReTwid.l < 0)&&(ReTwid.i[1]))
////            TempReA2 = ReTwid.i[0] + 1;
////          else TempReA2 = ReTwid.i[0];
//          if ((ReTwid.l < 0)&&(ReTwid.i[0]))  // Swapped Kyrre
//            TempReA2 = ReTwid.i[1] + 1;       // Swapped Kyrre
//          else TempReA2 = ReTwid.i[1];        // Swapped Kyrre
//          // Calculate new value for ReArray[indexB]
//          TempL.l = TempL.l << 1;
//          TempL.l -= ReTwid.l;
////          if ((TempL.l < 0)&&(TempL.i[1]))
////            TempReB2 = TempL.i[0] + 1;
////          else TempReB2 = TempL.i[0];
//          if ((TempL.l < 0)&&(TempL.i[0]))    // Swapped Kyrre
//            TempReB2 = TempL.i[1] + 1;        // Swapped Kyrre
//          else TempReB2 = TempL.i[1];         // Swapped Kyrre
//#ifndef NO_IMAG
//          // Calculate new value for ImArray[indexA]
////          TempL.i[1] = 0;
////          TempL.i[0] = TempImA;
//          TempL.i[0] = 0;                     // Swapped Kyrre
//          TempL.i[1] = TempImA;               // Swapped Kyrre
//          TempL.l = TempL.l >> 1;
//          ImTwid.l += TempL.l;
////          if ((ImTwid.l < 0)&&(ImTwid.i[1]))
////            TempImA = ImTwid.i[0] + 1;
////          else TempImA = ImTwid.i[0];
//          if ((ImTwid.l < 0)&&(ImTwid.i[0]))  // Swapped Kyrre
//            TempImA = ImTwid.i[1] + 1;        // Swapped Kyrre
//          else TempImA = ImTwid.i[1];         // Swapped Kyrre
//          // Calculate new value for ImArray[indexB]
//          TempL.l = TempL.l << 1;
//          TempL.l -= ImTwid.l;
////          if ((TempL.l < 0)&&(TempL.i[1]))
////            TempImB = TempL.i[0] + 1;
////          else TempImB = TempL.i[0];
//          if ((TempL.l < 0)&&(TempL.i[0]))    // Swapped Kyrre
//            TempImB = TempL.i[1] + 1;         // Swapped Kyrre
//          else TempImB = TempL.i[1];          // Swapped Kyrre
//#endif
//        }
//        ReArray[indexA] = TempReA2;
//        ReArray[indexB] = TempReB2;
//#ifndef NO_IMAG
//        ImArray[indexA] = TempImA;
//        ImArray[indexB] = TempImB;
//#endif
//        indexA++;
//        sin_index += group;
//      } // END of stage FOR loop (s_cnt)
//      indexA = indexB + 1;
//      sin_index = 0;
//    } // END of group FOR loop (g_cnt)
//    group /= 2;
//    stage *= 2;
//    
//  } // END of While loop
//} // END Int_FFT

//-----------------------------------------------------------------------------    
// Int_FFT    
//-----------------------------------------------------------------------------    
//    
// Performs a Radix-2 Decimation-In-Time FFT on the input array ReArray[]    
//    
// During each stage of the FFT, the values are calculated using a set of    
// "Butterfly" equations, as listed below:    
//    
// Re1 = Re1 + (Cos(x)*Re2 + Sin(x)*Im2)    
// Re2 = Re1 - (Cos(x)*Re2 + Sin(x)*Im2)    
// Im1 = Im1 + (Cos(x)*Im2 - Sin(x)*Re2)    
// Im2 = Im1 - (Cos(x)*Im2 - Sin(x)*Re2)    
//    
// The routine implements this calculation using the following values:    
//    
// Re1 = ReArray[indexA], Re2 = ReArray[indexB]    
// Im1 = ImArray[indexA], Im2 = ImArray[indexB]    
// x = the angle: 2*pi*(sin_index/NUM_FFT), in radians.  The necessary values    
//    are stored in code space in the SinTable[] array.    
//    
//    
// Key Points for using this FFT routine:    
//    
// 1) It expects REAL data (in ReArray[]), in 2's complement, 16-bit binary    
//    format and assumes a value of 0 for all imaginary locations    
//    (in ImArray[]).    
//    
// 2) It expects the REAL input data to be sorted in bit-reversed index order.    
//    
// 3) SIN and COS values are retrieved and calculated from a table consisting    
//    of 1/4 of a period of a SIN function.    
//    
// 4) It is optimized to use integer math only (no floating-point operations),    
//    and for storage space.  The input, all intermediate stages, and the    
//    output of the FFT are stored as 16-bit INTEGER values. This limits the    
//    precision of the routine.  When using input data of less than 16-bits,    
//    the best results are produced by left-justifying the data prior to    
//    windowing and performing the FFT.    
//    
// 5) The algorithm is a Radix-2 type, meaning that the number of samples must    
//    be 2^N, where N is an integer.  The minimum number of samples to process    
//    is 4.  The constant NUM_FFT contains the number of samples to process.    
//    
//    
   
//int xdata Real[NUM_FFT] _at_ DATA_BEGIN;   
//int xdata Imag[NUM_FFT] _at_ (DATA_BEGIN + (NUM_FFT * 2));   

void fft_Int_FFT(int16_t ReArray[], int16_t ImArray[])
{   
   
#if (NUM_FFT >= 512)    
unsigned int sin_index, g_cnt, s_cnt;        // Keeps track of the proper index    
unsigned int indexA, indexB;                 // locations for each calculation    
#endif    
   
#if (NUM_FFT <= 256)    
unsigned char sin_index, g_cnt, s_cnt;       // Keeps track of the proper index    
unsigned char indexA, indexB;                // locations for each calculation    
#endif    
   
unsigned int group = NUM_FFT/4, stage = 2;   
long CosVal, SinVal;   
long TempImA, TempImB, TempReA, TempReB, TempReA2, TempReB2;   
IBALONG ReTwid, ImTwid, TempL;   
   
// FIRST STAGE - optimized for REAL input data only.  This will set all    
// Imaginary locations to zero.    
//    
// Shortcuts have been taken to remove unnecessary multiplications during this    
// stage. The angle "x" is 0 radians for all calculations at this point, so    
// the SIN value is equal to 0.0 and the COS value is equal to 1.0.    
// Additionally, all Imaginary locations are assumed to be '0' in this stage of    
// the algorithm, and are set to '0'.    
   
   indexA = 0;   
   for (g_cnt = 0; g_cnt < NUM_FFT/2; g_cnt++)   
   {   
      indexB = indexA + 1;   
   
      TempReA = ReArray[indexA];   
      TempReB = ReArray[indexB];   
   
      // Calculate new value for ReArray[indexA]    
      TempL.l = (long)TempReA + TempReB;   
      if ((TempL.l < 0)&&(0x01 & TempL.b[3]))   
         TempReA2 = (TempL.l >> 1) + 1;   
      else TempReA2 = TempL.l >> 1;   
   
      // Calculate new value for ReArray[indexB]    
      TempL.l = (long)TempReA - TempReB;   
      if ((TempL.l < 0)&&(0x01 & TempL.b[3]))   
         ReArray[indexB] = (TempL.l >> 1) + 1;   
      else ReArray[indexB] = TempL.l >> 1;   
   
      ReArray[indexA] = TempReA2;   
   
      ImArray[indexA] = 0;                   // set Imaginary locations to '0'    
      ImArray[indexB] = 0;   
   
      indexA = indexB + 1;   
   }   
   
// END OF FIRST STAGE    
   
   
while (stage <= NUM_FFT/2)   
{   
   indexA = 0;   
   sin_index = 0;   
   
   for (g_cnt = 0; g_cnt < group; g_cnt++)   
   {   
      for (s_cnt = 0; s_cnt < stage; s_cnt++)   
      {   
         indexB = indexA + stage;   
   
         TempReA = ReArray[indexA];   
         TempReB = ReArray[indexB];   
         TempImA = ImArray[indexA];   
         TempImB = ImArray[indexB];   
   
// The following first checks for the special cases when the angle "x" is    
// equal to either 0 or pi/2 radians.  In these cases, unnecessary    
// multiplications have been removed to improve the processing speed.    
   
         if (sin_index == 0)  // corresponds to "x" = 0 radians    
         {   
   
            // Calculate new value for ReArray[indexA]    
            TempL.l = (long)TempReA + TempReB;   
            if ((TempL.l < 0)&&(0x01 & TempL.b[3]))   
               TempReA2 = (TempL.l >> 1) + 1;   
            else TempReA2 = TempL.l >> 1;   
   
            // Calculate new value for ReArray[indexB]    
            TempL.l = (long)TempReA - TempReB;   
            if ((TempL.l < 0)&&(0x01 & TempL.b[3]))   
               TempReB2 = (TempL.l >> 1) + 1;   
            else TempReB2 = TempL.l >> 1;   
   
            // Calculate new value for ImArray[indexB]    
           TempL.l = (long)TempImA - TempImB;   
            if ((TempL.l < 0)&&(0x01 & TempL.b[3]))   
               TempImB = (TempL.l >> 1) + 1;   
            else TempImB = TempL.l >> 1;   
   
            // Calculate new value for ImArray[indexA]    
            TempL.l = (long)TempImA + TempImB;   
            if ((TempL.l < 0)&&(0x01 & TempL.b[3]))   
               TempImA = (TempL.l >> 1) + 1;   
            else TempImA = TempL.l >> 1;   
   
         }   
         else if (sin_index == NUM_FFT/4) // corresponds to "x" = pi/2 radians    
         {   
   
            // Calculate new value for ReArray[indexB]    
            TempL.l = (long)TempReA - TempImB;   
            if ((TempL.l < 0)&&(0x01 & TempL.b[3]))   
            TempReB2 = (TempL.l >> 1) + 1;   
            else TempReB2 = TempL.l >> 1;   
   
            // Calculate new value for ReArray[indexA]    
            TempL.l = (long)TempReA + TempImB;   
            if ((TempL.l < 0)&&(0x01 & TempL.b[3]))   
               TempReA2 = (TempL.l >> 1) + 1;   
            else TempReA2 = TempL.l >> 1;   
   
            // Calculate new value for ImArray[indexB]    
            TempL.l = (long)TempImA + TempReB;   
            if ((TempL.l < 0)&&(0x01 & TempL.b[3]))   
               TempImB = (TempL.l >> 1) + 1;   
            else TempImB = TempL.l >> 1;   
   
            // Calculate new value for ImArray[indexA]    
            TempL.l = (long)TempImA - TempReB;   
            if ((TempL.l < 0)&&(0x01 & TempL.b[3]))   
               TempImA = (TempL.l >> 1) + 1;   
            else TempImA = TempL.l >> 1;   
   
         }   
   
         else   
         {   
            // If no multiplication shortcuts can be taken, the SIN and COS    
            // values for the Butterfly calculation are fetched from the    
            // SinTable[] array.    
   
            if (sin_index > NUM_FFT/4)   
            {   
               SinVal = SinTable[(NUM_FFT/2) - sin_index];   
               CosVal = -SinTable[sin_index - (NUM_FFT/4)];   
            }   
            else   
            {   
               SinVal = SinTable[sin_index];   
               CosVal = SinTable[(NUM_FFT/4) - sin_index];   
            }   
   
            // The SIN and COS values are used here to calculate part of the    
            // Butterfly equation    
            ReTwid.l = ((long)TempReB * CosVal) +   
                  ((long)TempImB * SinVal);   
   
            ImTwid.l = ((long)TempImB * CosVal) -   
                  ((long)TempReB * SinVal);   
   
            // Using the values calculated above, the new variables    
            // are computed    
   
            // Calculate new value for ReArray[indexA]    
            TempL.i[1] = 0;   
            TempL.i[0] = TempReA;   
            TempL.l = TempL.l >> 1;   
            ReTwid.l += TempL.l;   
            if ((ReTwid.l < 0)&&(ReTwid.i[1]))   
               TempReA2 = ReTwid.i[0] + 1;   
            else TempReA2 = ReTwid.i[0];   
           // Calculate new value for ReArray[indexB]    
            TempL.l = TempL.l << 1;   
            TempL.l -= ReTwid.l;   
            if ((TempL.l < 0)&&(TempL.i[1]))   
               TempReB2 = TempL.i[0] + 1;   
            else TempReB2 = TempL.i[0];   
   
            // Calculate new value for ImArray[indexA]    
            TempL.i[1] = 0;   
            TempL.i[0] = TempImA;   
            TempL.l = TempL.l >> 1;   
            ImTwid.l += TempL.l;   
            if ((ImTwid.l < 0)&&(ImTwid.i[1]))   
               TempImA = ImTwid.i[0] + 1;   
            else TempImA = ImTwid.i[0];   
   
            // Calculate new value for ImArray[indexB]    
            TempL.l = TempL.l << 1;   
            TempL.l -= ImTwid.l;   
            if ((TempL.l < 0)&&(TempL.i[1]))   
               TempImB = TempL.i[0] + 1;   
            else TempImB = TempL.i[0];   
   
         }   
   
         ReArray[indexA] = TempReA2;   
         ReArray[indexB] = TempReB2;   
         ImArray[indexA] = TempImA;   
         ImArray[indexB] = TempImB;   
   
         indexA++;   
         sin_index += group;   
      }                                      // END of stage FOR loop (s_cnt)    
      indexA = indexB + 1;   
      sin_index = 0;   
   }                                         // END of group FOR loop (g_cnt)    
   
   group /= 2;   
   stage *= 2;   
}                                            // END of While loop    
   
}  // END Int_FFT    

// Kyrre:
// Calculated the magnitude of each complex vector stored in Real[] and Imag[]
// from 0 to L and stores in result.
void fft_mag(uint16_t L, int16_t * result)
{
  if (L>255) {
    for (uint16_t n=0;n<L;n++) {
      result[n] = (int16_t)((int32_t)sqrt(Real[n]*Real[n] + (int32_t)Imag[n]*Imag[n]));
    }
  } else {
    for (uint8_t n=0;n<L;n++) {
      if (Real[n]!=0) {  // Is real bigger than 0?
        if (Imag[n]!=0) {  // Ok, Then do full calc if Imag bigger than zero too
          result[n] = (int16_t)(sqrt((int32_t)Real[n]*Real[n] + (int32_t)Imag[n]*Imag[n]));
        } else result[n] = abs(Real[n]); // No sence in calculating if Imag was zero
      } else result[n] = abs(Imag[n]);  // No sence in calculating if Real was zero
    }
  }
}



// Kyrre:
// Calculated the magnitude of each complex vector stored in Real[] and Imag[]
// from 0 to L and stores in result.
void fft_phase(uint16_t L, int16_t * result)
{
  uint8_t L8 = L;
  if (L>255) {
    for (uint16_t n=0;n<L;n++) {
      result[n] = (int16_t)((int32_t)sqrt(Real[n]*Real[n] + (int32_t)Imag[n]*Imag[n]));
    }
  } else {
    for (uint8_t n=0;n<L8;n++) {
      result[n] = (int16_t)((atan(((float)Imag[n])/Real[n])*180)/3.14);
    }
  }
}


// Kyrre:
// Calculated the magnitude of each complex vector stored in Real[] and Imag[]
// from 0 to L and stores in result.
void fft_magphase(uint16_t L, int16_t * mag_out, int16_t * phase_out)
{
  int8_t L8=L;
  int16_t mag;
  int16_t phase;
  if (L>255) {
    for (uint16_t n=0;n<L;n++) {
	  phase = (int16_t)((atan(((float)Imag[n])/Real[n])*180)/3.14);
      if (Real[n]!=0) {  // Is real bigger than 0?
        if (Imag[n]!=0) {  // Ok, Then do full calc if Imag bigger than zero too
          mag = (int16_t)(sqrt((int32_t)Real[n]*Real[n] + (int32_t)Imag[n]*Imag[n]));
        } else mag = abs(Real[n]); // No sence in calculating if Imag was zero
      } else mag = abs(Imag[n]);  // No sence in calculating if Real was zero
	  mag_out[n] = mag;
	  phase_out[n] = phase;
    }
  } else {
    for (uint8_t n=0;n<L8;n++) {
	  phase = (int16_t)((atan(((float)Imag[n])/Real[n])*180)/3.14);
      if (Real[n]!=0) {  // Is real bigger than 0?
        if (Imag[n]!=0) {  // Ok, Then do full calc if Imag bigger than zero too
          mag = (int16_t)(sqrt((int32_t)Real[n]*Real[n] + (int32_t)Imag[n]*Imag[n]));
        } else mag = abs(Real[n]); // No sence in calculating if Imag was zero
      } else mag = abs(Imag[n]);  // No sence in calculating if Real was zero
      mag_out[n] = mag;
	  phase_out[n] = phase;
    }
  }
}
