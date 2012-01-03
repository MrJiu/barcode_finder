/*
 * fft.h
 *
 * Created: 06.12.2011 17:47:46
 *  Author: kya
 */ 


#ifndef FFT_H_
#define FFT_H_

#ifndef SWAP
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#endif

void four1(double data[], unsigned long nn, int isign);
void realft(double data[], unsigned long n, int isign);
void twofft(double data1[], double data2[], double fft1[], double fft2[], unsigned long n);
void convlv(double data[], unsigned long n, double respns[], unsigned long m, int isign, double ans[]);

#endif /* FFT_H_ */
