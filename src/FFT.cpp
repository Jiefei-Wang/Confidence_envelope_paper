//The code is from :https://rosettacode.org/wiki/Fast_Fourier_transform#C.2B.2B with minor changes

#include <complex>
#include <iostream>
#include <valarray>
#include <Rcpp.h>

#define PI 3.14159265358979323846

#define complex_mul_real(x_real,x_img,y_real,y_img) (x_real*y_real - x_img*y_img)
#define complex_mul_img(x_real,x_img,y_real,y_img) (x_real*y_img + x_img*y_real)

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

void fft_swap(R_xlen_t N, double* x_real, double* x_img) {
	// Decimate
	unsigned int m = (unsigned int)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		unsigned int b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			//Complex t = x[a];
			//x[a] = x[b];
			//x[b] = t;
			std::swap(x_real[a], x_real[b]);
			std::swap(x_img[a], x_img[b]);
		}
	}
}
void fft(R_xlen_t N,double* x_real, double* x_img)
{
	// DFT
	R_xlen_t k = N, n;
	double thetaT = PI / N;
	Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				//Complex t = x[a] - x[b];
				Complex t = Complex(x_real[a] - x_real[b], x_img[a] - x_img[b]);
				//x[a] += x[b];
				//Rprintf("before:%u,%u,%u,%u,%f,%f,%f,%f\n", n, k, a, b, x_real[a], x_img[a], x_real[b], x_img[b]);
				x_real[a] += x_real[b];
				x_img[a] += x_img[b];
				//x[b] = t * T;
				Complex tmp = t * T;
				x_real[b] = tmp.real();
				x_img[b] = tmp.imag();
				//Rprintf("after:%u,%u,%u,%u,%f,%f,%f,%f\n", n, k, a, b, x_real[a], x_img[a], x_real[b], x_img[b]);
			}
			T *= phiT;
		}
	}
	fft_swap(N, x_real, x_img);
}



void fft_no_swap(R_xlen_t N, double* x_real, double* x_img)
{
	// DFT
	R_xlen_t k = N, n;
	double thetaT = PI / N;
	Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				//Complex t = x[a] - x[b];
				Complex t = Complex(x_real[a] - x_real[b], x_img[a] - x_img[b]);
				//x[a] += x[b];
				//Rprintf("before:%u,%u,%u,%u,%f,%f,%f,%f\n", n, k, a, b, x_real[a], x_img[a], x_real[b], x_img[b]);
				x_real[a] += x_real[b];
				x_img[a] += x_img[b];
				//x[b] = t * T;
				Complex tmp = t * T;
				x_real[b] = tmp.real();
				x_img[b] = tmp.imag();
				//Rprintf("after:%u,%u,%u,%u,%f,%f,%f,%f\n", n, k, a, b, x_real[a], x_img[a], x_real[b], x_img[b]);
			}
			T *= phiT;
		}
	}
}


void fft_negative_img(R_xlen_t N,double* x_real, double* x_img)
{
	// DFT
	R_xlen_t k = N, n;
	double thetaT = PI / N;
	Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				//Complex t = x[a] - x[b];
				Complex t = Complex(x_real[a] - x_real[b], -x_img[a] + x_img[b]);
				//x[a] += x[b];
				x_real[a] += x_real[b];
				x_img[a] += x_img[b];
				//x[b] = t * T;
				Complex tmp = t * T;
				x_real[b] = tmp.real();
				x_img[b] = -tmp.imag();
			}
			T *= phiT;
		}
	}

	fft_swap(N, x_real, x_img);
}

// inverse fft (in-place)
void ifft(R_xlen_t N, double* x_real, double* x_img)
{
	for (R_xlen_t i = 0; i < N; ++i) {
		x_img[i] = -x_img[i];
	}
    // conjugate the complex numbers
    //x = x.apply(std::conj);
    

    // forward fft
	fft(N,x_real, x_img);
    
    // conjugate the complex numbers again
    //x = x.apply(std::conj);
	for (R_xlen_t i = 0; i < N; ++i) {
		x_img[i] = -x_img[i]/N;
		x_real[i] /= N;
	}
    // scale the numbers
    //x /= x.size();
}

// inverse fft (in-place)
void ifft_no_img(R_xlen_t N, double* x_real, double* x_img)
{
	// forward fft
	fft_negative_img(N,x_real, x_img);

	// conjugate the complex numbers again
	//x = x.apply(std::conj);
	for (R_xlen_t i = 0; i < N; ++i) {
		x_real[i] /= N;
	}
	// scale the numbers
	//x /= x.size();
}




//Results will be in x
void convolution(R_xlen_t N, 
	double* x_real, double* x_img,
	double* y_real, double* y_img) {
	fft_no_swap(N, x_real, x_img);
	fft_no_swap(N, y_real, y_img);
	for (R_xlen_t i = 0; i < N; ++i) {
		double tmp = complex_mul_real(x_real[i], x_img[i], y_real[i], y_img[i]);
		x_img[i] = complex_mul_img(x_real[i], x_img[i], y_real[i], y_img[i]);
		x_real[i] = tmp;
	}
	fft_swap(N, x_real, x_img);
	ifft_no_img(N, x_real, x_img);
}
