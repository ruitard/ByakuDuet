// clang-format off

#include <cstdint>
#include <stdexcept>
#include <vector>
#include <memory>

#include "complex.hpp"
#include "autobuffer.hpp"

#define CV_SWAP(a,b,t) ((t) = (a), (a) = (b), (b) = (t))
#define CV_PI   3.1415926535897932384626433832795

namespace cv {

typedef std::uint32_t uint;
typedef unsigned char uchar;

static unsigned char bitrevTab[] =
{
  0x00,0x80,0x40,0xc0,0x20,0xa0,0x60,0xe0,0x10,0x90,0x50,0xd0,0x30,0xb0,0x70,0xf0,
  0x08,0x88,0x48,0xc8,0x28,0xa8,0x68,0xe8,0x18,0x98,0x58,0xd8,0x38,0xb8,0x78,0xf8,
  0x04,0x84,0x44,0xc4,0x24,0xa4,0x64,0xe4,0x14,0x94,0x54,0xd4,0x34,0xb4,0x74,0xf4,
  0x0c,0x8c,0x4c,0xcc,0x2c,0xac,0x6c,0xec,0x1c,0x9c,0x5c,0xdc,0x3c,0xbc,0x7c,0xfc,
  0x02,0x82,0x42,0xc2,0x22,0xa2,0x62,0xe2,0x12,0x92,0x52,0xd2,0x32,0xb2,0x72,0xf2,
  0x0a,0x8a,0x4a,0xca,0x2a,0xaa,0x6a,0xea,0x1a,0x9a,0x5a,0xda,0x3a,0xba,0x7a,0xfa,
  0x06,0x86,0x46,0xc6,0x26,0xa6,0x66,0xe6,0x16,0x96,0x56,0xd6,0x36,0xb6,0x76,0xf6,
  0x0e,0x8e,0x4e,0xce,0x2e,0xae,0x6e,0xee,0x1e,0x9e,0x5e,0xde,0x3e,0xbe,0x7e,0xfe,
  0x01,0x81,0x41,0xc1,0x21,0xa1,0x61,0xe1,0x11,0x91,0x51,0xd1,0x31,0xb1,0x71,0xf1,
  0x09,0x89,0x49,0xc9,0x29,0xa9,0x69,0xe9,0x19,0x99,0x59,0xd9,0x39,0xb9,0x79,0xf9,
  0x05,0x85,0x45,0xc5,0x25,0xa5,0x65,0xe5,0x15,0x95,0x55,0xd5,0x35,0xb5,0x75,0xf5,
  0x0d,0x8d,0x4d,0xcd,0x2d,0xad,0x6d,0xed,0x1d,0x9d,0x5d,0xdd,0x3d,0xbd,0x7d,0xfd,
  0x03,0x83,0x43,0xc3,0x23,0xa3,0x63,0xe3,0x13,0x93,0x53,0xd3,0x33,0xb3,0x73,0xf3,
  0x0b,0x8b,0x4b,0xcb,0x2b,0xab,0x6b,0xeb,0x1b,0x9b,0x5b,0xdb,0x3b,0xbb,0x7b,0xfb,
  0x07,0x87,0x47,0xc7,0x27,0xa7,0x67,0xe7,0x17,0x97,0x57,0xd7,0x37,0xb7,0x77,0xf7,
  0x0f,0x8f,0x4f,0xcf,0x2f,0xaf,0x6f,0xef,0x1f,0x9f,0x5f,0xdf,0x3f,0xbf,0x7f,0xff
};

static const double DFTTab[][2] =
{
{ 1.00000000000000000, 0.00000000000000000 },
{-1.00000000000000000, 0.00000000000000000 },
{ 0.00000000000000000, 1.00000000000000000 },
{ 0.70710678118654757, 0.70710678118654746 },
{ 0.92387953251128674, 0.38268343236508978 },
{ 0.98078528040323043, 0.19509032201612825 },
{ 0.99518472667219693, 0.09801714032956060 },
{ 0.99879545620517241, 0.04906767432741802 },
{ 0.99969881869620425, 0.02454122852291229 },
{ 0.99992470183914450, 0.01227153828571993 },
{ 0.99998117528260111, 0.00613588464915448 },
{ 0.99999529380957619, 0.00306795676296598 },
{ 0.99999882345170188, 0.00153398018628477 },
{ 0.99999970586288223, 0.00076699031874270 },
{ 0.99999992646571789, 0.00038349518757140 },
{ 0.99999998161642933, 0.00019174759731070 },
{ 0.99999999540410733, 0.00009587379909598 },
{ 0.99999999885102686, 0.00004793689960307 },
{ 0.99999999971275666, 0.00002396844980842 },
{ 0.99999999992818922, 0.00001198422490507 },
{ 0.99999999998204725, 0.00000599211245264 },
{ 0.99999999999551181, 0.00000299605622633 },
{ 0.99999999999887801, 0.00000149802811317 },
{ 0.99999999999971945, 0.00000074901405658 },
{ 0.99999999999992983, 0.00000037450702829 },
{ 0.99999999999998246, 0.00000018725351415 },
{ 0.99999999999999567, 0.00000009362675707 },
{ 0.99999999999999889, 0.00000004681337854 },
{ 0.99999999999999978, 0.00000002340668927 },
{ 0.99999999999999989, 0.00000001170334463 },
{ 1.00000000000000000, 0.00000000585167232 },
{ 1.00000000000000000, 0.00000000292583616 }
};

namespace {
template <typename T>
struct Constants {
    static const T sin_120;
    static const T fft5_2;
    static const T fft5_3;
    static const T fft5_4;
    static const T fft5_5;
};

template <typename T>
const T Constants<T>::sin_120 = (T)0.86602540378443864676372317075294;

template <typename T>
const T Constants<T>::fft5_2 = (T)0.559016994374947424102293417182819;

template <typename T>
const T Constants<T>::fft5_3 = (T)-0.951056516295153572116439333379382;

template <typename T>
const T Constants<T>::fft5_4 = (T)-1.538841768587626701285145288018455;

template <typename T>
const T Constants<T>::fft5_5 = (T)0.363271264002680442947733378740309;

}  //namespace

struct OcvDftOptions;

typedef void (*DFTFunc)(const OcvDftOptions & c, const void* src, void* dst);

struct OcvDftOptions {
    int nf;
    int *factors;
    double scale;

    int* itab;
    void* wave;
    int tab_size;
    int n;

    bool isInverse;
    bool noPermute;
    bool isComplex;

    bool haveSSE3;

    DFTFunc dft_func;
    bool useIpp;

    OcvDftOptions()
    {
        nf = 0;
        factors = 0;
        scale = 0;
        itab = 0;
        wave = 0;
        tab_size = 0;
        n = 0;
        isInverse = false;
        noPermute = false;
        isComplex = false;
        useIpp = false;

        dft_func = 0;
        // haveSSE3 = checkHardwareSupport(CV_CPU_SSE3);
    }
};

// Reference radix-2 implementation.
template<typename T> struct DFT_R2
{
    void operator()(Complex<T>* dst, const int c_n, const int n, const int dw0, const Complex<T>* wave) const {
        const int nx = n/2;
        for(int i = 0 ; i < c_n; i += n)
        {
            Complex<T>* v = dst + i;
            T r0 = v[0].re + v[nx].re;
            T i0 = v[0].im + v[nx].im;
            T r1 = v[0].re - v[nx].re;
            T i1 = v[0].im - v[nx].im;
            v[0].re = r0; v[0].im = i0;
            v[nx].re = r1; v[nx].im = i1;

            for( int j = 1, dw = dw0; j < nx; j++, dw += dw0 )
            {
                v = dst + i + j;
                r1 = v[nx].re*wave[dw].re - v[nx].im*wave[dw].im;
                i1 = v[nx].im*wave[dw].re + v[nx].re*wave[dw].im;
                r0 = v[0].re; i0 = v[0].im;

                v[0].re = r0 + r1; v[0].im = i0 + i1;
                v[nx].re = r0 - r1; v[nx].im = i0 - i1;
            }
        }
    }
};

// Reference radix-3 implementation.
template<typename T> struct DFT_R3
{
    void operator()(Complex<T>* dst, const int c_n, const int n, const int dw0, const Complex<T>* wave) const {
        const int nx = n / 3;
        for(int i = 0; i < c_n; i += n )
        {
            {
                Complex<T>* v = dst + i;
                T r1 = v[nx].re + v[nx*2].re;
                T i1 = v[nx].im + v[nx*2].im;
                T r0 = v[0].re;
                T i0 = v[0].im;
                T r2 = Constants<T>::sin_120*(v[nx].im - v[nx*2].im);
                T i2 = Constants<T>::sin_120*(v[nx*2].re - v[nx].re);
                v[0].re = r0 + r1; v[0].im = i0 + i1;
                r0 -= (T)0.5*r1; i0 -= (T)0.5*i1;
                v[nx].re = r0 + r2; v[nx].im = i0 + i2;
                v[nx*2].re = r0 - r2; v[nx*2].im = i0 - i2;
            }

            for(int j = 1, dw = dw0; j < nx; j++, dw += dw0 )
            {
                Complex<T>* v = dst + i + j;
                T r0 = v[nx].re*wave[dw].re - v[nx].im*wave[dw].im;
                T i0 = v[nx].re*wave[dw].im + v[nx].im*wave[dw].re;
                T i2 = v[nx*2].re*wave[dw*2].re - v[nx*2].im*wave[dw*2].im;
                T r2 = v[nx*2].re*wave[dw*2].im + v[nx*2].im*wave[dw*2].re;
                T r1 = r0 + i2; T i1 = i0 + r2;

                r2 = Constants<T>::sin_120*(i0 - r2); i2 = Constants<T>::sin_120*(i2 - r0);
                r0 = v[0].re; i0 = v[0].im;
                v[0].re = r0 + r1; v[0].im = i0 + i1;
                r0 -= (T)0.5*r1; i0 -= (T)0.5*i1;
                v[nx].re = r0 + r2; v[nx].im = i0 + i2;
                v[nx*2].re = r0 - r2; v[nx*2].im = i0 - i2;
            }
        }
    }
};

// Reference radix-5 implementation.
template<typename T> struct DFT_R5
{
    void operator()(Complex<T>* dst, const int c_n, const int n, const int dw0, const Complex<T>* wave) const {
        const int nx = n / 5;
        for(int i = 0; i < c_n; i += n )
        {
            for(int j = 0, dw = 0; j < nx; j++, dw += dw0 )
            {
                Complex<T>* v0 = dst + i + j;
                Complex<T>* v1 = v0 + nx*2;
                Complex<T>* v2 = v1 + nx*2;

                T r0, i0, r1, i1, r2, i2, r3, i3, r4, i4, r5, i5;

                r3 = v0[nx].re*wave[dw].re - v0[nx].im*wave[dw].im;
                i3 = v0[nx].re*wave[dw].im + v0[nx].im*wave[dw].re;
                r2 = v2[0].re*wave[dw*4].re - v2[0].im*wave[dw*4].im;
                i2 = v2[0].re*wave[dw*4].im + v2[0].im*wave[dw*4].re;

                r1 = r3 + r2; i1 = i3 + i2;
                r3 -= r2; i3 -= i2;

                r4 = v1[nx].re*wave[dw*3].re - v1[nx].im*wave[dw*3].im;
                i4 = v1[nx].re*wave[dw*3].im + v1[nx].im*wave[dw*3].re;
                r0 = v1[0].re*wave[dw*2].re - v1[0].im*wave[dw*2].im;
                i0 = v1[0].re*wave[dw*2].im + v1[0].im*wave[dw*2].re;

                r2 = r4 + r0; i2 = i4 + i0;
                r4 -= r0; i4 -= i0;

                r0 = v0[0].re; i0 = v0[0].im;
                r5 = r1 + r2; i5 = i1 + i2;

                v0[0].re = r0 + r5; v0[0].im = i0 + i5;

                r0 -= (T)0.25*r5; i0 -= (T)0.25*i5;
                r1 = Constants<T>::fft5_2*(r1 - r2); i1 = Constants<T>::fft5_2*(i1 - i2);
                r2 = -Constants<T>::fft5_3*(i3 + i4); i2 = Constants<T>::fft5_3*(r3 + r4);

                i3 *= -Constants<T>::fft5_5; r3 *= Constants<T>::fft5_5;
                i4 *= -Constants<T>::fft5_4; r4 *= Constants<T>::fft5_4;

                r5 = r2 + i3; i5 = i2 + r3;
                r2 -= i4; i2 -= r4;

                r3 = r0 + r1; i3 = i0 + i1;
                r0 -= r1; i0 -= i1;

                v0[nx].re = r3 + r2; v0[nx].im = i3 + i2;
                v2[0].re = r3 - r2; v2[0].im = i3 - i2;

                v1[0].re = r0 + r5; v1[0].im = i0 + i5;
                v1[nx].re = r0 - r5; v1[nx].im = i0 - i5;
            }
        }
    }
};

template<typename T> struct DFT_VecR2
{
    void operator()(Complex<T>* dst, const int c_n, const int n, const int dw0, const Complex<T>* wave) const {
        DFT_R2<T>()(dst, c_n, n, dw0, wave);
    }
};

template<typename T> struct DFT_VecR3
{
    void operator()(Complex<T>* dst, const int c_n, const int n, const int dw0, const Complex<T>* wave) const {
        DFT_R3<T>()(dst, c_n, n, dw0, wave);
    }
};

template<typename T> struct DFT_VecR4
{
    int operator()(Complex<T>*, int, int, int&, const Complex<T>*) const { return 1; }
};

#define BitRev(i,shift) \
   ((int)((((unsigned)bitrevTab[(i)&255] << 24)+ \
           ((unsigned)bitrevTab[((i)>> 8)&255] << 16)+ \
           ((unsigned)bitrevTab[((i)>>16)&255] <<  8)+ \
           ((unsigned)bitrevTab[((i)>>24)])) >> (shift)))

static int
DFTFactorize( int n, int* factors )
{
    int nf = 0, f, i, j;

    if( n <= 5 )
    {
        factors[0] = n;
        return 1;
    }

    f = (((n - 1)^n)+1) >> 1;
    if( f > 1 )
    {
        factors[nf++] = f;
        n = f == n ? 1 : n/f;
    }

    for( f = 3; n > 1; )
    {
        int d = n/f;
        if( d*f == n )
        {
            factors[nf++] = f;
            n = d;
        }
        else
        {
            f += 2;
            if( f*f > n )
                break;
        }
    }

    if( n > 1 )
        factors[nf++] = n;

    f = (factors[0] & 1) == 0;
    for( i = f; i < (nf+f)/2; i++ )
        CV_SWAP( factors[i], factors[nf-i-1+f], j );

    return nf;
}

static void
DFTInit( int n0, int nf, const int* factors, int* itab, int elem_size, void* _wave, int inv_itab )
{
    int digits[34], radix[34];
    int n = factors[0], m = 0;
    int* itab0 = itab;
    int i, j, k;
    Complex<double> w, w1;
    double t;

    if( n0 <= 5 )
    {
        itab[0] = 0;
        itab[n0-1] = n0-1;

        if( n0 != 4 )
        {
            for( i = 1; i < n0-1; i++ )
                itab[i] = i;
        }
        else
        {
            itab[1] = 2;
            itab[2] = 1;
        }
        if( n0 == 5 )
        {
            if( elem_size == sizeof(Complex<double>) )
                ((Complex<double>*)_wave)[0] = Complex<double>(1.,0.);
            else
                ((Complex<float>*)_wave)[0] = Complex<float>(1.f,0.f);
        }
        if( n0 != 4 )
            return;
        m = 2;
    }
    else
    {
        // radix[] is initialized from index 'nf' down to zero
        radix[nf] = 1;
        digits[nf] = 0;
        for( i = 0; i < nf; i++ )
        {
            digits[i] = 0;
            radix[nf-i-1] = radix[nf-i]*factors[nf-i-1];
        }

        if( inv_itab && factors[0] != factors[nf-1] )
            itab = (int*)_wave;

        if( (n & 1) == 0 )
        {
            int a = radix[1], na2 = n*a>>1, na4 = na2 >> 1;
            for( m = 0; (unsigned)(1 << m) < (unsigned)n; m++ )
                ;
            if( n <= 2  )
            {
                itab[0] = 0;
                itab[1] = na2;
            }
            else if( n <= 256 )
            {
                int shift = 10 - m;
                for( i = 0; i <= n - 4; i += 4 )
                {
                    j = (bitrevTab[i>>2]>>shift)*a;
                    itab[i] = j;
                    itab[i+1] = j + na2;
                    itab[i+2] = j + na4;
                    itab[i+3] = j + na2 + na4;
                }
            }
            else
            {
                int shift = 34 - m;
                for( i = 0; i < n; i += 4 )
                {
                    int i4 = i >> 2;
                    j = BitRev(i4,shift)*a;
                    itab[i] = j;
                    itab[i+1] = j + na2;
                    itab[i+2] = j + na4;
                    itab[i+3] = j + na2 + na4;
                }
            }

            digits[1]++;

            if( nf >= 2 )
            {
                for( i = n, j = radix[2]; i < n0; )
                {
                    for( k = 0; k < n; k++ )
                        itab[i+k] = itab[k] + j;
                    if( (i += n) >= n0 )
                        break;
                    j += radix[2];
                    for( k = 1; ++digits[k] >= factors[k]; k++ )
                    {
                        digits[k] = 0;
                        j += radix[k+2] - radix[k];
                    }
                }
            }
        }
        else
        {
            for( i = 0, j = 0;; )
            {
                itab[i] = j;
                if( ++i >= n0 )
                    break;
                j += radix[1];
                for( k = 0; ++digits[k] >= factors[k]; k++ )
                {
                    digits[k] = 0;
                    j += radix[k+2] - radix[k];
                }
            }
        }

        if( itab != itab0 )
        {
            itab0[0] = 0;
            for( i = n0 & 1; i < n0; i += 2 )
            {
                int k0 = itab[i];
                int k1 = itab[i+1];
                itab0[k0] = i;
                itab0[k1] = i+1;
            }
        }
    }

    if( (n0 & (n0-1)) == 0 )
    {
        w.re = w1.re = DFTTab[m][0];
        w.im = w1.im = -DFTTab[m][1];
    }
    else
    {
        t = -CV_PI*2/n0;
        w.im = w1.im = sin(t);
        w.re = w1.re = std::sqrt(1. - w1.im*w1.im);
    }
    n = (n0+1)/2;

    if( elem_size == sizeof(Complex<double>) )
    {
        Complex<double>* wave = (Complex<double>*)_wave;

        wave[0].re = 1.;
        wave[0].im = 0.;

        if( (n0 & 1) == 0 )
        {
            wave[n].re = -1.;
            wave[n].im = 0;
        }

        for( i = 1; i < n; i++ )
        {
            wave[i] = w;
            wave[n0-i].re = w.re;
            wave[n0-i].im = -w.im;

            t = w.re*w1.re - w.im*w1.im;
            w.im = w.re*w1.im + w.im*w1.re;
            w.re = t;
        }
    }
    else
    {
        Complex<float>* wave = (Complex<float>*)_wave;

        wave[0].re = 1.f;
        wave[0].im = 0.f;

        if( (n0 & 1) == 0 )
        {
            wave[n].re = -1.f;
            wave[n].im = 0.f;
        }

        for( i = 1; i < n; i++ )
        {
            wave[i].re = (float)w.re;
            wave[i].im = (float)w.im;
            wave[n0-i].re = (float)w.re;
            wave[n0-i].im = (float)-w.im;

            t = w.re*w1.re - w.im*w1.im;
            w.im = w.re*w1.im + w.im*w1.re;
            w.re = t;
        }
    }
}

// mixed-radix complex discrete Fourier transform: double-precision version
template<typename T> static void
DFT(const OcvDftOptions & c, const Complex<T>* src, Complex<T>* dst)
{
    const Complex<T>* wave = (Complex<T>*)c.wave;
    const int * itab = c.itab;

    int n = c.n;
    int f_idx, nx;
    int inv = c.isInverse;
    int dw0 = c.tab_size, dw;
    int i, j, k;
    Complex<T> t;
    T scale = (T)c.scale;

    int tab_step = c.tab_size == n ? 1 : c.tab_size == n*2 ? 2 : c.tab_size/n;

    // 0. shuffle data
    if( dst != src )
    {
        if( !inv )
        {
            for( i = 0; i <= n - 2; i += 2, itab += 2*tab_step )
            {
                int k0 = itab[0], k1 = itab[tab_step];
                dst[i] = src[k0]; dst[i+1] = src[k1];
            }

            if( i < n )
                dst[n-1] = src[n-1];
        }
        else
        {
            for( i = 0; i <= n - 2; i += 2, itab += 2*tab_step )
            {
                int k0 = itab[0], k1 = itab[tab_step];
                t.re = src[k0].re; t.im = -src[k0].im;
                dst[i] = t;
                t.re = src[k1].re; t.im = -src[k1].im;
                dst[i+1] = t;
            }

            if( i < n )
            {
                t.re = src[n-1].re; t.im = -src[n-1].im;
                dst[i] = t;
            }
        }
    }
    else
    {
        if( !c.noPermute )
        {
            if( c.nf == 1 )
            {
                if( (n & 3) == 0 )
                {
                    int n2 = n/2;
                    Complex<T>* dsth = dst + n2;

                    for( i = 0; i < n2; i += 2, itab += tab_step*2 )
                    {
                        j = itab[0];

                        CV_SWAP(dst[i+1], dsth[j], t);
                        if( j > i )
                        {
                            CV_SWAP(dst[i], dst[j], t);
                            CV_SWAP(dsth[i+1], dsth[j+1], t);
                        }
                    }
                }
                // else do nothing
            }
            else
            {
                for( i = 0; i < n; i++, itab += tab_step )
                {
                    j = itab[0];
                    if( j > i )
                        CV_SWAP(dst[i], dst[j], t);
                }
            }
        }

        if( inv )
        {
            for( i = 0; i <= n - 2; i += 2 )
            {
                T t0 = -dst[i].im;
                T t1 = -dst[i+1].im;
                dst[i].im = t0; dst[i+1].im = t1;
            }

            if( i < n )
                dst[n-1].im = -dst[n-1].im;
        }
    }

    n = 1;
    // 1. power-2 transforms
    if( (c.factors[0] & 1) == 0 )
    {
        if( c.factors[0] >= 4 && c.haveSSE3)
        {
            DFT_VecR4<T> vr4;
            n = vr4(dst, c.factors[0], c.n, dw0, wave);
        }

        // radix-4 transform
        for( ; n*4 <= c.factors[0]; )
        {
            nx = n;
            n *= 4;
            dw0 /= 4;

            for( i = 0; i < c.n; i += n )
            {
                Complex<T> *v0, *v1;
                T r0, i0, r1, i1, r2, i2, r3, i3, r4, i4;

                v0 = dst + i;
                v1 = v0 + nx*2;

                r0 = v1[0].re; i0 = v1[0].im;
                r4 = v1[nx].re; i4 = v1[nx].im;

                r1 = r0 + r4; i1 = i0 + i4;
                r3 = i0 - i4; i3 = r4 - r0;

                r2 = v0[0].re; i2 = v0[0].im;
                r4 = v0[nx].re; i4 = v0[nx].im;

                r0 = r2 + r4; i0 = i2 + i4;
                r2 -= r4; i2 -= i4;

                v0[0].re = r0 + r1; v0[0].im = i0 + i1;
                v1[0].re = r0 - r1; v1[0].im = i0 - i1;
                v0[nx].re = r2 + r3; v0[nx].im = i2 + i3;
                v1[nx].re = r2 - r3; v1[nx].im = i2 - i3;

                for( j = 1, dw = dw0; j < nx; j++, dw += dw0 )
                {
                    v0 = dst + i + j;
                    v1 = v0 + nx*2;

                    r2 = v0[nx].re*wave[dw*2].re - v0[nx].im*wave[dw*2].im;
                    i2 = v0[nx].re*wave[dw*2].im + v0[nx].im*wave[dw*2].re;
                    r0 = v1[0].re*wave[dw].im + v1[0].im*wave[dw].re;
                    i0 = v1[0].re*wave[dw].re - v1[0].im*wave[dw].im;
                    r3 = v1[nx].re*wave[dw*3].im + v1[nx].im*wave[dw*3].re;
                    i3 = v1[nx].re*wave[dw*3].re - v1[nx].im*wave[dw*3].im;

                    r1 = i0 + i3; i1 = r0 + r3;
                    r3 = r0 - r3; i3 = i3 - i0;
                    r4 = v0[0].re; i4 = v0[0].im;

                    r0 = r4 + r2; i0 = i4 + i2;
                    r2 = r4 - r2; i2 = i4 - i2;

                    v0[0].re = r0 + r1; v0[0].im = i0 + i1;
                    v1[0].re = r0 - r1; v1[0].im = i0 - i1;
                    v0[nx].re = r2 + r3; v0[nx].im = i2 + i3;
                    v1[nx].re = r2 - r3; v1[nx].im = i2 - i3;
                }
            }
        }

        for( ; n < c.factors[0]; )
        {
            // do the remaining radix-2 transform
            n *= 2;
            dw0 /= 2;

            if(c.haveSSE3)
            {
                DFT_VecR2<T> vr2;
                vr2(dst, c.n, n, dw0, wave);
            }
            else
            {
                DFT_R2<T> vr2;
                vr2(dst, c.n, n, dw0, wave);
            }
        }
    }

    // 2. all the other transforms
    for( f_idx = (c.factors[0]&1) ? 0 : 1; f_idx < c.nf; f_idx++ )
    {
        int factor = c.factors[f_idx];
        nx = n;
        n *= factor;
        dw0 /= factor;

        if( factor == 3 )
        {
            if(c.haveSSE3)
            {
                DFT_VecR3<T> vr3;
                vr3(dst, c.n, n, dw0, wave);
            }
            else
            {
                DFT_R3<T> vr3;
                vr3(dst, c.n, n, dw0, wave);
            }
        }
        else if( factor == 5 )
        {
            DFT_R5<T> vr5;
            vr5(dst, c.n, n, dw0, wave);
        }
        else
        {
            // radix-"factor" - an odd number
            int p, q, factor2 = (factor - 1)/2;
            int d, dd, dw_f = c.tab_size/factor;
            AutoBuffer<Complex<T> > buf(factor2 * 2);
            Complex<T>* a = buf.data();
            Complex<T>* b = a + factor2;

            for( i = 0; i < c.n; i += n )
            {
                for( j = 0, dw = 0; j < nx; j++, dw += dw0 )
                {
                    Complex<T>* v = dst + i + j;
                    Complex<T> v_0 = v[0];
                    Complex<T> vn_0 = v_0;

                    if( j == 0 )
                    {
                        for( p = 1, k = nx; p <= factor2; p++, k += nx )
                        {
                            T r0 = v[k].re + v[n-k].re;
                            T i0 = v[k].im - v[n-k].im;
                            T r1 = v[k].re - v[n-k].re;
                            T i1 = v[k].im + v[n-k].im;

                            vn_0.re += r0; vn_0.im += i1;
                            a[p-1].re = r0; a[p-1].im = i0;
                            b[p-1].re = r1; b[p-1].im = i1;
                        }
                    }
                    else
                    {
                        const Complex<T>* wave_ = wave + dw*factor;
                        d = dw;

                        for( p = 1, k = nx; p <= factor2; p++, k += nx, d += dw )
                        {
                            T r2 = v[k].re*wave[d].re - v[k].im*wave[d].im;
                            T i2 = v[k].re*wave[d].im + v[k].im*wave[d].re;

                            T r1 = v[n-k].re*wave_[-d].re - v[n-k].im*wave_[-d].im;
                            T i1 = v[n-k].re*wave_[-d].im + v[n-k].im*wave_[-d].re;

                            T r0 = r2 + r1;
                            T i0 = i2 - i1;
                            r1 = r2 - r1;
                            i1 = i2 + i1;

                            vn_0.re += r0; vn_0.im += i1;
                            a[p-1].re = r0; a[p-1].im = i0;
                            b[p-1].re = r1; b[p-1].im = i1;
                        }
                    }

                    v[0] = vn_0;

                    for( p = 1, k = nx; p <= factor2; p++, k += nx )
                    {
                        Complex<T> s0 = v_0, s1 = v_0;
                        d = dd = dw_f*p;

                        for( q = 0; q < factor2; q++ )
                        {
                            T r0 = wave[d].re * a[q].re;
                            T i0 = wave[d].im * a[q].im;
                            T r1 = wave[d].re * b[q].im;
                            T i1 = wave[d].im * b[q].re;

                            s1.re += r0 + i0; s0.re += r0 - i0;
                            s1.im += r1 - i1; s0.im += r1 + i1;

                            d += dd;
                            d -= -(d >= c.tab_size) & c.tab_size;
                        }

                        v[k] = s0;
                        v[n-k] = s1;
                    }
                }
            }
        }
    }

    if( scale != 1 )
    {
        T re_scale = scale, im_scale = scale;
        if( inv )
            im_scale = -im_scale;

        for( i = 0; i < c.n; i++ )
        {
            T t0 = dst[i].re*re_scale;
            T t1 = dst[i].im*im_scale;
            dst[i].re = t0;
            dst[i].im = t1;
        }
    }
    else if( inv )
    {
        for( i = 0; i <= c.n - 2; i += 2 )
        {
            T t0 = -dst[i].im;
            T t1 = -dst[i+1].im;
            dst[i].im = t0;
            dst[i+1].im = t1;
        }

        if( i < c.n )
            dst[c.n-1].im = -dst[c.n-1].im;
    }
}

/* FFT of real vector
   output vector format:
     re(0), re(1), im(1), ... , re(n/2-1), im((n+1)/2-1) [, re((n+1)/2)] OR ...
     re(0), 0, re(1), im(1), ..., re(n/2-1), im((n+1)/2-1) [, re((n+1)/2), 0] */
template<typename T> static void
RealDFT(const OcvDftOptions & c, const T* src, T* dst)
{
    int n = c.n;
    int complex_output = c.isComplex;
    T scale = (T)c.scale;
    int j;
    dst += complex_output;

    if( n == 1 )
    {
        dst[0] = src[0]*scale;
    }
    else if( n == 2 )
    {
        T t = (src[0] + src[1])*scale;
        dst[1] = (src[0] - src[1])*scale;
        dst[0] = t;
    }
    else if( n & 1 )
    {
        dst -= complex_output;
        Complex<T>* _dst = (Complex<T>*)dst;
        _dst[0].re = src[0]*scale;
        _dst[0].im = 0;
        for( j = 1; j < n; j += 2 )
        {
            T t0 = src[c.itab[j]]*scale;
            T t1 = src[c.itab[j+1]]*scale;
            _dst[j].re = t0;
            _dst[j].im = 0;
            _dst[j+1].re = t1;
            _dst[j+1].im = 0;
        }
        OcvDftOptions sub_c = c;
        sub_c.isComplex = false;
        sub_c.isInverse = false;
        sub_c.noPermute = true;
        sub_c.scale = 1.;
        DFT(sub_c, _dst, _dst);
        if( !complex_output )
            dst[1] = dst[0];
    }
    else
    {
        T t0, t;
        T h1_re, h1_im, h2_re, h2_im;
        T scale2 = scale*(T)0.5;
        int n2 = n >> 1;

        c.factors[0] >>= 1;

        OcvDftOptions sub_c = c;
        sub_c.factors += (c.factors[0] == 1);
        sub_c.nf -= (c.factors[0] == 1);
        sub_c.isComplex = false;
        sub_c.isInverse = false;
        sub_c.noPermute = false;
        sub_c.scale = 1.;
        sub_c.n = n2;

        DFT(sub_c, (Complex<T>*)src, (Complex<T>*)dst);

        c.factors[0] <<= 1;

        t = dst[0] - dst[1];
        dst[0] = (dst[0] + dst[1])*scale;
        dst[1] = t*scale;

        t0 = dst[n2];
        t = dst[n-1];
        dst[n-1] = dst[1];

        const Complex<T> *wave = (const Complex<T>*)c.wave;

        for( j = 2, wave++; j < n2; j += 2, wave++ )
        {
            /* calc odd */
            h2_re = scale2*(dst[j+1] + t);
            h2_im = scale2*(dst[n-j] - dst[j]);

            /* calc even */
            h1_re = scale2*(dst[j] + dst[n-j]);
            h1_im = scale2*(dst[j+1] - t);

            /* rotate */
            t = h2_re*wave->re - h2_im*wave->im;
            h2_im = h2_re*wave->im + h2_im*wave->re;
            h2_re = t;
            t = dst[n-j-1];

            dst[j-1] = h1_re + h2_re;
            dst[n-j-1] = h1_re - h2_re;
            dst[j] = h1_im + h2_im;
            dst[n-j] = h2_im - h1_im;
        }

        if( j <= n2 )
        {
            dst[n2-1] = t0*scale;
            dst[n2] = -t*scale;
        }
    }

    if( complex_output && ((n & 1) == 0 || n == 1))
    {
        dst[-1] = dst[0];
        dst[0] = 0;
        if( n > 1 )
            dst[n] = 0;
    }
}

} // namespace cv

namespace cv
{

/* DCT is calculated using DFT, as described here:
   http://www.ece.utexas.edu/~bevans/courses/ee381k/lectures/09_DCT/lecture9/:
*/
template<typename T> static void
DCT( const OcvDftOptions & c, const T* src, size_t src_step, T* dft_src, T* dft_dst, T* dst, size_t dst_step,
     const Complex<T>* dct_wave )
{
    static const T sin_45 = (T)0.70710678118654752440084436210485;

    int n = c.n;
    int j, n2 = n >> 1;

    src_step /= sizeof(src[0]);
    dst_step /= sizeof(dst[0]);
    T* dst1 = dst + (n-1)*dst_step;

    if( n == 1 )
    {
        dst[0] = src[0];
        return;
    }

    for( j = 0; j < n2; j++, src += src_step*2 )
    {
        dft_src[j] = src[0];
        dft_src[n-j-1] = src[src_step];
    }

    RealDFT(c, dft_src, dft_dst);
    src = dft_dst;

    dst[0] = (T)(src[0]*dct_wave->re*sin_45);
    dst += dst_step;
    for( j = 1, dct_wave++; j < n2; j++, dct_wave++,
                                    dst += dst_step, dst1 -= dst_step )
    {
        T t0 = dct_wave->re*src[j*2-1] - dct_wave->im*src[j*2];
        T t1 = -dct_wave->im*src[j*2-1] - dct_wave->re*src[j*2];
        dst[0] = t0;
        dst1[0] = t1;
    }

    dst[0] = src[n-1]*dct_wave->re;
}

static void
DCTInit( int n, int elem_size, void* _wave, int inv )
{
    static const double DctScale[] =
    {
    0.707106781186547570, 0.500000000000000000, 0.353553390593273790,
    0.250000000000000000, 0.176776695296636890, 0.125000000000000000,
    0.088388347648318447, 0.062500000000000000, 0.044194173824159223,
    0.031250000000000000, 0.022097086912079612, 0.015625000000000000,
    0.011048543456039806, 0.007812500000000000, 0.005524271728019903,
    0.003906250000000000, 0.002762135864009952, 0.001953125000000000,
    0.001381067932004976, 0.000976562500000000, 0.000690533966002488,
    0.000488281250000000, 0.000345266983001244, 0.000244140625000000,
    0.000172633491500622, 0.000122070312500000, 0.000086316745750311,
    0.000061035156250000, 0.000043158372875155, 0.000030517578125000
    };

    int i;
    Complex<double> w, w1;
    double t, scale;

    if( n == 1 )
        return;

    if( (n & (n - 1)) == 0 )
    {
        int m;
        for( m = 0; (unsigned)(1 << m) < (unsigned)n; m++ )
            ;
        scale = (!inv ? 2 : 1)*DctScale[m];
        w1.re = DFTTab[m+2][0];
        w1.im = -DFTTab[m+2][1];
    }
    else
    {
        t = 1./(2*n);
        scale = (!inv ? 2 : 1)*std::sqrt(t);
        w1.im = sin(-CV_PI*t);
        w1.re = std::sqrt(1. - w1.im*w1.im);
    }
    n >>= 1;

    if( elem_size == sizeof(Complex<double>) )
    {
        Complex<double>* wave = (Complex<double>*)_wave;

        w.re = scale;
        w.im = 0.;

        for( i = 0; i <= n; i++ )
        {
            wave[i] = w;
            t = w.re*w1.re - w.im*w1.im;
            w.im = w.re*w1.im + w.im*w1.re;
            w.re = t;
        }
    }
    else
    {
        Complex<float>* wave = (Complex<float>*)_wave;

        w.re = (float)scale;
        w.im = 0.f;

        for( i = 0; i <= n; i++ )
        {
            wave[i].re = (float)w.re;
            wave[i].im = (float)w.im;
            t = w.re*w1.re - w.im*w1.im;
            w.im = w.re*w1.im + w.im*w1.re;
            w.re = t;
        }
    }
}

typedef void (*DCTFunc)(const OcvDftOptions & c, const void* src, size_t src_step, void* dft_src,
                        void* dft_dst, void* dst, size_t dst_step, const void* dct_wave);

static void DCT_32f(const OcvDftOptions & c, const float* src, size_t src_step, float* dft_src, float* dft_dst,
                    float* dst, size_t dst_step, const Complexf* dct_wave)
{
    DCT(c, src, src_step, dft_src, dft_dst, dst, dst_step, dct_wave);
}

static void DCT_64f(const OcvDftOptions & c, const double* src, size_t src_step, double* dft_src, double* dft_dst,
                    double* dst, size_t dst_step, const Complexd* dct_wave)
{
    DCT(c, src, src_step, dft_src, dft_dst, dst, dst_step, dct_wave);
}



//! @name DFT flags
//! @{
#define CV_HAL_DFT_INVERSE        1
#define CV_HAL_DFT_SCALE          2
#define CV_HAL_DFT_ROWS           4
#define CV_HAL_DFT_COMPLEX_OUTPUT 16
#define CV_HAL_DFT_REAL_OUTPUT    32
#define CV_HAL_DFT_TWO_STAGE      64
#define CV_HAL_DFT_STAGE_COLS    128
#define CV_HAL_DFT_IS_CONTINUOUS 512
#define CV_HAL_DFT_IS_INPLACE 1024
//! @}

#define CV_8U   0
#define CV_8S   1
#define CV_16U  2
#define CV_16S  3
#define CV_32S  4
#define CV_32F  5
#define CV_64F  6
#define CV_16F  7

class OcvDctImpl
{
public:
    OcvDftOptions opt;

    int _factors[34];
    AutoBuffer<uint> wave_buf;
    AutoBuffer<int> itab_buf;

    DCTFunc dct_func;
    bool isRowTransform;
    bool isInverse;
    bool isContinuous;
    int start_stage;
    int end_stage;
    int width;
    int height;
    int depth;

    void init(int _width, int _height, int _depth, int flags)
    {
        width = _width;
        height = _height;
        depth = _depth;
        isInverse = (flags & CV_HAL_DFT_INVERSE) != 0;
        isRowTransform = (flags & CV_HAL_DFT_ROWS) != 0;
        isContinuous = (flags & CV_HAL_DFT_IS_CONTINUOUS) != 0;
        static DCTFunc dct_tbl[4] =
        {
            (DCTFunc)DCT_32f,
            (DCTFunc)NULL,
            (DCTFunc)DCT_64f,
            (DCTFunc)NULL
        };
        dct_func = dct_tbl[(int)isInverse + (depth == CV_64F)*2];
        opt.nf = 0;
        opt.isComplex = false;
        opt.isInverse = false;
        opt.noPermute = false;
        opt.scale = 1.;
        opt.factors = _factors;

        if (isRowTransform || height == 1 || (width == 1 && isContinuous))
        {
            start_stage = end_stage = 0;
        }
        else
        {
            start_stage = (width == 1);
            end_stage = 1;
        }
    }
    void apply(const uchar *src, size_t src_step, uchar *dst, size_t dst_step)
    {
        // CV_IPP_RUN(IPP_VERSION_X100 >= 700 && depth == CV_32F, ippi_DCT_32f(src, src_step, dst, dst_step, width, height, isInverse, isRowTransform))

        AutoBuffer<uchar> dct_wave;
        AutoBuffer<uchar> src_buf, dst_buf;
        uchar *src_dft_buf = 0, *dst_dft_buf = 0;
        int prev_len = 0;
        int elem_size = (depth == CV_32F) ? sizeof(float) : sizeof(double);
        int complex_elem_size = elem_size*2;

        for(int stage = start_stage ; stage <= end_stage; stage++ )
        {
            const uchar* sptr = src;
            uchar* dptr = dst;
            size_t sstep0, sstep1, dstep0, dstep1;
            int len, count;

            if( stage == 0 )
            {
                len = width;
                count = height;
                if( len == 1 && !isRowTransform )
                {
                    len = height;
                    count = 1;
                }
                sstep0 = src_step;
                dstep0 = dst_step;
                sstep1 = dstep1 = elem_size;
            }
            else
            {
                len = height;
                count = width;
                sstep1 = src_step;
                dstep1 = dst_step;
                sstep0 = dstep0 = elem_size;
            }

            opt.n = len;
            opt.tab_size = len;

            if( len != prev_len )
            {
                if( len > 1 && (len & 1) )
                    throw std::runtime_error("Odd-size DCT\'s are not implemented");

                opt.nf = DFTFactorize( len, opt.factors );
                bool inplace_transform = opt.factors[0] == opt.factors[opt.nf-1];

                wave_buf.allocate(len*complex_elem_size);
                opt.wave = wave_buf.data();
                itab_buf.allocate(len);
                opt.itab = itab_buf.data();
                DFTInit( len, opt.nf, opt.factors, opt.itab, complex_elem_size, opt.wave, isInverse );

                dct_wave.allocate((len/2 + 1)*complex_elem_size);
                src_buf.allocate(len*elem_size);
                src_dft_buf = src_buf.data();
                if(!inplace_transform)
                {
                    dst_buf.allocate(len*elem_size);
                    dst_dft_buf = dst_buf.data();
                }
                else
                {
                    dst_dft_buf = src_buf.data();
                }
                DCTInit( len, complex_elem_size, dct_wave.data(), isInverse);
                prev_len = len;
            }
            // otherwise reuse the tables calculated on the previous stage
            for(unsigned i = 0; i < static_cast<unsigned>(count); i++ )
            {
                dct_func( opt, sptr + i*sstep0, sstep1, src_dft_buf, dst_dft_buf,
                          dptr + i*dstep0, dstep1, dct_wave.data());
            }
            src = dst;
            src_step = dst_step;
        }
    }
};

}

std::vector<float> opencv_dct(const std::vector<float> &src, int square_size, size_t step) {
    std::vector<float> dst(src.size());
    auto dct = std::make_unique<cv::OcvDctImpl>();
    dct->init(square_size, square_size, CV_32F, 0);
    dct->apply(reinterpret_cast<const cv::uchar*>(src.data()), step, reinterpret_cast<cv::uchar*>(dst.data()), step);
    return dst;
}