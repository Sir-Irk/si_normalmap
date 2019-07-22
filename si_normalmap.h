/* LICENSE AT END OF FILE */

/***************************************************************************
 * Sir Irk's normal map generator
 *
 * Basic use:
 *     #define SI_NORMALMAP_IMPLEMENTATION before including this file to get
        the implementation. Otherwise this acts as a regualr header file
        
 *     uint32_t *in = ...load pixels from image
 *     uint32_t *nm = sinm_normal_map(in, w, h, scale, blurRadius, greyscaleType);
 *     ...write normal map to a file
 *
 *  Other defines you can use(before including this file):
 *  #define SI_NORMALMAP_NO_SIMD to ignore SSE/AVX functions
 *  #define SI_NORMALMAP_STATIC for static defintions(no extern functions)
 ***************************************************************************/


#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef SINM_DEF
#ifdef SI_NORMALMAP_STATIC
#define SINM_DEF static
#else
#define SINM_DEF extern
#endif
#endif

#ifndef _MSV_VER
    #ifdef __cplusplus
    #define sinm_inline inline
    #else
    #define sinm_inline
    #endif
#else
    #define sinm_inline __forceinline
#endif

#define sinm__min(a, b) ((a) < (b) ? (a) : (b))
#define sinm__max(a, b) ((a) > (b) ? (a) : (b))

typedef enum
{
    sinm_greyscale_none,
    sinm_greyscale_lightness,
    sinm_greyscale_average,
    sinm_greyscale_luminance,
    sinm_greyscale_count,
} sinm_greyscale_type;

#ifndef SI_NORMALMAP_IMPLEMENTATION
SINM_DEF void sinm_greyscale(const uint32_t *in, uint32_t *out, int32_t w, int32_t h, sinm_greyscale_type type);
//Converts values in "buffer" to greyscale  using either the
//lightness, average or luminance methods
//Result can be produced in-place if "in" and "out" are the same buffers

SINM_DEF uint32_t *sinm_normal_map(const uint32_t *in, int32_t w, int32_t h, float scale, float blurRadius, sinm_greyscale_type greyscaleType);
//Converts input buffer to a normal map and returns a pointer to it. 
//  "scale" controls the intensity of the result
//  "blurRadius" controls the radius for gaussian blurring before generating normals
//  "greyscaleType" specifies the conversion method from color to greyscale before 
//   generating the normal map. This step is skipped when using sinm_greyscale_none.

#else //SI_NORMALMAP_IMPLEMENTATION

#ifndef SI_NORMALMAP_NO_SIMD

#include <intrin.h> 
#include <emmintrin.h> 

#ifdef __AVX__

#define SINM_SIMD_INCREMENT 8

#define simd__int __m256i
#define simd__float __m256
#define simd__set1_epi32(a)     _mm256_set1_epi32(a)
#define simd__setzero_ix()      _mm256_setzero_si256()
#define simd__setzero_rx()      _mm256_setzero_ps()
#define simd__and_ix(a, b)      _mm256_and_si256(a, b)
#define simd__or_ix(a, b)       _mm256_or_si256(a, b)
#define simd__add_epi32(a, b)   _mm256_add_epi32(a, b)
#define simd__sub_epi32(a, b)   _mm256_sub_epi32(a, b)
#define simd__max_epi32(a, b)   _mm256_max_epi32(a, b)
#define simd__min_epi32(a, b)   _mm256_min_epi32(a, b)
#define simd__loadu_ix(a)       _mm256_loadu_si256(a)
#define simd__storeu_ix(ptr, v) _mm256_storeu_si256(ptr, v)
#define simd__srli_epi32(a, i)  _mm256_srli_epi32(a, i)
#define simd__slli_epi32(a, i)  _mm256_slli_epi32(a, i)
#define simd__set1_ps(a)        _mm256_set1_ps(a)
#define simd__cvtepi32_ps(a)    _mm256_cvtepi32_ps(a)
#define simd__cvtps_epi32(a)    _mm256_cvtps_epi32(a)
#define simd__add_ps(a, b)      _mm256_add_ps(a, b)
#define simd__mul_ps(a, b)      _mm256_mul_ps(a, b)

#else

#define SINM_SIMD_INCREMENT 4

#define simd__int __m128i
#define simd__float __m128
#define simd__set1_epi32(a)     _mm_set1_epi32(a)
#define simd__setzero_ix()      _mm_setzero_si128()
#define simd__setzero_rx()      _mm_setzero_ps()
#define simd__and_ix(a, b)      _mm_and_si128(a, b)
#define simd__or_ix(a, b)       _mm_or_si128(a, b)
#define simd__add_epi32(a, b)   _mm_add_epi32(a, b)
#define simd__sub_epi32(a, b)   _mm_sub_epi32(a, b)
#define simd__max_epi32(a, b)   _mm_max_epi32(a, b)
#define simd__min_epi32(a, b)   _mm_min_epi32(a, b)
#define simd__loadu_ix(a)       _mm_loadu_si128(a)
#define simd__storeu_ix(ptr, v) _mm_storeu_si128(ptr, v)
#define simd__srli_epi32(a, i)  _mm_srli_epi32(a, i)
#define simd__slli_epi32(a, i)  _mm_slli_epi32(a, i)
#define simd__set1_ps(a)        _mm_set1_ps(a)
#define simd__cvtepi32_ps(a)    _mm_cvtepi32_ps(a)
#define simd__cvtps_epi32(a)    _mm_cvtps_epi32(a)
#define simd__add_ps(a, b)      _mm_add_ps(a, b)
#define simd__mul_ps(a, b)      _mm_mul_ps(a, b)

#endif //AVX_AVAILABLE
#endif //SI_NORMALMAP_NO_SIMD

typedef struct 
{
    int32_t x, y;
} sinm__v2i;

typedef struct 
{
    float x,y,z;
} sinm__v3;


sinm_inline static float 
sinm__length(float x, float y, float z) 
{ 
    return sqrtf(x*x + y*y + z*z);
}

sinm_inline static sinm__v3 
sinm__normalized(float x, float y, float z) 
{
    sinm__v3 result;
    float len = sinm__length(x, y, z);

    if(len > 1e-04f) {
        float invLen = 1.0f / len;
        result.x = x*invLen;
        result.y = y*invLen;
        result.z = z*invLen;
    } else {
        result.x = result.y = result.z = 0.0f;
    }

    return result;
}

sinm_inline static uint32_t 
sinm__lightness_average(uint32_t r, uint32_t g, uint32_t b)
{
    return (sinm__max(sinm__max(r, g), b)+sinm__min(sinm__min(r, g), b))/2;
}

sinm_inline static uint32_t 
sinm__average(uint32_t r, uint32_t g, uint32_t b)
{
    return (r+g+b)/3;
}

//NOTE: bias is based on human eye sensitivity 
sinm_inline static uint32_t 
sinm__luminance(uint32_t r, uint32_t g, uint32_t b)
{
    return (uint32_t)(0.21f*r+0.72f*g+0.07f*b);
}

sinm_inline static uint32_t
sinm__greyscale_from_byte(uint8_t c)
{
    return (c | c << 8u | c << 16u | 255u << 24u);  
}

SINM_DEF double *
sinm__generate_gaussian_box(uint32_t n, double sigma)
{
    double wIdeal = sqrt((12.0*sigma*sigma/n)+1);
    double wl = floor(wIdeal);
    if((int64_t)wl%2 == 0) --wl;
    double wu = wl+2;

    double mIdeal = (12.0*sigma*sigma - n*wl*wl - 4*n*wl - 3*n)/(-4*wl - 4);
    double m = round(mIdeal);

    double *boxes = (double *)malloc(n*sizeof(double));
    for(int i = 0; i < n; ++i) boxes[i] = (i < m) ? wl : wu;
    return boxes;
}

//NOTE: decently optimized box blur based on http://blog.ivank.net/fastest-gaussian-blur.html
SINM_DEF void 
sinm__box_blur_h(uint32_t *in, uint32_t *out, int32_t w, int32_t h, int32_t r)
{
    float invR = 1.0f/(r+r+1);
    for(int i = 0; i < h; ++i) {
        int32_t oi = i*w;
        int32_t li = oi;
        int32_t ri = oi+r;

        uint32_t fv  = in[oi] & 0xFFu;
        uint32_t lv  = in[oi+w-1] & 0xFFu;
        uint32_t sum = (r+1)*fv;

        for(int j = 0; j < r; ++j) {
            sum += in[oi+j] & 0xFFu;
        }
        for(int j = 0; j <= r; ++j) {
            sum += (in[ri++] & 0xFFu) - fv;
            out[oi++] = sinm__greyscale_from_byte((uint8_t)(sum*invR));
        }
        for(int j = r+1; j < w-r; ++j) {
            sum += (in[ri++] & 0xFFu) - (in[li++] & 0xFFu);
            out[oi++] = sinm__greyscale_from_byte((uint8_t)(sum*invR));
        }
        for(int j = w-r; j < w; ++j) {
            sum += lv - (in[li++] & 0xFFu);
            out[oi++] = sinm__greyscale_from_byte((uint8_t)(sum*invR));
        }
    }
}

SINM_DEF void 
sinm__box_blur_v(uint32_t *in, uint32_t *out, int32_t w, int32_t h, float r)
{
    float invR = 1.0f/(r+r+1);
    for(int i = 0; i < w; ++i) {
        int32_t oi = i;
        int32_t li = oi;
        int32_t ri = oi+r*w;

        uint32_t fv  = in[oi] & 0xFFu;
        uint32_t lv  = in[oi+w*(h-1)] & 0xFFu;
        uint32_t sum = (r+1)*fv;

        for(int j = 0; j < r; j++) {
            sum += in[oi+j*w] & 0xFFu;
        }
        for(int j = 0; j <= r; j++) {
            sum += (in[ri] & 0xFFu) - fv;
            out[oi] = sinm__greyscale_from_byte((uint8_t)(sum*invR));
            ri+=w; oi+=w;
        }
        for(int j = r+1; j < h-r; j++) {
            sum += (in[ri] & 0xFFu) - (in[li] & 0xFFu);
            out[oi] = sinm__greyscale_from_byte((uint8_t)(sum*invR));
            li += w; ri+=w; oi+=w;
        }
        for(int j = h-r; j < h; j++) {
            sum += lv - (in[li] & 0xFFu);
            out[oi] = sinm__greyscale_from_byte((uint8_t)(sum*invR));
            li += w; oi+=w;
        }
    }
}

SINM_DEF void 
sinm__gaussian_box(uint32_t *in, uint32_t *out, int32_t w, int32_t h, float r)
{
    double *boxes = sinm__generate_gaussian_box(3, r);
    if(boxes) {
        sinm__box_blur_h(in, out, w, h, (boxes[0]-1)/2);
        sinm__box_blur_v(out, in, w, h, (boxes[0]-1)/2);
        sinm__box_blur_h(in, out, w, h, (boxes[1]-1)/2);
        sinm__box_blur_v(out, in, w, h, (boxes[1]-1)/2);
        sinm__box_blur_h(in, out, w, h, (boxes[2]-1)/2);
        sinm__box_blur_v(out, in, w, h, (boxes[2]-1)/2);
        memcpy(out, in, w*h*sizeof(uint32_t));
        free(boxes);
    }
}


SINM_DEF void 
sinm__sobel3x3_normals(const uint32_t *in, uint32_t *out, int32_t w, int32_t h, float scale)
{
    const float xk[3][3] = {
        {-1,  0,  1},
        {-2,  0,  2},
        {-1,  0,  1},
    };

    const float yk[3][3] = {
        {-1, -2, -1},
        { 0,  0,  0},
        { 1,  2,  1},
    };

    //TODO: optimize
    for(int32_t y = 0; y < h; ++y) {
        for(int32_t x = 0; x < w; ++x) {
            float xmag = 0.0f;
            float ymag = 0.0f;

            for(int32_t a = 0; a < 3; ++a) {
                for(int32_t b = 0; b < 3; ++b) {
                    int32_t xIdx = sinm__min(w-1, sinm__max(1, x+b-1));
                    int32_t yIdx = sinm__min(h-1, sinm__max(1, y+a-1));
                    int32_t index = yIdx*w+xIdx;
                    uint32_t pixel = in[index] & 0xFFu;
                    xmag += pixel*xk[a][b];
                    ymag += pixel*yk[a][b];
                }
            }

            sinm__v3 color = sinm__normalized(xmag*scale, ymag*scale, 255.0f);
            uint32_t r = (uint32_t)((1.0f+color.x)*0.5f*255);
            uint32_t g = (uint32_t)((1.0f+color.y)*0.5f*255);
            uint32_t b = (uint32_t)((1.0f+color.z)*0.5f*255);
            out[y*w+x] = (r | g << 8u | b << 16u | 255 << 24u);
        }
    }
}

static void
sinm__greyscale(const uint32_t *in, uint32_t *out, int32_t w, int32_t h, sinm_greyscale_type type)
{
    int32_t count = w*h;
    switch(type) {
        case sinm_greyscale_lightness : {
            for(int32_t i = 0; i < count; ++i) {
                uint32_t c = in[i];
                uint32_t l = sinm__lightness_average(c & 0xFFu, (c >> 8) & 0xFFu, (c >> 16) & 0xFFu);
                out[i] = sinm__greyscale_from_byte(l);
            }
        } break;

        case sinm_greyscale_average : {
            for(int32_t i = 0; i < count; ++i) {
                uint32_t c = in[i];
                uint32_t l = sinm__average(c & 0xFFu, (c >> 8) & 0xFFu, (c >> 16) & 0xFFu);
                out[i] = sinm__greyscale_from_byte(l);
            }
        } break;

        case sinm_greyscale_luminance : {
            for(int32_t i = 0; i < count; ++i) {
                uint32_t c = in[i];
                uint32_t l = sinm__luminance(c & 0xFFu, (c >> 8) & 0xFFu, (c >> 16) & 0xFFu);
                out[i] = sinm__greyscale_from_byte(l);
            }
        } break;
    }
}

#ifndef SI_NORMALMAP_NO_SIMD
static void
sinm__simd_greyscale(const uint32_t *in, uint32_t *out, int32_t w, int32_t h, sinm_greyscale_type type)
{
    simd__int redMask   = simd__set1_epi32(0xFF);
    simd__int greenMask = simd__set1_epi32(0xFF00u);
    simd__int blueMask  = simd__set1_epi32(0xFF0000u);
    simd__int alpha     = simd__set1_epi32(0xFF000000u);

    int32_t count = w*h;

    switch(type) {
        case sinm_greyscale_lightness : {
            for(int32_t i = 0; i < count; i += SINM_SIMD_INCREMENT) {
                simd__int c = simd__loadu_ix((simd__int *)&in[i]);
                simd__int r = simd__and_ix(c, redMask);
                simd__int g = simd__srli_epi32(simd__and_ix(c, greenMask), 8);
                simd__int b = simd__srli_epi32(simd__and_ix(c, blueMask), 16);

                simd__int max = simd__max_epi32(simd__max_epi32(r, g), b);
                simd__int min = simd__min_epi32(simd__min_epi32(r, g), b);
                simd__int l   = simd__srli_epi32(simd__add_epi32(min, max), 1);

                l = simd__or_ix(simd__slli_epi32(l, 16), 
                    simd__or_ix(simd__slli_epi32(l, 8), 
                    simd__or_ix(l, alpha)));

                simd__storeu_ix((simd__int *)&out[i], l);
            }
        } break;

        case sinm_greyscale_average : {
            simd__float inverse3 = simd__set1_ps(1.0f/3.0f);
            for(int32_t i = 0; i < count; i += SINM_SIMD_INCREMENT) {
                simd__int c = simd__loadu_ix((simd__int *)&in[i]);
                simd__int r = simd__and_ix(c, redMask);
                simd__int g = simd__srli_epi32(simd__and_ix(c, greenMask), 8);
                simd__int b = simd__srli_epi32(simd__and_ix(c, blueMask), 16);

                //NOTE: integer division is only available in SVML(not sse or avx) which I'm 
                //not going to use. Agner Fog has an efficient division implementation but I am 
                //simply going to use an inverse multiply in float and convert back to int. 
                //This may be changed later
                simd__int s = simd__add_epi32(simd__add_epi32(r, g), b);
                s = simd__cvtps_epi32(simd__mul_ps(simd__cvtepi32_ps(s), inverse3));
                s = simd__or_ix(simd__slli_epi32(s, 16), 
                    simd__or_ix(simd__slli_epi32(s, 8), 
                    simd__or_ix(s, alpha)));

                simd__storeu_ix((simd__int *)&out[i], s);
            }
        } break;

        case sinm_greyscale_luminance : {
            simd__float rBias = simd__set1_ps(0.21f);
            simd__float gBias = simd__set1_ps(0.72f);
            simd__float bBias = simd__set1_ps(0.07f);

            for(int32_t i = 0; i < count; i += SINM_SIMD_INCREMENT) {
                simd__int c = simd__loadu_ix((simd__int *)&in[i]);
                simd__float r = simd__cvtepi32_ps(simd__and_ix(c, redMask));
                simd__float g = simd__cvtepi32_ps(simd__srli_epi32(simd__and_ix(c, greenMask), 8));
                simd__float b = simd__cvtepi32_ps(simd__srli_epi32(simd__and_ix(c, blueMask), 16));

                r = simd__mul_ps(r, rBias);
                g = simd__mul_ps(g, gBias); 
                b = simd__mul_ps(b, bBias);

                simd__int sum = simd__cvtps_epi32(simd__add_ps(r, simd__add_ps(g, b)));

                sum = simd__or_ix(simd__slli_epi32(sum, 16), 
                      simd__or_ix(simd__slli_epi32(sum, 8), 
                      simd__or_ix(sum, alpha)));

                simd__storeu_ix((simd__int *)&out[i], sum);
            }
        } break;
    }
}
#endif //SI_NORMALMAP_USE_SIMD

SINM_DEF int
sinm_greyscale(const uint32_t *in, uint32_t *out, int32_t w, int32_t h, sinm_greyscale_type type)
{
#ifndef SI_NORMALMAP_NO_SIMD
    int32_t count = w*h;
    if(count > 0 && count % SINM_SIMD_INCREMENT == 0) { 
        sinm__simd_greyscale(in, out, w, h, type);
    } else {
        sinm__greyscale(in, out, w, h, type);
    }
#else
    sinm__greyscale(in, out, w, h, type);
#endif
    return 1;
}

SINM_DEF uint32_t *
sinm_normal_map(const uint32_t *in, int32_t w, int32_t h, float scale, float blurRadius, sinm_greyscale_type greyscaleType)
{
    if(!in) return NULL;

    //Intermediate buffer for processing so we don't have to change the input buffer
    int shouldFreeIntermediate = 0;
    uint32_t *intermediate = (uint32_t *)malloc(w*h*sizeof(uint32_t));
    if(!intermediate) return NULL;

    uint32_t *result = (uint32_t *)malloc(w*h*sizeof(uint32_t));
    if(result) {
        if(greyscaleType != sinm_greyscale_none) {
            sinm_greyscale(in, result, w, h, greyscaleType);
        } else {
            memcpy(result, in, w*h*sizeof(uint32_t));
        }

        float radius = sinm__min(sinm__min(w,h), sinm__max(0, blurRadius));
        sinm__gaussian_box(result, intermediate, w, h, radius);
        sinm__sobel3x3_normals(intermediate, result, w, h, scale);
    }

    free(intermediate);
    return result;
}

