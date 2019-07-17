/* LICENSE AT END OF FILE */

/***************************************************************************
 * Sir Irk's normal map generator
 *
 * Basic use:
 *     ...load your pixels to a buffer and convert to greyscale
 *     sinm_greyscale(inBuffer, pixelCount, greyscaleConversionType);
 *
 *     ...make an output buffer for normal map gen from greyscale buffer
 *     sinm_normal_map(inBuffer, outBuffer, w, h, scale);
 *
 *     ...write output buffer to a file
 ***************************************************************************/


#include <stdint.h>

#ifndef SINM_DEF
#ifdef SINM_STATIC
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
    #define snmp_inline __forceinline
#endif

#ifndef SI_NORMALMAP_IMPLEMENTATION
enum sinm_greyscale_type;

//Converts values in "buffer" to greyscale  using either the
//lightness, average or luminance methods
SINM_DEF void sinm_greyscale(uint32_t *buffer, int32_t count, greyscale_type type);

//Converts input buffer to a normal map and writes it to the output buffer
//Colors should be converted to greyscale for proper results
//"scale" adjusts the intensity of the result
SINM_DEF void sinm_normal_map(const uint32_t *inBuffer, uint32_t *outBuffer, int32_t w, int32_t h, int32_t comp, float scale);

#else //SI_NORMALMAP_IMPLEMENTATION

#ifdef SI_NORMALMAP_USE_SIMD

#include <intrin.h> 
#include <emmintrin.h> 

#ifdef __AVX__
#define SINM_SIMD_INCREMENT 8

#define simd__int __m256i
#define simd__float __m256

#define simd__set1_epi32(a) _mm256_set1_epi32(a)
#define simd__setzero_ix()  _mm256_setzero_si256()
#define simd__setzero_rx()  _mm256_setzero_ps()

#define simd__and_ix(a, b) _mm256_and_si256(a, b)
#define simd__or_ix(a, b)  _mm256_or_si256(a, b)

#define simd__add_epi32(a, b) _mm256_add_epi32(a, b)
#define simd__sub_epi32(a, b) _mm256_sub_epi32(a, b)

#define simd__max_epi32(a, b) _mm256_max_epi32(a, b)
#define simd__min_epi32(a, b) _mm256_min_epi32(a, b)

#define simd__loadu_ix(a)  _mm256_loadu_si256(a)
#define simd__storeu_ix(ptr, v) _mm256_storeu_si256(ptr, v)

#define simd__srli_epi32(a, i) _mm256_srli_epi32(a, i)
#define simd__slli_epi32(a, i) _mm256_slli_epi32(a, i)

#define simd__set1_ps(a) _mm256_set1_ps(a)
#define simd__cvtepi32_ps(a) _mm256_cvtepi32_ps(a)
#define simd__cvtps_epi32(a) _mm256_cvtps_epi32(a)

#define simd__add_ps(a, b) _mm256_add_ps(a, b)
#define simd__mul_ps(a, b) _mm256_mul_ps(a, b)

#else

#define simd__int __m128i
#define simd__float __m128

#define SINM_SIMD_INCREMENT 4
#define simd__set1_epi32(a) _mm_set1_epi32(a)
#define simd__setzero_ix()  _mm_setzero_si128()
#define simd__setzero_rx()  _mm_setzero_ps()

#define simd__and_ix(a, b) _mm_and_si128(a, b)
#define simd__or_ix(a, b)  _mm_or_si128(a, b)

#define simd__add_epi32(a, b) _mm_add_epi32(a, b)
#define simd__sub_epi32(a, b) _mm_sub_epi32(a, b)

#define simd__max_epi32(a, b) _mm_max_epi32(a, b)
#define simd__min_epi32(a, b) _mm_min_epi32(a, b)

#define simd__loadu_ix(a)  _mm_loadu_si128(a)
#define simd__storeu_ix(ptr, v) _mm_storeu_si128(ptr, v)

#define simd__srli_epi32(a, i) _mm_srli_epi32(a, i)
#define simd__slli_epi32(a, i) _mm_slli_epi32(a, i)

#define simd__set1_ps(a) _mm_set1_ps(a)
#define simd__cvtepi32_ps(a) _mm_cvtepi32_ps(a)
#define simd__cvtps_epi32(a) _mm_cvtps_epi32(a)

#define simd__add_ps(a, b) _mm_add_ps(a, b)
#define simd__mul_ps(a, b) _mm_mul_ps(a, b)

#endif //AVX_AVAILABLE
#endif //SI_NORMALMAP_USE_SIMD

typedef enum
{
    sinm_greyscaleType_lightness,
    sinm_greyscaleType_average,
    sinm_greyscaleType_luminance,
} sinm_greyscale_type;

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

sinm_inline static float 
sinm__linearize_srgb(float value)
{
    return value*value;
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
    #define sinm__min(a, b) (((a) < (b)) ? (a) : (b))
    #define sinm__max(a, b) (((a) > (b)) ? (a) : (b))
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

static void
sinm__greyscale(uint32_t *buffer, int32_t count, sinm_greyscale_type type)
{
    switch(type) {
        case sinm_greyscaleType_lightness : {
            for(int32_t i = 0; i < count; ++i) {
                uint32_t c = buffer[i];
                uint32_t l = sinm__lightness_average(c & 0xFFu, (c >> 8) & 0xFFu, (c >> 16) & 0xFFu);
                buffer[i] = (255 << 24 | l << 16 | l << 8 | l);
            }
        } break;

        case sinm_greyscaleType_average : {
            for(int32_t i = 0; i < count; ++i) {
                uint32_t c = buffer[i];
                uint32_t l = sinm__average(c & 0xFFu, (c >> 8) & 0xFFu, (c >> 16) & 0xFFu);
                buffer[i] = (255 << 24 | l << 16 | l << 8 | l);
            }
        } break;

        case sinm_greyscaleType_luminance : {
            for(int32_t i = 0; i < count; ++i) {
                uint32_t c = buffer[i];
                uint32_t l = sinm__luminance(c & 0xFFu, (c >> 8) & 0xFFu, (c >> 16) & 0xFFu);
                buffer[i] = (255 << 24 | l << 16 | l << 8 | l);
            }
        } break;
    }
}

#ifdef SI_NORMALMAP_USE_SIMD
static void
sinm__simd_greyscale(uint32_t *buffer, int32_t count, sinm_greyscale_type type)
{
    simd__int redMask   = simd__set1_epi32(0xFF);
    simd__int greenMask = simd__set1_epi32(0xFF00u);
    simd__int blueMask  = simd__set1_epi32(0xFF0000u);
    simd__int alpha     = simd__set1_epi32(0xFF000000u);

    switch(type) {
        case sinm_greyscaleType_lightness : {
            for(int32_t i = 0; i < count; i += SINM_SIMD_INCREMENT) {
                simd__int c = simd__loadu_ix((simd__int *)&buffer[i]);
                simd__int r = simd__and_ix(c, redMask);
                simd__int g = simd__srli_epi32(simd__and_ix(c, greenMask), 8);
                simd__int b = simd__srli_epi32(simd__and_ix(c, blueMask), 16);

                simd__int max = simd__max_epi32(simd__max_epi32(r, g), b);
                simd__int min = simd__min_epi32(simd__min_epi32(r, g), b);
                simd__int l   = simd__srli_epi32(simd__add_epi32(min, max), 1);

                l = simd__or_ix(simd__slli_epi32(l, 16), simd__or_ix(simd__slli_epi32(l, 8), simd__or_ix(l, alpha)));

                simd__storeu_ix((simd__int *)&buffer[i], l);
            }
        } break;

        case sinm_greyscaleType_average : {
            simd__float inverse3 = simd__set1_ps(1.0f/3.0f);
            for(int32_t i = 0; i < count; i += SINM_SIMD_INCREMENT) {
                simd__int c = simd__loadu_ix((simd__int *)&buffer[i]);
                simd__int r = simd__and_ix(c, redMask);
                simd__int g = simd__srli_epi32(simd__and_ix(c, greenMask), 8);
                simd__int b = simd__srli_epi32(simd__and_ix(c, blueMask), 16);

                //NOTE: integer division is only available in SVML(not sse or avx) which I'm 
                //not going to use. Agner Fog has an efficient division implementation but I am 
                //simply going to use an inverse multiply in float and convert back to int. 
                //This may be changed later
                simd__int a = simd__add_epi32(simd__add_epi32(r, g), b);
                a = simd__cvtps_epi32(simd__mul_ps(simd__cvtepi32_ps(a), inverse3));
                simd__storeu_ix((simd__int *)&buffer[i], a);
            }
        } break;

        case sinm_greyscaleType_luminance : {
            simd__float rBias = simd__set1_ps(0.21f);
            simd__float gBias = simd__set1_ps(0.72f);
            simd__float bBias = simd__set1_ps(0.07f);

            for(int32_t i = 0; i < count; i += SINM_SIMD_INCREMENT) {
                simd__int c = simd__loadu_ix((simd__int *)&buffer[i]);
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

                simd__storeu_ix((simd__int *)&buffer[i], sum);
            }
        } break;
    }
}
#endif //SI_NORMALMAP_USE_SIMD

SINM_DEF void
sinm_greyscale(uint32_t *buffer, int32_t count, sinm_greyscale_type type)
{
#ifdef SI_NORMALMAP_USE_SIMD
    sinm__simd_greyscale(buffer, count, type);
#else
    sinm__greyscale(buffer, count, type);
#endif
}

SINM_DEF void 
sinm_normal_map(const uint32_t *inBuffer, uint32_t *outBuffer, int32_t w, int32_t h, float scale)
{
    const float to01space = 1.0f / 255.0f;

    for(int32_t y = 0; y < h; ++y) {
        for(int32_t x = 0; x < w; ++x) {

            int32_t index = y * w + x;
            uint32_t curP = inBuffer[index];

            //NOTE: If curP is at an edge then just use curP as the "adjacent" pixel.
            //This gives us a clean edge as opposed to interpolating to some default value.
            uint32_t prevP  = (x > 0)      ? inBuffer[index-1] : curP;
            uint32_t aboveP = (y > 0)      ? inBuffer[index-w] : curP;
            uint32_t nextP  = (x != (w-1)) ? inBuffer[index+1] : curP;
            uint32_t underP = (y != (h-1)) ? inBuffer[index+w] : curP;

            //TODO: Support other channels? or just convert input to the expected format(using alpha only)?
            float h0 = (curP   & 0xFFu)*to01space;
            float h1 = (nextP  & 0xFFu)*to01space;
            float h2 = (underP & 0xFFu)*to01space;
            float h3 = (aboveP & 0xFFu)*to01space;
            float h4 = (prevP  & 0xFFu)*to01space;

            //NOTE: Interpolate half-way between adjacent height values
            //for a "ramp" effect between pixels as opposed to a "flat wall"
            //if we use exact pixel values.
            float rDiff = ((h0-h1)+(h4-h0))*0.5f;
            float gDiff = ((h0-h2)+(h3-h0))*0.5f;

            //NOTE: Scale here will increase/decrease the angle between the x/y axis and the z axis.
            //In other words(for x axis) for dot(V3(rDiff*scale, 0.0f, 0.0f), V3(1.0f, 0.0f, 0.0f)) 
            //the higher the scale the closer the result is to 1.0f;
            sinm__v3 color = sinm__normalized(rDiff*scale, gDiff*scale, 1.0f);

            //NOTE: convert from fp -1.0f/1.0f space to 0.0f/1.0f space then to 0u/255u space
            uint32_t r = (uint32_t)((1.0f+color.x)*0.5f*255);
            uint32_t g = (uint32_t)((1.0f+color.y)*0.5f*255);
            uint32_t b = (uint32_t)((1.0f+color.z)*0.5f*255);

            outBuffer[index] = (r | g << 8u | b << 16u | 255 << 24u);
        }
    }
}
#endif //ifndef SI_NORMALMAP_IMPLEMENTATION
/*
Copyright (c) 2019 Jeremy Montgomery
Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to do 
so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
SOFTWARE.
*/
