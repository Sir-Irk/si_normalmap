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

#ifndef _MSC_VER
    #ifdef __cplusplus
    #define sinm__inline inline
    #else
    #define sinm__inline
    #endif
#else
    #define sinm__inline __forceinline
#endif

#ifdef _MSC_VER
    #define sinm__aligned_var(type, bytes) __declspec(align(bytes)) type
#else
    #define sinm__aligned_var(type, bytes) type __attribute__((aligned(bytes)))
#endif

#ifndef SINM_TYPES
#define SINM_TYPES
typedef enum
{
    sinm_greyscale_none,
    sinm_greyscale_lightness,
    sinm_greyscale_average,
    sinm_greyscale_luminance,
    sinm_greyscale_count,
} sinm_greyscale_type;

typedef int8_t  i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef float f32;
typedef double f64;

typedef i32 b32;

#endif

#ifndef SI_NORMALMAP_IMPLEMENTATION

SINM_DEF void sinm_greyscale(const uint32_t *in, uint32_t *out, int32_t w, int32_t h, sinm_greyscale_type type);
//Converts values in "buffer" to greyscale  using either the
//lightness, average or luminance methods
//Result can be produced in-place if "in" and "out" are the same buffers

SINM_DEF uint32_t *sinm_normal_map(const uint32_t *in, int32_t w, int32_t h, float scale, float blurRadius, sinm_greyscale_type greyscaleType, int flipY);
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
#define simd__setzero_ps()      _mm256_setzero_ps()
#define simd__and_ix(a, b)      _mm256_and_si256(a, b)
#define simd__or_ix(a, b)       _mm256_or_si256(a, b)
#define simd__andnot_ps(a, b)   _mm256_andnot_ps(a, b)
#define simd__add_epi32(a, b)   _mm256_add_epi32(a, b)
#define simd__sub_epi32(a, b)   _mm256_sub_epi32(a, b)
#define simd__max_epi32(a, b)   _mm256_max_epi32(a, b)
#define simd__min_epi32(a, b)   _mm256_min_epi32(a, b)
#define simd__loadu_ix(a)       _mm256_loadu_si256(a)
#define simd__loadu_ps(a)       _mm256_loadu_ps(a)
#define simd__storeu_ix(ptr, v) _mm256_storeu_si256(ptr, v)
#define simd__srli_epi32(a, i)  _mm256_srli_epi32(a, i)
#define simd__slli_epi32(a, i)  _mm256_slli_epi32(a, i)
#define simd__set1_ps(a)        _mm256_set1_ps(a)
#define simd__cvtepi32_ps(a)    _mm256_cvtepi32_ps(a)
#define simd__cvtps_epi32(a)    _mm256_cvtps_epi32(a)
#define simd__add_ps(a, b)      _mm256_add_ps(a, b)
#define simd__mul_ps(a, b)      _mm256_mul_ps(a, b)
#define simd__sqrt_ps(a)        _mm256_sqrt_ps(a)
#define simd__cmp_ps(a, b, c)   _mm256_cmp_ps(a, b, c)
#define simd__div_ps(a, b)      _mm256_div_ps(a, b)
#define simd__hadd_ps(a, b)     _mm256_hadd_ps(a, b)
#define simd__cvtss_f32(a)      _mm256_cvtss_f32(a)

#else

#define SINM_SIMD_INCREMENT 4

#define simd__int __m128i
#define simd__float __m128
#define simd__set1_epi32(a)     _mm_set1_epi32(a)
#define simd__setzero_ix()      _mm_setzero_si128()
#define simd__setzero_ps()      _mm_setzero_ps()
#define simd__and_ix(a, b)      _mm_and_si128(a, b)
#define simd__or_ix(a, b)       _mm_or_si128(a, b)
#define simd__andnot_ps(a, b)   _mm_andnot_ps(a, b)
#define simd__add_epi32(a, b)   _mm_add_epi32(a, b)
#define simd__sub_epi32(a, b)   _mm_sub_epi32(a, b)
#define simd__max_epi32(a, b)   _mm_max_epi32(a, b)
#define simd__min_epi32(a, b)   _mm_min_epi32(a, b)
#define simd__loadu_ix(a)       _mm_loadu_si128(a)
#define simd__loadu_ps(a)       _mm_loadu_ps(a)
#define simd__storeu_ix(ptr, v) _mm_storeu_si128(ptr, v)
#define simd__srli_epi32(a, i)  _mm_srli_epi32(a, i)
#define simd__slli_epi32(a, i)  _mm_slli_epi32(a, i)
#define simd__set1_ps(a)        _mm_set1_ps(a)
#define simd__cvtepi32_ps(a)    _mm_cvtepi32_ps(a)
#define simd__cvtps_epi32(a)    _mm_cvtps_epi32(a)
#define simd__add_ps(a, b)      _mm_add_ps(a, b)
#define simd__mul_ps(a, b)      _mm_mul_ps(a, b)
#define simd__sqrt_ps(a)        _mm_sqrt_ps(a)
#define simd__cmp_ps(a, b, c)   _mm_cmp_ps(a, b, c)
#define simd__div_ps(a, b)      _mm_div_ps(a, b)
#define simd__hadd_ps(a, b)     _mm_hadd_ps(a, b)
#define simd__cvtss_f32(a)      _mm_cvtss_f32(a)

#endif //AVX_AVAILABLE
#endif //SI_NORMALMAP_NO_SIMD

#define sinm__min(a, b) ((a) < (b) ? (a) : (b))
#define sinm__max(a, b) ((a) > (b) ? (a) : (b))

typedef struct 
{
    int32_t x, y;
} sinm__v2i;

typedef struct 
{
    float x,y,z;
} sinm__v3;

sinm__inline static float 
sinm__length(float x, float y, float z) 
{ 
    return sqrtf(x*x + y*y + z*z);
}

sinm__inline static simd__float 
sinm__length_simd(simd__float x, simd__float y, simd__float z) 
{
    return simd__sqrt_ps(simd__add_ps(simd__add_ps(simd__mul_ps(x, x), simd__mul_ps(y, y)), simd__mul_ps(z, z)));
}

sinm__inline static sinm__v3 
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

sinm__inline static uint32_t 
sinm__lightness_average(uint32_t r, uint32_t g, uint32_t b)
{
    return (sinm__max(sinm__max(r, g), b)+sinm__min(sinm__min(r, g), b))/2;
}

sinm__inline static uint32_t 
sinm__average(uint32_t r, uint32_t g, uint32_t b)
{
    return (r+g+b)/3;
}

//NOTE: bias is based on human eye sensitivity 
sinm__inline static uint32_t 
sinm__luminance(uint32_t r, uint32_t g, uint32_t b)
{
    return (uint32_t)(0.21f*r+0.72f*g+0.07f*b);
}

sinm__inline static uint32_t
sinm__greyscale_from_byte(uint8_t c)
{
    return (c | c << 8u | c << 16u | 255u << 24u);  
}

static sinm__inline sinm__v3
sinm__rgba_to_v3(uint32_t c)
{
    sinm__v3 result = {
        (float)((c >>  0) & 0xFFu)-127.0f,
        (float)((c >>  8) & 0xFFu)-127.0f,
        (float)((c >> 16) & 0xFFu)-127.0f
    };
    
    return result;
}

static sinm__inline void
sinm__rgba_to_v3_simd(simd__int c, simd__float *x, simd__float *y, simd__float *z)
{
    simd__int ff   = simd__set1_epi32(0xFF);
    simd__int v127 = simd__set1_epi32(127);
    *x = simd__cvtepi32_ps(simd__sub_epi32(simd__and_ix(simd__srli_epi32(c, 0) , ff), v127));
    *y = simd__cvtepi32_ps(simd__sub_epi32(simd__and_ix(simd__srli_epi32(c, 8) , ff), v127));
    *z = simd__cvtepi32_ps(simd__sub_epi32(simd__and_ix(simd__srli_epi32(c, 16), ff), v127));
}

static sinm__inline uint32_t
sinm__unit_vector_to_rgba(sinm__v3 v)
{
    uint32_t r = (uint32_t)((1.0f+v.x)*127.0f);
    uint32_t g = (uint32_t)((1.0f+v.y)*127.0f);
    uint32_t b = (uint32_t)((1.0f+v.z)*127.0f);
    return r | g << 8u | b << 16u | 255u << 24u;
}

static sinm__inline simd__int
sinm__v3_to_rgba_simd(simd__float x, simd__float y, simd__float z)
{
    simd__float one = simd__set1_ps(1.0f);
    simd__float v127 = simd__set1_ps(127.0f);
    simd__int a = simd__set1_epi32(255u << 24u);
    simd__int r = simd__cvtps_epi32(simd__mul_ps(simd__add_ps(one, x), v127));
    simd__int g = simd__cvtps_epi32(simd__mul_ps(simd__add_ps(one, y), v127));
    simd__int b = simd__cvtps_epi32(simd__mul_ps(simd__add_ps(one, z), v127));
    simd__int c = simd__or_ix(simd__or_ix(simd__or_ix(r, simd__slli_epi32(g, 8)), simd__slli_epi32(b, 16)), a);
    return c;
}

SINM_DEF void
sinm__generate_gaussian_box(float *outBoxes, int32_t n, float sigma)
{
    float wIdeal = sqrtf((12.0f*sigma*sigma/n)+1.0f);
    int32_t wl = floorf(wIdeal);
    if(wl%2 == 0) --wl;
    int32_t wu = wl+2;

    float mIdeal = (12.0f*sigma*sigma - n*wl*wl - 4.0f*n*wl - 3.0f*n)/(-4.0f*wl - 4.0f);
    int32_t m = roundf(mIdeal);

    for(int i = 0; i < n; ++i) {
        outBoxes[i] = (i < m) ? wl : wu;
    }
}

//NOTE: decently optimized box blur based on http://blog.ivank.net/fastest-gaussian-blur.html
SINM_DEF void 
sinm__box_blur_h(uint32_t *in, uint32_t *out, int32_t w, int32_t h, float r)
{
    float invR = 1.0f/(r+r+1);
    for(int i = 0; i < h; ++i) {
        int32_t  oi  = i*w;
        int32_t  li  = oi;
        int32_t  ri  = oi+r;
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
        int32_t  oi  = i;
        int32_t  li  = oi;
        int32_t  ri  = oi+r*w;
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
            li+=w; ri+=w; oi+=w;
        }
        for(int j = h-r; j < h; j++) {
            sum += lv - (in[li] & 0xFFu);
            out[oi] = sinm__greyscale_from_byte((uint8_t)(sum*invR));
            li+=w; oi+=w;
        }
    }
}

SINM_DEF void 
sinm__gaussian_box(uint32_t *in, uint32_t *out, int32_t w, int32_t h, float r)
{
    float boxes[3];
    sinm__generate_gaussian_box(boxes, sizeof(boxes)/sizeof(boxes[0]), r);

    sinm__box_blur_h(in, out, w, h, (boxes[0]-1)/2);
    sinm__box_blur_v(out, in, w, h, (boxes[0]-1)/2);
    sinm__box_blur_h(in, out, w, h, (boxes[1]-1)/2);
    sinm__box_blur_v(out, in, w, h, (boxes[1]-1)/2);
    sinm__box_blur_h(in, out, w, h, (boxes[2]-1)/2);
    sinm__box_blur_v(out, in, w, h, (boxes[2]-1)/2);

    memcpy(out, in, w*h*sizeof(uint32_t));
}

SINM_DEF void 
sinm__sobel3x3_normals_row_range(const uint32_t *in, uint32_t *out, int32_t xs, int32_t xe, int32_t w, int32_t h, float scale, int flipY)
{
    const float xk[3][3] = { {-1,  0,  1}, {-2,  0,  2}, {-1,  0,  1}, };
    const float yk[3][3] = { {-1, -2, -1}, { 0,  0,  0}, { 1,  2,  1}, };

    for(int32_t y = 0; y < h; ++y) {
        for(int32_t x = xs; x < xe; ++x) {
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
            out[y*w+x] = sinm__unit_vector_to_rgba(color);
        }
    }
}

static sinm__inline void 
sinm__sobel3x3_normals(const uint32_t *in, uint32_t *out, int32_t w, int32_t h, float scale, int flipY)
{
    sinm__sobel3x3_normals_row_range(in, out, 0, w, w, h, scale, flipY);
}

static void
sinm__sobel3x3_normals_simd(const uint32_t *in, uint32_t *out, int32_t w, int32_t h, float scale, int flipY)
{
    const float xk[3][4] = { {-1,  0,  1,  0}, {-2,  0,  2,  0}, {-1,  0,  1,  0}, };
    const float yk[3][4] = { {-1, -2, -1,  0}, { 0,  0,  0,  0}, { 1,  2,  1,  0}, };

    simd__float simdScale = simd__set1_ps(scale);
    simd__float simdFlipY = simd__set1_ps((flipY) ? -1.0f : 1.0);
    simd__float simd1     = simd__set1_ps(1.0f);
    simd__float simd127   = simd__set1_ps(127.0f);

    int batchCount = 0;
    sinm__aligned_var(float, SINM_SIMD_INCREMENT) xBatch[SINM_SIMD_INCREMENT];
    sinm__aligned_var(float, SINM_SIMD_INCREMENT) yBatch[SINM_SIMD_INCREMENT];

    for(int32_t yIter = 0; yIter < h; ++yIter) {
        for(int32_t xIter = SINM_SIMD_INCREMENT; xIter < w-SINM_SIMD_INCREMENT; ++xIter) {
            __m128 xmag = _mm_set1_ps(0.0f);
            __m128 ymag = _mm_set1_ps(0.0f);

            for(int32_t a = 0; a < 3; ++a) {
                int32_t xIdx = sinm__min(w-1, sinm__max(1, xIter-1));
                int32_t yIdx = sinm__min(h-1, sinm__max(1, yIter+a-1));
                int32_t index = yIdx*w+xIdx;

                __m128i pixel = _mm_loadu_si128((__m128i *)&in[index]);
                pixel  = _mm_and_si128(pixel, _mm_set1_epi32(0xFFu));
                __m128 pixelf = _mm_cvtepi32_ps(pixel);
                __m128 kx = _mm_loadu_ps((float *)&xk[a]);
                __m128 ky = _mm_loadu_ps((float *)&yk[a]);
                xmag = _mm_add_ps(_mm_mul_ps(pixelf, kx), xmag);
                ymag = _mm_add_ps(_mm_mul_ps(pixelf, ky), ymag);
            }

            __m128 xSum = _mm_hadd_ps(xmag, xmag);
            __m128 ySum = _mm_hadd_ps(ymag, ymag);
            float xn = _mm_cvtss_f32(_mm_hadd_ps(xSum, xSum));
            float yn = _mm_cvtss_f32(_mm_hadd_ps(ySum, ySum));

            xBatch[batchCount] = xn;
            yBatch[batchCount++] = yn;
            if(batchCount == SINM_SIMD_INCREMENT) {
                batchCount = 0;
                simd__float x = simd__loadu_ps(xBatch);
                simd__float y = simd__loadu_ps(yBatch);
                simd__float z = simd__set1_ps(255.0f);
                simd__float len = sinm__length_simd(x, y, z);
                simd__float invLen = simd__div_ps(simd__set1_ps(1.0f), len);
                simd__float mask = simd__cmp_ps(simd__setzero_ps(), len, _CMP_EQ_OQ);                                      
                x = simd__andnot_ps(mask, simd__mul_ps(x, invLen));                                                                
                y = simd__andnot_ps(mask, simd__mul_ps(y, invLen));                                                                
                z = simd__andnot_ps(mask, simd__mul_ps(z, invLen));                        
                int index = yIter*w+(xIter-(SINM_SIMD_INCREMENT-1));
                simd__storeu_ix((simd__int *)&out[index], sinm__v3_to_rgba_simd(x, y, z));
            }
        }
    }

    sinm__sobel3x3_normals_row_range(in, out, 0, SINM_SIMD_INCREMENT, w, h, scale, flipY);
    sinm__sobel3x3_normals_row_range(in, out, w-SINM_SIMD_INCREMENT, w, w, h, scale, flipY);
}

SINM_DEF void
sinm_normalize(uint32_t *in, int32_t w, int32_t h, float scale, int flipY)
{
    float invScale = 1.0f/scale;
    float yDir = (flipY) ? -1.0f : 1.0f;
    for(int32_t i = 0; i < w*h; ++i) {
        sinm__v3 v = sinm__rgba_to_v3(in[i]);
        in[i] = sinm__unit_vector_to_rgba(sinm__normalized(v.x, v.y*yDir, v.z*invScale));
    }
}

SINM_DEF void
sinm__normalize_simd(uint32_t *in, int32_t w, int32_t h, float scale, int flipY)
{
    assert(w % SINM_SIMD_INCREMENT == 0);
    for(int32_t i = 0; i < w*h; i += SINM_SIMD_INCREMENT) {
        simd__int pixel = simd__loadu_ix((simd__int *)&in[i]);
        simd__float x, y, z;
        sinm__rgba_to_v3_simd(pixel, &x, &y, &z);
        simd__float len = sinm__length_simd(x, y, z);
        simd__float invLen = simd__div_ps(simd__set1_ps(1.0f), len);
        simd__float mask = simd__cmp_ps(simd__setzero_ps(), len, _CMP_EQ_OQ);                                      
        x = simd__andnot_ps(mask, simd__mul_ps(x, invLen));                                                                
        y = simd__andnot_ps(mask, simd__mul_ps(y, invLen));                                                                
        z = simd__andnot_ps(mask, simd__mul_ps(z, invLen));                        
        simd__storeu_ix((simd__int *)&in[i], sinm__v3_to_rgba_simd(x, y, z));
    }
}

SINM_DEF void
sinm_composite(const uint32_t *in, uint32_t *out, int32_t w, int32_t h)
{
    for(int32_t i = 0; i < w*h; ++i) {
        uint32_t ic = in[i], oc = out[i];
        uint32_t r = (uint32_t)((((ic >> 0 ) & 0xFFu) + ((oc >> 0 ) & 0xFFu)) * 0.5f);
        uint32_t g = (uint32_t)((((ic >> 8 ) & 0xFFu) + ((oc >> 8 ) & 0xFFu)) * 0.5f); 
        uint32_t b = (uint32_t)((((ic >> 16) & 0xFFu) + ((oc >> 16) & 0xFFu)) * 0.5f); 
        out[i] = (r | g << 8u | b << 16u | 255u << 24u);
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
#endif //SI_NORMALMAP_NO_SIMD

SINM_DEF void
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
}

SINM_DEF int
sinm_normal_map_buffer(const uint32_t *in, uint32_t *out, int32_t w, int32_t h, float scale, float blurRadius, sinm_greyscale_type greyscaleType, int flipY, int useSimd)
{
    uint32_t *intermediate = (uint32_t *)malloc(w*h*sizeof(uint32_t));

    if(intermediate) {
        if(greyscaleType != sinm_greyscale_none) {
            sinm_greyscale(in, out, w, h, greyscaleType);
        } else {
            memcpy(out, in, w*h*sizeof(uint32_t));
        }

        float radius = sinm__min(sinm__min(w,h), sinm__max(0, blurRadius));
        if(radius >= 1.0f) {
            sinm__gaussian_box(out, intermediate, w, h, radius);
        } else {
            memcpy(intermediate, out, w*h*sizeof(uint32_t));
        }
#ifdef SI_NORMALMAP_NO_SIMD
            sinm__sobel3x3_normals(intermediate, out, w, h, scale, flipY);
#else
            sinm__sobel3x3_normals_simd(intermediate, out, w, h, scale, flipY);
#endif
        }
        free(intermediate);
        return 1;
    }
    return 0;
}

SINM_DEF uint32_t *
sinm_normal_map(const uint32_t *in, int32_t w, int32_t h, float scale, float blurRadius, sinm_greyscale_type greyscaleType, int flipY, int useSimd)
{
    uint32_t *result = (uint32_t *)malloc(w*h*sizeof(uint32_t));
    if(result) sinm_normal_map_buffer(in, result, w, h, scale, blurRadius, greyscaleType, flipY, useSimd);
    return result;
}

#if 0
//NOTE: returns RGB or RGBA depending on input
SINM_DEF uint8_t *
sinm__normal_map_u8(const uint8_t *in, uint32_t w, int32_t h, float scale, float blurRadius, 
                    sinm_greyscale_type greyscaleType, int flipY, uint32_t numChannels)
{
    if(numChannels > 4 || numChannels == 0) return NULL;

    if(numChannels == 4) {
        return (uint8_t *)sinm__normal_map((uint32_t *)in, w, h, scale, blurRadius, greyscaleType, flipY);
    }

    size_t sizeInBytes = w*h*sizeof(uint32_t);
    uint32_t *intermediate = (uint32_t *)malloc(w*h*sizeof(uint32_t));
    if(!intermediate) return NULL;
    uint32_t *result = (uint32_t *)malloc(w*h*sizeof(uint32_t));

    switch(numChannels) {
        case 3 : { 
            for(int i = 0 i < w*h; ++i) {
                int idx = numChannels*i;
                result[i] = (in[idx] | in[idx+1] << 8u | in[idx+2] << 16u | 255u << 24u);
            }
        } break;
        case 1 : { 
            for(int i = 0 i < w*h; ++i) {
                int idx = numChannels*i;
                result[i] = sinm__greyscale_from_byte(in[idx]);
            }
        }
    }

    if(result) {
        if(numChannels > 1 && greyscaleType != sinm_greyscale_none) {
            //sinm_greyscale_u8(in, result, w, h, greyscaleType, numChannels);
        } else {
            //memcpy(result, in, sizeInBytes);
        }

        float radius = sinm__min(sinm__min(w,h), sinm__max(0, blurRadius));
        sinm__gaussian_box(result, intermediate, w, h, radius);
        sinm__sobel3x3_normals(intermediate, result, w, h, scale, flipY);
    }

    free(intermediate);
    return result;
}
#endif

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

