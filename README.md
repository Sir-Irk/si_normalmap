# si_normalmap
A very simple normal map generator written as a single header library. 

Features:
 - Convert color buffer to greyscale using either lightness, average or luminance methods
 - Convert greyscale buffer to a normal map with gaussian blur pre-filtering
 - Uses AVX/SSE for images with power of 2 dimensions(TODO: make this more flexible for other dimensions)

TODO: add gpu support using Opengl(this is a work in progress).
 
### Interface
```C
typedef enum
{
    sinm_greyscale_none,
    sinm_greyscale_lightness,
    sinm_greyscale_average,
    sinm_greyscale_luminance,
    sinm_greyscale_count, //For iteration only. Not a valid option
} sinm_greyscale_type;

//run greyscale conversion. Can be done in-place if in and out are the same buffer
void sinm_greyscale(const uint32_t *in, uint32_t *out, int32_t w, int32_t h, sinm_greyscale_type type);

//generate normal map. Allocates a new buffer for its output
uint32_t *sinm_normal_map(const uint32_t *in, int32_t w, int32_t h, float scale, float blurRadius, sinm_greyscale_type greyscaleType, int flipY);

//generates normal map. Takes a pre-allocated buffer for its output
int sinm_normal_map_buffer(const uint32_t* in, uint32_t* out, int32_t w, int32_t h, float scale, float blurRadius, sinm_greyscale_type greyscaleType, int flipY)

//normalize the values of "in"
void sinm_normalize(uint32_t* in, int32_t w, int32_t h, float scale, int flipY)

//composite two buffers together. Useful for layering normal maps
void sinm_composite(const uint32_t* in1, const uint32_t* in2, uint32_t* out, int32_t w, int32_t h)

//same as sinm_composite but allocates a new buffer for its output
uint32_t* sinm_composite_alloc(const uint32_t* in1, const uint32_t* in2, int32_t w, int32_t h)
```

### Basic Usage:
```C
#define SI_NORMALMAP_IMPLEMENTATION
#include "si_normalmap.h"

int main()
{
    uint32_t *image = //...load pixels from some image;
    uint32_t *nm = sinm_normal_map(image, imageWidth, imageHeight, 1.0f, 2.0f, sinm_greyscale_average); 

    //... write normalmap to a file
    
    return 0;
}
```

### More specific Example using more features(error checking excluded for brevity. Using stb libraries to load and write https://github.com/nothings/stb)
```C
#define SI_NORMALMAP_IMPLEMENTATION
#define SI_NORMALMAP_STATIC
#include "si_normalmap.h"
#include "inttypes.h"
#include "stdio.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int main(void)
{
    int x, y;
    uint32_t* pixels = (uint32_t*)stbi_load("albedo.png", &x, &y, NULL, 4);

    // Make 2 normal maps at different blur radii and strength and composite them together.
    uint32_t* nm0 = sinm_normal_map(pixels, x, y, 2.0f, 1.0f, sinm_greyscale_luminance, 0);
    uint32_t* nm1 = sinm_normal_map(pixels, x, y, 3.0f, 8.0f, sinm_greyscale_luminance, 0);
    uint32_t* composite = sinm_composite_alloc(nm0, nm1, x, y);

    stbi_write_png("normal.png", x, y, 4, composite, 0);

    free(nm0);
    free(nm1);
    free(composite);
    free(pixels);
}
```

### Example Output(texture from opengameart.org):

Input:

![input](https://github.com/Sir-Irk/si_normalmap/blob/master/Examples/albedo.png) 

Output:

![output](https://github.com/Sir-Irk/si_normalmap/blob/master/Examples/normal_composite_example.png)

Comparison of lighting without and with normal map(and a second "detail" normal map)

![lighting](https://github.com/Sir-Irk/si_normalmap/blob/master/Examples/lighting.png)


