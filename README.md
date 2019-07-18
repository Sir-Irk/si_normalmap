# si_normalmap
A very simple normal map generator written as a single header library. 

## WORK IN PROGRESS. USE AT OWN RISK

Features:
 - Convert color buffer to greyscale using either lightness, average or luminance methods(with SSE/AVX versions available)
 - Convert greyscale buffer to a normal map
 
Possible features to add
 - SIMD version of normal map generation
 - Filtering options for normal map generation (ie: 3x3 sobel, etc)
 - Blurring functionality. But most likely this will be made in a seperate library. I'd prefer to keep this simple.

### Example(texture from opengameart.org):

Input:

![input](https://imgur.com/Grx9Uvs.png) 

Output:

![output](https://imgur.com/SWFhlh7.png)

### Available functions
```C
void sinm_greyscale(uint32_t *buffer, int32_t count, sinm_greyscale_type type);
void sinm_normal_map(const uint32_t *inBuffer, uint32_t *outBuffer, int32_t w, int32_t h, float scale);
```

### Basic Usage:
```C
#define SI_NORMALMAP_IMPLEMENTATION
#include "si_normalmap.h"

...

uint32_t *image = ...load pixels from some image;
sinm_greyscale(image, imageWidth*imageHeight, sinm_greyscaleType_average);

uint32_t *normalmap = ...allocate a second buffer of the same dimensions
sinm_normal_map(image, normalmap, imageWidth, imageHeight, 1.0f);
  
... write normalmap to a file

```
 
 Currently if you want to enable simd(#define SI_NORMALMAP_USE_SIMD) for the greyscale function
 you need to ensure the pixel count of the image is a multiple of either 4(for SSE) or 8(for AVX).
