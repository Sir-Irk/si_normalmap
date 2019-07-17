# si_normalmap
Simple Single Header Normal Map Generator

## WORK IN PROGRESS. DO NOT USE(YET)

A very simple normal map generator written as a single header library. 

Features:
 - Convert color buffer to greyscale using either lightness, average or luminance methods(with SSE/AVX versions available)
 - Convert greyscale buffer to a normal map
 
Possible features to add
 - SIMD version of normal map generation
 - Filtering options for normal map generation (ie: 3x3 sobel, etc)
 - Blurring functionality. But most likely this will be made in a seperate library. I'd prefer to keep this simple.

Basic example:
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
