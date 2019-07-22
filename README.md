# si_normalmap
A very simple normal map generator written as a single header library. 

## WORK IN PROGRESS. USE AT OWN RISK

Features:
 - Convert color buffer to greyscale using either lightness, average or luminance methods(with SSE/AVX versions available)
 - Convert greyscale buffer to a normal map with gaussian blur pre-filtering
 
Possible Todo(s):
 - SIMD optimization of normal map generation and gaussian blur

### Example(texture from opengameart.org):

Input:

![input](https://imgur.com/Grx9Uvs.png) 

Output:

![output](https://i.imgur.com/m64imlB.png)

Comparison of lighting without and with normal map(and a second "detail" normal map)

![lighting](https://imgur.com/CIw2oFB.png)

### Available functions
```C
void sinm_greyscale(uint32_t *buffer, int32_t count, sinm_greyscale_type type);
uint32_t *sinm_normal_map(const uint32_t *in, int32_t w, int32_t h, float scale, float blurRadius, sinm_greyscale_type greyscaleType);
```

### Basic Usage:
```C
#define SI_NORMALMAP_IMPLEMENTATION
#include "si_normalmap.h"

...

uint32_t *image = ...load pixels from some image;
uint32_t *nm = sinm_normal_map(image, imageWidth, imageHeight, 1.0f, 2.0f, sinm_greyscale_average); 

... write normalmap to a file

```
