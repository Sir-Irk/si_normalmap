#define SI_NORMALMAP_IMPLEMENTATION
#define SI_NORMALMAP_STATIC
#include "../si_normalmap.h"

// stb_image and stb_image_write not included in this repo
// get them from here https://github.com/nothings/stb
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int
main(void)
{
    int x, y;
    uint32_t *pixels = (uint32_t *)stbi_load("plasma_214.png", &x, &y, NULL, 4);
    assert(pixels);

    // Make 2 normal maps at different blur radii and strength and composite them
    // together.
    uint32_t *nm0 = sinm_normal_map(pixels, x, y, 2.0f, 1.0f, sinm_greyscale_luminance, 1);
    uint32_t *nm1 = sinm_normal_map(pixels, x, y, 2.0f, 8.0f, sinm_greyscale_luminance, 1);
    uint32_t *composite = sinm_composite_alloc(nm0, nm1, x, y);
    sinm_normalize(composite, x, y, 1.0f, 0);

    stbi_write_png("normal.png", x, y, 4, composite, 0);

    free(nm0);
    free(nm1);
    free(composite);
    free(pixels);
}
