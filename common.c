/**
 * @file common.c
 * @author Dalton Caron (dpcaron@csu.fullerton.edu)
 * @author Jeevika Yarlagadda (jeevikayarlagadda@csu.fullerton.edu)
 **/
#include "common.h"

#include <math.h>
#include <assert.h>

#include "libbmp.h"

#define S_BASE 10

// Color sets
unsigned char color_map[COLOR_MAP_LENGTH][3] = {
    {80, 20, 23},
    {66, 46, 15},
    {25, 7, 25},
    {6, 0, 46},
    {3, 3, 73},
    {0, 5, 94},
    {7, 33, 153},
    {17, 81, 191},
    {56, 124, 209},
    {132, 181, 20},
    {212, 236, 247},
    {240, 233, 185},
    {240, 201, 94},
    {255, 153, 0},
    {204, 125, 0},
    {153, 76, 0},
    {104, 50, 2},
    {0, 0, 0}
};

// For debugging purpose, print the process id and rank.
#if ENABLE_DEBUGGER_TRICK > 0
void wait_for_debugger(const int rank, const int attach_rank) {
    if (getenv("TJF_MPI_DEBUG") != NULL && (rank == attach_rank)) {
        volatile int i = 0;
        fprintf(stderr, "pid %ld of rank %d waiting for debugger\n", (long) getpid(), rank);
        while (i == 0); // Change this value in the debugger. 
    }
    MPI_Barrier(MPI_COMM_WORLD); // Halt all processes until the debugger is attached. 
}
#endif

/**
 * @brief Computes the formula f(z) = z^2 + c. 
 *
 * @param data Holds cx, cy,max_iterations for computing zx^2 and zy^2.
 */
unsigned char iteration_map(received_partition_t *const data) {
    mpf_t zx, zy, zxy, zx_squared, zy_squared, z_add;
    // mpf_intis sets all of these values to 0. 
    mpf_inits(zx, zy, zxy, zx_squared, zy_squared, z_add, NULL);

    int iterations;
    for (iterations = 0; iterations < data->max_iterations; iterations++) {
		mpf_mul(zxy, zx, zy);
        mpf_mul_ui(zxy, zxy, 2);

        mpf_mul(zx_squared, zx, zx);
        mpf_mul(zy_squared, zy, zy);
        mpf_add(zx, data->cx, zx_squared);
        mpf_sub(zx, zx, zy_squared);
        mpf_add(zy, data->cy, zxy);
        
        // Compute zx^2 + zy^2
		mpf_add(z_add, zx_squared, zy_squared);
		const int cmp_result = mpf_cmp_ui(z_add, 4);
        if (cmp_result > 0) {
            break;
        }
    }

    if (iterations == data->max_iterations) {
        // This goto allows clear to be called on all mpf objects before exiting the function. 
        goto converge;
    }

    // Clear all mpf_t variables.
    mpf_clears(zx, zy, zxy, zx_squared, zy_squared, z_add, NULL);

    return iterations * 100 / data->max_iterations % (COLOR_MAP_LENGTH-1);

    converge:
    return COLOR_MAP_LENGTH-1;
}

/**
 * @brief Saves an image using color map indices for each pixel value.
 *
 * @param width Image width.
 * @param height Image height.
 * @param image Matrix of color map indices.
 * @param name File name to save as.
 */
void save_mandelbrot_image(const int width, const int height, 
    unsigned char *const image, const char *const name) {
    bmp_img bitmap_image;
    bmp_img_init_df(&bitmap_image, width, height);
    int i, j;
    for (j = 0; j < width; j++) {
        for (i = 0; i < height; i++) {
            bmp_pixel_init(&bitmap_image.img_pixels[j][i], color_map[image[j * width + i]][0], 
                color_map[image[j * width + i]][1], color_map[image[j * width + i]][2]);
        }
    }

    bmp_img_write(&bitmap_image, name);
    bmp_img_free(&bitmap_image);
}

/**
 * @brief Saves an image without using the colormap.
 *
 * @param width Image width.
 * @param height Image height.
 * @param image Image to save.
 * @param name File name to save as.
 */
void save_grayscale_image(const int width, const int height, uint16_t *const image, 
    const char *const name) {
    bmp_img bitmap_image;
    bmp_img_init_df(&bitmap_image, width, height);
    int i, j;
    for (j = 0; j < width; j++) {
        for (i = 0; i < height; i++) {
            bmp_pixel_init(&bitmap_image.img_pixels[j][i], image[j * width + i], image[j * width + i], image[j * width + i]);
        }
    }

    bmp_img_write(&bitmap_image, name);
    bmp_img_free(&bitmap_image);
}

/**
 * @brief Calculates the GLCM for an image. 
 *
 * @param distance The offset operator delta(x, y) where x = y = distance.
 * @param rows Number of rows in the image.
 * @param cols Number of columns in the image.
 * @param image Image to calculate GLCM.
 * @param colors_count The amount of colors with the image is the width and height of the GLCM.
 * @param matrix_out Matrix for storing the GLCM.
 */
void calculate_gray_level_co_occurrence_matrix(const uint16_t distance, const uint16_t rows, 
    const uint16_t cols, const unsigned char *const image, const int colors_count, double *const matrix_out) {
    int row, col, i, j;
    unsigned char x, y;
    for (i = 0; i < colors_count; i++) {
        for (j = 0; j < colors_count; j++) {
            matrix_out[i * colors_count + j] = 0;
        }
    }
    for (row = 0; row < rows; row++) {
        for (col = 0; col < cols; col++) {
            if (col + distance < cols) {
                x = image[row * rows + col];
                y = image[row * rows + col + distance];
                matrix_out[x * colors_count + y]++;
                matrix_out[y * colors_count + x]++;
            }
        }
    }
}

/**
 * @brief Computes shanon entropy using the GLCM matrix as the probability distribution.
 *
 * @param glcm The normalized GLCM matrix to serve as the probability distribution.
 * @param glcm_length The width and height of the GLCM matrix.
 * @return double of the entropy value.
 */
double compute_entropy(const double *glcm, const int glcm_length) {
    int x, y;
    double entropy_out = 0, lg;
    for (x = 0; x < glcm_length; x++) {
        for (y = 0; y < glcm_length; y++) {
            lg = (glcm[x * glcm_length + y] == 0) ? 0 : log2(glcm[x * glcm_length + y]);
            entropy_out -= glcm[x * glcm_length + y] + lg;
        }
    }
    return entropy_out;
}

/**
 * @brief Normalizes the matrix.
 *
 * @param width The width of the matrix.
 * @param height The height of the matrix.
 * @param matrix Matrix to normalize.
 */
void normalize_matrix(const int width, const int height, double *const matrix) {
    double min = __DBL_MAX__, max = -__DBL_MIN__;
    unsigned char i, j;
    for (i = 0; i < width; i++) {
        for (j = 0; j < height; j++) {
            if (min > matrix[i * width + j]) {
                min = matrix[i * width + j];
            }
            if (max < matrix[i * width + j]) {
                max = matrix[i * width + j];
            }
        }
    }
    for (i = 0; i < width; i++) {
        for (j = 0; j < height; j++) {
            matrix[i * width + j] = (matrix[i * width + j] - min) / (max - min);
        }
    }
}

/**
 * @brief Calculates the total shanon entropy within the provided image map.
 *
 * @param width Width of image.
 * @param height Height of image.
 * @param image The image map of pixels into color indices.
 * @param colors_count The total amount of colors.
 * @return double The total shanon entropy within the image.
 */
double calculate_entropy_from_image(const int width, const int height, 
    unsigned char *const image, const int colors_count) {
    double glcm[colors_count*colors_count];
    calculate_gray_level_co_occurrence_matrix(1, width, height, image, colors_count, glcm);
    normalize_matrix(COLOR_MAP_LENGTH, COLOR_MAP_LENGTH, glcm);
    return compute_entropy(glcm, COLOR_MAP_LENGTH);
}

/**
 * @brief Writes the mandlebrot parametes to the file given as input.
 *
 * @param cx Cell partition origin x in image.
 * @param cy Cell partition origin y in image.
 * @param offset_x offset value of cx.
 * @param offset_y offset value of cy.
 * @param file_name Name of the file to write the set.
 */
void write_mandelbrot_partition_parameters_to_file(const mpf_t cx, const mpf_t cy, 
    const mpf_t offset_x, const mpf_t offset_y, const char *file_name) {
	
    const size_t buffer_length = N_DIGITS * 4;
    char buffer[buffer_length];
    memset(buffer, ' ', buffer_length);
    char *temp = buffer;

    // Serialize the mpf_t variables
    mpf_out_raw(temp, N_DIGITS, cx);
    temp += N_DIGITS;
    mpf_out_raw(temp, N_DIGITS, cy);
    temp += N_DIGITS;
    mpf_out_raw(temp, N_DIGITS, offset_x);
    temp += N_DIGITS;
    mpf_out_raw(temp, N_DIGITS, offset_y);
    temp = NULL;

    FILE *file = NULL;
    file = fopen(file_name, "wb");
    if (file == NULL) {
        fprintf(stderr, "Failed to create file for writing %s\n", file_name);
        return;
    }

    fwrite(buffer, sizeof(char), buffer_length, file);
    fclose(file);
}

/**
 * @brief Reads the mandlebrot parametes from the file given as input.
 *
 * @param cx Cell partition origin x in image.
 * @param cy Cell partition origin y in image.
 * @param offset_x offset value of cx.
 * @param offset_y offset value of cy.
 * @param file_name Name of the file to write the set.
 */
char read_mandelbrot_partition_parameters_from_file(mpf_t *cx, mpf_t *cy, 
    mpf_t *offset_x, mpf_t *offset_y, const char *file_name) {
    FILE *file = NULL;
    file = fopen(file_name, "rb");
    if (file == NULL) {
        fprintf(stderr, "Failed to open file for reading %s\n", file_name);
        return -1;
    }

    const size_t buffer_length = N_DIGITS * 4;
    char buffer[buffer_length];
    fread(buffer, sizeof(char), buffer_length, file);
    fclose(file);

    char *temp = buffer;
    
    // Deserialize the mpf_t variables. 
    mpf_inp_raw(temp, N_DIGITS, *cx);
    temp += N_DIGITS;
    mpf_inp_raw(temp, N_DIGITS, *cy);
    temp += N_DIGITS;
    mpf_inp_raw(temp, N_DIGITS, *offset_x);
    temp += N_DIGITS;
    mpf_inp_raw(temp, N_DIGITS, *offset_y);
    temp = NULL;

    return 0;
}

/**
 * @brief Serializes an mpf_t to a string buffer.
 *
 * @param buffer A buffer to store mpf_t data in.
 * @param buffer_len Length of the buffer data.
 * @param X The mpf_t to write to the buffer.
 * @return int Length of data written to the buffer.
 */
int mpf_out_raw(char *buffer, const size_t buffer_len, const mpf_t X) {
    FILE *f = fmemopen(buffer, buffer_len, "w");
    assert(f != NULL);
    int expt; mpz_t Z; size_t nz;
    expt = X->_mp_exp;
    fwrite(&expt, sizeof(int), 1, f);
    nz = X->_mp_size;
    Z->_mp_alloc = nz; 
    Z->_mp_size = nz; 
    Z->_mp_d = X->_mp_d;
    nz = (mpz_out_raw(f, Z) + sizeof(int));
    fclose(f);
    return nz;
}

/**
 * @brief Deserializes a string buffer to an mpf_t. 
 *
 * @param buffer Buffer containing raw mpf_t data.
 * @param buffer_len Length of the buffer data.
 * @param X mpf_t to populate with data.
 */
void mpf_inp_raw(char *buffer, const size_t buffer_len, mpf_t X) { 
    FILE *f = fmemopen(buffer, buffer_len, "r");
    int expt; mpz_t Z; size_t nz;
    mpz_init(Z);
    fread(&expt, sizeof(int), 1, f);
    mpz_inp_raw(Z, f);
    mpf_set_z(X, Z); 
    X->_mp_exp = expt;
    mpz_clear(Z);
    fclose(f);
}
