/**
 * @file common.h
 * @author Dalton Caron (dpcaron@csu.fullerton.edu)
 * @author Jeevika Yarlagadda (jeevikayarlagadda@csu.fullerton.edu)
 **/
#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <stdio.h>

// Enables the wait_for_debugger function. 
#define ENABLE_DEBUGGER_TRICK 0

#if ENABLE_DEBUGGER_TRICK > 0
    #include <unistd.h> // Only for getpid().
#endif

#define COLOR_MAP_LENGTH 18 // Length of the color map array

// Calculate offset
#define CALCULATE_C_OFFSET(c, image_size) (c / ((double)(image_size - 1)))

#define size_cx_cy 2^20 // We are using this to set maximum precision. 
#define N_DIGITS 100 

// DO NOT SEND OVER MPI. 
typedef struct mandelbrot_settings {
    mpf_t z_start_x;
    mpf_t z_start_y;
    mpf_t offset_x;
    mpf_t offset_y;
    int32_t max_iterations;
} mandelbrot_settings_t;

typedef enum follower_request {PARTITION_NO_IMAGE, PARTITION_SEND_IMAGE, TERMINATE} follower_request_t;

// DO NOT SEND OVER MPI. 
typedef struct mandelbrot_partition {
    mpf_t cx;
    mpf_t cy;
    int16_t i;
    int16_t j;
    int32_t clength;
} mandelbrot_partition_t;

typedef struct actual_partition {
    int32_t clength;
    int32_t max_iterations;
    int64_t cx_s_len; 
    int64_t cy_s_len; 
    int64_t offset_x_s_len;
    int64_t offset_y_s_len;
} partition_t;

// There must be no variable width items as the struct is sent through MPI. 
typedef struct follower_message {
    follower_request_t request_type;
    union {
        partition_t partition;
    } payload;
} follower_message_t;

// Received partition structure.
typedef struct {
    int32_t clength;
    int32_t max_iterations;
    mpf_t cx;
    mpf_t cy;
    mpf_t offset_x;
    mpf_t offset_y;
} received_partition_t;

#if ENABLE_DEBUGGER_TRICK > 0
    /**
     * @brief Traps process with rank attach_rank so a debugger may attach 
     * to it before the MPI program begins execution. All processes are 
     * held behind a barrier until the process of rank attach_rank is 
     * modified by a debugger to manually leave an infinite loop. See 
     * @link{http://www.sci.utah.edu/~tfogal/academic/Fogal-ParallelDebugging.pdf} 
     * for more information about this trick. Set ENABLE_DEBUGGER_TRICK to enable 
     * this feature during development. 
     * 
     * @param rank The rank of the current process executing this function.
     * @param attach_rank The rank of the process we want to trap and attach 
     * debugger to. 
     */
    void wait_for_debugger(const int rank, const int attach_rank);
#endif

/**
 * @brief Computes the iterative formula f(z) = z^2 + c. 
 * 
 * @return unsigned char The color map index of the iteration. 18 implies no divergence (black). 
 */
unsigned char iteration_map(received_partition_t *const data);

/**
 * @brief Saves an image using color map indices for each pixel value. 
 * 
 * @param width Image width. 
 * @param height Image height. 
 * @param image_map Matrix of color map indices. 
 * @param name File name to save as. 
 */
void save_mandelbrot_image(const int width, const int height, 
    unsigned char *const image_map, const char *const name);

/**
 * @brief Saves an image without using the colormap. 
 * 
 * @param width Image width.
 * @param height Image height.
 * @param image Image to save.
 * @param name File name to save as. 
 */
void save_grayscale_image(const int width, const int height, uint16_t *const image, 
    const char *const name);

/**
 * @brief Calculates the GLCM for an image. The GLCM is normalized. There 
 * are as many pixels values as COLOR_MAP_LENGTH. 
 * See @link{https://en.wikipedia.org/wiki/Co-occurrence_matrix} for more info on the GLCM. 
 * 
 * @param distance The offset operator delta(x, y) where x = y = distance. 
 * @param rows Number of rows in the image. 
 * @param cols Number of columns in the image. 
 * @param image Image to calculate GLCM for. 
 * @param colors_count The amount of colors with the image is the width and height of the GLCM. 
 * @param matrix_out Matrix for storing the GLCM. 
 */
void calculate_gray_level_co_occurrence_matrix(const uint16_t distance, const uint16_t rows, 
    const uint16_t cols, const unsigned char *const image, const int colors_count, double *const matrix_out);

/**
 * @brief Computes shanon entropy using the GLCM matrix as the probability distribution. 
 * See @link{https://en.wikipedia.org/wiki/Entropy_(information_theory)} for more info on entropy. 
 * 
 * @param glcm The normalized GLCM matrix to serve as the probability distribution.
 * @param glcm_length The width and height of the GLCM matrix.  
 * @return double of the entropy value.
 */
double compute_entropy(const double *glcm, const int glcm_length);

/**
 * @brief Normalizes the GLCM matrix. 
 * 
 * @param width The width of the matrix.
 * @param height The height of the matrix. 
 * @param matrix Matrix to normalize. 
 */
void normalize_matrix(const int width, const int height, double *const matrix);

/**
 * @brief Calculates the total shanon entropy within the provided image map. 
 * 
 * @param width Width of image.
 * @param height Height of image.
 * @param image The image map of pixels into color indices.
 * @param colors_count The total amount of colors. 
 * @return double The total shanon entropy within the image. 
 */
double calculate_entropy_from_image(const int width, const int height, unsigned char *const image, const int colors_count);

/**
 * @brief Writes the mandlebrot parametes to the file given as input.
 *
 * @param cx Cell partition origin x in image.
 * @param cy Cell partition origin y in image.
 * @param offset_x offset value of cx.
 * @param offset_y offset value of cy.
 * @param file_name Name of the file to write the set.
 */
void write_mandelbrot_partition_parameters_to_file(const mpf_t cx, const mpf_t cy, const mpf_t offset_x, const mpf_t offset_y, const char *file_name);

/**
 * @brief Reads the mandlebrot parametes from the file given as input.
 *
 * @param cx Cell partition origin x in image.
 * @param cy Cell partition origin y in image.
 * @param offset_x offset value of cx.
 * @param offset_y offset value of cy.
 * @param file_name Name of the file to write the set.
 */
char read_mandelbrot_partition_parameters_from_file(mpf_t *cx, mpf_t *cy, mpf_t *offset_x, mpf_t *offset_y, const char *file_name);

/**
 * @brief Serializes an mpf_t to a string buffer. Designed by 
 * Mr. John Ruble, but modified to take a buffer rather than 
 * a file stream. 
 * 
 * @param buffer A buffer to store mpf_t data in.
 * @param buffer_len Length of the buffer data. 
 * @param X The mpf_t to write to the buffer. 
 * @return int Length of data written to the buffer. 
 */
int mpf_out_raw(char *buffer, const size_t buffer_len, const mpf_t X);

/**
 * @brief Deserializes a string buffer to an mpf_t. Designed by 
 * Mr. John Ruble, but modified to take a buffer rather than 
 * a file stream. 
 * 
 * @param buffer Buffer containing raw mpf_t data. 
 * @param buffer_len Length of the buffer data.
 * @param X mpf_t to populate with data. 
 */
void mpf_inp_raw(char *buffer, const size_t buffer_len, mpf_t X);

#endif

