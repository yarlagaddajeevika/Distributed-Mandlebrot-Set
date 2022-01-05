/**
 * @file leader.h
 * @author Dalton Caron (dpcaron@csu.fullerton.edu)
 * @author Jeevika Yarlagadda (jeevikayarlagadda@csu.fullerton.edu)
 **/
#ifndef LEADER_H
#define LEADER_H

#include "common.h"

 /**
  * @brief Populates each mandelbrot_partition_t with data related to the mandelbrot
  * partition assigned to process p at index p-1 of partitions. Note that there are
  * size-1 partitions due to there being a leader node.
  *
  * @param receive_data 1 or 0 if the follower should send image data along with the entropy.
  * @param settings The settings associated with the mandelbrot calculation.
  * @param partitions The mandelbrot partition information for each process. Initialized but empty as input.
  * @param partition_size The length of a partition cell; symmetrical on the x and y axises.
  * @param image_size The size of the image.
  * @param cx Cell partition origin x in image.
  * @param cy Cell partition origin y in image.
  * @param temporary_float temporary variable to hold the float values.
  */
void populate_and_send_partitions(const unsigned char receive_data, mandelbrot_settings_t *settings, 
    mandelbrot_partition_t *partitions, const int partition_size, const int image_size, 
    mpf_t cx, mpf_t cy, mpf_t temporary_float);

/**
 * @brief Receives the calculated partition and entropy from the follower. The partition 
 * is combined with the solution and the process rank-1 with the most entropy is returned. 
 * 
 * @param image Complete solution image. 
 * @param partition_size The length of a partition cell; symmetrical on the x and y axises. 
 * @param size The total amount of processes in this communicatior. 
 * @param partition_map A one-to-one map of processes IDs of rank-1 to mandelbrot_partition_t. 
 * @param successor_process_out The process with the maximum entropy set by this function. 
 */
void receive_and_populate_solution(const int image_size, unsigned char image[image_size][image_size], 
    const int partition_size, const int size, const mandelbrot_partition_t *const partition_map,  
    int *const successor_process_out);

/**
 * @brief Receives the entropy from the follower.
 *
 * @param size The total amount of processes in this communicatior.
 * @param partition_map A one-to-one map of processes IDs of rank-1 to mandelbrot_partition_t.
 * @param successor_process The process with the maximum entropy.
 */
void receive_entropies(const int size, int *const successor_process);

/**
 * @brief Calculates a the mandelbrot set based upon the current mandelbrot settings. 
 * 
 * @param settings The settings associated with the mandelbrot calculation.
 * @param partitions The mandelbrot partition information for each process. Initialized but empty as input.
 * @param size The size of the MPI communicator.
 * @param image_size The size of the image.
 * @param save_image Integer either 0 or 1 to save image, used at the last itereation.
 * @param cx Cell partition origin x in image.
 * @param cy Cell partition origin y in image.
 * @param temporary_float temporary variable to hold the float values.
 */
void perform_mandelbrot_iteration(mandelbrot_settings_t *const settings, mandelbrot_partition_t *partition_map, 
    const int size, const int image_size, const unsigned char save_image, mpf_t cx, mpf_t cy, mpf_t temporary_float);

/**
 * @brief Starts initiator node work.
 * Based on the iterations provided perform the mandelbrot iteration and save the final image, terminate all the followers.
 *
 * @param size The size of the MPI communicator.
 * @param settings The settings associated with the mandelbrot calculation.
 * @param image_size The size of the image.
 * @param iterations The number of zoom iterations to be performed.
 * @param save_each_iteration Integer with 0 or 1, to save the iteration or not.
 */
void start_initiator(const int size, mandelbrot_settings_t *settings, const int image_size, 
    const int iterations, const int save_each_iteration);

#endif