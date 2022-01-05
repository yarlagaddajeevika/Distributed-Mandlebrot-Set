/**
 * @file follower.h
 * @author Dalton Caron (dpcaron@csu.fullerton.edu)
 * @author Jeevika Yarlagadda (jeevikayarlagadda@csu.fullerton.edu)
 **/
#ifndef FOLLOWER_H
#define FOLLOWER_H

#include "common.h"

/**
 * @brief Attempts to receive a partition description from the leader node of rank 0. 
 */
void follower_receive_partition(received_partition_t *partition_out);

void compute_image_partition(unsigned char *const image, received_partition_t *const partition);

// Starts worker process work. 
void start_follower(const int rank);

#endif