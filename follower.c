/**
 * @file follower.c
 * @author Dalton Caron (dpcaron@csu.fullerton.edu)
 * @author Jeevika Yarlagadda (jeevikayarlagadda@csu.fullerton.edu)
 **/
#include "follower.h"

#include <assert.h>

/**
 * @brief receives message from sender i.e, leader.
 *
 * @param message_received Consists of request_type and payload
 */
void follower_receive_message(follower_message_t *message_received) {
    MPI_Status status;
    MPI_Request request = MPI_REQUEST_NULL;

    MPI_Irecv(message_received, sizeof(follower_message_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);
}

/**
 * @brief receives partition data and de-serializes the data.
 *
 * @param partition The partition of the messageReceived.
 * @param data_out Pointer which holds clength(Cell partition width and height) and max_iterations.
 * @param rank Rank of a processor.
 */
void receive_partition_data(const partition_t *const partition, received_partition_t *data_out, const int rank) {
    MPI_Status status;
    MPI_Request request = MPI_REQUEST_NULL;

    data_out->clength = partition->clength;
    data_out->max_iterations = partition->max_iterations;

    size_t buffer_len = partition->cx_s_len + partition->cy_s_len + partition->offset_x_s_len + partition->offset_y_s_len;
    char buffer[buffer_len];
    char *temp = buffer;

    MPI_Irecv(buffer, buffer_len, MPI_BYTE, 0, 0, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);

    // Set the double values to data_out.
    mpf_set_d(data_out->cx, 0);
    mpf_set_d(data_out->cy, 0);
    mpf_set_d(data_out->offset_x, 0);
    mpf_set_d(data_out->offset_y, 0);

    // Deserialize data received.
    mpf_inp_raw(temp, partition->cx_s_len, data_out->cx);
    temp += partition->cx_s_len;
    mpf_inp_raw(temp, partition->cy_s_len, data_out->cy);
    temp += partition->cy_s_len;
    mpf_inp_raw(temp, partition->offset_x_s_len, data_out->offset_x);
    temp += partition->offset_x_s_len;
    mpf_inp_raw(temp, partition->offset_y_s_len, data_out->offset_y);
    temp = NULL;
}

/**
 * @brief computes the image partitions by using iteration_map function.
 *
 * @param image Image of the complete solution.
 * @param partition The partition of the messageReceived.
 */
void compute_image_partition(unsigned char *const image, received_partition_t *const partition) {
	mpf_t cx_copy;
    // Set cx_copy with partition->cx value.
    mpf_init_set(cx_copy, partition->cx);
    for (int i = 0; i < partition->clength; i++) {
        for (int j = 0; j < partition->clength; j++) {
            image[i * partition->clength + j] = iteration_map(partition);
			mpf_add(partition->cx, partition->cx, partition->offset_x);
        }
		mpf_set(partition->cx, cx_copy);
		mpf_sub(partition->cy, partition->cy, partition->offset_y);
    }
    // Clear mpf_t variable defined.
    mpf_clear(cx_copy);
}

/**
 * @brief compute and send image partition to receiver.
 *
 * @param rank Rank of a processor.
 * @param partition The partition of the messageReceived.
 */
void compute_and_send_image_partition(const int rank, received_partition_t partition) {
    const int partition_power = partition.clength*partition.clength;
    const uint64_t data_buffer_length = partition_power + sizeof(double);
    
    unsigned char data[data_buffer_length];
    unsigned char *image_partition = &(data[0]);
    double *entropy = (double *)(&(data[partition_power]));

    compute_image_partition(image_partition, &partition);

    *entropy = calculate_entropy_from_image(partition.clength, partition.clength, 
        image_partition, COLOR_MAP_LENGTH);

    MPI_Send(&(data[0]), data_buffer_length, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
}

/**
 * @brief calculate the entropy by using calculate_entropy_from_image function and send the entropy values.
 *
 * @param rank Rank of a processor.
 * @param partition The partition of the messageReceived.
 */
void compute_and_send_entropy_only(const int rank, received_partition_t partition) {
    const int partition_power = partition.clength*partition.clength;
    unsigned char image_partition[partition_power];
    double entropy;

    compute_image_partition(image_partition, &partition);

    entropy = calculate_entropy_from_image(partition.clength, partition.clength, 
        image_partition, COLOR_MAP_LENGTH);

    MPI_Send(&entropy, sizeof(double), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
}

/**
 * @brief starts the follower code.
 * based on the messageReceived request_type determine whether it is a send image or no image or terminate the process.
 *
 * @param rank Rank of a processor.
 */
void start_follower(const int rank) {

    unsigned char running = 1;
    follower_message_t messageReceived;
    received_partition_t partition;
    
    // Initialize a NULL-terminated list of mpf_t variables, and set their values to 0
    mpf_inits(partition.cx, partition.cy, partition.offset_x, partition.offset_y, NULL);

    while (running) {
        follower_receive_message(&messageReceived);
        switch (messageReceived.request_type) {
            case PARTITION_SEND_IMAGE:
                receive_partition_data(&messageReceived.payload.partition, &partition, rank);
                compute_and_send_image_partition(rank, partition);
                break;
            case PARTITION_NO_IMAGE:
                receive_partition_data(&messageReceived.payload.partition, &partition, rank);
                compute_and_send_entropy_only(rank, partition);
                break;
            case TERMINATE:
                running = 0;
                break;
            default:
                fprintf(stderr, "Received an unknown message. Terminating process %d.\n", rank);
                running = 0;
                break;
        }
    }

    // Free the space occupied by a NULL-terminated list of mpf_t variables.
    mpf_clears(partition.cx, partition.cy, partition.offset_x, partition.offset_y, NULL);

}
