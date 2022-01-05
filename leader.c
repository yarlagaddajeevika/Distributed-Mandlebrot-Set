/**
 * @file leader.c
 * @author Dalton Caron (dpcaron@csu.fullerton.edu)
 * @author Jeevika Yarlagadda (jeevikayarlagadda@csu.fullerton.edu)
 **/
#include "leader.h"

#include <math.h>
#include <assert.h>
#include <time.h>

/**
 * @brief Copies the sub-solution of image partition cell into 
 * image into the rectangle (px, py, px + clength, py + clength).
 * 
 * @param px Cell partition origin x in image. 
 * @param py Cell partition origin y in image. 
 * @param clength Cell partition width and height. 
 * @param image_size The size of the image. 
 * @param image Image of the complete solution. 
 * @param cell Image partition to merge into image. 
 */
void copy_solution_cell_into_solution(int px, int py, int clength, int image_size, 
    unsigned char image[image_size][image_size], unsigned char *cell) {
    int x = 0, y = 0;
    for (int i = px; i < px + clength; i++) {
        for (int j = py; j < py + clength; j++) {
            image[i][j] = cell[y * clength + x];
            x++;
        }
        x = 0;
        y++;
    }
}
 
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
    mpf_t cx, mpf_t cy, mpf_t temporary_float) {
    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Status status;
    uint16_t i, j, p = 1;

    // Set cx and cy variables
	mpf_set(cx, settings->z_start_x);
	mpf_set(cy, settings->z_start_y);
	
    follower_message_t messageToSend;

    for (i = 0; i < image_size; i += partition_size) {
        for (j = 0; j < image_size; j += partition_size) {

            // Populate partition information to send to the followers. 
            mandelbrot_partition_t *partition = (mandelbrot_partition_t *)(&(partitions[p-1]));
            partition->i = i;
            partition->j = j;
			// Set the new values of cx and cy.
			mpf_set(partition->cx, cx);
			mpf_set(partition->cy, cy);
            partition->clength = partition_size;

            size_t cx_s_len, cy_s_len, offset_x_s_len, offset_y_s_len;
            size_t buffer_len = N_DIGITS * 4;
            size_t length_left = buffer_len;
            size_t written_length = 0;

            char buffer[buffer_len];
            char *temp = buffer;

            // Serialize the cx,cy,offset_x,offset_y lengths
            cx_s_len = mpf_out_raw(temp, buffer_len, partition->cx);
            written_length += cx_s_len;
            temp += cx_s_len;
            length_left -= cx_s_len;

            cy_s_len = mpf_out_raw(temp, length_left, partition->cy);
            written_length += cy_s_len;
            temp += cy_s_len;
            length_left -= cy_s_len;

            offset_x_s_len = mpf_out_raw(temp, length_left, settings->offset_x);
            written_length += offset_x_s_len;
            temp += offset_x_s_len;
            length_left -= offset_x_s_len;

            offset_y_s_len = mpf_out_raw(temp, length_left, settings->offset_y);
            written_length += offset_y_s_len;

            // Populate message to send with data. 
            messageToSend.request_type = (receive_data) ? PARTITION_SEND_IMAGE : PARTITION_NO_IMAGE;
            messageToSend.payload.partition.clength = partition->clength;
            messageToSend.payload.partition.max_iterations = settings->max_iterations;
            messageToSend.payload.partition.cx_s_len = cx_s_len;
            messageToSend.payload.partition.cy_s_len = cy_s_len;
            messageToSend.payload.partition.offset_x_s_len = offset_x_s_len;
            messageToSend.payload.partition.offset_y_s_len = offset_y_s_len;

            // Send partition information to next worker. 
            MPI_Isend(&messageToSend, sizeof(follower_message_t), MPI_BYTE, p, 0, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);

            MPI_Request request = MPI_REQUEST_NULL;
            MPI_Status status;

            MPI_Isend(buffer, written_length, MPI_BYTE, p, 0, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);

            // Save coordinates assigned to process p for combining solutions later. 
            p++;

			mpf_mul_ui(temporary_float, settings->offset_x, (unsigned long) partition_size);
			mpf_add(cx, cx, temporary_float);
        }

		mpf_set(cx, settings->z_start_x);
		mpf_mul_ui(temporary_float, settings->offset_y, (unsigned long) partition_size);
		mpf_sub(cy, cy, temporary_float);
    }

}

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
    int *const successor_process) {
    MPI_Status status;
    const uint64_t receive_buffer_length = partition_size*partition_size + sizeof(double);

    unsigned char receive_buffer[receive_buffer_length];
    unsigned char *image_partition = &(receive_buffer[0]);
    double *entropy = (double *)(&(receive_buffer[partition_size*partition_size]));
    
    double max_entropy = -__DBL_MAX__; 
    uint16_t process_with_max_entropy;

    for (int p = 1; p < size; p++) {
        MPI_Recv(&(receive_buffer[0]), receive_buffer_length, MPI_BYTE, p, 0, MPI_COMM_WORLD, &status);
        copy_solution_cell_into_solution(partition_map[p-1].i, partition_map[p-1].j, partition_size, image_size, image, image_partition);
        if (max_entropy < *entropy) {
            max_entropy = *entropy;
            process_with_max_entropy = p;
        }
    }

    *successor_process = process_with_max_entropy;
}

/**
 * @brief Receives the entropy from the follower. 
 *
 * @param size The total amount of processes in this communicatior.
 * @param partition_map A one-to-one map of processes IDs of rank-1 to mandelbrot_partition_t.
 * @param successor_process The process with the maximum entropy.
 */
void receive_entropies(const int size, int *const successor_process) {
    MPI_Status status;
    double entropy;
    double max_entropy = -__DBL_MAX__; 
    uint16_t process_with_max_entropy;

    for (int p = 1; p < size; p++) {
        MPI_Recv(&entropy, sizeof(double), MPI_BYTE, p, 0, MPI_COMM_WORLD, &status);
        if (max_entropy < entropy) {
            max_entropy = entropy;
            process_with_max_entropy = p;
        }
    }

    *successor_process = process_with_max_entropy;
}

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
void perform_mandelbrot_iteration(mandelbrot_settings_t *const settings, mandelbrot_partition_t * partition_map,
    const int size, const int image_size, const unsigned char save_image, mpf_t cx, mpf_t cy, mpf_t temporary_float) {
    const int cell_count = floor(sqrt(size));
    const int clength = image_size / cell_count;

    populate_and_send_partitions(save_image, settings, partition_map, clength, 
        image_size, cx, cy, temporary_float);

    // Technically, we do not need to assemble the full image until the final iteration.
    // Removing this and performing the save on the last iteration would be faster, but I 
    // do enjoy watching the zoom progress. I will make an saving optional. 
    int successor_process;
    if (save_image) {
        time_t now = time(NULL);
		
        struct tm tm = *localtime(&now);
        char time_buffer[20];
        sprintf(time_buffer, "%04d-%02d-%02d_%02d-%02d-%02d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
        char result_buffer[40];
        sprintf(result_buffer, "result%s.bmp", time_buffer);
        char param_buffer[40];
        sprintf(param_buffer, "params%s.bin", time_buffer);

        unsigned char image[image_size][image_size];
        receive_and_populate_solution(image_size, image, clength, size, partition_map, &successor_process);
        save_mandelbrot_image(image_size, image_size, &(image[0][0]), result_buffer);
        write_mandelbrot_partition_parameters_to_file(settings->z_start_x, settings->z_start_y, settings->offset_x, settings->offset_y, param_buffer);
        printf("Solution image written to %s and partition parameters written to %s on last iteration.\n", result_buffer, param_buffer);
    } else {
        receive_entropies(size, &successor_process);
    }

    // Select a subpartition of the image to calculate next by changing the settings. 
	mpf_set(settings->z_start_x, partition_map[successor_process-1].cx);
	mpf_set(settings->z_start_y, partition_map[successor_process-1].cy);
    mpf_div_ui(settings->offset_x, settings->offset_x, (unsigned long) sqrt(size-1));
    mpf_div_ui(settings->offset_y, settings->offset_y, (unsigned long) sqrt(size-1));
}

/**
 * @brief Terminate the followers by sending the MPI message with request type as TERMINATE.
 * 
 * @param size The size of the MPI communicator.
 */
void terminate_all_followers(const int size) {
    MPI_Request request;
    follower_message_t terminate_message;
    terminate_message.request_type = TERMINATE;
	// Send the terminate message to each follower.
    for (int p = 1; p < size; p++) {
        MPI_Isend(&terminate_message, sizeof(follower_message_t), MPI_BYTE, p, 0, MPI_COMM_WORLD, &request);
    }
}

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
    const int iterations, const int save_each_iteration) {

    mpf_t cx, cy, temporary_float;
    mpf_inits(cx, cy, temporary_float, NULL);

    // Mandelbrot partition storage is initialized here to make less mpf_init and clear calls. 
    // Processes may be received in a different order, so we save coordinates of a partition 
    // for each process in the partition map. It is indexed to exclude the leader process.  
    mandelbrot_partition_t partition_map[size-1];
    for (int p = 0; p < size-1; p++) {
        mpf_inits(partition_map[p].cx, partition_map[p].cy, NULL);
    }

    int i;
    for (i = 1; i < iterations; i++) {
        perform_mandelbrot_iteration(settings, partition_map, size, image_size, save_each_iteration, cx, cy, temporary_float);
        printf("Iteration %d completed.\n", i);
    }
	//to save the final image send the save_image parameter as 1
    perform_mandelbrot_iteration(settings, partition_map, size, image_size, 1, cx, cy, temporary_float);
    printf("Iteration %d completed.\n", i);

    for (int p = 0; p < size-1; p++) {
        mpf_clears(partition_map[p].cx, partition_map[p].cy, NULL);
    }

    mpf_clears(cx, cy, temporary_float, NULL);

	// Once all the iterations are completed, terminate the followers. 
    terminate_all_followers(size);
}