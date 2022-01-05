/**
 * @file main.c
 * @author Dalton Caron (dpcaron@csu.fullerton.edu)
 * @author Jeevika Yarlagadda (jeevikayarlagadda@csu.fullerton.edu)
 * @brief This file defines the arguments, serializes, and deserializes the data.
 * and starts the initiator or follower based on the rank. 
 **/
#include <argp.h>
#include <assert.h>

#include "leader.h"
#include "follower.h"

// A typical argp setup. 
const char *argp_program_version = "distributed_mandelbrot 1.0.0"; // Version of the program
const char *argp_program_bug_address = "dpcaron@csu.fullerton.edu jeevikayarlagadda@csu.fullerton.edu"; // Report error to email provided
static char doc[] = "A program for performing a distributed mandelbrot zoom through MPI."; 
static char args_doc[] = "[-f filename] [-z zoom_iterations] [-s image_size_pixels] [-i mandelbrot_iterations] [-e]";
static struct argp_option options[] = {
    {"filename", 'f', "FILE", 0, "Read in mandelbrot coordinates and offsets as a text file."},
    {"zoom_iterations", 'z', "ZOOM_ITRS", 0, "How many times to zoom and recompute the mandelbrot set."},
    {"image_size_pixels", 's', "IMG_PIXELS", 0, "The size of the mandelbrot image (must be a power of 2)."},
    {"mandelbrot_iterations", 'i', "MANDL_ITRS", 0, "How many times to repeat f(z) before giving up. Must be higher the deeper the zoom."},
    {"save_each_iteration", 'e', 0, 0, "Boolean flag to save the computed mandelbrot each iteration."},
    { 0 }
};

// Define a struct for the arguments.
struct arguments {
    char *filename;
    int zoom_iterations;
    int image_size_pixels;
    int mandelbrot_iterations;
    char save_each_iteration;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = state->input;
    switch (key) {
        case 'f': arguments->filename = arg; break;
        case 'z': arguments->zoom_iterations = atoi(arg); break;
        case 's': arguments->image_size_pixels = atoi(arg); break;
        case 'i': arguments->mandelbrot_iterations = atoi(arg); break;
        case 'e': arguments->save_each_iteration = 1; break;
        case ARGP_KEY_ARG: return 0;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };

/**
 * @brief checks if the number x is power of two.
 *
 * @param x Size of a processor or image size in pixels.
 */
unsigned char is_power_of_two(unsigned long x) {
    return (x != 0) && ((x & (x - 1)) == 0);
}

int main(int argc, char *argv[]) {
    int rank, size;

    // Initiate MPI environment.
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // For debugging purpose it prints the rank and process id of a process.
    #if ENABLE_DEBUGGER_TRICK > 0
        printf("rank %d has process id %u\n", rank, getpid());
        wait_for_debugger(rank, 1);
    #endif

    // Transforms all mpf_init(X) calls into mpf_init2(X, size_cx_cy). 
    mpf_set_default_prec(size_cx_cy);

    if (size < 5 || !is_power_of_two(size-1)) {
        printf("Number of processes must be perfectly divisible by 4 when 1 is subtracted.\n");
    } else if (rank == 0) {
        // Parse command line arguments only in the leader. 
        struct arguments arguments;
        arguments.filename = NULL;
        arguments.zoom_iterations = 4;
        arguments.image_size_pixels = 512;
        arguments.mandelbrot_iterations = 2000;
        arguments.save_each_iteration = 0;

        argp_parse(&argp, argc, argv, 0, 0, &arguments);

        if (!is_power_of_two(arguments.image_size_pixels)) {
            printf("Image size must be a power of 2: got %d; using 512 default instead.\n", arguments.image_size_pixels);
            arguments.image_size_pixels = 512;
        }

        const double c_width_x = 3.0;
        const double c_height_y = 2.0;

        mandelbrot_settings_t settings;
		// Assigning the default mandelbrot structure values.
        mpf_inits(settings.z_start_x, settings.z_start_y, settings.offset_x, settings.offset_y, NULL);
        mpf_set_d(settings.z_start_x, -2.0);
        mpf_set_d(settings.z_start_y, 1.0);
        mpf_set_d(settings.offset_x, CALCULATE_C_OFFSET(c_width_x, arguments.image_size_pixels));
        mpf_set_d(settings.offset_y, CALCULATE_C_OFFSET(c_height_y, arguments.image_size_pixels));
        settings.max_iterations = arguments.mandelbrot_iterations;

		// If file name is not null, read the parameters from the file specified.
        if (arguments.filename != NULL) {
            mpf_t cx, cy, offset_x, offset_y;
            mpf_inits(cx, cy, offset_x, offset_y, NULL);

            char error = read_mandelbrot_partition_parameters_from_file(&cx, &cy, 
                &offset_x, &offset_y, arguments.filename);

            // Use the defaults if we cannot read the file. Notify the user of this. 
            if (error == 0) {
				mpf_set(settings.z_start_x, cx);
				mpf_set(settings.z_start_y, cy);
				mpf_set(settings.offset_x, offset_x);
				mpf_set(settings.offset_y, offset_y);
            } else {
                fprintf(stderr, "Failed to read parameter file. Using default settings.\n");
            }
            mpf_clears(cx, cy, offset_x, offset_y, NULL);
        }

        printf("Using settings cx ");
        mpf_out_str(stdout, 10, N_DIGITS, settings.z_start_x);
        printf(" cy ");
        mpf_out_str(stdout, 10, N_DIGITS, settings.z_start_y);
        printf(" offset_x ");
        mpf_out_str(stdout, 10, N_DIGITS, settings.offset_x);
        printf(" offset_y ");
        mpf_out_str(stdout, 10, N_DIGITS, settings.offset_y); 
        printf("\n");

		// Start the initiator as rank is zero.
        start_initiator(size, &settings, arguments.image_size_pixels, 
            arguments.zoom_iterations, arguments.save_each_iteration);

        // Deallocate all mpf types. 
        mpf_clears(settings.z_start_x, settings.z_start_y, settings.offset_x, 
            settings.offset_y, NULL);
    } else {
		// Start the follower with rank.
        start_follower(rank);
    }

    MPI_Finalize();
    return 0;
}
