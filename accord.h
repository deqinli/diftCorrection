#ifndef ACCORD_H
#define ACCORD_H


#include <stdlib.h>
#include <string.h>
#include <tiffio.h>
#include <fftw3.h>
#include <math.h>
#include <assert.h>

#define ACC_IM_TYPE double
#define ACC_DIM_TYPE unsigned long int
#define ACC_SIG_DIM_TYPE long int
#define IM_8BIT_TYPE unsigned char
#define IM_8BIT_MAX 255


enum {
    IM_INITIALIZED=0,
    IM_NOT_INITIALIZED
};


#define IM_SIZES_NO_MATCH 0
#define IM_SIZES_MATCH 1

#define IM_DESTRUCTIVE 1
#define IM_NON_DESTRUCTIVE 0

#define IM_SIZE (image->sizex * image->sizey)

extern void (*status_callback)(int i, char *message); // this has to be initialized


/**
 * \brief Image structure.
 *
 * This structure contains everything related to an image including pixel size, its Fourier
 * images, and FFTW plans.
 */
typedef struct im_struct{
    ACC_IM_TYPE *buffer; ///< pointer to the image buffer
    fftw_complex *fft_buffer; ///< pointer to the complex buffer for the Fourier image
    unsigned int *weight_buffer; ///< this buffer keeps track, how many contributions are in each pixel. NULL if not defined.
    ACC_DIM_TYPE sizex; ///< width of the image in pixels
    ACC_DIM_TYPE sizey; ///< height of the image in pixels
    fftw_plan fft_plan; ///< FFTW plan handle for the forward FFT
    fftw_plan ifft_plan; ///< FFT plan handle for the inverse FFT
    ACC_SIG_DIM_TYPE maskx; ///< x-coordinate of the mask (needed because of invalid data due to FFT-based shift)
    ACC_SIG_DIM_TYPE masky; ///< y-coordinate of the mask
    char fourier_ok; ///< indicator whether the Fourier image is valid
} t_im_struct;


/**
 * Error-numbers
 */
enum accord_errors {
    ACC_ERR_NO_ERROR=0,///< This means that everything is OK.
    ACC_ERR_TOO_MANY_ITERATIONS, ///< Newton maximum search hit the maximum iterations limit.
    ACC_ERR_TIFF_OPEN, ///< Cannot open TIFF file.
    ACC_ERR_TIFF_WRITE, ///< Writing to TIFF file failed.
    ACC_ERR_SIZES_MISMATCH, ///< Image sizes do not match.
    ACC_ERR_NO_FFT_BUFFER ///< FFT buffer is not initialized

};

/**
 * \brief This function inputs two images and finds their displacement vector.
 *
 * \param image1 pointer to the first image structure
 * \param image2 pointer to the second image structure
 * \param filter_radius radius of the frequency filter. (Try 30)
 * \param shx pointer to a double, where this function stores the x-component of the displacement vector
 * \param shy analogously for the y-coordinate
 * \return error number: #accord_errors (ACC_ERR_NO_ERROR = 0 when OK)
 */
char find_displacement(t_im_struct *image1, t_im_struct *image2, double filter_radius, double *shx, double *shy);

/**
 * \brief Adds image2 to image1 with correction of drift.
 *
 * \param image1 pointer to the first image structure
 * \param image2 pointer to the second image structure
 * \param destructive Provide IM_NON_DESTRUCTIVE or IM_DESTRUCTIVE. If
 * IM_DESTRUCTIVE, the image2 will be shifted and thus changed.
 * \param shx pointer to a double, where this function stores the x-component of the displacement vector. Set to NULL if not interested in displacement info.
 * \param shy analogously for the y-coordinate
 * \return error number: #accord_errors (ACC_ERR_NO_ERROR = 0 when OK)
 * \sa imageio
 */
char add_with_correction(t_im_struct *image1, t_im_struct *image2, char destructive, double *shx_out, double *shy_out);

void print_error(int err);
void print_status(int i, char *message);
// @}



// ****                             IMAGE_IO                            ****
// *************************************************************************

/** \defgroup imageio Image Operations
 * \brief Functions in this group provide basic operations with images.
@{
 */

/**
 * \brief This function initializes an empty image structure.
 *
 * The dimensions must be
 * provided. The image buffer is allocated here. Thus, the #destroy_image function
 * must be called, when the image is no more needed.
 * \param image pointer to an image structure.
 * \param sizex width of the image in pixels
 * \param sizey height of the image in pixels
 * \return error number: #accord_errors (ACC_ERR_NO_ERROR = 0 when OK)
 */
void initialize_image(t_im_struct *image, ACC_DIM_TYPE sizex, ACC_DIM_TYPE sizey);

/**
 * \brief Initializes weight buffer.
 *
 * When the image is shifted, part of it is invalid. This invalid part is not
 * added to the final image. However, in order to keep the gray-levels correct,
 * it's necessary to account for the number of contributions in each pixel.
 * \param image pointer to an image structure.
 */
void initialize_weight_counting(t_im_struct *image);

/**
 * \brief Loads image from an 8-bit TIFF file.
 *
 * \param image pointer to an image structure.
 * \param filename The file name.
 * \param initialized IM_INITIALIZED - overwrites an initialized image\n
 * IM_NOT_INITIALIZED - properly initializes the image (including the right size)
 * \return error number: #accord_errors (ACC_ERR_NO_ERROR = 0 when OK)
 */
char load_image(t_im_struct *image, char *filename, char initialized);

/**
 * \brief Saves image to an 8-bit TIFF file.
 *
 * \param image pointer to an image structure.
 * \param filename The file name.
 * \return error number: #accord_errors (ACC_ERR_NO_ERROR = 0 when OK)
 */
char save_image(t_im_struct *image, char *filename);

/**
 * \brief Destroys the image structure.
 *
 * Frees all buffers.
 * Call this function, when the image is no more needed.
 * \param image pointer to an image structure.
 */
void destroy_image(t_im_struct *image);

/**
 * \brief Multiplies each pixel of an image by a coefficient.
 *
 * \param image pointer to an image structure.
 * \param coef a real coefficient
 */
void multiply_image(t_im_struct *image, double coef);

/**
 * \brief Normalizes the image.
 *
 * This means that after this operation, the minimum value
 * is 0 and the maximum is 1.
 * \param image pointer to an image structure.
 */
void normalize_image(t_im_struct *image);


/**
 * \brief Normalizes the image with respect to information weight of pixels.
 *
 * Like with #normalize_image, the minimum value
 * is 0 and the maximum is 1 after normalization. This function only observes
 * pixels with full information weight.
 *
 * When the image is composed from frames with correction of drift, the frames are shifted
 * before summing. This invalidates part of the frame. For example, if the frame
 * is shifted by 3 pixels to the left, the right-most 3 pixels in all rows of
 * this frame do not contain valid data. This is accounted for in the
 * im_struct::weight_buffer
 * \param image pointer to an image structure.
 * \param full_weight_only 0 for standard normalization <br>
 *   1 for normalization, where only allways acquired pixels count (full weight)
 */
void normalize_image_w(t_im_struct *image, char full_weight_only);

/**
 * \brief Finds the maximum grey-level and the corresponding index.
 *
 * \param image pointer to an image structure.
 * \param i pointer to the index variable, which will be overwritten by this
 * function. If you're not interested in the index, set it to NULL
 * \param gl pointer to the grey-level variable. This function stores the found
 * maximum grey-level to it. If you're not interested, set this to NULL
 */
void find_maximum(t_im_struct *image, unsigned long int *i, ACC_IM_TYPE *gl);

/**
 * \brief Converts from two-dimensional coordinates to one-dimensional array index.
 *
 * See
 * #i_to_xy
 * \param x x-coordinate
 * \param y y-coordinate
 * \param image pointer to an image structure.
 * \return index i corresponding to the x and y coordinates
 */
unsigned long int xy_to_i(ACC_DIM_TYPE x, ACC_DIM_TYPE y, t_im_struct *image);

/**
 * \brief Converts from one-dimensional array index to two-dimensional coordinates.
 *
 * See #xy_to_i
 * \param image pointer to an image structure.
 * \param i the array index to be converted
 * \param x pointer to the x-coordinate (overwritten by this function)
 * \param y pointer to the y-coordinate  (overwritten by this function)
 */
void i_to_xy(unsigned long int i, ACC_DIM_TYPE *x, ACC_DIM_TYPE *y, t_im_struct *image);

/**
 * \brief Checks if the sizes of two images match.
 *
 * \param image1 pointer to the first image structure.
 * \param image2 pointer to the second image structure.
 * \return IM_SIZES_MATCH if the sizes match\n
 * IM_SIZES_NO_MATCH if the sizes do not match.
 */
char check_sizes(t_im_struct *image1, t_im_struct *image2);

/**
 * \brief Sums two images
 *
 * Adds (image_s + image_d) and puts the result to the first one
 * (image_d).
 * \param image_d pointer to the destination image structure. (This structure
 * will be overwritten.)
 * \param image_s pointer to the source image structure.
 * \return error number: #accord_errors (ACC_ERR_NO_ERROR = 0 when OK)
 */
char add_images(t_im_struct *image_d, t_im_struct *image_s);

/**
 * \brief Copies image_s to the image_d
 *
 * \param image_d pointer to the destination image structure. (This structure
 * will be overwritten by the result.)
 * \param image_s pointer to the source image structure.
 * \param initialized IM_INITIALIZED - overwrites an initialized image\n
 * IM_NOT_INITIALIZED - first properly initializes the image (including the right size)
 * \return error number: #accord_errors (ACC_ERR_NO_ERROR = 0 when OK)
 */
char copy_image(t_im_struct *image_d, t_im_struct *image_s, char initialized);

/**
 * \brief copies image including the weight buffer.
 *
 * see #copy_image
 */
char copy_image_with_weight(t_im_struct *image_d, t_im_struct *image_s, char initialized);

/**
 * \brief Crops the given image.
 *
 * \param image pointer to an image structure.
 * \param left number of pixels to crop off from left
 * \param right number of pixels to crop off from right
 * \param top number of pixels to crop off from top
 * \param bottom number of pixels to crop off from bottom
 */

char crop_image(t_im_struct *image, ACC_DIM_TYPE left, ACC_DIM_TYPE right, ACC_DIM_TYPE top, ACC_DIM_TYPE bottom);
/*@}*/


// ****                             FOURIER                             ****
// *************************************************************************

/* \defgroup fourier Fourier-Related Functions
 * Functions in this group provide the core functionality of this library.
 * Contains the functions for calculating forward and inverse Fourier transform,
 * frequency filtering, and cross-correlation.
@{
 */

/**
 * \brief Creates FFTW plans
 *
 * and stores them in the #im_struct::fft_plan and
 * #im_struct::ifft_plan. These plans
 * are needed for the FFTW library. For the first time, this function may take
 * longer time, nex time it should be quicker, because the FFTW library already
 * has its "wisdom".
 * \param image pointer to an image structure.
 * \return error number: #accord_errors (ACC_ERR_NO_ERROR = 0 when OK)
 * \sa imageio
 */
char create_fft_plans(t_im_struct *image);

/**
 * \brief Destroys the FFT plans.
 *
 * This cleans forward ind inverse FFTW plans, frees the
 * fft_buffer and assigns its pointer to NULL.
 * \param image pointer to an image structure.
 * \sa im_struct
 */
void destroy_fft_plans(t_im_struct *image);

/**
 * \brief Calculates the forward FFT.
 *
 * The result is stored in the
 * #im_struct::fft_buffer. This function also sets the #im_struct::fourier_ok to
 * 1.\n
 * <b>Attention:</b> You are responsible for clearing it, when you write to
 * either #im_struct::fft_buffer or #im_struct::buffer.
 * \param image pointer to an image structure.
 */
void calculate_fft(t_im_struct *image);

/**
 * \brief Calculates the inverse FFT.
 *
 * Analogous to #calculate_fft.
 * \param image pointer to an image structure.
 * \return error number: #accord_errors (ACC_ERR_NO_ERROR = 0 when OK)
 */
char calculate_ifft(t_im_struct *image);

/**
 * \brief Conjugates the first complex value (x) and multiplies it with the second one
 * (y).
 *
 * The result is written to x.
 * \n What this function does can be written as:
 * \f$x = \bar{x} * y\f$, where the \f$\bar{x}\f$ denotes the conjugated value
 * of \f$x\f$.
 * \param x pointer to a complex value
 * \param y pointer to a complex value
 */
void conjugate_multiply(fftw_complex *x, fftw_complex *y);

/**
 * \brief Conjugates and multiplies two images with application of a frequency
 * filter.
 *
 * Calls #conjugate_multiply for every pixel in the image and multiplies it with
 * the number returned by <i>filter_callback</i> function. This callback
 * function must be provided. For an example of such a function, see
 * #cylinder_filter.
 * \param im1 pointer to an image structure.
 * \param im2 pointer to an image structure.
 * \param filter_callback a pointer to a call-back function. See #cylinder_filter or #no_filter functions for examples.
 * \param user_par a pointer to a user parameter. This is forwarded to the
 * filter_callback function.
 * \return error number: #accord_errors (ACC_ERR_NO_ERROR = 0 when OK)
 */
char conjugate_multiply_images_with_filter(t_im_struct *im1, t_im_struct *im2,
                                           double (*filter_callback)(t_im_struct *im, unsigned long int i, void *user_par),
                                           void *user_par);

/**
 * \brief The simplest frequency-filter call-back function
 * \param im pointer to an image structure.
 * \param i array index
 * \param user_par parameter forwarded from the
 * #conjugate_multiply_images_with_filter. Ignored.
 * \sa cylinder_filter
 */
double no_filter(t_im_struct *im, unsigned long int i, void *user_par);

/**
 * \brief Circular low-pass-frequency-filter call-back function
 * \param im pointer to an image structure.
 * \param i array index
 * \param user_par parameter forwarded from the
 * #conjugate_multiply_images_with_filter. In this case it is a pointer to the
 * radius if the low-pass filter. Its type is <i>float</i>.
 * \sa no_filter
 */
double cylinder_filter(t_im_struct *im, unsigned long int i, void *user_par);

double BW_filter(t_im_struct *im, unsigned long int i, void *user_par);

/**
 * \brief Conversion of x and y coordinate indexes to i with quadrant-shift.
 *
 * The x varies from
 * (- im_struct::sizex/2) to (im_struct::sizex/2). y varies analogously. This is
 * translated to a one-dimensional array index i varying from 0 to (im_struct::sizex * im_struct::sizey -1).
 * \param x the signed x-coordinate
 * \param y the signed y-coordinate
 * \param image pointer to an image structure.
 * \return array index i
 */
unsigned long int xy_to_i_fourier_shift(ACC_SIG_DIM_TYPE x, ACC_SIG_DIM_TYPE y, t_im_struct *image);

/**
 * \brief Conversion of array index i to x and y coordinate indexes with quadrant-shift.
 *
 * see #xy_to_i_fourier_shift.
 * \param i array index
 * \param x pointer to the signed x-coordinate
 * \param y pointer to the signed y-coordinate
 * \param image pointer to an image structure.
 */
void i_to_xy_fourier_shift(unsigned long int i, ACC_SIG_DIM_TYPE *x, ACC_SIG_DIM_TYPE *y, t_im_struct *image);

/**
 * \brief Shifts image with a sub-pixel accuracy
 *
 * \param image pointer to an image structure.
 * \param shiftx x-component of the shift vector
 * \param shifty y-component of the shift vector
 */
char shift_image(t_im_struct *image, double shiftx, double shifty);

// ****                   LINEAR EQUATION SOLVER                        ****
// *************************************************************************

void cholesky(double *a, double *l, int dim);
void lesolve(double *a, double *y, double *x, int dim);

// ****                        POLYNOMIAL FIT                           ****
// *************************************************************************

void grab_fitting_area(t_im_struct *image, unsigned long int i, ACC_DIM_TYPE area_size, ACC_SIG_DIM_TYPE *x, ACC_SIG_DIM_TYPE *y, ACC_IM_TYPE *z);
void newton_maximum_search(double *sq_delta, double *a, double *xin, double *yin);
char cubic_fit_maximum(ACC_SIG_DIM_TYPE *x, ACC_SIG_DIM_TYPE *y, ACC_IM_TYPE *z, unsigned int count, double *xs, double *ys);
char find_maximum_subpixel(t_im_struct *image, double *xs, double *ys);
/*@}*/


#endif // ACCORD_H
