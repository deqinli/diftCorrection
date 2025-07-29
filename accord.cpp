/*
    This file is a part of:
    Charged-Particle-Microscope Image Composition with Correction of Drift (ACCORD)

    2010  Petr Cizmar @ National Institute of Standards and Technology
    E-mail: petr.cizmar@nist.gov

    As this software was developed as part of work done by the United States
    Government, it is not subject to copyright, and is in the public domain.
    Note that according to GNU.org public domain is compatible with GPL.

 */

#define CDIM 10 // number of coefficients of the 2-dim cubic polynomial (do not change)

#include "accord.h"
#include <string>
#include <QtMath>
void (*status_callback)(int i, char *message) = print_status;

char find_displacement(t_im_struct *image1, t_im_struct *image2, double filter_radius, double *shx, double *shy){/*{{{*/
    // finds displacement of the image2 with respect to image1 with s sub-pixel accuracy
    // result is output to shx and shy
    // *************
    // ATTENTION: destroys the image1 buffer - you can find the cross-correlation
    // function there instead
    // *************

    char err = 0;

    double radius = filter_radius;

    if ((err = conjugate_multiply_images_with_filter(image1, image2, BW_filter, (void *)&radius))) return err;// arg3 choices: cylinder_filter,BW_filter
    if ((err = calculate_ifft(image1))) return err;
    normalize_image(image1);
    if ((err = find_maximum_subpixel(image1, shx, shy))) return err;
    return ACC_ERR_NO_ERROR;
}/*}}}*/

char add_with_correction(t_im_struct *image1, t_im_struct *image2, char destructive, double *shx_out, double *shy_out){/*{{{*/
    // adds image2  to the image1 with correction of drift
    // shx and shy are outputs for the displacement coordinates. set to NULL if
    // not interested.

    char err = ACC_ERR_NO_ERROR;

    t_im_struct helper1, helper2;
    t_im_struct *h2;

    if ((err = copy_image(&helper1, image1, IM_NOT_INITIALIZED))) return err;
    normalize_image(&helper1);

    if (destructive == IM_NON_DESTRUCTIVE) {
        if ((err = copy_image(&helper2, image2, IM_NOT_INITIALIZED))) {
            destroy_image(&helper1);
            return err;
        }
        h2 = &helper2;
    } else
        h2 = image2;

    double shx, shy;
    err = find_displacement(&helper1, h2, 30, &shx, &shy); // fixed radius.
    if (err) return err;

    if (shx_out) *shx_out = shx;
    if (shy_out) *shy_out = shy;

    shift_image(h2, shx, shy);

    add_images(image1, h2);

    destroy_image(&helper1);
    if (destructive == IM_NON_DESTRUCTIVE) destroy_image(&helper2);
    return ACC_ERR_NO_ERROR;
}/*}}}*/

void print_error(int err){/*{{{*/
    std::string err_string;
    switch (err){
    case ACC_ERR_NO_ERROR:
        err_string="All OK";
        break;
    case ACC_ERR_TOO_MANY_ITERATIONS:
        err_string="Too many iterations in the Newton solver";
        break;
    case ACC_ERR_TIFF_WRITE:
        err_string="Cannot write TIFF file";
        break;
    case ACC_ERR_TIFF_OPEN:
        err_string="Cannot open TIFF file";
        break;
    case ACC_ERR_SIZES_MISMATCH:
        err_string="Image sizes mismatch";
        break;
    case ACC_ERR_NO_FFT_BUFFER:
        err_string="FFT buffer not initialized";
        break;
    default:
        err_string="Unknown error - This is probably a bug. Please report it.";
    }
    fprintf(stderr, "Error: %s.\n", err_string.c_str());
    return;
}/*}}}*/

void print_status(int i, char *message){/*{{{*/
    if (i) printf("#%d ", i);
    if (message) printf("%s\n", message);
}/*}}}*/

// vim: cindent


void initialize_image(t_im_struct *image, ACC_DIM_TYPE sizex, ACC_DIM_TYPE sizey){/*{{{*/
    memset(image, 0, sizeof(t_im_struct)); // zeroes everything
    image->buffer = (ACC_IM_TYPE *) malloc(sizex * sizey * sizeof(ACC_IM_TYPE));
    assert(image->buffer); // Fatal error - malloc failed
    image->sizex = sizex;
    image->sizey = sizey;
}/*}}}*/

char load_image(t_im_struct *image, char *filename, char initialized){/*{{{*/
    // loads the image from an 8-bit TIFF file. if the image is not initialized,
    // the last parameter should be set to 1
    uint32 i, j, sizex_h, sizey_h;
    int a = 0;
    a = a / 2;

    image->fourier_ok = a; // loading a TIFF file invalidates the possibly existing data in fft_buffer
    unsigned char *helper_buffer;
    assert(filename!=nullptr&&"load_image::filename is nullptr");
    TIFF *tif = TIFFOpen(filename, "r");
    if (tif) {
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &sizex_h);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &sizey_h);

        if (initialized == IM_NOT_INITIALIZED) initialize_image(image, sizex_h, sizey_h); // initialize if not initialized
        assert(helper_buffer = (unsigned char*) malloc(sizeof(uint32) * sizex_h * sizey_h)); // 8-bit images only

        TIFFReadRGBAImage(tif, sizex_h, sizey_h, (uint32 *)helper_buffer, 0);
        for (j=0; j < sizey_h; j++)
            for (i=0; i < sizex_h; i++) {
                uint32 v_source;
                v_source = 4 * (sizex_h * j + i);
                uint32 v_dest;
                v_dest = sizex_h * (sizey_h - 1 - j) + i;
                image->buffer[v_dest]= (ACC_IM_TYPE) helper_buffer[v_source] / IM_8BIT_MAX; //copy every 4. value to the grey-scale buffer
            }
        free(helper_buffer);

        //copy the sizes to the image structure
        image->sizex = sizex_h;
        image->sizey = sizey_h;

    } else {
        return ACC_ERR_TIFF_OPEN;
    }
    TIFFClose(tif);
    return ACC_ERR_NO_ERROR;
}/*}}}*/

char save_image(t_im_struct *image, char *filename){/*{{{*/

    IM_8BIT_TYPE *outbuf = (IM_8BIT_TYPE *) malloc(sizeof(IM_8BIT_TYPE) * IM_SIZE);
    assert(outbuf);
    uint32 i;
    for (i=0; i < IM_SIZE; i++) outbuf[i] = image->buffer[i] * IM_8BIT_MAX;
    TIFF *output;
    if((output = TIFFOpen(filename, "w")) == NULL) return ACC_ERR_TIFF_OPEN;
    // Write the tiff tags to the file
    TIFFSetField(output, TIFFTAG_IMAGEWIDTH, image->sizex);
    TIFFSetField(output, TIFFTAG_IMAGELENGTH, image->sizey);

    //TIFFSetField(output, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE); --
    //compression disabled ^
    TIFFSetField(output, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(output, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(output, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(output, TIFFTAG_BITSPERSAMPLE, sizeof(IM_8BIT_TYPE) * 8); // 8-bit images for now
    TIFFSetField(output, TIFFTAG_SAMPLESPERPIXEL, 1);

    // Actually write the image
    if(TIFFWriteEncodedStrip(output, 0, outbuf, IM_SIZE * sizeof(IM_8BIT_TYPE)) == 0)
        return ACC_ERR_TIFF_WRITE;
    TIFFClose(output);
    free(outbuf);

    return ACC_ERR_NO_ERROR;
}/*}}}*/


void destroy_image(t_im_struct *image){/*{{{*/
    free(image->buffer);
    destroy_fft_plans(image);
    image->fourier_ok = 0;
    if (image->weight_buffer) free (image->weight_buffer);
}/*}}}*/

#define sx image->sizex
#define sy image->sizey
void initialize_weight_counting(t_im_struct *image){/*{{{*/
    if (!image->weight_buffer) image->weight_buffer = (unsigned int *) malloc( sx * sy * sizeof(unsigned int));
    unsigned long int i;
    for (i=0; i < sx*sy; i++) image->weight_buffer[i] = 1;
}/*}}}*/

void multiply_image(t_im_struct *image, double coef){/*{{{*/
    unsigned long int i;
    for (i=0; i<IM_SIZE; i++) image->buffer[i] *= coef;
    image->fourier_ok = 0;
}/*}}}*/

void normalize_image(t_im_struct *image){/*{{{*/
    normalize_image_w(image, 0);
}/*}}}*/

void normalize_image_w(t_im_struct *image, char full_weight_only){/*{{{*/
    unsigned long int i;
    unsigned int max_weight;

    ACC_IM_TYPE max_value = image->buffer[0]; // initialization
    ACC_IM_TYPE min_value = image->buffer[0]; // initialization

    if (image->weight_buffer) {
        max_weight = image->weight_buffer[0]; // initialization as well
        for (i=0; i<IM_SIZE; i++) if (image->weight_buffer[i] > max_weight) max_weight = image->weight_buffer[i];
        for (i=0; i<IM_SIZE; i++) image->buffer[i] *= (ACC_IM_TYPE) max_weight / image->weight_buffer[i];
    }


    for (i=0; i<IM_SIZE; i++)
        if ((!image->weight_buffer) || (!full_weight_only) || (image->weight_buffer[i] == max_weight)) {
            if (image->buffer[i] > max_value) max_value = image->buffer[i];
            if (image->buffer[i] < min_value) min_value = image->buffer[i];
        }

    if (image->weight_buffer) {
        free(image->weight_buffer);
        image->weight_buffer = NULL;
    }

    if (max_value - min_value)
        for (i=0; i<IM_SIZE; i++) {
            image->buffer[i] = (image->buffer[i] - min_value) / (max_value - min_value);
            if (image->buffer[i] > 1) image->buffer[i] = 1; // overflow protection
            if (image->buffer[i] < 0) image->buffer[i] = 0; // underflow protection
        }

    image->fourier_ok = 0;
}/*}}}*/

unsigned long int xy_to_i(ACC_DIM_TYPE x, ACC_DIM_TYPE y, t_im_struct *image){/*{{{*/
    return y*sx + x;
}/*}}}*/


void i_to_xy(unsigned long int i, ACC_DIM_TYPE *x, ACC_DIM_TYPE *y, t_im_struct *image){/*{{{*/
    (*y) = i / sx;
    (*x) = i - sx * (*y);
    return;
}/*}}}*/

void find_maximum(t_im_struct *image, unsigned long int *i, ACC_IM_TYPE *gl){/*{{{*/
    // finds maximum value of the image and fills the index i corresponding to
    // the position of the maximum, and the maximum grey level gl
    // if you're not interested in i or gl, set them to NULL
    if (i) *i = 0; //IMPORTANT!!!!
    unsigned long int j;
    ACC_DIM_TYPE max_val;
    max_val = image->buffer[0];
    for (j=0; j< IM_SIZE; j++) if (max_val < image->buffer[j]) {
            max_val = image->buffer[j];
            if (i) (*i) = j;
        }
    if (gl) (*gl = max_val);
    return;
}/*}}}*/

char check_sizes(t_im_struct *image1, t_im_struct *image2){/*{{{*/
    if ((image1->sizex == image2->sizex) && (image1->sizey == image2->sizey)) return IM_SIZES_MATCH;
    return IM_SIZES_NO_MATCH;
}/*}}}*/

char add_images(t_im_struct *image_d, t_im_struct *image_s){/*{{{*/
    // adds image_d and image_s and stores the result in image_d

    if (check_sizes(image_d, image_s) == IM_SIZES_NO_MATCH) return ACC_ERR_SIZES_MISMATCH;
    unsigned long int i;
    for (i=0; i< image_d->sizex * image_d->sizey; i++){
        ACC_DIM_TYPE x,y;
        i_to_xy(i, &x, &y, image_s);
        if ((image_s->maskx < 0) && (x <= -image_s->maskx)) continue;
        if ((image_s->maskx > 0) && (x >= image_s->sizex - image_s->maskx)) continue;
        if ((image_s->masky < 0) && (y <= -image_s->masky)) continue;
        if ((image_s->masky > 0) && (y >= image_s->sizey - image_s->masky)) continue;
        image_d->buffer[i] += image_s->buffer[i];
        if (image_d->weight_buffer) image_d->weight_buffer[i]++;
    }
    return ACC_ERR_NO_ERROR;

}/*}}}*/


char copy_image_with_weight(t_im_struct *image_d, t_im_struct *image_s, char initialized){/*{{{*/
    char err;
    if ((err = copy_image(image_d, image_s, initialized))) return err;
    if (image_s->weight_buffer) {
        if (!image_d->weight_buffer) initialize_weight_counting(image_d);
        memcpy(image_d->weight_buffer, image_s->weight_buffer, image_s->sizex * image_s->sizey * sizeof(unsigned int));
    }
    return ACC_ERR_NO_ERROR;
}/*}}}*/

char copy_image(t_im_struct *image_d, t_im_struct *image_s, char initialized){/*{{{*/
    // copies the image_s to image_d
    char err = ACC_ERR_NO_ERROR;
    if (initialized == IM_NOT_INITIALIZED)
        initialize_image(image_d, image_s->sizex, image_s->sizey);
    else {
        if ((err = check_sizes(image_d, image_s)) == IM_SIZES_NO_MATCH) return ACC_ERR_SIZES_MISMATCH;
    }
    memcpy(image_d->buffer, image_s->buffer, image_s->sizex * image_s->sizey * sizeof(ACC_IM_TYPE));
    image_d->maskx = image_s->maskx;
    image_d->masky = image_s->masky;
    return ACC_ERR_NO_ERROR;
}/*}}}*/

char crop_image(t_im_struct *image, ACC_DIM_TYPE left, ACC_DIM_TYPE right, ACC_DIM_TYPE top, ACC_DIM_TYPE bottom){/*{{{*/
    ACC_DIM_TYPE newsizex = image->sizex - left - right;
    ACC_DIM_TYPE newsizey = image->sizey - top - bottom;

    ACC_IM_TYPE *newbuffer = (ACC_IM_TYPE *) malloc (newsizex * newsizey * sizeof(ACC_IM_TYPE));

    ACC_DIM_TYPE x,y;
    unsigned long int i = 0;
    for (y=top; y<image->sizey-bottom; y++)
        for (x=left; x<image->sizex-right; x++){
            newbuffer[i++] = image->buffer[y*image->sizex + x];
        }
    free(image->buffer);
    image->buffer = newbuffer;
    image->sizex = newsizex;
    image->sizey = newsizey;
    image->fourier_ok = 0;
    if (image->weight_buffer) {
        free(image->weight_buffer);
        image->weight_buffer = NULL;
    }
    destroy_fft_plans(image);
    return ACC_ERR_NO_ERROR;
}/*}}}*/

// ****                             FOURIER                             ****
char create_fft_plans(t_im_struct *image){/*{{{*/
    image->fft_buffer = (fftw_complex *) fftw_malloc(IM_SIZE * sizeof(fftw_complex));
    assert(image->fft_buffer);

    image->fft_plan = fftw_plan_dft_2d(image->sizey, image->sizex, image->fft_buffer, image->fft_buffer, FFTW_FORWARD, FFTW_MEASURE);
    image->ifft_plan = fftw_plan_dft_2d(image->sizey, image->sizex, image->fft_buffer, image->fft_buffer, FFTW_BACKWARD, FFTW_MEASURE);
    return ACC_ERR_NO_ERROR;
}/*}}}*/

void destroy_fft_plans(t_im_struct *image){/*{{{*/
    if (image->fft_buffer){
        fftw_destroy_plan(image->fft_plan);
        fftw_destroy_plan(image->ifft_plan);
        fftw_free(image->fft_buffer);
        image->fft_buffer = NULL;
        image->fourier_ok = 0;
    }
}/*}}}*/

void calculate_fft(t_im_struct *image){/*{{{*/
    if (!image->fft_buffer) create_fft_plans(image);

    ACC_DIM_TYPE i;
    for (i=0; i < IM_SIZE; i++) {
        image->fft_buffer[i][0] = image->buffer[i];
        image->fft_buffer[i][1] = 0;
    }
    fftw_execute(image->fft_plan);
}/*}}}*/

char calculate_ifft(t_im_struct *image){/*{{{*/
    if (!image->fft_buffer){
        perror("Trying to calculate IFFT, while FFT buffer is not even initilaized.");
        return ACC_ERR_NO_FFT_BUFFER;
    }
    ACC_DIM_TYPE i;
    fftw_execute(image->ifft_plan);
    for (i=0; i < IM_SIZE; i++) image->buffer[i] = image->fft_buffer[i][0] / (IM_SIZE);

    if (image->weight_buffer) {
        free(image->weight_buffer);
        image->weight_buffer = NULL;
    }


    return ACC_ERR_NO_ERROR;
}/*}}}*/

void conjugate_multiply(fftw_complex *x, fftw_complex *y){/*{{{*/
    // Conjugates x and multiplies it with y. The result goes to x.
    double q0, q1;
    q0 = (*x)[0] * (*y)[0] + (*x)[1] * (*y)[1];
    q1 = (*x)[0] * (*y)[1] - (*x)[1] * (*y)[0];
    (*x)[0] = q0;
    (*x)[1] = q1;
}/*}}}*/

char conjugate_multiply_images_with_filter(t_im_struct *im1, t_im_struct *im2, /*{{{*/
                                           double (*filter_callback)(t_im_struct *im, unsigned long int i, void *user_par),
                                           void *user_par){
    // Conjugates the first image and multiplies with the second images. All this
    // happens only where filter outputs a nonzero value (to save time).
    // the user_par is a user parameter, which is forwarded to the filter.

    if (check_sizes(im1, im2) == IM_SIZES_NO_MATCH) {
        return ACC_ERR_SIZES_MISMATCH;
    }

    if (!im1->fourier_ok) calculate_fft(im1);
    if (!im2->fourier_ok) calculate_fft(im2);

    unsigned long int i;
    double filter_coef;
    for (i=0; i< im1->sizex * im1->sizey; i++){
        filter_coef=filter_callback(im1, i, user_par);
        if (filter_coef){
            conjugate_multiply(&(im1->fft_buffer[i]), &(im2->fft_buffer[i]));
            im1->fft_buffer[i][0] *= filter_coef; // apply filter
            im1->fft_buffer[i][1] *= filter_coef;
        } else {
            im1->fft_buffer[i][0] = 0; // filter zeroes this value
            im1->fft_buffer[i][1] = 0;

        }
    }
    return ACC_ERR_NO_ERROR;
}/*}}}*/

double no_filter(t_im_struct *im, unsigned long int i, void *user_par){/*{{{*/
    return 1;
}/*}}}*/

#define sx image->sizex
#define sy image->sizey
void i_to_xy_fourier_shift(unsigned long int i, ACC_SIG_DIM_TYPE *x, ACC_SIG_DIM_TYPE *y, t_im_struct *image){/*{{{*/
    // converts index from i to x and y coordinates, plus does the FFT shift
    // (like in Matlab http://www.mathworks.com/access/helpdesk/help/techdoc/ref/fftshift.html)
    i_to_xy(i, (ACC_DIM_TYPE *) x, (ACC_DIM_TYPE *) y, image);
    if ((*x) > sx/2) (*x) -= sx;
    if ((*y) > sy/2) (*y) -= sy;
}/*}}}*/

unsigned long int xy_to_i_fourier_shift(ACC_SIG_DIM_TYPE x, ACC_SIG_DIM_TYPE y, t_im_struct *image){/*{{{*/
    if (x < 0) x += (ACC_SIG_DIM_TYPE) sx;
    if (y < 0) y += (ACC_SIG_DIM_TYPE) sy;
    return xy_to_i((ACC_DIM_TYPE) x, (ACC_DIM_TYPE) y, image);
}/*}}}*/

double cylinder_filter(t_im_struct *im, unsigned long int i, void *user_par){/*{{{*/
    // this filter cuts off all frequencies hihger than corresponding to a given
    // radius. Radius is obtained as the user_par.
    double radius = *(double *)user_par;
    ACC_SIG_DIM_TYPE x,y;
    i_to_xy_fourier_shift(i, &x, &y, im);
    if (x*x+y*y < radius*radius) return 1; else return 0;
}/*}}}*/

double BW_filter(t_im_struct *im, unsigned long int i, void *user_par){/*{{{*/
    // this filter smoothly depresses high frequencies using Butterworth Filter
    // Radius is obtained as the user_par. parameter N must be hard coded below
    double radius = *(double *)user_par;
    int N = 4;// N = {1,2,3,..} , high N makes filter aproach ideal
    ACC_SIG_DIM_TYPE x,y;
    i_to_xy_fourier_shift(i, &x, &y, im);
    if (x*x+y*y > radius*radius*pow(97/3,1/N)) return 0;
    else return (1/(1+pow((y*y+x*x)/(radius*radius),N)));
}/*}}}*/

char shift_image(t_im_struct *image, double shiftx, double shifty){/*{{{*/
    // this function shifts the image by the given vector. The shiftx and shifty
    // don't have to be inegers.

    int err = 0;

    if (!image->fourier_ok) calculate_fft(image);

    unsigned long int i;
    double reg;

    image->maskx = shiftx;
    image->masky = shifty;

    ACC_SIG_DIM_TYPE x,y;
    for (i=0; i<IM_SIZE; i++) {
        i_to_xy_fourier_shift(i, &x, &y, image);
        double l = 2*M_PI*(x*shiftx/(double)sx + y*shifty/(double)sy);
        reg = image->fft_buffer[i][0] * cos(l) - image->fft_buffer[i][1] * sin(l);
        image->fft_buffer[i][1] = image->fft_buffer[i][1] * cos(l) + image->fft_buffer[i][0] * sin(l);
        image->fft_buffer[i][0] = reg;
    }

    err = calculate_ifft(image);
    if (err) return err;
    return ACC_ERR_NO_ERROR;
}/*}}}*/

#define AREA_SIZE 6
char find_maximum_subpixel(t_im_struct *image, double *xs, double *ys){/*{{{*/
    // this function finds the maximum with the sub-pixel resolution

    ACC_SIG_DIM_TYPE xf, yf; //position of the full-pixel maximum
    unsigned long int i;
    find_maximum(image, &i, NULL);
    i_to_xy_fourier_shift(i, &xf, &yf, image);

    ACC_SIG_DIM_TYPE *x = (ACC_SIG_DIM_TYPE *) malloc(sizeof(ACC_SIG_DIM_TYPE) * (2*AREA_SIZE) * (2*AREA_SIZE));
    assert(x);
    ACC_SIG_DIM_TYPE *y = (ACC_SIG_DIM_TYPE *) malloc(sizeof(ACC_SIG_DIM_TYPE) * (2*AREA_SIZE) * (2*AREA_SIZE));
    assert(y);
    ACC_IM_TYPE *z = (ACC_IM_TYPE *) malloc(sizeof(ACC_IM_TYPE) * (2*AREA_SIZE) * (2*AREA_SIZE));
    assert(z);

    grab_fitting_area(image, i, (ACC_DIM_TYPE) AREA_SIZE, x, y, z);

    *xs = xf; // set the initial values. IMPORTANT
    *ys = yf;

    char err = cubic_fit_maximum(x, y, z, 4*AREA_SIZE*AREA_SIZE, xs, ys);

    free(x);
    free(y);
    free(z);

    if (err) return err;
    return ACC_ERR_NO_ERROR;
}/*}}}*/

void grab_fitting_area(t_im_struct *image, unsigned long int i, ACC_DIM_TYPE area_size, ACC_SIG_DIM_TYPE *x, ACC_SIG_DIM_TYPE *y, ACC_IM_TYPE *z){/*{{{*/
    // this function creates three vectors of samples that will be used for the
    // two-dimensional polynomial least-squares fitting. x - the x coordinate,
    // y - the y coordinate, z - the corresponding grey level. The area is
    // defined by the central index i and the size.  x and y varies from  i-size to y+size
    // ATTENTION: x,y,z must be pre-allocated

    ACC_SIG_DIM_TYPE cx, cy; //the center coordinates
    i_to_xy_fourier_shift(i, &cx, &cy, image);
    ACC_SIG_DIM_TYPE p, q;
    unsigned long int r=0;
    unsigned long int j;
    for (q = cy - (signed) area_size; q < cy + (signed) area_size; q++)
        for (p = cx - (signed) area_size; p < cx + (signed) area_size; p++){
            j = xy_to_i_fourier_shift(p,q,image);
            x[r] = p;
            y[r] = q;
            z[r++] = image->buffer[j];
        }
}/*}}}*/

char cubic_fit_maximum(ACC_SIG_DIM_TYPE *x, ACC_SIG_DIM_TYPE *y, ACC_IM_TYPE *z, unsigned int count, double *xs, double *ys){/*{{{*/
    // fits the cubic polynomial on the data, which provides the coefficients.
    // These coefficients are then used to analytically find the maximum.

    // The least-squares problem can be converted into a set of linear equations,
    // which can be solved using Cholesky decomposition. See Wikipedia for
    // details.
    //
    // The matrix representation of the problem is: M * A = B.
    // M is the matrix of the sums of different combinations of different powers
    // of x and y (see
    // https://secure.wikimedia.org/wikipedia/de/wiki/Methode_der_kleinsten_Quadrate )

    unsigned int i, j, k; // general use indexes

    double *m = (double *) malloc (sizeof(double) * CDIM * CDIM); // allocate the matrix M
    assert(m);
    memset(m, 0, CDIM*CDIM*sizeof(double));
    double *b = (double *) malloc(CDIM * sizeof(double)); // the right side
    assert(b);
    memset(b, 0, CDIM*sizeof(double));
    double *a = (double *) malloc(CDIM * sizeof(double)); // the searched coefficients
    assert(a);

    // construction of the matrix M and the vector B
    double *r = (double *) malloc(CDIM * sizeof(double));
    for (i=0; i<count; i++){
        r[0] = x[i] * x[i] * x[i];
        r[1] = y[i] * y[i] * y[i];
        r[2] = x[i] * x[i] * y[i];
        r[3] = x[i] * y[i] * y[i];
        r[4] = x[i] * x[i];
        r[5] = y[i] * y[i];
        r[6] = x[i] * y[i];
        r[7] = x[i];
        r[8] = y[i];
        r[9] = 1;
        for (j=0; j< CDIM; j++) {
            for (k=j; k<CDIM; k++) m[k+CDIM*j] += r[k]*r[j]; // matrix M
            b[j] += r[j]*z[i]; // vector B
        }
    }
    free(r);

    lesolve(m, b, a, CDIM); // solve the linear-equation problem using Cholesky decomposition.
    free(m);
    free(b);

// maximum search using the Newton method. Maximum number of iterations is
// limited by MAX_ITER
#define MAX_ITER 10000
    i = 0;
    double sq_delta;
    do{
        newton_maximum_search(&sq_delta, a, xs, ys);
    } while ((sq_delta > 1e-8) && (i++ < MAX_ITER));
    free(a);
    if (i >= MAX_ITER) {
        return ACC_ERR_TOO_MANY_ITERATIONS;
    }
    return ACC_ERR_NO_ERROR;
}/*}}}*/

void newton_maximum_search(double *sq_delta, double *a, double *xin, double *yin){/*{{{*/
    double derx, dery;
    double x, y;
    int i;

    double *dxcoe = (double *) malloc(CDIM * sizeof(double));
    assert(dxcoe);
    double *dycoe = (double *) malloc(CDIM * sizeof(double));
    assert(dycoe);
    double *fcoe = (double *) malloc(CDIM * sizeof(double));
    assert(fcoe);

    x = *xin;
    y = *yin;

    fcoe[0] = x * x * x;
    fcoe[1] = y * y * y;
    fcoe[2] = x * x * y;
    fcoe[3] = x * y * y;
    fcoe[4] = x * x;
    fcoe[5] = y * y;
    fcoe[6] = x * y;
    fcoe[7] = x;
    fcoe[8] = y;
    fcoe[9] = 1;

    dxcoe[0] = 3 * x * x;
    dxcoe[1] = 0;
    dxcoe[2] = 2 * x * y;
    dxcoe[3] = y * y;
    dxcoe[4] = 2 * x;
    dxcoe[5] = 0;
    dxcoe[6] = y;
    dxcoe[7] = 1;
    dxcoe[8] = 0;
    dxcoe[9] = 0;

    dycoe[0] = 0;
    dycoe[1] = 3 * y * y;
    dycoe[2] = x * x;
    dycoe[3] = 2 * x * y;
    dycoe[4] = 0;
    dycoe[5] = 2 * y;
    dycoe[6] = x;
    dycoe[7] = 0;
    dycoe[8] = 1;
    dycoe[9] = 0;

    derx = 0;
    dery = 0;
    for (i=0; i<CDIM; i++) {
        derx += a[i]*dxcoe[i];
        dery += a[i]*dycoe[i];
    }

    *xin += derx;
    *yin += dery;

    free(fcoe);
    free(dxcoe);
    free(dycoe);

    *sq_delta = (derx * derx + dery * dery);
}/*}}}*/


void cholesky(double *a, double *l, int dim){/*{{{*/
    // Cholesky decomposition
    // a must be positive definite. l must be allocated. dim is dimension of
    // the matrix

    int i, j, k;
    double sm;

    memset(l,0,dim*dim*sizeof(double)); // zero the l-matrix

    for (j=0; j<dim; j++)
        for (i=0; i<dim; i++){
            if (i==j){
                sm=0;
                for (k=0; k<i; k++) sm += l[i+dim*k] * l[i+dim*k];
                l[i+dim*i] = sqrt(a[i+dim*i]-sm);
            }
            if (i < j){
                sm = 0;
                for (k=0; k<i; k++) sm += l[j+dim*k] * l[i+dim*k];
                l[j+dim*i] = (a[j+dim*i] - sm) / l[i+dim*i];
            }
        }
}/*}}}*/

void lesolve(double *a, double *y, double *x, int dim){/*{{{*/
    // Linear equation solver.
    // a = (dimXdim) matrix, y = right side vector , x <- the solution
    // vector
    double *b,*z;
    int i,j;

    assert(b = (double *) malloc(dim*dim*sizeof(double)));
    assert(z = (double *) malloc(dim*sizeof(double)));
//#define CDIM dim

    cholesky(a, b, dim);
    for (j=0; j<dim; j++){
        z[j] = y[j];
        for (i=j-1; i>=0; i--){
            z[j] -= b[j+dim*i]*z[i];
        }
        z[j] /= b[j+dim*j];
    }
    for (j=dim-1; j>=0; j--){
        x[j] = z[j];
        for (i=dim-1; i>j; i--){
            x[j] -= b[i+dim*j]*x[i];
        }
        x[j] /= b[j+dim*j];
    }
    free(b);
    free(z);
}/*}}}*/

