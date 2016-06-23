#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

#ifdef NETCDF
#include "netcdf.h"
#endif

#ifdef OMP
#include <omp.h>
#endif

#define ABS(x)		((x) < 0 ? -(x) : (x))
#define MAX(a,b)    ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b)    ( ((a) < (b)) ? (a) : (b) )

// radius of the Earth in metres
//#define earthradius (6356750.52)
#define R	(6378100.0)
#define earthradius (6378100.0)
//#define R	(earthradius)

#define PI	(3.14159265359)
#define	d2r(x)	(x*PI/180.0)
#define r2d(x)  (x*180.0/PI)


#define TOLERANCE	1e-6

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\nFunction: %s\nFile: %s\nLine %d\n", nc_strerror(e),__func__,__FILE__,__LINE__); exit(ERRCODE);}

#define fail(...) my_fail(__LINE__,__func__,__FILE__,__VA_ARGS__)



typedef struct{
	double	node_coord[4][2];
	double	node_value[4];
	double	interp_weights[4];
}element;

typedef struct{
	double	*x;
	double	*y;
}mesh;

typedef struct{
	// netcdf variables
    int     ncid;
    int     varid;
    int     retval;

	// index variables
    size_t      x;
    size_t      y;
	size_t		z;
	size_t		t;

}netcdf;

typedef struct{

	int		nElements;
	int		nodesPerEl;
	int		nx;
	int		ny;

	// anti-clockwise element numbering
	double	xi[4];
	double	eta[4];

	double	pos[2];

	element	*el;
	mesh	msh;


	// fields data read from netcdf
	netcdf	nc;
	double ***field;
	double ***u;
	double ***v;
	double *lon;
	double *lat;
	double *time;

}e;

// prototypes
void interpolate_point(element *el, double *interp_value);

void calculate_interpolation_weights(element *el, double *xi, double *eta, double *pos);

double	evaluate_linear_quad_shape_function( double *xi, double *eta, double *pos_to_interp_at, int node );

void do_test(e *E);

void init_xi_eta(e *E);

void generate_mesh(e *E, int nx, int ny);
void print_elements(e *E);


// stuff in main
double relative_difference(double a, double b);

double spheriq_dist(double lon1, double lat1, double lon2, double lat2, int debug);
double distance(double lat1, double lon1, double lat2, double lon2);
double bearing(double lat1, double lon1, double lat2, double lon2);


double ***malloc3d_double(int dim1, int dim2, int dim3);
void my_fail( const int line, const char *func, const char *file, const char *msg, ... );

float RandomFloat(float min, float max);


// netcdf read stuff
void get_field(e *E, char *field_name, void* field);
void get_attribute(e *E, char *var, char *att_name, double *att);
void get_dimension(e *E, char *dim_name, size_t *dim);
