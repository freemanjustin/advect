#include "advect.h"

double relative_difference(double a, double b){
	double c = ABS(a);
	double d = ABS(b);

	d = MAX(c, d);

	return d == 0.0 ? 0.0 : ABS(a - b) / d;	// wtf is this?
}

double analytic_function(double x, double y){

	//return(sin(x)*cos(y));
	return(sin(x*y)+cos(y*x));
}

double ***malloc3d_double(int dim1, int dim2, int dim3)
{

	size_t		layer1_count = dim1;
	size_t		layer2_count = dim1 * dim2;
	size_t		layer1_size = sizeof(double ***) * layer1_count;
	size_t		layer2_size = sizeof(double **) * layer2_count;

	size_t 	layers_size = layer1_size + layer2_size ;

	size_t	 	data_count = dim1 * dim2 * dim3 ;
	size_t 	data_size = sizeof(double) * data_count;

	void 	*raw_bytes = (void*)malloc(layers_size + data_size);

	double ***layer1 = (double ***)(raw_bytes);
	double **layer2 = (double **)(raw_bytes + layer1_size);
	double *double_data = (double *)(raw_bytes + layers_size);

	int i, j;
	double **this_layer2;
	//double *this_layer3;

	for (i = 0; i < dim1; i++) {
		layer1[i] = layer2 + (i * dim2);
		this_layer2 = layer1[i];

		for (j = 0; j < dim2; j++) {
			this_layer2[j] = double_data + (i * dim2 * dim3 + j * dim3);
		}
	}
	return layer1;

}

void my_fail( const int line, const char *func, const char *file, const char *msg, ... )
{
    va_list args;

    fprintf( stdout, "mdbcat error: file %s in function %s at line %d \n", file, func, line );

    va_start(args, msg);
    vprintf(msg, args);
    va_end(args);

    fflush( stdout );

    exit(1);

}

float RandomFloat(float min, float max)
{
    // this  function assumes max > min, you may want
    // more robust error checking for a non-debug build
    //assert(max > min);
    float random = ((float) rand()) / (float) RAND_MAX;

    // generate (in your case) a float between 0 and (4.5-.78)
    // then add .78, giving you a float between .78 and 4.5
    float range = max - min;
    return (random*range) + min;
}

void get_field(e *E, char *field_name, void* field){

    if(( E->nc.retval = nc_inq_varid(E->nc.ncid, field_name, &E->nc.varid)))
        ERR(E->nc.retval);

    if ((E->nc.retval = nc_get_var_double(E->nc.ncid, E->nc.varid, (double*)field)))
        ERR(E->nc.retval);

}

void get_attribute(e *E, char *var, char *att_name, double *att){

    if(( E->nc.retval = nc_inq_varid(E->nc.ncid, var, &E->nc.varid)))
        ERR(E->nc.retval);

    if((E->nc.retval = nc_get_att_double(E->nc.ncid, E->nc.varid, att_name, att)))
        ERR(E->nc.retval);

}

void get_dimension(e *E, char *dim_name, size_t *dim){


    if((E->nc.retval = nc_inq_dimid(E->nc.ncid, dim_name, &E->nc.varid)))
        ERR(E->nc.retval);
    if((E->nc.retval = nc_inq_dimlen(E->nc.ncid,E->nc.varid,dim)))
        ERR(E->nc.retval);
}



void do_test(e *E){

	double	interp_value;

	// construct a test element
	// anticlockwise numbering
	E->el[0].node_value[0] = 100.0;
	E->el[0].node_value[1] = 60.0;
	E->el[0].node_value[2] = 50.0;
	E->el[0].node_value[3] = 90.0;

	printf("node values: 0 = %f, 1 = %f, 2 = %f, 3 = %f\n", E->el[0].node_value[0], E->el[0].node_value[1], E->el[0].node_value[2], E->el[0].node_value[3]);

	// this is the element coordinates
	E->el[0].node_coord[0][0] = 2.0;	E->el[0].node_coord[0][1] = 2.0; // 0
	E->el[0].node_coord[1][0] = 4.0;	E->el[0].node_coord[1][1] = 2.0; // 1
	E->el[0].node_coord[2][0] = 4.0;	E->el[0].node_coord[2][1] = 3.0; // 2
	E->el[0].node_coord[3][0] = 2.0;	E->el[0].node_coord[3][1] = 3.0; // 3

	// this is the interpolation point
	E->pos[0] = 2.5;
	E->pos[1] = 2.5;

	calculate_interpolation_weights(&E->el[0], E->xi, E->eta, E->pos);

	interpolate_point(&E->el[0], &interp_value);
	printf("Interp value from precomput weights = %f\n", interp_value);

}

void init_xi_eta(e *E){

	// anti-clockwise element numbering
	// master element nodal positions
    //
    //  E->xi is the x coordinate
    //  E->eta is the y coordinate
    //
    //            ^
    // (-1,1)     |      (1,1)
    //    o----------------o
    //    |       |        |
    //    |     (0,0)------|- >
    //    |                |
    //    o----------------o
    // (-1,-1)           (1,-1)
    //

	E->xi[0] =	-1.0;
	E->xi[1] =	 1.0;
	E->xi[2] =	 1.0;
	E->xi[3] =	-1.0;

	E->eta[0] = -1.0;
	E->eta[1] = -1.0;
	E->eta[2] =  1.0;
	E->eta[3] =  1.0;

}




void generate_mesh(e *E, int nx, int ny){

	int i,j, el;
	double x,y,xstart, ystart, xend,yend,dx, dy;

	dx = 2.0/(double)(nx);
	dy = 2.0/(double)(ny);

	//printf("# dx = %f, dy = %f\n", dx, dy);

	// these are the starting positions for x and y
	xstart = -1.0;
	xend = 1.0;

	ystart = -1.0;
	yend = 1.0;

	x = xstart;
	y = ystart;
	el = 0;
	for(i=0;i<ny;i++){	// rows
		for(j=0;j<nx;j++){	// cols

			E->el[el].node_coord[0][0] = x;		E->el[el].node_coord[0][1] = y; // 0
			E->el[el].node_coord[1][0] = x+dx;	E->el[el].node_coord[1][1] = y; // 1
			E->el[el].node_coord[2][0] = x+dx;	E->el[el].node_coord[2][1] = y+dy; // 2
			E->el[el].node_coord[3][0] = x;		E->el[el].node_coord[3][1] = y+dy; // 3

			E->el[el].node_value[0] = analytic_function(x,y);
			E->el[el].node_value[1] = analytic_function(x+dx,y);
			E->el[el].node_value[2] = analytic_function(x+dx,y+dy);
			E->el[el].node_value[3] = analytic_function(x,y+dy);

			el++;
			// update lower left coordinte for next element calc

			x = x + (dx);
		}
		y = y + (dy);
		x = xstart;

	}

}

void print_elements(e *E){

	int i;

	for(i=0;i<E->nElements;i++){
		printf("Element %d:\n", i);
		printf("\tnode pos: 0 = (%f,%f) 1 = (%f,%f) 2 = (%f,%f) 3 = (%f,%f)\n", E->el[i].node_coord[0][0], E->el[i].node_coord[0][1],E->el[i].node_coord[1][0], E->el[i].node_coord[1][1],
			                                                                    E->el[i].node_coord[2][0], E->el[i].node_coord[2][1],E->el[i].node_coord[3][0], E->el[i].node_coord[3][1]);
		printf("\tnode values: 0 = %f, 1 = %f, 2 = %f, 3 = %f\n", E->el[i].node_value[0], E->el[i].node_value[1], E->el[i].node_value[2], E->el[i].node_value[3]);

	}

}


/* this function prints out the node coordinates and values for the mesh
 that was generated by generate_mesh.
 the output format is suitable for plotting in gnuplot

 ./app > out
 gnuplot
 set pm3d
 unset surface
 splot "out" u 1:2:3

*/
void print_mesh(e *E, int nx, int ny){

	int i,j,el;
	FILE	*out;

	out = fopen("generated.dat","w");

	el=0;
	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			fprintf(out, "%f %f %f\n", E->el[el].node_coord[0][0], E->el[el].node_coord[0][1], E->el[el].node_value[0]);
			el++;
		}
		fprintf(out,"%f %f %f\n", E->el[el-1].node_coord[1][0], E->el[el-1].node_coord[1][1], E->el[el-1].node_value[1]);
		fprintf(out,"\n");
	}
	// print the top row
	el = el-nx;
	for(j=0;j<nx;j++){
		fprintf(out,"%f %f %f\n", E->el[el].node_coord[3][0], E->el[el].node_coord[3][1], E->el[el].node_value[3]);
		el++;
	}
	fprintf(out,"%f %f %f\n", E->el[el-1].node_coord[2][0], E->el[el-1].node_coord[2][1], E->el[el-1].node_value[2]);
	fclose(out);
}

void store_mesh(e *E, int nx, int ny){

	int i,j,el;

	el=0;

	//printf("x node coordinates:\n");
	for(i=0;i<nx;i++){
		E->msh.x[i] = E->el[i].node_coord[0][0];
		//printf("i = %d: %f\n",i, E->msh.x[i]);
	}
	E->msh.x[nx] = E->el[nx-1].node_coord[1][0];
	//printf("i = %d: %f\n",nx, E->msh.x[nx]);

	//printf("y node coordinates:\n");
	for(i=0;i<ny;i++){
		E->msh.y[i] = E->el[i*nx].node_coord[0][1];
		//printf("i = %d: %f\n",i, E->msh.y[i]);
	}
	E->msh.y[ny] = E->el[(ny-1)*nx].node_coord[3][1];
	//printf("i = %d: %f\n",ny, E->msh.y[ny]);
}


void generate_target_nodes(e *E, int nx, int ny){

	int i,j,el;

	el=0;
	for(i=0;i<ny;i++){
		for(j=0;j<nx;j++){
			//printf("%f %f %f\n", E->el[el].node_coord[0][0], E->el[el].node_coord[0][1], E->el[el].node_value[0]);
			el++;
		}
		//printf("%f %f %f\n", E->el[el-1].node_coord[1][0], E->el[el-1].node_coord[1][1], E->el[el-1].node_value[1]);
		//printf("\n");
	}
	// print the top row
	el = el-nx;
	for(j=0;j<nx;j++){
		//printf("%f %f %f\n", E->el[el].node_coord[3][0], E->el[el].node_coord[3][1], E->el[el].node_value[3]);
		el++;
	}
	//printf("%f %f %f\n", E->el[el-1].node_coord[2][0], E->el[el-1].node_coord[2][1], E->el[el-1].node_value[2]);

}

int	get_owner_element(e *E, double *pos){

	int		i;
	int		max_its;
	int		element;
	int		element_new;
	int		x_element;
	int		y_element;
	double	point;
	double	interval_size;
	int		check = 1;

	max_its = 100;

	// search over elements to determine which element contains this position

	// apply a bisection method on the elements

	// search in x
	i=0;
	element = (double)E->nx/2.0;
	element_new = 0;
	interval_size = E->nx/2.0;

	//printf("begin: element = %d, element_new = %d\n", element, element_new);
	do{
		//printf("it = %d:   element = %d, element_new = %d, interval_size = %f, E->msh.x[%d] = %f, pos[0] = %f\n",i, element, element_new,interval_size, element, E->msh.x[element], pos[0]);
		//printf("E->msh.x[%d] = %f, E->msh.x[%d] = %f\n", element, E->msh.x[element], element+1, E->msh.x[element+1]);
		//if( (pos[0] <= E->msh.x[element+1]) && (pos[0] >= E->msh.x[element])  ){
		if( (relative_difference(pos[0],E->msh.x[element+1])<=TOLERANCE) || (relative_difference(pos[0],E->msh.x[element])<=TOLERANCE)
		   ||
		   (pos[0] <= E->msh.x[element+1]) && (pos[0] >= E->msh.x[element]) ){
			// found it
			element = element;
			//printf("found it: element is %d\n", element);
			break;
		}
		else if( E->msh.x[element] >= pos[0] ){ // means that the pos point is on the left of our current index
			//printf("pre >: element = %d, (double)element/2.0 = %f\n", element, (double)element/2.0);
			//printf("pre >: element - element/2.0 = %f\n", (double)element - (double)element/2.0);
			element_new = ( (double)element - ceil(interval_size/2.0) );
			if(element_new < 0) element_new=0;
			//printf("is greater than: element = %d, element_new = %d\n", element, element_new);
		}
		else{	// means that the pos point is on the right of our current index

			element_new = ceil( (double)element + (interval_size/2.0) );

			if(element_new >= E->nx){
				//printf("SET element_new: element_new = %d!!!\n", element_new);
				element_new=E->nx-1;
				//printf("SET element_new to E->nx-1: element_new = %d!!!\n", element_new);
			}

			//printf("is less than: element = %d, element_new = %d, interval_size = %f\n", element, element_new, interval_size);
		}
		interval_size = fabs(element_new - element);
		element = element_new;
		i++;
	}while(i<max_its);

	if(i>=max_its){
		printf("x get owner failed!\n");
		printf("pos = %f, %f\n", pos[0], pos[1]);
		printf("its = %d\n", i);

		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n", E->el[element].node_coord[3][0], E->el[element].node_coord[3][1],
			   E->el[element].node_coord[2][0], E->el[element].node_coord[2][1]);
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n\n", E->el[element].node_coord[0][0], E->el[element].node_coord[0][1],
			   E->el[element].node_coord[1][0], E->el[element].node_coord[1][1]);

		exit(1);
	}
	//printf("x: get_owner converged in %d iterations: pos[0] = %f is in element %d:\n", i, pos[0], element);
	//printf("E->msh.x[%d] =  %f, E->msh.x[%d+1] =  %f\n", element, E->msh.x[element],element, E->msh.x[element+1]);


	// x_element stores the column that this point will be in
	// lets use this information to constrain our search in y

	x_element = element;
	element = ceil((double)E->ny/2.0);
	element_new = 0;//(double)E->ny/2.0;
	interval_size = E->ny/2.0;
	// now search in y
	i=0;
	do{
		//printf("  y: element = %d, E->el[element].node_coord[0][1] = %f,E->el[element].node_coord[3][1] = %f, pos[1] = %f\n",element,  E->el[element].node_coord[0][1],E->el[element].node_coord[3][1], pos[1]);
		//if( ((pos[1] >= E->msh.y[element]) && (pos[1] <= E->msh.y[element+1])) ){
		if( (relative_difference(pos[1],E->msh.y[element+1])<=TOLERANCE) || (relative_difference(pos[1],E->msh.y[element])<=TOLERANCE)
		   ||
		   (pos[1] <= E->msh.y[element+1]) && (pos[1] >= E->msh.y[element]) ){
			//printf("found it\n");
			// found it
			break;
		}

		if( (E->msh.y[element] >= pos[1]) ){
			element_new = ( (double)element - ceil(interval_size/2.0) ) ;
			if(element_new < 0) element_new=0;
			//printf("y: is greater than: element = %d, element_new = %d\n", element, element_new);
		}
		else{
			element_new = ( (double)element + ceil(interval_size/2.0) );

			if(element_new >= E->ny) {
				element_new=E->ny-1;
				//printf("SET element_new to E->ny-1!!!\n");
			}

			//printf("y: is less than: element = %d, element_new = %d\n", element, element_new);
		}

		interval_size = fabs(element_new - element);
		element = element_new;


		i++;
	}while(i<max_its);

	if(i>=max_its){
		printf("y get owner failed!\n");
		printf("pos = %f, %f\n", pos[0], pos[1]);
		printf("its = %d\n", i);

		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n", E->el[element].node_coord[3][0], E->el[element].node_coord[3][1],
			   E->el[element].node_coord[2][0], E->el[element].node_coord[2][1]);
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n\n", E->el[element].node_coord[0][0], E->el[element].node_coord[0][1],
			   E->el[element].node_coord[1][0], E->el[element].node_coord[1][1]);

		exit(1);
	}

	//printf("y: get_owner converged in %d iterations: pos[1] = %f is in element %d:\n", i, pos[1], element);
	//printf("y coords are: %f, %f\n", E->msh.y[element], E->msh.y[element+1]);

	// the owner element for this point is element * x_element
	element = (element*(E->nx)) + x_element;

	/*
	if( ( pos[0] < E->el[element].node_coord[0][0] ) || ( pos[0] > E->el[element].node_coord[1][0] ) ){

		printf("get_owner failed: x location is not in element:\n");
		if( pos[0] < E->el[element].node_coord[0][0] ){
			printf("i reckon pos[0] = %.6g is less than coord [0][0] = %.6g, the diff is %.12g\n", pos[0],E->el[element].node_coord[0][0], pos[0] - E->el[element].node_coord[0][0] );
			printf("the relative diff is = %.6g\n", relative_difference(pos[0],E->el[element].node_coord[0][0]) );
		}
		if( pos[0] > E->el[element].node_coord[1][0] )
			printf("i reckon pos[0] is greater than [1][0]\n");

		printf("\ninterp point is %f, %f\n", pos[0],pos[1]);
		printf("owner element is %d\n", element);
		printf("element coords:\n \n");
		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n", E->el[element].node_coord[3][0], E->el[element].node_coord[3][1],
												 E->el[element].node_coord[2][0], E->el[element].node_coord[2][1]);
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n\n", E->el[element].node_coord[0][0], E->el[element].node_coord[0][1],
			   E->el[element].node_coord[1][0], E->el[element].node_coord[1][1]);

		exit(1);

	}

   if( ( pos[1] < E->el[element].node_coord[0][1]) || ( pos[1] > E->el[element].node_coord[3][1]) ){
	//if( (relative_difference(pos[1],E->el[element].node_coord[0][1]) < 0.0) || (relative_difference(pos[1],E->el[element].node_coord[3][1]) < 0.0) ){
	   printf("get_owner failed: y location is not in element:\n");

	  printf("\ninterp point is %f, %f\n", pos[0],pos[1]);
	  printf("owner element is %d\n", element);
	  printf("element coords:\n \n");
	  printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n", E->el[element].node_coord[3][0], E->el[element].node_coord[3][1],
			 E->el[element].node_coord[2][0], E->el[element].node_coord[2][1]);
	  printf("   |                              |\n");
	  printf("   |                              |\n");
	  printf("   |                              |\n");
	  printf("   |                              |\n");
	  printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n\n", E->el[element].node_coord[0][0], E->el[element].node_coord[0][1],
			 E->el[element].node_coord[1][0], E->el[element].node_coord[1][1]);
	  exit(1);
	}
	*/

	/*
	if(check){
		printf("\n\n\npos = %f, %f\n", pos[0], pos[1]);
		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n", E->el[element].node_coord[3][0], E->el[element].node_coord[3][1],
			   E->el[element].node_coord[2][0], E->el[element].node_coord[2][1]);
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n\n", E->el[element].node_coord[0][0], E->el[element].node_coord[0][1],
			   E->el[element].node_coord[1][0], E->el[element].node_coord[1][1]);
	}
	*/

	return element;
}







// robust method but very slow (N^2)
int	get_owner_elementII(e *E, double *pos){

	int	i;
	int	element, element_x;
	double	*xdiff;
	double	*ydiff;
	double	min;
	double value;

	xdiff = malloc((E->nx+1)*sizeof(double));
	ydiff = malloc((E->ny+1)*sizeof(double));


	// find the minimun diff between left hand mesh nodes and our x position
	min = 999999;
	for(i=0;i<E->nx;i++){
		xdiff[i] = E->msh.x[i] - pos[0];
		if( xdiff[i] <= 0.0){
			value = fabs(xdiff[i]);
			if( (value < min) || (relative_difference(value,min) <= TOLERANCE) ){
				min = value;
				element_x = i;
				//printf("setting min: i = %d, min = %.12f, xdiff = %.12g\n",i,min, xdiff[i]);
			}
		}
		else
			break;
	}



	//printf("element_x = %d\n", element_x);


	// find the minimun diff between left hand mesh nodes and our x position
	min = 999999;
	for(i=0;i<=E->ny;i++){
		ydiff[i] = E->msh.y[i] - pos[1];
		if( ydiff[i] <= 0.0){
			value = fabs(ydiff[i]);
			if( (value < min) || (relative_difference(value,min) <= TOLERANCE) ){

				min = value;
				element = i;
				//printf("setting min: i = %d, min = %.12f, xdiff = %.12g\n",i,min, xdiff[i]);
			}
		}
		else
			break;

	}
	//printf("element = %d\n", element);
	//printf("owner element = %d\n", (element*(E->nx)) + element_x);

	element = (element*(E->nx)) + element_x;


	if( ( pos[0] < E->el[element].node_coord[0][0] ) || ( pos[0] > E->el[element].node_coord[1][0] ) ){
		printf("get_ownerII failed: x location is not in element:\n");

		printf("\ninterp point is %f, %f\n", pos[0],pos[1]);
		printf("xdiff array is:\n");
		for(i=0;i<=E->nx;i++){
			printf("xdiff[%d] = %f\n", i, xdiff[i]);
		}
		printf("owner element is %d\n", element);
		printf("element coords:\n \n");
		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n", E->el[element].node_coord[3][0], E->el[element].node_coord[3][1],
			   E->el[element].node_coord[2][0], E->el[element].node_coord[2][1]);
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n\n", E->el[element].node_coord[0][0], E->el[element].node_coord[0][1],
			   E->el[element].node_coord[1][0], E->el[element].node_coord[1][1]);
		exit(1);

	}

	if( ( pos[1] < E->el[element].node_coord[0][1]) || ( pos[1] > E->el[element].node_coord[3][1]) ){
		printf("get_owner failed: y location is not in element:\n");

		printf("\ninterp point is %f, %f\n", pos[0],pos[1]);
		printf("owner element is %d\n", element);
		printf("element coords:\n \n");
		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n", E->el[element].node_coord[3][0], E->el[element].node_coord[3][1],
			   E->el[element].node_coord[2][0], E->el[element].node_coord[2][1]);
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf("   |                              |\n");
		printf(" (%.5f, %.5f) --------- (%.5f, %.5f)\n\n", E->el[element].node_coord[0][0], E->el[element].node_coord[0][1],
			   E->el[element].node_coord[1][0], E->el[element].node_coord[1][1]);
		exit(1);
	}


	return element;
}

void test_interp(e *E){

	int	i,j;
	int	nx, ny;
	double x,y,xstart, ystart, dx, dy;
	double xend, yend;
	double pos[2];
	int	element;
	double	interp_value;
	double	analytic_value;
	FILE	*out;
	FILE	*error;



	// set up the source mesh
	E->nx = 128;
	E->ny = 128;

	E->nElements = E->nx*E->ny ;

	E->el = malloc(E->nElements*sizeof(element));
	E->nodesPerEl = 4;

	E->msh.x = malloc((E->nx+1) * sizeof(double));
	E->msh.y = malloc((E->ny+1) * sizeof(double));

	init_xi_eta(E);

	generate_mesh(E, E->nx, E->ny);
	//print_elements(E);
	print_mesh(E,E->nx,E->ny);
	store_mesh(E, E->nx, E->ny);



	out = fopen("interp.dat","w");
	error = fopen("error.dat","w");


	// set up the target mesh
	// generate a list of nodal coordinates to interpolate onto
	nx = 257;//E->nx;
	ny = 317;//E->ny;

	xstart = -1.0;
	xend = 1.0;
	ystart = -1.0;
	yend= 1.0;

	x = xstart;
	y = ystart;

	dx = (xend-xstart)/(double)(nx);
	dy = (yend-ystart)/(double)(ny);


	//#pragma omp parallel for private(i,j,pos, element, x, y, interp_value, analytic_value)
	for(i=0;i<=ny;i++){
		for(j=0;j<=nx;j++){

			if(j==0) pos[0] = xstart;
			else if(j==nx) pos[0] = xend;
			else if(relative_difference(x,0.0)<TOLERANCE) x = 0.0;
			else pos[0] = x;

			if(i==0) pos[1] = ystart;
			else if(i==ny) pos[1] = yend;
			else if(relative_difference(y,0.0)<TOLERANCE) y = 0.0;
			else pos[1] = y;

			//printf("\n\n\npos = ( %f, %f)\n", pos[0], pos[1]);

			// find out which element this lies within
			element = get_owner_element(E, pos);
			//element = get_owner_elementII(E, pos);

			calculate_interpolation_weights(&E->el[element], E->xi, E->eta, pos);
			interpolate_point(&E->el[element], &interp_value);

			analytic_value = analytic_function(pos[0],pos[1]);

			fprintf(out,"%f %f %f\n", pos[0], pos[1], interp_value);
			fprintf(error,"%f %f %.4g\n", pos[0], pos[1], interp_value-analytic_value);//100.0*(fabs(interp_value - analytic_value))/fabs(analytic_value)  );

			x = x + (dx);
		}
		fprintf(out,"\n");
		fprintf(error,"\n");
		y = y + (dy);
		x = xstart;
	}

	fclose(out);
	fclose(error);

}


void update_particle_position_euler(double *pos, double *vel, double dt)
{
	int d;

	for( d=0; d<2; d++ ) {
		//printf("pre: pos[%d] = %f, vel[%d] = %f, dt = %f\n", d, pos[d], d, vel[d], dt);
		pos[d] = pos[d] + dt * vel[d];
		//printf("post: pos[%d] = %f, vel[%d] = %f, dt = %f\n", d, pos[d], d, vel[d], dt);
	}

}


// calculate the timestep value
// returns dt (seconds)
// input is:
//		vMax = maximum velocity in mesh (m/s)
//		sMin = minimum mesh size (metres)
void get_dt( double* dt, double	max_vel, double min_mesh ) {

				int	dim = 2;	// only 2d mesh

				// Set advective time step.
        *dt = min_mesh / ( sqrt(dim) * max_vel );;

}



// to build with netcdf support:
// gcc -DNETCDF -O3 -o test `/Users/jfreeman/necdf401/bin/nc-config --cflags` main.c /Users/jfreeman/necdf401/lib/libnetcdf.a

void test_interp_netcdf(e *E, char *file1, char *file2){

	int	i,j,k,t,el;

	int	nx, ny;

	double x,y,xstart, ystart, dx, dy;
	double xend, yend;
	double dt;

	// element dimensions in metres
	double	cell_width;
	double	cell_height;

	double	dlon;
	double	dlat;

	double pos[2];
	double cart_pos[3];
	double vel[2];

	double z;

	double	interp_value;
	char	fname[256];
	FILE	*out;

	double length;
	double	min_lon, max_lon, min_lat, max_lat;
	double	min_x, min_y;
	double	min_mesh, max_vel;
	double	current_time;

	printf("file = %s\n", file1);

	if ((E->nc.retval = nc_open(file1, NC_NOWRITE, &E->nc.ncid)))
		ERR(E->nc.retval);

	// get dims from netcdf file
	get_dimension(E, "xu_ocean", &E->nc.x);
	get_dimension(E, "yu_ocean", &E->nc.y);
	get_dimension(E, "Time", &E->nc.t);

	printf("dimensions:\n");
	printf("xu_ocean = %d\n", E->nc.x);
	printf("yu_ocean = %d\n", E->nc.y);
	printf("Time = %d\n", E->nc.t);

	// malloc memeory for the dimensions
	E->lon = malloc(E->nc.x*sizeof(double));
	E->lat = malloc(E->nc.y*sizeof(double));
	E->time = malloc(E->nc.t*sizeof(double));

	// malloc room for the field
	//E->field = malloc3d_double(E->nc.t+1, E->nc.y, E->nc.x); // cause the field is time, lat, lon
	E->u = malloc3d_double(E->nc.t+1, E->nc.y, E->nc.x);
	E->v = malloc3d_double(E->nc.t+1, E->nc.y, E->nc.x);

	// read the field
	get_field(E,"xu_ocean", &E->lon[0]);
	get_field(E,"yu_ocean", &E->lat[0]);
	get_field(E,"Time", &E->time[0]);
	get_field(E,"usurf", &E->u[0][0][0]);

	nc_close(E->nc.ncid);

	if ((E->nc.retval = nc_open(file2, NC_NOWRITE, &E->nc.ncid)))
		ERR(E->nc.retval);

	get_field(E,"vsurf", &E->v[0][0][0]);

	nc_close(E->nc.ncid);

	printf("u[0][0][0] = %f, v[0][0][0] = %f\n", E->u[0][0][0], E->v[0][0][0]);
	// unpack the netcdf data
	for(k=0;k<E->nc.t;k++){
		for(i=0;i<E->nc.y;i++){
			for(j=0;j<E->nc.x;j++){
				//E->u[k][i][j] *= 0.0003051851;
				//E->v[k][i][j] *= 0.0003051851;

				if(E->u[k][i][j] < -10.0){
					E->u[k][i][j] = 0.0;
					E->v[k][i][j] = 0.0;
				}
			}
		}
	}
	printf("u[0][0][0] = %f, v[0][0][0] = %f\n", E->u[0][0][0], E->v[0][0][0]);

	// setup the interpolation source grid data structures
	E->nx = E->nc.x-1;
	E->ny = E->nc.y-1;
	E->nElements = E->nx*E->ny ;
	E->el = malloc(E->nElements*sizeof(element));
	E->nodesPerEl = 4;

	E->msh.x = malloc((E->nx+1) * sizeof(double));
	E->msh.y = malloc((E->ny+1) * sizeof(double));

	init_xi_eta(E);

	printf("E->nx = %d, E->ny = %d\n", E->nx, E->ny);

	//generate_mesh(E, E->nx, E->ny);
	el = 0;
	t = 0;
	min_mesh = 9e12;	// is this too small for some meshes?
	for(i=0;i<E->ny;i++){
		for(j=0;j<E->nx;j++){

			// set positions for this element
			// element node numbering is anti-clockwise
			E->el[el].node_coord[0][0] = E->lon[j];
			E->el[el].node_coord[0][1] = E->lat[i]; // 0

			E->el[el].node_coord[1][0] = E->lon[j+1];
			E->el[el].node_coord[1][1] = E->lat[i]; // 1

			E->el[el].node_coord[2][0] = E->lon[j+1];
			E->el[el].node_coord[2][1] = E->lat[i+1]; // 2

			E->el[el].node_coord[3][0] = E->lon[j];
			E->el[el].node_coord[3][1] = E->lat[i+1]; // 3

			// determine the minimum mesh distance in metres
			// we only check the bottom and right side of the element
			// as we assume the mesh is orthogonal
			for(k=0;k<2;k++){
				length = spheriq_dist(E->el[el].node_coord[k][0], E->el[el].node_coord[k][1],
					 													E->el[el].node_coord[k+1][0], E->el[el].node_coord[k+1][1], 0);

				if(length < min_mesh ) {
					min_mesh = length;
					printf("\t\tcurrent min x distance is %f\n", min_mesh);
					printf("\t\tthis is at: element = %d, i = %d, j = %d\n",el, i, j);
				}

			}

			// now set the nodal value for this element
			E->el[el].node_value[0] = E->u[t][i][j];
			E->el[el].node_value[1] = E->u[t][i][j+1];
			E->el[el].node_value[2] = E->u[t][i+1][j+1];
			E->el[el].node_value[3] = E->u[t][i+1][j];

			el++;
		}
	}

	// for plotting only
	print_mesh(E,E->nx,E->ny);
	// the store_mesh call is required for the function get_owner_element
	store_mesh(E, E->nx, E->ny);

	// set up the target mesh
	// generate a list of nodal coordinates to interpolate onto
	/*
	nx = 1000;//E->nx;
	ny = 800;//E->ny;

	xstart = E->lon[0];//0.0;
	xend = E->lon[E->nc.x-1];//358.0;
	ystart = E->lat[0];//-90.0;
	yend = E->lat[E->nc.y-1];//90.0;



	dx = (xend-xstart)/(double)(nx);
	dy = (yend-ystart)/(double)(ny);
	*/

	// test passive particle advection

	xstart = E->lon[0];//0.0;
	xend = E->lon[E->nc.x-1];//358.0;
	ystart = E->lat[0];//-90.0;
	yend = E->lat[E->nc.y-1];//90.0;



	sprintf(fname,"pos_vel.dat");
	out = fopen(fname,"w");

	// multiple particles

	srand ( time(NULL) );
	double center_lat =  -47.336 ;
  double center_lon =  144.8669;
	double box_width = 0.0;
  min_lon = center_lon - box_width;     max_lon = center_lon + box_width;
  min_lat = center_lat - box_width;     max_lat = center_lat + box_width;

	int	nParticles = 1;
	fprintf(out,"#%d %d\n", nParticles, E->nc.t);
	for(k=0;k<nParticles;k++){


		pos[0] = RandomFloat(min_lon, max_lon) ;
		pos[1] = RandomFloat(min_lat, max_lat) ;

		printf("particle %d: pos = %f, %f\n", k, pos[0], pos[1]);

		// print the original particle position to the output file
		fprintf(out,"%f %f 0.0 0.0\n", pos[0], pos[1]);

		current_time = 0.0;
		t = 0;
		//for(t=0;t<E->nc.t;t++){
		do{

			// find maximum velocity for this time level
			max_vel = -10000.0;
			for(i=0;i<E->ny;i++){
				for(j=0;j<E->nx;j++){

					if(E->u[t][i][j] > max_vel )
						max_vel = E->u[t][i][j];

					if(E->v[t][i][j] > max_vel )
						max_vel = E->v[t][i][j];

				}
			}

			printf("t = %d, max_velocity is = %f\n", t, max_vel);

			// convert pos from lon lat to cartesian
			//x = R * cos(lat) * cos(lon)
			//y = R * cos(lat) * sin(lon)
			//z = R * sin(lat)

			// convert back
			// lon = atan2(y, x)
			// lat = asin(z / R)

			//cart_pos[0] = R*cos(d2r(pos[1]))*cos(d2r(pos[0]));	// x in conversion formula
			//cart_pos[1] = R*cos(d2r(pos[1]))*sin(d2r(pos[0]));	// y in conversion formula
			//cart_pos[2] = R*sin(d2r(pos[1]));										// z in conversion formula

			// simple 2d mapping to equidistant rectangular map projection
			// x = R * lon *  cos(lat)
			// y = R * lat
			cart_pos[0] = R * d2r(pos[0]) * cos(d2r(pos[1]));
			cart_pos[1] = R * d2r(pos[1]);

			printf("## pre: pos: %f, %f\n", pos[0],pos[1]);
			printf("## becomes: cart: %f, %f\n", cart_pos[0], cart_pos[1]);

			// convert cart_pos back to lon lat
			//pos[0] = r2d(atan2(cart_pos[1], cart_pos[0]));	// calculate longitude first
			//pos[1] = r2d(asin(cart_pos[2]/R));
			// convert back:
			// lat = y/R;
			// lon = x / (R*cos(lat))
			pos[1] = r2d(cart_pos[1]/R);
			pos[0] = r2d(cart_pos[0] / (R*cos(d2r(pos[1]))));

			printf("## conv back: pos: %f, %f\n", pos[0],pos[1]);
			//printf("## becomes: cart: %f, %f\n\n", cart_pos[0], cart_pos[1]);
			//exit(1);

			// set up the nodal value for u
			el = 0;
			for(i=0;i<E->ny;i++){
				for(j=0;j<E->nx;j++){

					E->el[el].node_value[0] = E->u[t][i][j];
					E->el[el].node_value[1] = E->u[t][i][j+1];
					E->el[el].node_value[2] = E->u[t][i+1][j+1];
					E->el[el].node_value[3] = E->u[t][i+1][j];

					el++;
				}
			}

			// find out which element this lies within
			el = get_owner_element(E, pos);
			calculate_interpolation_weights(&E->el[el], E->xi, E->eta, pos);
			// interpolate the nodal velocity to the current point
			interpolate_point(&E->el[el], &interp_value);
			vel[0] = interp_value;	// velocity is in metres per second

			// set up nodal value for v
			el = 0;
			for(i=0;i<E->ny;i++){
				for(j=0;j<E->nx;j++){

					E->el[el].node_value[0] = E->v[t][i][j];
					E->el[el].node_value[1] = E->v[t][i][j+1];
					E->el[el].node_value[2] = E->v[t][i+1][j+1];
					E->el[el].node_value[3] = E->v[t][i+1][j];

					el++;
				}
			}
			el = get_owner_element(E, pos);
			// interpolate velocity to this point
			interpolate_point(&E->el[el], &interp_value);
			vel[1] = interp_value;	// velocity is in metres per second

			// advect the partcle
			//dt = 10800.0; // ofam surface fields are 3 hourly
			get_dt( &dt, max_vel, min_mesh );
			current_time += dt;
			printf("dt = %f\nu = %f, v = %f\nx = %f, y = %f\n", dt, vel[0], vel[1], cart_pos[0], cart_pos[1]);
			update_particle_position_euler(cart_pos,vel,dt);
			printf("new position after advection:\n\tx = %f, y = %f\n\n", cart_pos[0], cart_pos[1]);
			// convert cart_pos back to lon lat
			pos[1] = r2d(cart_pos[1]/R);	// calc lat first
			pos[0] = r2d(cart_pos[0] / (R*cos(d2r(pos[1]))));	// now lon
			printf("new position after advection:\n\tlon = %f, lat = %f\n\n", pos[0], pos[1]);


			fprintf(out,"%f %f %f %f\n", pos[0], pos[1], vel[0], vel[1]);
			// figure out what time level we are at
			t = floor(current_time/10800.0);
			printf("**** current_time = %f, dt = %f, t = %d\n",current_time,dt,t);
		}while(t<E->nc.t);
		fprintf(out,"\n");
	}

	fprintf(out,"\n");
	fclose(out);


}

int main( int argc, char *argv[] )
{
	e	*E;

	E = malloc(sizeof(e));


	// single element interpolation test
	//do_test(E);

	// analytic funtion interpolation test
	//test_interp(E);

	// netcdf accessg interpolation tests
	if(argc >= 2){
		test_interp_netcdf(E, argv[1], argv[2]);
	}
	else{
		printf("no input netcdf file specified\n");
		printf("doing nothing\n");
	}



	//
	// calculate_interpolation_weights(&E->el[0], E->xi, E->eta, E->pos);
	//
	// interpolate_point(&E->el[0], &interp_value);

	return 0;
}




void interpolate_point(element *el, double *interp_value){

	int i;

	*interp_value = 0.0;
	for(i=0; i < 4; i++)
		*interp_value += el->node_value[i] * el->interp_weights[i];

}


void calculate_interpolation_weights(element *el, double *xi, double *eta, double *pos){

	int i;
	double	pos_local[2];

	// transform from global interpolation position to local element coordinates
	pos_local[0] = (pos[0] - 0.5*(el->node_coord[1][0] + el->node_coord[0][0]))/( 0.5 * (el->node_coord[1][0]-el->node_coord[0][0]) );
	pos_local[1] = (pos[1] - 0.5*(el->node_coord[3][1] + el->node_coord[0][1]))/( 0.5 * (el->node_coord[3][1]-el->node_coord[0][1]) );

	for(i=0; i < 4; i++)
		el->interp_weights[i] = evaluate_linear_quad_shape_function( xi, eta, pos_local, i );

}

double evaluate_linear_quad_shape_function( double *xi, double *eta, double *pos_to_interp_at, int node ){

	return ( 0.25 * (1.0 + pos_to_interp_at[0] * xi[node]) * (1.0 + pos_to_interp_at[1] * eta[node]) );

}
