//Author: Marta Razza
//This code allows to compute the Lyapunov Error and the reversibility error method for the 4D Henon map.
//As for the 2D Henon map, the aim is to study the stability of the the orbits for a particle in a circular accelerator.

//The modulation is due to unavoidable external effects 
//such as the coupling between synchrotron and betatron motion 
//(respectively motion on longitudinal and transversal plane) 
//or power supply ripple and can be modelled by a set of non-linear oscillations

//For further details, see the report in the repository.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 


//Defining the functions
void henon(double *x, double *p_x, double *y, double *p_y, double nu_x, double nu_y, int i, double eps);
void tangent(double *eta_x, double *eta_px, double *eta_y, double *eta_py, double x, double p_x, double y, double p_y, double nu_x, double nu_y, int i, double eps);
void henonInv(double *x, double *p_x, double *y, double *p_y, double nu_x, double nu_y, int i, double eps);
double calculateREM(double x, double y, double p_x, double p_y, double x0, double y0, double px_0, double py_0);
double calculateLE(double eta1x, double eta1px, double eta1y, double eta1py, double eta2x, double eta2px, double eta2y, double eta2py, double eta3x, double eta3px, double eta3y, double eta3py, double eta4x, double eta4px, double eta4y, double eta4py);

int main(int argc, char *argv[])
{

//Checking the number of arguments
	if (argc != 7)
	{
		printf("Usage: %s nu_x nu_y eps N_steps lattice_length N_turns\n", argv[0]);
		return 1;
	}

//Start measuring time
	clock_t start_time, end_time;
	start_time = clock();

//Opening the files
char filename [50];
	sprintf(filename, "lyapunov_nux%s_nuy%s_eps%s_n%s.txt", argv[1], argv[2],argv[3], argv[6]);
	FILE *lyap = fopen(filename, "w");
	if (lyap == NULL) {
    	perror("Error opening file");
    	return 1;
	}
char filename2 [50];
	sprintf(filename2, "reversibility_nux%s_nuy%s_eps%s_n%s.txt", argv[1], argv[2],argv[3], argv[6]);
	FILE *REM = fopen(filename2, "w");
	if (REM == NULL) {
		perror("Error opening file");
		return 1;
	}
//Defining the variables
	double nu_x; //tune along x
	double nu_y; //tune along y
	double eps;  //modulation parameter
	nu_x = atof(argv[1]); 
	nu_y = atof(argv[2]);
	eps = atof(argv[3]);
	int N_steps = atoi(argv[4]); //number of steps to explore the entire lattice
	double lattice_length = atof(argv[5]); //length of the lattice
	int N_turns = atoi(argv[6]); //number of turns


	double x0, px_0, y0, py_0, x, p_x, y, p_y, eta1x, eta1px, eta1y, eta1py, eta2x, eta2px, eta2y, eta2py, eta3x, eta3px, eta3y, eta3py, eta4x, eta4px, eta4y, eta4py;
	int k, i, j, h, f, l;


	//x and y lattice
	for (i = 0; i < N_steps; i++) //cycle for x
	{
		for (j = 0; j < N_steps; j++) //cycle for y
		{
			x0 = (double)i * (lattice_length/(double)N_steps);
			y0 = (double)j * (lattice_length/(double)N_steps);
			px_0 = 0.;
			py_0 = 0.;
			eta1x = 1;
			eta1px = 0;
			eta1y = 0;
			eta1py = 0;
			eta2x = 0;
			eta2px = 1;
			eta2y = 0;
			eta2py = 0;
			eta3x = 0;
			eta3px = 0;
			eta3y = 1;
			eta3py = 0;
			eta4x = 0;
			eta4px = 0;
			eta4y = 0;
			eta4py = 1;
			x = x0;
			y = y0;
			p_x = px_0;
			p_y = py_0;
//Cycles for the number of turns
//Evolving the map forward and backward to obtain the reversibility error method
			for (k = 0; k < N_turns; k++)
			{
				henon(&x, &p_x, &y, &p_y, nu_x, nu_y, k, eps);
			}
			for (h = 0; h < N_turns; h++)
			{
				henonInv(&x, &p_x, &y, &p_y, nu_x, nu_y, h, eps);
			}

			fprintf(REM, "%lf\t%lf\t%lf\n", x0, y0, log10(calculateREM(x, y, p_x, p_y, x0, y0, px_0, py_0)));
//Reset to initial conditions before computing Lyapunov Error
			x = x0;
			y = y0;
			p_x = px_0;
			p_y = py_0;
//Evolving the tangent map to obtain the Lyapunov Error
			for (f = 0; f < N_turns; f++)
			{
				tangent(&eta1x, &eta1px, &eta1y, &eta1py, x, p_x, y, p_y, nu_x, nu_y, f, eps);
				tangent(&eta2x, &eta2px, &eta2y, &eta2py, x, p_x, y, p_y, nu_x, nu_y, f, eps);
				tangent(&eta3x, &eta3px, &eta3y, &eta3py, x, p_x, y, p_y, nu_x, nu_y, f, eps);
				tangent(&eta4x, &eta4px, &eta4y, &eta4py, x, p_x, y, p_y, nu_x, nu_y, f, eps);
				henon(&x, &p_x, &y, &p_y, nu_x, nu_y, f, eps);
			}

			fprintf(lyap, "%lf\t%lf\t%lf\n", x0, y0, log10(log10(calculateLE(eta1x, eta1px, eta1y, eta1py, eta2x, eta2px, eta2y, eta2py, eta3x, eta3px, eta3y, eta3py, eta4x, eta4px, eta4y, eta4py))));
			printf("Processing lattice positions: x0=%d\t y0=%d\n", i, j);
			
		}
	}
//Closing the files
	printf("Processing terminated successfully\n");
	fclose(lyap);
	fclose(REM);

//End measuring time
	end_time = clock();
	double time_spent = (double)(end_time-start_time)/CLOCKS_PER_SEC;
	printf("Time taken to run the program: %f seconds\n", time_spent);

	return 0;
	
}
//function to compute the 4D Henon map
void henon(double *x, double *p_x, double *y, double *p_y, double nu_x, double nu_y, int i, double eps)
{
	double omega_x;
	double omega_y;
	double omega1 = 2. * M_PI / 868.12;
	double tk[7] = {omega1, 2. * omega1, 3. * omega1, 6. * omega1, 7. * omega1, 10. * omega1, 12. * omega1}; // omega k
	double epsilon[7] = {1.E-4, 0.218E-4, 0.708E-4, 0.254E-4, 0.100E-4, 0.078E-4, 0.218E-4};				 // parametro epsilon k
	double t0_x = 2. * M_PI * nu_x;
	double t0_y = 2. * M_PI * nu_y;

	omega_x = t0_x * (1. + eps * (epsilon[0] * cos(tk[0] * i) + epsilon[1] * cos(tk[1] * i) + epsilon[2] * cos(tk[2] * i) + epsilon[3] * cos(tk[3] * i) + epsilon[4] * cos(tk[4] * i) + epsilon[5] * cos(tk[5] * i) + epsilon[6] * cos(tk[6] * i)));
	omega_y = t0_y * (1. + eps * (epsilon[0] * cos(tk[0] * i) + epsilon[1] * cos(tk[1] * i) + epsilon[2] * cos(tk[2] * i) + epsilon[3] * cos(tk[3] * i) + epsilon[4] * cos(tk[4] * i) + epsilon[5] * cos(tk[5] * i) + epsilon[6] * cos(tk[6] * i)));

	double x_old;
	double y_old;
	double px_old;
	double py_old;

	x_old = *x;
	px_old = *p_x;
	y_old = *y;
	py_old = *p_y;

	*x = x_old * cos(omega_x) + (px_old + (x_old * x_old - y_old * y_old)) * sin(omega_x);
	*p_x = -sin(omega_x) * x_old + cos(omega_x) * (px_old + (x_old * x_old - y_old * y_old));
	*y = y_old * cos(omega_y) + (py_old - (2. * x_old * y_old)) * sin(omega_y);
	*p_y = -sin(omega_y) * y_old + cos(omega_y) * (py_old - (2. * x_old * y_old));
}
//function to compute the tangent map
void tangent(double *eta_x, double *eta_px, double *eta_y, double *eta_py, double x, double p_x, double y, double p_y, double nu_x, double nu_y, int i, double eps)
{
	double eta_x_old = *eta_x;
	double eta_px_old = *eta_px;
	double eta_y_old = *eta_y;
	double eta_py_old = *eta_py;

	double omega_x;
	double omega_y;
	double omega1 = 2. * M_PI / 868.12;
	double tk[7] = {omega1, 2. * omega1, 3. * omega1, 6. * omega1, 7. * omega1, 10. * omega1, 12. * omega1}; // omega k
	double epsilon[7] = {1.E-4, 0.218E-4, 0.708E-4, 0.254E-4, 0.100E-4, 0.078E-4, 0.218E-4};				 // parametro epsilon k
	double t0_x = 2. * M_PI * nu_x;
	double t0_y = 2. * M_PI * nu_y;

	omega_x = t0_x * (1. + eps * (epsilon[0] * cos(tk[0] * i) + epsilon[1] * cos(tk[1] * i) + epsilon[2] * cos(tk[2] * i) + epsilon[3] * cos(tk[3] * i) + epsilon[4] * cos(tk[4] * i) + epsilon[5] * cos(tk[5] * i) + epsilon[6] * cos(tk[6] * i)));
	omega_y = t0_y * (1. + eps * (epsilon[0] * cos(tk[0] * i) + epsilon[1] * cos(tk[1] * i) + epsilon[2] * cos(tk[2] * i) + epsilon[3] * cos(tk[3] * i) + epsilon[4] * cos(tk[4] * i) + epsilon[5] * cos(tk[5] * i) + epsilon[6] * cos(tk[6] * i)));

	*eta_x = (cos(omega_x) + 2. * x * sin(omega_x)) * eta_x_old + sin(omega_x) * eta_px_old - 2. * y * sin(omega_x) * eta_y_old;
	*eta_px = (2. * x * cos(omega_x) - sin(omega_x)) * eta_x_old + cos(omega_x) * eta_px_old - 2. * y * cos(omega_x) * eta_y_old;
	*eta_y = -2. * y * sin(omega_y) * eta_x_old + (cos(omega_y) - 2. * x * sin(omega_y)) * eta_y_old + sin(omega_y) * eta_py_old;
	*eta_py = (-2. * x * cos(omega_y) - sin(omega_y)) * eta_y_old - 2. * y * cos(omega_y) * eta_x_old + cos(omega_y) * eta_py_old;
}
//function to compute the inverse map
void henonInv(double *x, double *p_x, double *y, double *p_y, double nu_x, double nu_y, int i, double eps)
{
	double x_old;
	double y_old;
	double px_old;
	double py_old;
	double omega_x;
	double omega_y;
	double omega1 = 2. * M_PI / 868.12;
	double tk[7] = {omega1, 2. * omega1, 3. * omega1, 6. * omega1, 7. * omega1, 10. * omega1, 12. * omega1}; // omega k
	double epsilon[7] = {1.E-4, 0.218E-4, 0.708E-4, 0.254E-4, 0.100E-4, 0.078E-4, 0.218E-4};				 // parametro epsilon k
	double t0_x = 2. * M_PI * nu_x;
	double t0_y = 2. * M_PI * nu_y;

	omega_x = t0_x * (1. + eps * (epsilon[0] * cos(tk[0] * i) + epsilon[1] * cos(tk[1] * i) + epsilon[2] * cos(tk[2] * i) + epsilon[3] * cos(tk[3] * i) + epsilon[4] * cos(tk[4] * i) + epsilon[5] * cos(tk[5] * i) + epsilon[6] * cos(tk[6] * i)));
	omega_y = t0_y * (1. + eps * (epsilon[0] * cos(tk[0] * i) + epsilon[1] * cos(tk[1] * i) + epsilon[2] * cos(tk[2] * i) + epsilon[3] * cos(tk[3] * i) + epsilon[4] * cos(tk[4] * i) + epsilon[5] * cos(tk[5] * i) + epsilon[6] * cos(tk[6] * i)));

	x_old = *x;
	y_old = *y;
	px_old = *p_x;
	py_old = *p_y;

	*x = x_old * cos(omega_x) - px_old * sin(omega_x);
	*y = y_old * cos(omega_y) - py_old * sin(omega_y);
	*p_x = (*y) * (*y) - (*x) * (*x) + x_old * sin(omega_x) + px_old * cos(omega_x);
	*p_y = 2. * (*x) * (*y) + sin(omega_y) * y_old + py_old * cos(omega_y);
}
//function to calculate the reversibility error method
double calculateREM(double x, double y, double p_x, double p_y, double x0, double y0, double px_0, double py_0)
{
	double sq2 = 1.4142135623;
	double EPS = 1E-16;
	double erREM;
	erREM = (sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0) + (p_x - px_0) * (p_x - px_0) + (p_y - py_0) * (p_y - py_0)) / (sq2 * EPS));
	return erREM;
}
//function to calculate the Lyapunov Error
double calculateLE(double eta1x, double eta1px, double eta1y, double eta1py, double eta2x, double eta2px, double eta2y, double eta2py, double eta3x, double eta3px, double eta3y, double eta3py, double eta4x, double eta4px, double eta4y, double eta4py)
{
	double LE;
	double sq2 = 1.4142135623;
	LE = (sqrt(eta1x * eta1x + eta1px * eta1px + eta2x * eta2x + eta2px * eta2px + eta3x * eta3x + eta3px * eta3px + eta4x * eta4x + eta4px * eta4px + eta1y * eta1y + eta1py * eta1py + eta2y * eta2y + eta2py * eta2py + eta3y * eta3y + eta3py * eta3py + eta4y * eta4y + eta4py * eta4py));
	return LE;
}
