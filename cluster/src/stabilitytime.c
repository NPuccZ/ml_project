#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void henon(double *x, double *p_x, double *y, double *p_y, double nu_x, double nu_y, long i, double eps);
double stabilityTime(double x0, double y0, double nu_x, double nu_y, double eps, double r_c, long N_iterations);

int main(int argc, char *argv[])
{
    if (argc != 7)
    {
        printf("Usage: %s nu_x nu_y eps N_steps lattice_length N_iterations\n", argv[0]);
        return 1;
    }

    clock_t start_time = clock();

    double nu_x          = atof(argv[1]);
    double nu_y          = atof(argv[2]);
    double eps           = atof(argv[3]);
    int    N_steps       = atoi(argv[4]);
    double lattice_length = atof(argv[5]);
    long   N_iterations  = atol(argv[6]);
    double r_c           = 10000.;

    // Filename encodes all parameters to avoid collisions between jobs
    char filename[256];
    sprintf(filename, "stabtime_nux%s_nuy%s_eps%s_N%s.txt", argv[1], argv[2], argv[3], argv[6]);
    FILE *out = fopen(filename, "w");
    if (out == NULL) {
        perror("Error opening output file");
        return 1;
    }

    for (int i = 0; i < N_steps; i++)
    {
        for (int j = 0; j < N_steps; j++)
        {
            double x0 = (double)i * (lattice_length / (double)N_steps);
            double y0 = (double)j * (lattice_length / (double)N_steps);
            fprintf(out, "%lf\t%lf\t%lf\n", x0, y0,
                    stabilityTime(x0, y0, nu_x, nu_y, eps, r_c, N_iterations));
        }
    }

    fclose(out);

    double time_spent = (double)(clock() - start_time) / CLOCKS_PER_SEC;
    printf("Done. Output: %s  (%.1f s)\n", filename, time_spent);

    return 0;
}

void henon(double *x, double *p_x, double *y, double *p_y, double nu_x, double nu_y, long i, double eps)
{
    double omega1 = 2. * M_PI / 868.12;
    double tk[7]      = {omega1, 2.*omega1, 3.*omega1, 6.*omega1, 7.*omega1, 10.*omega1, 12.*omega1};
    double epsilon[7] = {1.E-4, 0.218E-4, 0.708E-4, 0.254E-4, 0.100E-4, 0.078E-4, 0.218E-4};

    double mod = epsilon[0]*cos(tk[0]*i) + epsilon[1]*cos(tk[1]*i) + epsilon[2]*cos(tk[2]*i)
               + epsilon[3]*cos(tk[3]*i) + epsilon[4]*cos(tk[4]*i) + epsilon[5]*cos(tk[5]*i)
               + epsilon[6]*cos(tk[6]*i);

    double omega_x = 2. * M_PI * nu_x * (1. + eps * mod);
    double omega_y = 2. * M_PI * nu_y * (1. + eps * mod);

    double x_old  = *x,  px_old = *p_x;
    double y_old  = *y,  py_old = *p_y;

    *x   =  cos(omega_x) * x_old + sin(omega_x) * (px_old + x_old*x_old - y_old*y_old);
    *p_x = -sin(omega_x) * x_old + cos(omega_x) * (px_old + x_old*x_old - y_old*y_old);
    *y   =  cos(omega_y) * y_old + sin(omega_y) * (py_old - 2.*x_old*y_old);
    *p_y = -sin(omega_y) * y_old + cos(omega_y) * (py_old - 2.*x_old*y_old);
}

double stabilityTime(double x0, double y0, double nu_x, double nu_y, double eps, double r_c, long N_iterations)
{
    double x = x0, p_x = 0.;
    double y = y0, p_y = 0.;

    for (long k = 1; k <= N_iterations; k++)
    {
        henon(&x, &p_x, &y, &p_y, nu_x, nu_y, k, eps);
        if (x*x + p_x*p_x + y*y + p_y*p_y >= r_c)
            return log10((double)k);
    }

    return log10((double)N_iterations);
}
