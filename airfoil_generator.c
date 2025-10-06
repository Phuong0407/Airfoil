#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
  double* x_c;
  double* y_c;
  double* x_u;
  double* y_u;
  double* x_l;
  double* y_l;
  double c;
  int N;
} Airfoil_t;

void generate_airfoil_4digit(int m, int p, int t, double c, int N, Airfoil_t *af) {
    double m_real = m/100.0;
    double p_real = p/10.0;
    double t_real = t/100.0;

    af->x_c = (double*)malloc((N+1)*sizeof(double));
    af->y_c = (double*)malloc((N+1)*sizeof(double));
    af->x_u = (double*)malloc((N+1)*sizeof(double));
    af->y_u = (double*)malloc((N+1)*sizeof(double));
    af->x_l = (double*)malloc((N+1)*sizeof(double));
    af->y_l = (double*)malloc((N+1)*sizeof(double));
    af->c = c;
    af->N = N;

    for (int i = 0; i <= N; ++i) {
      double xp = (double)i/(double)N;
      double x = xp * c;
      af->x_c[i] = x;

      double A = 0.0;
      double dyc = 0.0;
      if (xp <= p_real) {
        A = m_real/(p_real*p_real);
        af->y_c[i] = c * A * (2.0*p_real*xp - xp*xp);
        dyc = 2.0 * A * (p_real - xp);
      } else {
        A = m_real/((1-p_real)*(1-p_real));
        af->y_c[i] = c * A * (1.0 - 2.0*p_real + 2.0*p_real*xp - xp*xp);
        dyc = 2.0 * A * (p_real - xp);
      }

      double theta = atan(dyc);
      double y_t = 5.0*t_real*c*(0.2969*sqrt(xp) - 0.1260*xp - 0.3516*xp*xp + 0.2843*xp*xp*xp - 0.1015*xp*xp*xp*xp);

      double sintheta = sin(theta);
      double costheta = cos(theta);

      af->x_u[i] = af->x_c[i] - y_t*sintheta;
      af->y_u[i] = af->y_c[i] + y_t*costheta;
      af->x_l[i] = af->x_c[i] + y_t*sintheta;
      af->y_l[i] = af->y_c[i] - y_t*costheta;
  }
}

void plot_airfoil(const char* filename) {
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (!gnuplotPipe) {
        perror("Cannot open gnuplot");
        return;
    }

    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set grid\n");
    fprintf(gnuplotPipe, "plot '%s' with lines lw 2 title 'Airfoil'\n", filename);
    
    fflush(gnuplotPipe);
}

int main() {
    int N = 100;
    int m = 2, p = 4, t = 15;
    double c = 1.0;
    Airfoil_t af;

    generate_airfoil_4digit(m, p, t, c, N, &af);

    FILE *fp = fopen("naca2412.dat", "w");
    if (!fp) {
        perror("Cannot open file");
        return 1;
    }

    for (int i = 0; i <= N; i++) {
        fprintf(fp, "%f %f\n", af.x_u[i], af.y_u[i]);
    }
    for (int i = N - 1; i >= 0; i--) {
        fprintf(fp, "%f %f\n", af.x_l[i], af.y_l[i]);
    }

    fclose(fp);

    printf("Data exported to file.dat\n");

    plot_airfoil("naca2412.dat");

    free(af.x_c); free(af.y_c);
    free(af.x_u); free(af.y_u);
    free(af.x_l); free(af.y_l);

    return 0;
}