#include "airfoil_generator.c"

int main() {
    int N = 100;
    int m = 2, p = 4, t = 15;
    double c = 1.0;
    Airfoil_t af;

    generate_airfoil_4digit(2412, c, N, &af);

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