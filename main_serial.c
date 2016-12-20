#include <stdio.h>
#include <math.h>

#define J 1000
#define alpha 1.
#define TIME_INTERVAL 10.
#define dt (alpha / J) // Courant Condition
#define N (TIME_INTERVAL / dt)


#define TWO_DOT(y) (alpha * alpha * J * J * (y[j + 1] + y[j - 1] - 2. * y[j]))

int main()
{
    int n;
    int j;
    double y[J + 1], v[J + 1], a[J + 1], da[J + 1], dda[J + 1];
    double y_next[J + 1], v_next[J + 1];
    FILE *fp;
    
    // Record of parameters
    fp = fopen("parameters.txt", "w");
    fprintf(fp, "%d %d %e", (int)N, J, TIME_INTERVAL);
    fclose(fp);

    
    // Initialization of functions
    y[0] = y[J] = v[0] = v[J] = a[0] = a[J] = 0.;
    for (j = 1; j < J; j++) {
        if (j < J/2) {
            y[j] = pow(sin(2. * M_PI * j / J), 2.);
            v[j] = - 4. * M_PI * alpha * cos(2. * M_PI * j / J) * sin(2. * M_PI * j / J);
        } else {
            y[j] = 0.;
            v[j] = 0.;
        }
    }

    // Evolution of functions
    fp = fopen("output", "wb");
    n = 0;
    do {
        for (j = 1; j < J; j++) {
              a[j] = TWO_DOT(y);
        }
        for (j = 1; j < J; j++) {
             da[j] = TWO_DOT(v);
            dda[j] = TWO_DOT(a);
        }
        
        fwrite(y, sizeof(double), J + 1, fp);
        
        if (n++ >= N) break;
        
        for (j = 1; j < J; j++) {
            y_next[j] = y[j] + v[j] * dt +  a[j] * dt * dt / 2. +  da[j] * dt * dt * dt / 4.;
            v_next[j] = v[j] + a[j] * dt + da[j] * dt * dt / 2. + dda[j] * dt * dt * dt / 4.;
        }
        for (j = 1; j < J; j++) {
            y[j] = y_next[j];
            v[j] = v_next[j];
        }
    } while (1);
    fclose(fp);
    
    return 0;
}
