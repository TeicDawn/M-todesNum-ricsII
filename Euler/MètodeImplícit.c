#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Definimos las mismas variables que en el caso del método explícito.
// La única nueva es el número de iteraciones, K.

#define Nx 100
#define K 100

const double x_max = 1.0;
double delta_x = x_max/Nx;

const double t_max = 2/(4*pow(M_PI,2));
double delta_t = 0.5*delta_x;
double G = delta_t/pow(delta_x,2);
int Nt = t_max/delta_t;

double Tf = 10.0;
double Ta = 300.0;
double r=1.0;

double anal (double x){
	return(Ta + Tf*(0.75*cos(x)*exp(-t_max*4*pow(M_PI,2)) + 0.25*cos(3*x)*exp(-9*t_max*4*pow(M_PI,2))));
}

int main(){
	int i;
	int j;
	int k;
	
	double T[Nx][Nt];
	
	// Las condiciones iniciales también son las mismas.
	
	for(i = 0; i < Nx; i++){
		T[i][0] = pow(cos(2*M_PI*i*delta_x),3);
	}
	
	// Ahora, en el for sobre i dentro del for sobre j, los valores de la siguiente
	// columna de la matriz (el siguiente valor de j) se obtienen de forma asintótica
	// y no directa. Así pues, introducimos un for más con una variable muda k
	// que itere el for de las i unas K veces para cada j.
	
	for(j = 1; j < Nt; j++){
		
		for (i = 0; i < Nx; i++){
			T[i][j]=T[i][j-1];
		}
		
		for (k = 0; k < K; k++){
			T[0][j] = (1/(1+(2*G)))*(T[0][j-1]+(G*(T[1][j]+T[Nx-1][j])));
			
			for(i = 1; i < Nx-1; i++){
				T[i][j] = (1/(1+(2*G)))*(T[i][j-1]+(G*(T[i+1][j]+T[i-1][j])));
			}
			
			T[Nx-1][j] = (1/(1+(2*G)))*(T[Nx-1][j-1]+(G*(T[0][j]+T[Nx-2][j])));
		}	
	}
	
	// Debemos recordar que el error es ahora sensible tanto a G como a K.
	
	for(i = 0; i < Nx; i++){
		printf("%lf %lf %lf %lf %lf\n", (i*delta_x)*2*M_PI*r, (delta_t*Nt)*pow(2*M_PI,2), Ta+Tf*T[i][Nt-1], anal((i*delta_x)*2*M_PI), fabs(anal((i*delta_x)*2*M_PI)-Ta-Tf*T[i][Nt-1]));
    }
}
