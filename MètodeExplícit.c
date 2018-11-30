#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Definimos las constantes que nos dan del mallado. El enunciado nos da Nx y gamma,
// y calculamos con ellas la parte entera de Nt buscando t=2s. El tiempo m�ximo real ser� Nt para delta_t.
// Tendremos una discretizaci�n de 100 puntos para el espacio que ir� de (0,1) y 993
// para el tiempo (0, 2/(2pi^2)) en las variables normalizadas.


#define Nx 100
#define Nt 1034
#define gamma 0.49

double x_max = 1.0;
double delta_x = x_max/Nx;

double delta_t = gamma*(delta_x*delta_x);
double t_max = delta_t*Nt;

double Tf = 10.0;
double Ta = 300.0;
double R = 1;

// Definimos la funci�n anal(x) como la soluci� anal�tica, que usaremos
// para calcular el error absoluto real.

double anal (double x){
    return(Ta + Tf*(0.75*cos(2*M_PI*x)*exp(-t_max*4*pow(M_PI,2)) + 0.25*cos(6*M_PI*x)*exp(-9*t_max*4*pow(M_PI,2))));
}

// Ahora comenzamos el programa. Consistir� en unos bucles que
// describen la evoluci�n temporal en una matriz y un print
// para obtener los datos que nos interesan.

int main(){
	int i;
	int j;
	
	double T[Nx][Nt];
	
	// Llenamos la primera fila con las condiciones iniciales
	// con un for sobre i.
	
	for(i = 0; i < Nx; i++){
	T[i][0] = pow(cos(2*M_PI*i*delta_x),3);
	}
	
	// Ahora, con un for sobre j, hacemos correr el tiempo.
	
	for(j = 0; j < Nt-1; j++){
		
		// Para obtener la temperatura a cada posic�n sabiendo
		// las del tiempo anterioir, hemos de escribir Nx ecuaciones.
		// Dos de elles las "hacemos a mano" porque neesitamos imponer la condici�n de contorno
		// de periodicidad. El resto, las abreviamos con un for sobre i.
		
		T[0][j+1] = T[0][j]+(gamma*(T[1][j]-2*T[0][j]+T[Nx-1][j]));
		T[Nx-1][j+1] = T[Nx-1][j]+(gamma*(T[0][j]-2*T[Nx-1][j]+T[Nx-2][j]));
	    
		for(i = 1; i < Nx-1; i++){
			T[i][j+1] = T[i][j]+(gamma*(T[i+1][j]-2*T[i][j]+T[i-1][j]));
		}
	}
	
	// Ahora que tenemos la evoluci�n temporal entera guardada en la matriz T[i][j],
	// definimos las variables Tf, Ta, r para deshacer la normalizaci�n de variables 
	// y hacemos print de la posic�n x = i*delta_x y la correspondiente temperatura final T[i][Nt-1] all�
	// A�adimos el print para saber el tiempo m�s cercano a 2 seg que le corresponde a la soluci�n num�rica
	// encontrada pues es el valor que se sustituir� en la soluci�n anal�tca para comparar el resultado.
	
	for(i = 0; i < Nx; i++){
		printf("%lf %lf %lf %lf %lf\n", i*delta_x, i*delta_x*2*M_PI*R, anal(i*delta_x), Ta+Tf*T[i][Nt-1],
		fabs(Ta+Tf*T[i][Nt-1]-anal(i*delta_x)));
	}
}