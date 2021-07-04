//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//                                            ECUACIÓN DE SCHRÖDINGER
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"    //Librería para utilizar los números complejos 
#include "gsl_rng.h"    //Librería para generación de números aleatorios


#define N 1500
#define T 6000              //Valor máximo del tiempo (fotogramas)
#define nciclos 150          //Esta variable va de 1 a N/4 como máximo
#define nD 1500              //Tiempo que pasa hasta que el detector 
#define iteraciones 1000   //Número de iteraciones del programa
#define lambda 0.7
#define PI 3.14159265359



gsl_rng *tau;   //Definimos el puntero para poder utilizarlo en todas las funciones--> variable externa


double k0,s,norma;
double V[N]; 
double Pd,Pi;
double x,k,y;
int contadornd;
int mt;
double K;   //coeficiente de transmisión

fcomplex phi[N];  //N+1 para tener en cuenta las condiciones de contorno
fcomplex chi[N];
fcomplex alpha[N-1];
fcomplex beta[N-1];
fcomplex A0[N];

FILE *f1;


//Definición de funciones

int CalcAlpha();
int CalcBeta();
int CalcChi();
int PrintfCondIni();




/////////////////////////////////////////////////////// FUNCIÓN PRINCIPAL //////////////////////////////////////////////////////////////////

int main()
{
    //PARA LA GENERACIÓN DE NÚMEROS ALEATORIOS

    extern gsl_rng *tau; //Puntero al estado del número aleatorio
    int semilla=278349; //Semilla del generador de números aleatorios


    tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero
    gsl_rng_set(tau,semilla); //Inicializamos la semilla

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //En primer lugar calculamos las condiciones iniciales para t=0
    //
    //Calculamos---> Constantes k0, s
    //               El potencial V
    //               El vector phi para t=0, para representarlo haremos uso del módulo al cuadrado que es un vector de reales
    //               Calculamos el valor del vector de alfas, que no depende del tiempo
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Abrimos los ficheros 

    f1 = fopen("Potencial.txt","w");


    //Damos valores a los parámetros iniciales (según el guión)

    k0=2*PI*nciclos/N; 
    s=(0.25)/(pow(k0,2));


    //Fijamos las condiciones de contorno 
    phi[0]=phi[N-1]=Complex(0.0,0.0);
    chi[0]=chi[N-1]=Complex(0.0,0.0);

    

    for (int j=0;j<N;j++)
    {

        //Calculamos el potencial cuadrado

        if(j>=((2*N/5)-1) && j<=((3*N/5)-1))
        {
            V[j]=lambda*pow(k0,2);

        }else
        {
            V[j]=0;
        }

        A0[j] = Complex(-2-V[j],2/s);   //Calculamos el vector A0

    }
    
    for(int i=0;i<T;i++)
    {
        //Guardamos el valor del potencial para cada unidad de tiempo en el fichero

        for(int j=0;j<N;j++)
		{
			fprintf(f1,"%i\t%lf\n",j,V[j]);
		}

		fprintf(f1,"\n\n");
    }


    

    mt=0;

    for(int z=0;z<iteraciones;z++) //Para una buena estadística simulamos el sistema 10^3 veces---> m
    {

        //Calculamos el vector phi inicial

        for(int j=1;j<N-1;j++)
        {
            phi[j]=Cgauss(k0*j,exp(-8*pow(4*j-N,2)/pow(N,2)));
        }

        CalcAlpha();

        //PrintfCondIni();

    
    


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    //  EMPEZAMOS EL BUCLE DEL TIEMPO---> Calculamos toda la función phi para cada unidad de tiempo
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
        contadornd=0;

        for(int i=0;i<T;i++)    //Bucle del tiempo
        {
            CalcBeta(); //Vector de coeficientes beta
            CalcChi();  //Calculamos la chi--> para calcular phi


            //Calculamos la phi y la guardamos en el fichero

            for(int j=1;j<N-1;j++)  //De 1 a N-1 por las condiciones de contorno
            {
                phi[j]=Csub(chi[j],phi[j]); 
            }
        
            if(contadornd==nD)  //Cuando evolucione nD pasos calculamos la probabilidad 
            {

                //Calculamos k para normalizar la función phi

                k=0.0;

                for(int j=0;j<N;j++)
                {
                    k=k + pow(Cabs(phi[j]),2);
                }

                y=(1.0)/sqrt(k);

                //Variamos la función de onda dividiendola por raiz de k

                for(int j=0;j<N;j++)
                {
                    phi[j]= RCmul(y,phi[j]);
                }

                //Calculamos la probabilidad derecha Pd

                Pd=0.0;

                for(int i=(4.0*N/5.0);i<N;i++)
                {
                    Pd = Pd + pow(Cabs(phi[i]),2);

                }

                //Generamos un número aleatorio x entre [0:1] y lo comparamos con la probabilidad

                x=gsl_rng_uniform(tau);

                if(x<Pd)
                {
                    mt=mt+1;
                }
                else
                {
                    //Hacemos 0 la función de onda para el intervalo de 4N/5 - N

                    for(int j=(4.0*N/5.0);j<N;j++)
                    {
                        phi[j]=Complex(0.0,0.0);
                    }

                    //Calculamos k

                    k=0.0;

                    for(int j=0;j<N;j++)
                    {
                        k=k + pow(Cabs(phi[j]),2);
                    }

                    y=(1.0)/sqrt(k);

                    //Variamos la función de onda dividiendola por raiz de k

                    for(int j=0;j<N;j++)
                    {
                        phi[j]= RCmul(y,phi[j]);
                    }

                    //Calculamos la probabilidad izquierda Pi

                    Pi=0.0;

                    for(int i=0;i<(N/5);i++)
                    {
                        Pi=Pi + pow(Cabs(phi[i]),2);

                    }

                    //Generamos un número aleatorio x entre [0:1] y lo comparamos con la probabilidad

                    x=gsl_rng_uniform(tau);

                    if(x>Pi)
                    {
                        //Recalculamos la función de onda
                        //Hacemos 0 la función de onda para el intervalo de 0 - N/5

                        for(int j=0;j<(N/5);j++)
                        {
                            phi[j]=Complex(0.0,0.0);
                        }

                        //Calculamos k

                        k=0.0;

                        for(int j=0;j<N;j++)
                        {
                            k=k + pow(Cabs(phi[j]),2);
                        }

                        y=(1.0)/sqrt(k);

                        //Variamos la función de onda dividiendola por raiz de k

                        for(int j=0;j<N;j++)
                        {
                            phi[j]= RCmul(y,phi[j]);
                        }
                        contadornd=(-1);
                    }

                }

            }

            contadornd=contadornd+1;
        }

    }
    
    //Calculamos el coeficiente de transmisión

    K=(1.0*mt)/(1.0*iteraciones);

    printf("El valor de mt es %i y el de m es %i por lo que\n\n",mt,iteraciones);
    printf("El valor del coeficiente de transmisión es: %E \n",K);

    //Cerramos los ficheros

    fclose(f1);

    return 0;
}











//////////////////////////////////////////////////////////// FUNCIONES /////////////////////////////////////////////////////////////////////////

int CalcAlpha()     //Función que calcula el valor del vector de alfas
{
    alpha[N-2]=Complex(0.0,0.0);

    for (int k=N-2;k>0;k--)
    {
        alpha[k-1]=Cdiv(Complex(-1.0,0.0),Cadd(A0[k],alpha[k]));

    }

    return 0;
}





int CalcBeta()      //Función que calcula el valor del vector de betas para cada unidad de tiempo
{
    beta[N-2]=Complex(0.0,0.0);

    for (int k=N-2;k>0;k--)
    {
        beta[k-1]=Cdiv(Csub(Cmul(Complex(0.0,4.0/s),phi[k]),beta[k]),Cadd(A0[k],alpha[k]));

    }

    return 0;
}




int CalcChi()   //Función que calcula el vector de chi 
{
    for(int j=1;j<N-1;j++)
    {
        chi[j] = Cadd(Cmul(alpha[j-1],chi[j-1]),beta[j-1]);
    }
    return 0;
}




int PrintfCondIni()     //Función que imprime por pantalla las condiciones iniciales 
{
    //VECTOR POTENCIAL 

    printf("VECTOR DE POTENCIAL\n\n");          

    for(int i=0; i<N; i++)               
    {
        printf("V[%i]= %lf\n",i, V[i]);
    }
    

    printf("\n\n\n");  

    //VECTOR ALFA

    printf("VECTOR ALFA\n\n");          

    for(int i=0; i<N-1; i++)               
    {
        printf("Alfa[%i]= %lf+%lf i\n",i, alpha[i].r,alpha[i].i);
    }


    printf("\n\n\n");  
    

    //VECTOR PHI
    
    printf("VECTOR PHI\n\n");          

    for(int i=0; i<N; i++)               
    {
        printf("Phi[%i]= %lf+%lf i\n",i, phi[i].r,phi[i].i);
    }

    printf("\n\n\n");  

    


    //NORMA TOTAL, k0 y s

    
    //printf("Norma total: %lf \n\n",norma);  
    printf("k0: %lf \n\n",k0);  
    printf("s: %lf \n\n",s);

    printf("\n\n\n");
    

    return 0;
}
