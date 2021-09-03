#include <stdio.h>
#include <math.h>
#include "omp.h"
#define M_PI 3.14159265358979323846


   FILE *fptr5;
   FILE *fptr6;
   FILE *fptr7;
   FILE *fptr8;

#define NUM_THREADS 8

void funcE();
void funcF();
void funcG();
void funcH();

void main()
{
    omp_set_num_threads(NUM_THREADS);
    const double startTime = omp_get_wtime();

    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            (void)funcE();
            #pragma omp section
            (void)funcF();
            #pragma omp section
            (void)funcG();
            #pragma omp section
            (void)funcH();
        }
    }
    const double endTime = omp_get_wtime();
    printf("tomo (%lf) segundos\n",(endTime - startTime));

}

void funcE(){
    long N = 50000;
    fptr5=fopen("RungeInt1.txt","w");
    //printf("Numero de pasos:%d Atendido por thread:%d\n", N,omp_get_thread_num());
    fprintf(fptr5, "Datos que encuentra el metodo de Euler(variable ind.\t variable dep.\t numero de thread)\n");
    double h,t,w;
    double k1= 0.0, k2= 0.0, k3=0.0, k4=0.0;
    double w0=M_PI/4, t0=0,a=0.0,b=M_PI;
    double ab=0.0;
    int i;
    w=w0;
    fprintf(fptr5, "%f\t %f\n", a, w);
    for(i=0;i<N;i++){
        h=(b-a)/N;
        t=a+(h*i);
        ab=t*t;
        k1=h*(t*exp(3*t)-(2*w));
        k2=h*(((t+(0.5*h))*exp(3*((t+(0.5*h))))-(2*(w+(0.5*k1)))));  
        k3=h*(((t+(0.5*h))*exp(3*((t+(0.5*h))))-(2*(w+(0.5*k2)))));  
        k4=h*(((t+h)*exp(3*((t+h)))-(2*(w+k3))));
        w=w+((0.1666666667)*(k1+(2*k2)+(2*k3)+k4));
        fprintf(fptr5, "%f\t %f \t numero de thread:%d\n", t+h, w,omp_get_thread_num());
        }

   fclose(fptr5);
}

void funcF(){
    long N = 50000;
    fptr6=fopen("RungeInt2.txt","w");
    //printf("Numero de pasos:%d Atendido por thread:%d\n", N,omp_get_thread_num());
    fprintf(fptr6, "Datos que encuentra el metodo de Euler(variable ind.\t variable dep.\t numero de thread)\n");
    double h,t,w;
    double k1= 0.0, k2= 0.0, k3=0.0, k4=0.0;
    double w0=M_PI/4, t0=0,a=0.0,b=M_PI;
    double ab=0.0;
    int i;
    w=w0;
    fprintf(fptr6, "%f\t %f\n", a, w);
    for(i=0;i<N;i++){
        h=(b-a)/N;
        t=a+(h*i);
        ab=t*t;
        k1=h*(pow((t-w),2)+1.0);
        k2=h*(pow(((t+(0.5*h))-(w+(0.5*k1))),2)+1.0); 
        k3=h*(pow(((t+(0.5*h))-(w+(0.5*k2))),2)+1.0); 
        k4=h*(pow(((t+h)-(w+k3)),2)+1.0);
        w=w+((0.1666666667)*(k1+(2*k2)+(2*k3)+k4));
        fprintf(fptr6, "%f\t %f \t numero de thread:%d\n", t+h, w,omp_get_thread_num());
        }

   fclose(fptr6);
}

void funcG(){
    long N = 50000;
    fptr7=fopen("RungeInt3.txt","w");
    //printf("Numero de pasos:%d Atendido por thread:%d\n", N,omp_get_thread_num());
    fprintf(fptr7, "Datos que encuentra el metodo de Euler(variable ind.\t variable dep.\t numero de thread)\n");
    double h,t,w;
    double k1= 0.0, k2= 0.0, k3=0.0, k4=0.0;
    double w0=M_PI/4, t0=0,a=0.0,b=M_PI;
    double ab=0.0;
    int i;
    w=w0;
    fprintf(fptr7, "%f\t %f\n", a, w);
    for(i=0;i<N;i++){
        h=(b-a)/N;
        t=a+(h*i);
        ab=t*t;
        k1=h*(1.0 + (w/t));
        k2=h*(1.0 + ((w+(0.5*k1))/(t+(0.5*h)))); 
        k3=h*(1.0 + ((w+(0.5*k2))/(t+(0.5*h))));
        k4=h*h*(1.0 + ((w+k3)/(t+h)));
        w=w+((0.1666666667)*(k1+(2*k2)+(2*k3)+k4));
        fprintf(fptr7, "%f\t %f \t numero de thread:%d\n", t+h, w,omp_get_thread_num());
        }

   fclose(fptr7);
}

void funcH(){
    long N = 50000;
    fptr8=fopen("RungeInt4.txt","w");
    //printf("Numero de pasos:%d Atendido por thread:%d\n", N,omp_get_thread_num());
    fprintf(fptr8, "Datos que encuentra el metodo de Euler(variable ind.\t variable dep.\t numero de thread)\n");
    double h,t,w;
    double k1= 0.0, k2= 0.0, k3=0.0, k4=0.0;
    double w0=M_PI/4, t0=0,a=0.0,b=M_PI;
    double ab=0.0;
    int i;
    w=w0;
    fprintf(fptr8, "%f\t %f\n", a, w);
    for(i=0;i<N;i++){
        h=(b-a)/N;
        t=a+(h*i);
        ab=t*t;
        k1=h*(cos(2*t*w)+sin(3*t*w));
        k2=h*(cos(2*(t+(0.5*h))*(w+(0.5*k1)))+sin(3*(t+(0.5*h))*(w+(0.5*k1)))); 
        k3=h*(cos(2*(t+(0.5*h))*(w+(0.5*k2)))+sin(3*(t+(0.5*h))*(w+(0.5*k2))));
        k4=h*(cos(2*(t+h)*(w+k3))+sin(3*(t+h)*(w+k3))); 
        w=w+((0.1666666667)*(k1+(2*k2)+(2*k3)+k4));
        fprintf(fptr8, "%f\t %f \t numero de thread:%d\n", t+h, w,omp_get_thread_num());
        }

   fclose(fptr8);
}
