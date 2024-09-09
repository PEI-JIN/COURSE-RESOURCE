//  Project: Error Control Coding - uncoded BPSK	(電通所 Q36114221 蘇沛錦)
//  Written by SU PEI-JIN on June 13, 2023.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>    

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

int u_reg[6];       // Register for Information Bits Generator
int u_reg_reset;    // Reset signal for Information Bits Generator

/* Information Bits Generator */
int *Info_Bits_Gen(int L)
{
    int i;       // for loop counter
    int *u = (int *)calloc(L, sizeof(int));   // Information bits

    for (i = 0; i < L; i++)
    {
        if (i < 6 && u_reg_reset == 1)
        {
            u[i] = u_reg[i];
            if (i == 5)
            {
                u_reg_reset = 0;
            }
        }
        else
        {
            u[i] = u_reg[0]^u_reg[1];
            u_reg[0] = u_reg[1];
            u_reg[1] = u_reg[2];
            u_reg[2] = u_reg[3];
            u_reg[3] = u_reg[4];
            u_reg[4] = u_reg[5];
            u_reg[5] = u[i];
        }
    }

    return u;
}

/* Modulator */
int Modulator(int c)
{
    int x;   // Modulated Symbols

    if (c == 0)
    {
        x = 1;
    }
    else if (c == 1)
    {
        x = -1;
    }

    return x;
}

/* AWGN channel */
double ran1(long *idum)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    double temp;
    if (*idum <= 0 || !iy){
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        for (j=NTAB+7;j>=0;j--){
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

void normal(double *n1, double* n2, double sigma, long *idum)
{
    double x1, x2, s;

    do{
        x1 = ran1(idum);
        x2 = ran1(idum);
        x1 = 2*x1 - 1;
        x2 = 2*x2 - 1;
        s = x1*x1 + x2*x2;
    } while (s >= 1.0);
    *n1 = sigma*x1*sqrt(-2*log(s)/s);
    *n2 = sigma*x2*sqrt(-2*log(s)/s);
}

/* Hard Limter */
int Hard_Limter(double y)
{
    return (y >= 0) ? 0:1;
}

int main(void)
{
    int i;    // for loop counter

    /* 讀取模擬參數 */
    int N;                  // Number of information bits
    double max_SNR_dB;      // Bit signal-to-noise ratio (dB) [SNR(dB) 模擬最大上限]
    long SEED;              // Seed: a negative integer

    FILE *fp;
    char filename[30] = "Sim.txt";
    char buffer[50];
    int buffer_flag;

    // Input Sim.txt:
    if((fp=fopen(filename,"r"))==NULL){
        printf("!!! Input File Does Not Exist !!!\n");
        system("pause");
        exit(0);
    }
    else{
        printf("Input %s\n",filename);

        fscanf(fp, "%d", &N);
        printf(" %d", N);
        buffer_flag = fscanf(fp, "%[^\n]", buffer);
        if(buffer_flag != 0){
            printf("%s\n", buffer);
        }
        else{
            putchar('\n');
        }

        fscanf(fp, "%lf", &max_SNR_dB);
        printf(" %.1lf", max_SNR_dB);
        buffer_flag = fscanf(fp, "%[^\n]", buffer);
        if(buffer_flag != 0){
            printf("%s\n", buffer);
        }
        else{
            putchar('\n');
        }

        fscanf(fp, "%ld", &SEED);
        printf(" %ld", SEED);
        buffer_flag = fscanf(fp, "%[^\n]", buffer);
        if(buffer_flag != 0){
            printf("%s\n", buffer);
        }
        else{
            putchar('\n');
        }

        fclose(fp);
    }

    int *u;         // Point to the memory of Information Bits
    int *x;         // Point to the memory of Modulated Symbols
    double *z;      // Point to the temporary register of Noise
    double *y;      // Point to the temporary register of Received Symbols
    int *u_est;     // Point to the memory of Estimated Information Bits

    clock_t start, finish;

    // 配置好所需的記憶體空間
    x = (int *)calloc(2, sizeof(int));          // Modulated Symbols
    z = (double *)calloc(2, sizeof(double));    // Noise
    y = (double *)calloc(2, sizeof(double));    // Noise
    u_est = (int *)calloc(N, sizeof(int));      // Estimated Information Bits

    /* 擷取程式開始執行的時間 */
    start = clock();

    // Prerequisite for AWGN channel
    long *idum;
    idum = (long *)malloc(sizeof(long));
    *idum = SEED;   // SEED must be a negative integer
    double SNR_dB;  // SNR: bit signal-to-noise ratio (dB)
    double SNR;     // Bit signal-to-noise ratio 
    double sigma;   // standard deviation of the noise
    int SNR_dB_Idx;

    int num_bits_error = 0;     // No. rx bit errors

    // Reset Register for Information Bits Generator
    u_reg[0] = 1;
    for (i = 1; i < 6; i++)
    {
        u_reg[i] = 0;            
    }     
    u_reg_reset = 1;

    // Information Bits Generator
    u = Info_Bits_Gen(N);

    for (SNR_dB_Idx = 0; SNR_dB_Idx <= ceil(max_SNR_dB/0.5); SNR_dB_Idx++)
    {
        SNR_dB = SNR_dB_Idx*0.5;
        SNR = pow(10, SNR_dB/10.0);
        sigma = sqrt(1/(2*SNR));
        
        num_bits_error = 0;     // No. received bit errors

        if (N%2 == 0)
        {
            for (i = 0; i < N/2; i++)
            {
                // Modulator
                x[0] = Modulator(u[2*i]);
                x[1] = Modulator(u[2*i + 1]);

                // AWGN channel
                normal(&z[0], &z[1], sigma, idum);
                y[0] = x[0] + z[0];
                y[1] = x[1] + z[1];

                // Hard Limter
                u_est[2*i] = Hard_Limter(y[0]);
                u_est[2*i + 1] = Hard_Limter(y[1]);
            }
        }
        else
        {
            for (i = 0; i < N/2; i++)
            {
                // Modulator
                x[0] = Modulator(u[2*i]);
                x[1] = Modulator(u[2*i + 1]);

                // AWGN channel
                normal(&z[0], &z[1], sigma, idum);
                y[0] = x[0] + z[0];
                y[1] = x[1] + z[1];

                // Hard Limter
                u_est[2*i] = Hard_Limter(y[0]);
                u_est[2*i + 1] = Hard_Limter(y[1]);
            }

            // Modulator
            x[0] = Modulator(u[N-1]);
            // AWGN channel
            normal(&z[0], &z[1], sigma, idum);
            y[0] = x[0] + z[0];

            // Hard Limter
            u_est[N-1] = Hard_Limter(y[0]);
        }

        // BER Measurement
        for (i = 0; i < N; i++)
        {
            if (u_est[i] != u[i])
            {
                num_bits_error = num_bits_error + 1;
            }            
        }

        // Printf
        if (SNR_dB_Idx == 0)
        {
            fp = fopen("Output.txt","w");
        }
        else
        {
            fp = fopen("Output.txt","a");
        }

        fprintf(fp,"\nSNR(dB): %.1lf\n", SNR_dB);
        fprintf(fp,"No. tx bits: %d, No. rx bit errors: %d\n", N, num_bits_error);
        fclose(fp);     

        printf("\n!!! Complete Simulation !!!\n");
        printf("SNR(dB): %.1lf\n", SNR_dB);
        printf("No. tx bits: %d, No. rx bit errors: %d\n", N, num_bits_error);       
    }

    /* 擷取程式結束執行的時間 */
    finish = clock();

    /* 計算程式執行時間 */
    double runtime = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("Runtime = %.3lf seconds\n", runtime);

    fp = fopen("Output.txt","a");
    fprintf(fp,"\nRuntime = %.3lf seconds\n", runtime);
    fclose(fp);

    // Deallocate the memory previously allocated by a call to calloc
    free(u);
    free(x);
    free(z);
    free(y);
    free(u_est);

    system("pause");
    return 0;
}