//  Project: Error Control Coding - Low-Density Parity-Check Codes	(電通所 Q36114221 蘇沛錦)
//  Written by SU PEI-JIN on June 13, 2023.
//  Consider the (1023, 781) Low-Density Parity-Check (LDCP) code
//  with block length n = 1023, dimension k = 781, row weight ⍴ = 32,
//  and column weight γ = 32.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Define Variables */
int n = 1023;           // n: block length
int k = 781;            // k: dimension
int num_weight = 32;    // row weight ⍴ = 32, and column weight γ = 32      

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

void Parser(int **G_matrix, int **H_info)
{
    int i, j;       // for loop counter
    FILE *fp;
    char *string;
    string = (char *)malloc( _MAX_PATH );   // 配置記憶體區塊

    // Input ldpc_G_1023.txt:
    if((fp=fopen("ldpc_G_1023.txt","r"))==NULL)
    {
        printf("!!! Input File Does Not Exist !!!\n");
        system("pause");
        exit(0);
    }
    else
    {
        for (i = 0; i < k; i++)
        {
            for (j = 0; j < n; j++)
            {
                fscanf(fp, "%1s", string);         // 讀取下一個非空白單一位元組字元
                G_matrix[i][j] = atoi(string);     // 使用 atoi 函式將儲存為字串的數字轉換成數值
                /* 另一種可行的方法 */
                // fscanf(fp, "%1d", &G_matrix[i][j]);  
            }
        }
        fclose(fp);
        free(string);  // 取消配置或釋放記憶體區塊  
    }

    // Input ldpc_H_1023.txt:
    if((fp=fopen("ldpc_H_1023.txt","r"))==NULL)
    {
        printf("!!! Input File Does Not Exist !!!\n");
        system("pause");
        exit(0);
    }
    else
    {
        for (i = 0; i < (2*n); i++)
        {
            for (j = 0; j < num_weight; j++)
            {
                fscanf(fp, "%d", &H_info[i][j]);
            }
        }
        fclose(fp);       
    }
}

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

/* Encoder */
int *Encoder(int *u, int **G_matrix)
{
    int i, j;       // for loop counter
    int *c = (int *)calloc(n, sizeof(int));   // Codeword

    for (j = 0; j < n; j++)
    {
        if (j < k)
        {
            c[j] = u[j];
        }
        else
        {
            for (i = 0; i < k; i++)
            {
                c[j] = c[j]^(u[i]*G_matrix[i][j]);
            }
        }   
    }

    return c;
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

void normal(double* n1, double* n2, double sigma, long* idum)
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

/* Bit-Flipping (BF) Decoder */
void BF_Decoder(int **H_info, double *y, int *u_est, int max_iter)
{
    int i, j, l;    // for loop counter
    int check_sum_signal = 0;
    int delta = 0;      // threshold

    // 配置所需的記憶體空間
    int *x = (int *)calloc(n, sizeof(int));   // Received symbols through the hard limter
    int *s;   // Parity-check sums
    int *f;   // No. unsatisfied parity-check equations

    /* Hard Limter */
    for (i = 0; i < n; i++)
    {
        x[i] = (y[i] >= 0) ? (0):(1);
        if (i < k)
        {
            u_est[i] = x[i];
        }
    }
    
    for (int l = 0; l < max_iter; l++)
    {
        check_sum_signal = 0;
        s = (int *)calloc(n, sizeof(int));   // Parity-check sums
        f = (int *)calloc(n, sizeof(int));   // No. unsatisfied parity-check equations
        for (i = 0; i < (2*n); i++)
        {
            /* Step 1: Compute the parity-check sums (syndrome bits) */
            if (i < n)
            {
                for (j = 0; j < num_weight; j++)
                {
                    s[i] = s[i]^x[H_info[i][j]-1];
                }
                if (s[i] != 0)
                {
                    check_sum_signal = 1;
                }                
            }

            /* Step 2: Find No. unsatisfied parity-check equations for each bit */
            else
            {
                // If S = 0, stop decoding
                if (check_sum_signal == 0)
                {
                    break;
                }
                
                for (j = 0; j < num_weight; j++)
                {
                    f[i - n] = f[i - n] + s[H_info[i][j]-1];
                }
            }
        }

        // If S = 0, stop decoding
        if (check_sum_signal == 0)
        {
            break;
        }

        /* Step 3: Identify the set S of bits for which f_i is the largest */
        delta = f[0];
        for (i = 1; i < n; i++)
        {
            if (f[i] >= delta)
            {
                delta = f[i];
            }
        }

        /* Step 4: Flip the bits in the set S */
        for (i = 0; i < n; i++)
        {
            if (f[i] >= delta)
            {
                x[i] = x[i]^(1);
                if (i < k)
                {
                    u_est[i] = x[i];
                }
            }        
        }
        // Deallocate the memory previously allocated by a call to calloc
        free(s);
        free(f);
    }    
    
    // Deallocate the memory previously allocated by a call to calloc
    free(x);
}

int main(void)
{
    int i, j;    // for loop counter
    double R = (double)k/(double)n;     // Code rate R = 781/1023

    int **G_matrix; // the generator matrix (G)
    int **H_info;   // information about the parity-check matrix (H)

    int *u;         // Point to the memory of Information Bits
    int *c;         // Point to the memory of Codeword
    int *x;         // Point to the memory of Modulated Symbols
    double *z;      // Point to the temporary register of Noise
    double *y;      // Point to the memory of Received Symbols
    int *u_est;     // Estimated information bits

    clock_t start, finish;

    // 配置好所需的記憶體空間
    G_matrix = (int **)malloc(sizeof(int *)*k);
    H_info = (int **)malloc(sizeof(int *)*(2*n));
    for (i = 0; i < (2*n); i++)
    {
        if (i < k)
        {
            G_matrix[i] = (int *)malloc(sizeof(int)*n);
        }        
        H_info[i] = (int *)malloc(sizeof(int)*num_weight);
    }
    // Initializing registers of Modulated Symbols, and Noise, resp.
    x = (int *)calloc(2, sizeof(int));   // Modulated Symbols
    z = (double *)calloc(2, sizeof(double));   // Noise
    // Initializing memory of Received Symbols
    y = (double *)calloc(n, sizeof(double));   // Received Symbols
    // Initializing memory of Estimated information bits
    u_est = (int *)calloc(k, sizeof(int));   // Estimated information bits

    /* 讀取模擬參數 */
    int max_block_errors;   // Number of decoded bits
    double max_SNR_dB;      // Bit signal-to-noise ratio (dB) [SNR(dB) 模擬最大上限]
    long SEED;              // Seed: a negative integer
    int max_iter;           // number of maximum iterations

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

        fscanf(fp, "%d", &max_block_errors);
        printf(" %d", max_block_errors);
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

        fscanf(fp, "%d", &max_iter);
        printf(" %d", max_iter);
        buffer_flag = fscanf(fp, "%[^\n]", buffer);
        if(buffer_flag != 0){
            printf("%s\n", buffer);
        }
        else{
            putchar('\n');
        }

        fclose(fp);
    }

    /* 擷取程式開始執行的時間 */
    start = clock();

    /* Parser */
    Parser(G_matrix, H_info);

    // Prerequisite for AWGN channel
    long *idum;
    idum = (long *)malloc(sizeof(long));
    *idum = SEED;   // SEED must be a negative integer
    double SNR_dB;  // SNR: bit signal-to-noise ratio (dB)
    double SNR;     // Bit signal-to-noise ratio 
    double sigma;   // standard deviation of the noise
    int SNR_dB_Idx;

    /* 配置統計不同 SNR(dB) 條件下的 BER 所需之記憶體空間 */
    int num_tx_bits = 0;        // No. decoded bits
    int num_bits_error = 0;     // No. decoded bit errors
    int num_tx_blocks = 0;      // No. decoded blocks
    int num_blocks_error = 0;   // No. decoded block errors
    int block_error_flag = 0;   // decoded block errors alarm

    for (SNR_dB_Idx = 0; SNR_dB_Idx <= ceil(max_SNR_dB/0.2); SNR_dB_Idx++)
    {
        // Reset Register for Information Bits Generator
        u_reg[0] = 1;
        for (i = 1; i < 6; i++)
        {
            u_reg[i] = 0;            
        }     
        u_reg_reset = 1;

        SNR_dB = SNR_dB_Idx*0.2;
        SNR = pow(10, SNR_dB/10.0);
        sigma = sqrt(1/(2*R*SNR));
        
        num_tx_bits = 0;        // No. decoded bits
        num_bits_error = 0;     // No. decoded bit errors
        num_tx_blocks = 0;      // No. decoded blocks
        num_blocks_error = 0;   // No. decoded block errors
        block_error_flag = 0;   // decoded block errors alarm

        while (num_blocks_error < max_block_errors)
        {
            num_tx_blocks = num_tx_blocks + 1;

            /* Information Bits Generator */
            u = Info_Bits_Gen(k);

            /* Encoder */
            c = Encoder(u, G_matrix);

            for (j = 0; j <= n/2; j++)
            {
                if (j == n/2)
                {
                    /* Modulator */
                    x[0] = Modulator(c[2*j]);

                    /* AWGN channel */
                    normal(&z[0], &z[1], sigma, idum);
                    y[2*j] = x[0] + z[0];
                }
                else
                {
                    /* Modulator */
                    x[0] = Modulator(c[2*j]);
                    x[1] = Modulator(c[2*j + 1]);

                    /* AWGN channel */
                    normal(&z[0], &z[1], sigma, idum);
                    y[2*j] = x[0] + z[0];
                    y[2*j + 1] = x[1] + z[1];
                }            
            }

            /* Bit-Flipping (BF) Decoder */
            BF_Decoder(H_info, y, u_est, max_iter);

            printf("\nSNR(dB): %.1lf, No. decoded blocks: %d\n", SNR_dB, num_tx_blocks);

            /* BER measurement */
            block_error_flag = 0;
            for (j = 0; j < k; j++)
            {
                if (u[j] != u_est[j])
                {
                    num_bits_error = num_bits_error + 1;
                    block_error_flag = 1;
                }        
            }

            num_tx_bits = num_tx_bits + k;

            if (block_error_flag == 1)
            {
                num_blocks_error = num_blocks_error + 1;
                printf("!!! Attention!!!\n");
                printf("No. decoded block errors: %d\nNo. decoded bits: %d, No. decoded bit errors: %d\n", num_blocks_error, num_tx_bits, num_bits_error);
            }       
        }

        if (SNR_dB_Idx == 0)
        {
            fp = fopen("Output_BF.txt","w");
        }
        else
        {
            fp = fopen("Output_BF.txt","a");
        }

        fprintf(fp,"\nSNR(dB): %.1lf\n", SNR_dB);
        fprintf(fp,"No. decoded bits: %d, No. decoded bit errors: %d\n", num_tx_bits, num_bits_error);
        fprintf(fp,"No. decoded blocks: %d, No. decoded block errors: %d\n", num_tx_blocks, num_blocks_error);
        fclose(fp);     

        printf("\n!!! Complete Bit-Flipping (BF) Decoding !!!\n");
        printf("SNR(dB): %.1lf\n", SNR_dB);
        printf("No. decoded bits: %d, No. decoded bit errors: %d\n", num_tx_bits, num_bits_error);
        printf("No. decoded blocks: %d, No. decoded block errors: %d\n", num_tx_blocks, num_blocks_error);

        // Deallocate the memory previously allocated by a call to calloc
        free(u);
        free(c);
    }

    /* 擷取程式結束執行的時間 */
    finish = clock();

    /* 計算程式執行時間 */
    double runtime = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("Runtime = %.3lf seconds\n", runtime);

    fp = fopen("Output_BF.txt","a");
    fprintf(fp,"\nRuntime = %.3lf seconds\n", runtime);
    fclose(fp);

    // Deallocate the memory previously allocated by a call to calloc
    for (i = 0; i < (2*n); i++)
    {
        if (i < k)
        {
            free(G_matrix[i]);
        }
        free(H_info[i]);
    }    
    free(G_matrix);
    free(H_info);
    free(x);
    free(z);
    free(y);
    free(u_est);

    system("pause");
    return 0;
}