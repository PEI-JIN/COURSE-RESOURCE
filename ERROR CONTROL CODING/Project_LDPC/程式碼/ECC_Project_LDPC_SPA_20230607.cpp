//  Project: Error Control Coding - Low-Density Parity-Check Codes	(電通所 Q36114221 蘇沛錦)
//  Written by SU PEI-JIN on June 7, 2023.
//  Consider the (1023, 781) Low-Density Parity-Check (LDCP) code
//  with block length n = 1023, dimension k = 781, row weight ⍴ = 32,
//  and column weight γ = 32.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define min(a,b) (((a) < (b)) ? (a) : (b))  // A preprocessor macro that returns the smaller of two values.

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

/* Sign function (sgn) */
double sgn(double value)
{
    return ((value >= 0.0) ? (1.0) : (-1.0));   // Return an integer indicating the sign of a number
}

/* The function delta() calculates the term △(L1, L2) for Sum-Product Algorithm (SPA) */
double delta(double L1, double L2)
{
    return log((1.0 + exp(-fabs(L1 + L2)))/(1.0 + exp(-fabs(L1 - L2))));
}

/* Check node operation (CHK) for Sum-Product Algorithm (SPA) */
double CHK(double L1, double L2)
{
    return sgn(L1)*sgn(L2)*min(fabs(L1), fabs(L2)) + delta(L1, L2);
}

/* Step 1: Bottom-up (horizontal) */
double Bottom_up(double **msg_mem, int m, int l)
{
    int i;       // for loop counter
    double CHK_temp = 0;    // 暫存函數CHK()計算完的訊息
    double *msg_reg = (double *)calloc((num_weight - 1), sizeof(double));   // 儲存 Variable node 所上傳的 q_m,l' 訊息 (扣除變數 l 以外)

    int count = 0;
    for (i = 0; i < num_weight; i++)
    {
        if (i != l)
        {
            msg_reg[count] = msg_mem[m][i];
            count = count + 1;
        }
    }

    CHK_temp = CHK(msg_reg[0], msg_reg[1]);
    for (i = 2; i < (num_weight - 1); i++)
    {
        CHK_temp = CHK(CHK_temp, msg_reg[i]);
    }
    
    free(msg_reg);  // Deallocate the memory previously allocated by a call to calloc
    return CHK_temp;
}

/* Variable node operation (VAR) */
double VAR(double **msg_mem, int m, int l, int n)
{
    if (n == 1 && l == 0)
    {
        return msg_mem[m][1];
    }
    else if (n == 1 && l == 1)
    {
        return msg_mem[m][0];
    }
    else if (n == 1)
    {
        return msg_mem[m][1] + msg_mem[m][0];
    }
    
    if (n == 2 && l == 3)
    {
        return msg_mem[m][2];
    }
    else if (n == 2  && l == 2)
    {
        return VAR( msg_mem, m, l, 1);
    }
    else if (n == 2)
    {
        return msg_mem[m][2] + VAR( msg_mem, m, l, 1);
    }

    if (n == l && n >= 3)
    {
        return msg_mem[m][n - 1] + VAR( msg_mem, m, l, n - 2);
    }
    else if (n >= 3)
    {
        return msg_mem[m][n] + VAR( msg_mem, m, l, n - 1);
    }
    return 0;
}

/* Decoder */
void Decoder(int **H_info, double *y, int *u_est, double sigma, int max_iter)
{
    int i, j, f;    // for loop counter
    int *msg_mem_count;    // 記憶體的儲存位置 (以儲存計算完的 q_m,l 和 r_m,l 訊息)
    int *x_est = (int *)calloc(n, sizeof(int));     // the estimated codeword
    int parity_check = 0;   // 用來檢查 H*x_est = 0 是否成立
    double Lc = 2.0/(sigma*sigma);

    // 配置好所需的記憶體空間
    double **msg_mem = (double **)calloc(2*n, sizeof(double *));    // 用來儲存 q_m,l 和 r_m,l 訊息的記憶體空間
    for (i = 0; i < (2*n); i++)
    {      
        msg_mem[i] = (double *)calloc(num_weight, sizeof(double));
    }
    double *q = (double *)calloc(n, sizeof(double));      // the log a posteriori probability for each variable node 'l'

    // Initialization
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < num_weight; j++)
        {
            msg_mem[i][j] = Lc*y[(H_info[i][j]) - 1];
        }        
    }

    /* Iterative Decoding */
    for (i = 0; i < max_iter; i++)
    {
        msg_mem_count = (int *)calloc(2*n, sizeof(int));
        
        for (j = 0; j < (2*n); j++)
        {
            if (j < n)
            {
                /* Step 1: Bottom-up (horizontal) */
                for (f = 0; f < num_weight; f++)
                {
                    msg_mem[(H_info[j][f]) - 1 + n][msg_mem_count[(H_info[j][f]) - 1 + n]] = Bottom_up(msg_mem, j, f);
                    msg_mem_count[(H_info[j][f]) - 1 + n] = msg_mem_count[(H_info[j][f]) - 1 + n] + 1;
                }
            }
            else
            {
                for (f = 0; f < num_weight; f++)
                {
                    /* Step 2: Top-down (Vertical) */
                    msg_mem[(H_info[j][f]) - 1][msg_mem_count[(H_info[j][f]) - 1]] = VAR(msg_mem, j, f, num_weight - 1) + (Lc*y[j - n]);
                    msg_mem_count[(H_info[j][f]) - 1] = msg_mem_count[(H_info[j][f]) - 1] + 1;

                    /* Step 3: Termination */
                    q[j - n] = msg_mem[j][f] + VAR(msg_mem, j, f, num_weight - 1) + (Lc*y[j - n]);
                }
            } 
        }

        /* Termination */
        for (j = 0; j < n; j++)
        {
            if (q[j] > 0)
            {
                x_est[j] = 0;
            }
            else
            {
                x_est[j] = 1;
            }

            if (j < k)
            {
                /* Output the estimated information bits */
                u_est[j] = x_est[j];
            }        
        }

        parity_check = 0;
        for (j = 0; j < n; j++)
        {
            parity_check = 0;
            for (f = 0; f < num_weight; f++)
            {
                parity_check = parity_check^x_est[H_info[j][f]-1];
            }
            
            if (parity_check != 0)  // 當 H*x_est = 0 不成立時，迭代解碼演算法 -> [繼續]
            {
                break; 
            }     
        }

        if (parity_check == 0)  // 當 H*x_est = 0 成立時，迭代解碼演算法 -> [停止]
        {
            break;
        }    

        free(msg_mem_count);    // Deallocate the memory previously allocated by a call to calloc        
    }
 
    // Deallocate the memory previously allocated by a call to calloc
    for (i = 0; i < (2*n); i++)
    {
        free(msg_mem[i]);
    }    
    free(msg_mem);
    free(q);
    free(x_est);
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

            /* Decoder */
            Decoder(H_info, y, u_est, sigma, max_iter);

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
            fp = fopen("Output_SPA.txt","w");
        }
        else
        {
            fp = fopen("Output_SPA.txt","a");
        }

        fprintf(fp,"\nSNR(dB): %.1lf\n", SNR_dB);
        fprintf(fp,"No. decoded bits: %d, No. decoded bit errors: %d\n", num_tx_bits, num_bits_error);
        fprintf(fp,"No. decoded blocks: %d, No. decoded block errors: %d\n", num_tx_blocks, num_blocks_error);
        fclose(fp);     

        printf("\n!!! Complete Soft-Decision Decoding (SPA) !!!\n");
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

    fp = fopen("Output_SPA.txt","a");
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