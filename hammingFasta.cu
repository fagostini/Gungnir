#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <vector>

#define MAXLEN 50

/*
Compile:
   nvcc -arch=sm_61 -o hammingFastaGPU hammingFasta.cu

Test: 
   ./hammingFastaGPU query.fa

Test (memcheck)
   nvprof --print-gpu-trace --device-buffer-size on ./hammingFastaGPU test1k.fa
   nvprof --print-gpu-trace --device-buffer-size on ./hammingFastaGPU test10k.fa
   nvprof --print-gpu-trace --device-buffer-size on ./hammingFastaGPU test100k.fa
*/

inline
cudaError_t checkCuda(cudaError_t result){
#if defined(DEBUG) || defined(_DEBUG)
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %sn", cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
#endif
  return result;
}

/* Code from https://www.geeksforgeeks.org/hamming-distance-two-strings */
__global__
void distance_hamming( char *str, size_t p_s, int *dist, size_t p_d, int N){
    int index =  blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    int i, ii;
    if (index < N ){
        for(i=index; i<N; i+=stride){
            int* row_dis = (int*)((char*)dist + i*p_d);
            for(ii=0; ii<N; ii++){
                int count = 0;
                for(int pos=0; pos<p_s; pos++){
                    if( str[i*p_s+pos] == '\0' )
                        break;
                    // printf("%c", str[i*p_s+pos]);
                    if( str[i*p_s+pos] != str[ii*p_s+pos] )
                        count++; 
                }
                // printf(" %d\n", count);
                row_dis[ii] = count;
            }
            // printf("\n");
        }
    }
    __syncthreads();
}

/* Code from https://www.geeksforgeeks.org/edit-distance-dp-5/ */
#define min(x,y) ((x) < (y) ? (x) : (y)) //calculate minimum between two values

int distance_levenshtein(std::string str1, std::string str2){ 
    int len1 = str1.length(); 
    int len2 = str2.length(); 
  
    // Create a DP array to memoize result of previous computations 
    int DP[2][len1 + 1]; 
  
    // To fill the DP array with 0 
    memset(DP, 0, sizeof DP); 
  
    // Base condition when second string is empty then we remove all characters 
    for (int i = 0; i <= len1; i++) 
        DP[0][i] = i; 
  
    // Start filling the DP 
    // This loop run for every character in second string
    for (int i = 1; i <= len2; i++) { 
        // This loop compares the char from second string with first string
        // characters 
        for (int j = 0; j <= len1; j++) { 
            // if first string is empty then we have to perform add character
            // operation to get second string 
            if (j == 0) 
                DP[i % 2][j] = i; 
  
            // if character from both string is same then we do not perform any
            // operation . here i % 2 is for bound the row number
            else if (str1[j - 1] == str2[i - 1]) { 
                DP[i % 2][j] = DP[(i - 1) % 2][j - 1]; 
            } 
  
            // if character from both string is not same then we take the minimum 
            // from three specified operation 
            else { 
                DP[i % 2][j] = 1 + min(DP[(i - 1) % 2][j], 
                                       min(DP[i % 2][j - 1], 
                                           DP[(i - 1) % 2][j - 1])); 
            } 
        } 
    } 
  
    // after complete fill the DP array if the len2 is even then we end
    // up in the 0th row else we end up in the 1th row so we take len2 % 2
    // to get row cout << DP[len2 % 2][len1] << endl;
    return DP[len2 % 2][len1];
}

/* Code from http://rosettacode.org/wiki/FASTA_format#C.2B.2B */
int main( int argc, char **argv ){
    if( argc <= 1 ){
        std::cerr << "Usage: "<<argv[0]<<" [infile]" << std::endl;
        return -1;
    }
 
    std::ifstream input(argv[1]);
    if( !input.good() ){
        std::cerr << "Error opening '" << argv[1] << "'. Bailing out." << std::endl;
        return -1;
    }

    // Getting the total number of entries from the input file
    int N = 0;
    std::string tmp_line;
    printf("Counting entries in the input file... ");	
    while( std::getline( input, tmp_line ).good() ){
        if( !tmp_line.empty() && tmp_line[0] == '>' ){
            N++;
        }
    }
    tmp_line.clear();
    printf("OK\n");	

    // Declaring Vector of String type 
    printf("Allocating 2D arrays on host... ");
    int *dis;
    char *ids, *seq;
    ids = (char*) malloc(sizeof(char)*(MAXLEN)*N);
    seq = (char*) malloc(sizeof(char)*(MAXLEN)*N);
    dis = (int*) calloc(sizeof(int), N*N);
    printf("OK\n");	

    int i = 0;
    char c;
    input.clear();
    input.seekg(0, std::ios::beg);
    input.get(c);
    printf("Reading input file and filling the arrays... ");	
    do{
        // std::cout << c << std::endl;
        if( c == '>' ){
            int p = 0;
            char name[(MAXLEN)];
            do{
                name[p++] = c;
                input.get(c);
            } while( c != '\n' || p > MAXLEN );
            name[p] = '\0';
            memcpy( &ids[i*MAXLEN], name, sizeof(char)*MAXLEN );
        } else {
            int p = 0;
            char content[(MAXLEN)];
            do{
                content[p++] = c;
                input.get(c);
            } while( c != '\n' || p > MAXLEN );
            content[p] = '\0';
            memcpy( &seq[i*MAXLEN], content, sizeof(char)*MAXLEN ); 
            i++;
        }
        input.get(c);
    } while( !input.eof() );
    input.close();
    printf("OK\n");	

    // printf("Array initialisation check:\n   ");
    // for (i=0; i<N; i++){
    //     for(int ii=0; ii<MAXLEN; ii++){
    //         if( seq[i*MAXLEN+ii] == '\0')
    //             break;
    //         std::cout << seq[i*MAXLEN+ii];
    //     }
    // printf("\n   ");
    // }
    // printf("\n");

    int *GPU_dis;
    char *GPU_seq;
    size_t pitch_seq, pitch_dis;

    printf("Allocating 2D arrays on device... ");	
    checkCuda( cudaMallocPitch((void**)&GPU_seq, &pitch_seq, sizeof(char)*(MAXLEN), N) );
    checkCuda( cudaMallocPitch((void**)&GPU_dis, &pitch_dis, sizeof(int)*N, N) );
    printf("OK\n");

    printf("Initialising 2D arrays to device... ");	
    checkCuda( cudaMemcpy2D(GPU_seq, pitch_seq, seq, sizeof(char)*(MAXLEN), sizeof(char)*(MAXLEN), N, cudaMemcpyHostToDevice) );
    checkCuda( cudaMemcpy2D(GPU_dis, pitch_dis, dis, sizeof(int)*N, sizeof(int)*N, N, cudaMemcpyHostToDevice) );
    printf("OK\n");

    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;

    printf("Running kernel operations... ");	
    // multi<<<1, 4>>>(GPU_seq, pitch_seq, GPU_dis, pitch_dis, N);
    distance_hamming<<<numBlocks, blockSize>>>(GPU_seq, pitch_seq, GPU_dis, pitch_dis, N);
    printf("OK\n");

    // Wait for GPU to finish before accessing on host
    printf("Synchronising host and device... ");	
    cudaDeviceSynchronize();
    printf("OK\n");

    printf("Copying 2D array from device to host... ");	
    cudaMemcpy2D(dis, sizeof(int)*N, GPU_dis, pitch_dis, sizeof(int)*N, N, cudaMemcpyDeviceToHost);
    printf("OK\n");

    // printf("Results check:\n   ");
    // for(i = 0; i < N*N; i++){
    //     printf("%d   ", dis[i]);
    //     if(i%N == N-1)
    //     printf("\n   ");
    // }
    // printf("\n");

    // release dynamically allocated memory
    printf("Releasing allocated memory:\n");	
    free(ids);
    free(seq);
    free(dis);
    printf("Host   --> OK\n");
    cudaFree(GPU_seq);
    cudaFree(GPU_dis);
    printf("Device --> OK\n");

    return 0;
}