// Include C++ header files
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
// CUDA libraries.
#include <cuda.h>
#include <cuda_runtime.h>

#define BLOCKSIZE 32

#define MAXREF 3000000000

#define MAXLEN 40
#define MISMATCH 9

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

__global__
void distance_hamming( char *str, size_t p_s, char *ref, size_t len, int *dist, size_t p_d, int N){
    int Xindex =  blockIdx.x * blockDim.x + threadIdx.x;
    int Xstride = blockDim.x * gridDim.x;
    int Yindex =  blockIdx.y * blockDim.y + threadIdx.y;
    int Ystride = blockDim.y * gridDim.y;
    int i, ii;
    // printf("\n%d %d", index, stride);
    if(Xindex < N && Xindex < (len-MAXLEN+1) ){
        for(i=Xindex; i<N; i+=Xstride){
            // printf("\n%d %.40s %s", i, &str[i*p_s], &ref[0]);
            int* row_dis = (int*)((char*)dist + i*p_d);
            for(ii=Yindex; ii<(len-MAXLEN+1); ii+=Ystride){
                // printf("\n%d %.40s %.40s", i, &str[i*p_s], &ref[ii]);
                int pos, count = 0;
                for(pos=0; pos<p_s; pos++){
                    if( str[i*p_s+pos] == '\0' )
                        break;
                    // printf("%c", str[i*p_s+pos]);
                    count += ( str[i*p_s+pos] != ref[ii+pos] );
                }
                // printf("\n%d %.40s %.40s %d", i, &str[i*p_s], &ref[ii], count);
                // printf("   %d", count);
                if( count <= MISMATCH ){
                    row_dis[count] += 1;
                    // printf("%d %d %d\n", i, count, row_dis[count]);
                }
            }
        // printf("\n");
        }
    }
    // __syncthreads();
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
        std::cerr << "Usage: " << argv[0] << " <QUERY_FASTA> [<REFERENCE_FASTA>]" << std::endl;
        return -1;
    }
 
    std::ifstream input(argv[1]);
    if( !input.good() ){
        std::cerr << "Error opening '" << argv[1] << "'. Bailing out." << std::endl;
        return -1;
    }

    std::ifstream ref(argv[2]);
    if( !ref.good() ){
        std::cerr << "Error opening '" << argv[2] << "'. Bailing out." << std::endl;
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
    ids = (char*) malloc(sizeof(char)*(MAXLEN+1)*N);
    seq = (char*) malloc(sizeof(char)*(MAXLEN+1)*N);
    dis = (int*) calloc(sizeof(int), (MISMATCH+1)*N);
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
            char name[(MAXLEN+1)];
            do{
                name[p++] = c;
                input.get(c);
            } while( c != '\n' || p > (MAXLEN+1) );
            name[p] = '\0';
            memcpy( &ids[i*(MAXLEN+1)], name, sizeof(char)*(MAXLEN+1) );
        } else {
            int p = 0;
            char content[(MAXLEN+1)];
            do{
                content[p++] = c;
                input.get(c);
            } while( c != '\n' || p > (MAXLEN+1) );
            content[p] = '\0';
            memcpy( &seq[i*(MAXLEN+1)], content, sizeof(char)*(MAXLEN+1) ); 
            i++;
        }
        input.get(c);
    } while( !input.eof() );
    input.close();
    printf("OK\n");	

    // printf("Array initialisation check:\n   ");
    // for (i=0; i<N; i++){
    //     for(int ii=0; ii<(MAXLEN+1); ii++){
    //         if( seq[i*(MAXLEN+1)+ii] == '\0')
    //             break;
    //         std::cout << seq[i*(MAXLEN+1)+ii];
    //     }
    // printf("\n   ");
    // }
    // printf("\n");

    int *GPU_dis;
    char *GPU_seq;
    size_t pitch_seq, pitch_dis;

    printf("Allocating 2D arrays on device... ");	
    checkCuda( cudaMallocPitch((void**)&GPU_seq, &pitch_seq, sizeof(char)*(MAXLEN+1), N) );
    checkCuda( cudaMallocPitch((void**)&GPU_dis, &pitch_dis, sizeof(int)*(MISMATCH+1), N) );
    printf("OK\n");

    printf("Initialising 2D arrays to device... ");	
    checkCuda( cudaMemcpy2D(GPU_seq, pitch_seq, seq, sizeof(char)*(MAXLEN+1), sizeof(char)*(MAXLEN+1), N, cudaMemcpyHostToDevice) );
    checkCuda( cudaMemcpy2D(GPU_dis, pitch_dis, dis, sizeof(int)*(MISMATCH+1), sizeof(int)*(MISMATCH+1), N, cudaMemcpyHostToDevice) );
    printf("OK\n");

    printf("Allocating reference array on host... ");
    char *reference;
    reference = (char*) malloc(sizeof(char)*(MAXREF+1));
    printf("OK\n");

    printf("Reading reference the input file... ");
    int ref_index = 0, ref_length = 0, n = 0;
    while( ref.get(c) && !ref.eof() ){
        if( c == '>' ){
            while( c != '\n' ){
                ref.get(c);
            }
        } else {
            do{
                if( c == 'N'){
                    n += 1;
                } else {
                    n = 0;
                }
                if( n < (MAXLEN+1) ){
                    reference[ref_length++] = toupper(c);
                }
                ref.get(c);
            } while( c != '\n' || ref_length > (MAXREF+1) );
        }
    }
    reference[ref_length] = '\0';
    ref.close();
    printf("OK\n");

    // printf("%s\n", reference);
    printf("Effective refence length: %d\n", ref_length);

    char *GPU_ref;
    size_t pitch_ref = sizeof(char)*ref_length;

    printf("Allocating reference array on device... ");
    // checkCuda( cudaMallocManaged((void**)&GPU_ref, sizeof(char)*ref_length) );
    checkCuda( cudaMallocPitch((void**)&GPU_ref, &pitch_ref, sizeof(char)*ref_length, 1) );
    printf("OK\n");

    printf("Initialising reference array on device... ");
    checkCuda( cudaMemcpy(GPU_ref, reference, sizeof(char)*ref_length, cudaMemcpyHostToDevice) );
    printf("OK\n");

    // int blockSize = N < BLOCKSIZE ? N : BLOCKSIZE;
    // int numBlocks = (N + blockSize - 1) / blockSize;

    dim3 block(BLOCKSIZE, BLOCKSIZE);
    dim3 grid ( (N + BLOCKSIZE - 1 ) / BLOCKSIZE, (N + BLOCKSIZE - 1) / BLOCKSIZE );

    printf("Running kernel operations... ");
    // distance_hamming<<<numBlocks, blockSize>>>(GPU_seq, pitch_seq, GPU_ref, pitch_ref, GPU_dis, pitch_dis, N);
    distance_hamming<<<grid, block>>>(GPU_seq, pitch_seq, GPU_ref, pitch_ref, GPU_dis, pitch_dis, N);
    printf("OK\n");

    // Wait for GPU to finish before accessing on host
    printf("Synchronising host and device... ");
    cudaDeviceSynchronize();
    printf("OK\n");

    printf("Copying 2D array from device to host... ");
    cudaMemcpy2D(dis, sizeof(int)*(MISMATCH+1), GPU_dis, pitch_dis, sizeof(int)*(MISMATCH+1), N, cudaMemcpyDeviceToHost);
    printf("OK\n");

    printf("Writing results...");
    std::ofstream output("results.tsv");
    for(i = 0; i < N; i++){
        for(int ii=0; ii<(MAXLEN+1); ii++){
            if( ids[i*(MAXLEN+1)+ii] == '\0'){
                output << ' ';
                break;
            }
            output << ids[i*(MAXLEN+1)+ii];
        }
        for(int ii=0; ii<(MAXLEN+1); ii++){
            if( seq[i*(MAXLEN+1)+ii] == '\0')
                break;
            output << seq[i*(MAXLEN+1)+ii];
        }
        for(int ii=0; ii<(MISMATCH+1); ii++){
            output << ' ' << dis[i*(MISMATCH+1)+ii];
        }
        output << '\n';
    }
    output.close();
    printf("OK\n");

    // release dynamically allocated memory
    printf("Releasing allocated memory:\n");	
    free(ids);
    free(seq);
    free(dis);
    free(reference);
    printf("Host   --> OK\n");
    cudaFree(GPU_seq);
    cudaFree(GPU_dis);
    cudaFree(GPU_ref);
    printf("Device --> OK\n");

    return 0;
}