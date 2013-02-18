#undef _GLIBCXX_ATOMIC_BUILTINS
#undef _GLIBCXX_USE_INT128

#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/system_error.h>



#ifdef	__cplusplus
extern "C" {
#endif

#include <cuda.h>   
#include "struct.h"

#define TILE_WIDTH 16
#define PITCH 64
#define MEM_SIZE_MAX 134217728L
#define CLOCK_PRECISION  1E9
    
struct greater_rmsd{
    
    __host__ __device__
    bool operator()(Triple x, Triple y) 
    {
        return x.rmsd < y.rmsd;
    }
            
};
    
    
__device__ int  sum_to_zero_gpu (float *a, float *b ) {
    int i;
    float sum,aux;

    sum = 0;
    for (i=0; i<3; i++ ) {
	aux = a[i] + b[i];
	sum += aux*aux;
    }
    return (sum < 0.001);
}


__device__ int  normalized_cross_gpu (float *x, float *y, float * v, float *norm_ptr) {

    /* v is the output */
    float norm = 0;
    float vec[3];
    int i;

    if (sum_to_zero_gpu (x, y) ) return 1;
    
    vec[0] = x[1]*y[2] -  x[2]*y[1];
    norm += vec[0]*vec[0];
    vec[1] = x[2]*y[0] -  x[0]*y[2]; 
    norm += vec[1]*vec[1];
    vec[2] = x[0]*y[1] -  x[1]*y[0];
    norm += vec[2]*vec[2];
    norm = sqrt (norm);

    for (i=0; i<3; i++ ) {
	v[i] = vec[i] /norm;
    }

    if ( norm_ptr) *norm_ptr = norm;
    

    return 0;
    
}

__device__ int unnorm_dot_gpu (float *x, float *y, float * dot) {


    float cosine = 0;
    int i;
    
    for (i=0; i<3; i++ ) {
	cosine += x[i]*y[i];
    }

    if (cosine > 1.0 ) 
	cosine = 1.0; /* this should be numerical */
   
    *dot = cosine;
   
    return 0;
    
}


__device__ int distance_of_nearest_approach_gpu(float *x_d, float *x_cm_d, float *y_d, float *y_cm_d, int *set_of_directions_x,  int *set_of_directions_y,
        int set_size, float * rmsd_ptr) {

    float cm_vector[3], distance_x, distance_y;
    float aux, rmsd;
    int i, ray_a, ray_b, a, b, prev_next, norm;
    float cross[3];

    //if (set_size <= 1) return 1;
       
    float * row_data_a, * row_data_b;

    
    rmsd = 0.0;
    norm = 0;
    /* the rmsd for the remaining vectors is ... */
    for (a = 0; a < set_size; a++) {

        for (prev_next = -1; prev_next <= +1; prev_next += 2) {
            b = (set_size + a + prev_next) % set_size;

            ray_a = set_of_directions_x[a];
            ray_b = set_of_directions_x[b];
            
            row_data_a = (float*)(((char*)x_cm_d) + (ray_a * PITCH));
            row_data_b = (float*)(((char*)x_cm_d) + (ray_b * PITCH));
    
            
            /* distance of nearest approach of ray b
               to the cm of a, in the set of directions x */
            
            for (i = 0; i < 3; ++i) {
                cm_vector[i] = row_data_b[i] - row_data_a[i];
            }
            
            
            row_data_a = (float*)(((char*)x_d) + (ray_a * PITCH));
            normalized_cross_gpu(cm_vector, row_data_a, cross, &distance_x);

            ray_a = set_of_directions_y[a];
            ray_b = set_of_directions_y[b];
            /* distance of nearest approach of ray b
               to the cm of a, in the set of directions y */
            
            row_data_a = (float*)(((char*)y_cm_d) + (ray_a * PITCH));
            row_data_b = (float*)(((char*)y_cm_d) + (ray_b * PITCH));
    
            for (i = 0; i < 3; ++i) {
                cm_vector[i] = row_data_b[i] - row_data_a[i];
            }
            
            row_data_a = (float*)(((char*)y_d) + (ray_a * PITCH));
            normalized_cross_gpu(cm_vector, row_data_a, cross, &distance_y);

            aux = distance_x - distance_y;
            rmsd += aux*aux;
            norm++;
        }

    }

    rmsd /= norm;
    rmsd = sqrt(rmsd);

    *rmsd_ptr = rmsd;

    return 0;
}

__device__ int  same_hand_triple_gpu(float *x_d, float *x_cm_d, float *y_d, float *y_cm_d, int *set_of_directions_x,
        int *set_of_directions_y, int set_size) {

    float cm_vector[3], avg_cm[3], cross[3], dx, dy;
    int i, ray_a, ray_b, ray_c;

    // if (set_size != 3) return 0;

    
    /*****************************/
    /*****************************/
    ray_a = set_of_directions_x[0];
    ray_b = set_of_directions_x[1];
    ray_c = set_of_directions_x[2];
    /* I better not use the cross prod here: a small diff
       int he angle makeing them pointing toward each other or
       away from each other changes the direction of cross prod;
       rather, use one vector, and distance between the cm's as the other */
    float * row_data_a, * row_data_b, * row_data_c;

    row_data_a = (float*)(((char*)x_d) + (ray_a * PITCH));
    row_data_b = (float*)(((char*)x_d) + (ray_b * PITCH));
   
    
    normalized_cross_gpu(row_data_a, row_data_b, cross, NULL);
    
    row_data_a = (float*)(((char*)x_cm_d) + (ray_a * PITCH));
    row_data_b = (float*)(((char*)x_cm_d) + (ray_b * PITCH));
    row_data_c = (float*)(((char*)x_cm_d) + (ray_c * PITCH)); 
    
    /* note I am making another cm vector here */
    for (i = 0; i < 3; i++) {
        avg_cm[i] = (row_data_b[i] - row_data_a[i])/2;
        cm_vector[i] = row_data_c[i] - avg_cm[i];
        //avg_cm[i] = (x_cm_d[ray_b][i] + x_cm_d[ray_a][i]) / 2;
        //cm_vector[i] = x_cm_d[ray_c][i] - avg_cm[i];
    }
    unnorm_dot_gpu(cm_vector, cross, &dx);

    /*****************************/
    /*****************************/
    ray_a = set_of_directions_y[0];
    ray_b = set_of_directions_y[1];
    ray_c = set_of_directions_y[2];
    
    row_data_a = (float*)(((char*)y_d) + (ray_a * PITCH));
    row_data_b = (float*)(((char*)y_d) + (ray_b * PITCH));
   
    normalized_cross_gpu(row_data_a, row_data_b, cross, NULL);
    
    row_data_a = (float*)(((char*)y_cm_d) + (ray_a * PITCH));
    row_data_b = (float*)(((char*)y_cm_d) + (ray_b * PITCH));
    row_data_c = (float*)(((char*)y_cm_d) + (ray_c * PITCH));    
    
    /* note I am making another cm vector here */
    for (i = 0; i < 3; i++) {
        avg_cm[i] = (row_data_b[i] - row_data_a[i])/2;
        cm_vector[i] = row_data_c[i]- avg_cm[i];
        
        //avg_cm[i] = (y_cm_d[ray_b][i] + y_cm_d[ray_a][i]) / 2;
        //cm_vector[i] = y_cm_d[ray_c][i] - avg_cm[i];
    }
    /*note: unnorm_dot thinks it is getting unit vectors,
      and evrything that is >1 will be "rounded" to 1
      (similarly for -1) - it doesn't do the normalization itself*/
    
    unnorm_dot_gpu(cm_vector, cross, &dy);

    if (dx * dy < 0) return 0;


    return 1; /* this isn't err value - the handedness is the same */
}


__global__ void find_triplets_gpu(float *x_d, float *x_cm_d, float *y_d, float *y_cm_d, int NX, int NY, int *x_triple_array_d, int *y_triple_array_d, 
        int cnt_x, int cnt_y, Triple *triple_array_d, float * rmsd_array_d){
    
    int row = blockIdx.x * TILE_WIDTH + threadIdx.x;
    int col = blockIdx.y * TILE_WIDTH + threadIdx.y;
    int i;

    
    int *row_data;
    float rmsd = -5; 
    
    int triple_x[3], triple_y[3];
    // double q_init[4] = {0.0}; // no change

    
    if (row < cnt_x && col < cnt_y) {
        triple_array_d[row*cnt_y + col].rmsd = BAD_RMSD + 1;
        rmsd_array_d[row*cnt_y + col] = BAD_RMSD + 1;
        row_data = (int*)(((char*)x_triple_array_d) + (row * PITCH));
        
        for (i = 0; i < 3; ++i) {
            triple_x[i] =  row_data[i];
        }
        
        row_data = (int*)(((char*)y_triple_array_d) + (col * PITCH));
        
        for (i = 0; i < 3; ++i) {
            triple_y[i] =  row_data[i];
        }  
        
        for (i = 0; i <3; ++i) {
            triple_array_d[row*cnt_y + col].triple_x[i] = triple_x[i];
            triple_array_d[row*cnt_y + col].triple_y[i] = triple_y[i];
        }

        //if (!same_hand_triple_gpu(x_d, x_cm_d, y_d, y_cm_d, triple_x,triple_y, 3) ) return;

        distance_of_nearest_approach_gpu(x_d, x_cm_d, y_d, y_cm_d, triple_x, triple_y, 3, &rmsd);
        if (rmsd > BAD_RMSD) {
            return; 
        }
                    // insert values to array of structs
            
        triple_array_d[row*cnt_y + col].rmsd = rmsd;
        rmsd_array_d[row*cnt_y + col] = rmsd;
        
    }
    

/*
   if (!same_hand_triple_gpu(x_d, x_cm_d, y_d, y_cm_d, triple_x,triple_y, 3) ) return;
*/
          // if (!same_hand_triple(X_rep, x_triple, Y_rep, y_triple, 3)) continue;

/*
    distance_of_nearest_approach_gpu(x_d, x_cm_d, y_d, y_cm_d, triple_x, triple_y, 3, &rmsd);
    if (rmsd > BAD_RMSD) {
        return; 
    }
*/



            // insert values to array of structs
            
//    triple_array_d[row*cnt_y + col].rmsd = rmsd;
            //printf("%d\n", myid);

    
}

extern int insert_triple_to_heap_gpu(Representation* X_rep, Representation* Y_rep,
				     int ** x_triple_array, int ** y_triple_array, int x_triple_cnt,
				     int y_triple_cnt, PriorityQueue * heap) {
    

    double **x_db = X_rep->full;
    double **x_cm_db = X_rep->cm;
    double **y_db = Y_rep->full;
    double **y_cm_db = Y_rep->cm;
    int NX = X_rep->N_full;
    int NY = Y_rep->N_full;
    
    float **x = fmatrix(NX, 3);
    float **x_cm = fmatrix(NX, 3);;
    float **y = fmatrix(NY, 3);;
    float **y_cm = fmatrix(NY, 3);;
    
    int i,j;
    
    for(i = 0; i < NX; ++i) {
        for(j = 0; j < 3; ++j) {
            x[i][j] = x_db[i][j];
            x_cm[i][j] = x_cm_db[i][j];
        }
    }
    
    for(i = 0; i < NY; ++i) {
        for(j = 0; j < 3; ++j) {
            y[i][j] = y_db[i][j];
            y_cm[i][j] = y_cm_db[i][j];
        }
    }
    
    float *x_d;
    float *x_cm_d;
    float *y_d;
    float *y_cm_d;
    
    int *x_triple_array_d;
    int *y_triple_array_d;
    
    
    int cnt_x = x_triple_cnt;
    int cnt_y = y_triple_cnt;
    
    struct timespec requestStart, requestEnd;
    clock_gettime(CLOCK_REALTIME, &requestStart);
        
    
    Triple * triple_array = (Triple *) malloc(TOP_RMSD * sizeof(Triple));
    
    
    if (cudaSuccess != cudaMalloc(&x_d, NX * 3 * PITCH * sizeof(float))) printf("CUDA allocation error!\n");
    cudaMemcpy2D(x_d,PITCH,x[0],sizeof(float)*3,sizeof(float)*3,NX,cudaMemcpyHostToDevice);   
    
    if (cudaSuccess != cudaMalloc(&x_cm_d, NX * 3 * PITCH * sizeof(float))) printf("CUDA allocation error!\n");
    cudaMemcpy2D(x_cm_d,PITCH,x_cm[0],sizeof(float)*3,sizeof(float)*3,NX,cudaMemcpyHostToDevice);
    
    if (cudaSuccess != cudaMalloc(&y_d, NY * 3 * PITCH * sizeof(float))) printf("CUDA allocation error!\n");
    cudaMemcpy2D(y_d,PITCH,y[0],sizeof(float)*3,sizeof(float)*3,NY,cudaMemcpyHostToDevice);
    
    if (cudaSuccess != cudaMalloc(&y_cm_d, NY * 3 * PITCH * sizeof(float))) printf("CUDA allocation error!\n");
    cudaMemcpy2D(y_cm_d,PITCH,y_cm[0],sizeof(float)*3,sizeof(float)*3,NY,cudaMemcpyHostToDevice);
    
    if (cudaSuccess != cudaMalloc(&x_triple_array_d, cnt_x * 3 * PITCH * sizeof(int))) printf("CUDA allocation error!\n");
    cudaMemcpy2D(x_triple_array_d,PITCH,x_triple_array[0],sizeof(int)*3,sizeof(int)*3,cnt_x,cudaMemcpyHostToDevice);
    
    if (cudaSuccess != cudaMalloc(&y_triple_array_d, cnt_y * 3 * PITCH * sizeof(int))) printf("CUDA allocation error!\n");
    cudaMemcpy2D(y_triple_array_d,PITCH,y_triple_array[0],sizeof(int)*3,sizeof(int)*3,cnt_y,cudaMemcpyHostToDevice);  
  
    
    clock_gettime(CLOCK_REALTIME, &requestEnd);
    
    // call gpu and find triplets
    
    Triple * triple_array_d;   
    float * rmsd_array_d;
    
    size_t size = cnt_x * cnt_y * sizeof(Triple);
    size_t size_ratio = sizeof(Triple)/sizeof(float);

    
    if (size < MEM_SIZE_MAX) {
        if (cudaSuccess != cudaMalloc((void **)&triple_array_d, size)) printf("CUDA allocation error!\n");
        if (cudaSuccess != cudaMalloc((void**)&rmsd_array_d, cnt_x * cnt_y * sizeof(float))) printf("CUDA allocation error!\n");
    } else {
        if (cudaSuccess != cudaMalloc((void **)&triple_array_d, MEM_SIZE_MAX)) printf("CUDA allocation error!\n");
        if (cudaSuccess != cudaMalloc((void**)&rmsd_array_d, MEM_SIZE_MAX / size_ratio)) printf("CUDA allocation error!\n");
    }
 
         
    size_t free_mem, total_mem;
    size_t size_curr = size;
    size_t  y_triple_array_pos = 0;
    
    
    int not_enough_memory = 1;
    
    size_t cnt_y_rest = cnt_y;
    
    
    while (not_enough_memory) {
    
        clock_gettime(CLOCK_REALTIME, &requestStart);
        cudaMemGetInfo(&free_mem, &total_mem);
        
        if (size_curr > MEM_SIZE_MAX){
        
            cnt_y = MEM_SIZE_MAX/(cnt_x * sizeof(Triple));
            cnt_y_rest -= cnt_y;
            
            size = MEM_SIZE_MAX;
            size_curr = cnt_x * cnt_y_rest * sizeof(Triple); 
        } else {
            cnt_y = cnt_y_rest;
            size = size_curr;
            not_enough_memory = 0;
        }
        
        int n_blocks_x = cnt_x/TILE_WIDTH + (cnt_x%TILE_WIDTH == 0 ? 0:1);
        int n_blocks_y = cnt_y/TILE_WIDTH + (cnt_y%TILE_WIDTH == 0 ? 0:1);
        dim3 numBlocks(n_blocks_x, n_blocks_y); 
        dim3 threadsPerBlock(TILE_WIDTH, TILE_WIDTH);
        
        int * y_triple_array_d_curr = (int*)(((char*)y_triple_array_d) + (y_triple_array_pos * PITCH));
        
        find_triplets_gpu <<<numBlocks, threadsPerBlock>>> (x_d, x_cm_d, y_d, y_cm_d, NX, NY, x_triple_array_d, y_triple_array_d_curr, 
        cnt_x, cnt_y, triple_array_d, rmsd_array_d);
 
        cudaThreadSynchronize();
       
        
	Triple * triple_array_output_d;
    
        size_t no_of_out_pairs = TOP_RMSD < (cnt_x * cnt_y)? TOP_RMSD : cnt_x * cnt_y; 
        if (cudaSuccess != cudaMalloc((void **)&triple_array_output_d, no_of_out_pairs*sizeof(Triple))) printf("CUDA allocation error!\n");    
        
        thrust::device_vector<int>  indices(cnt_x * cnt_y); 
        thrust::sequence(indices.begin(),indices.end());
        thrust::device_ptr<Triple> structures(triple_array_d);
        thrust::device_ptr<Triple> structures_out(triple_array_output_d);
        thrust::device_ptr<float>rmsd(rmsd_array_d);
         
        //thrust::sort_by_key(rmsd, rmsd + cnt_x * cnt_y, structures);
        thrust::sort_by_key(rmsd, rmsd + cnt_x * cnt_y, indices.begin());
        

        thrust::device_vector<int>::iterator iter = indices.begin() + no_of_out_pairs;
        
        thrust::gather(indices.begin(), iter, structures, structures_out);

        // thrust::sort(structures, structures+ cnt_x * cnt_y, greater_rmsd());
        
        
        
        //cudaMemcpy(triple_array, triple_array_d, no_of_out_pairs * sizeof(Triple), cudaMemcpyDeviceToHost);
        cudaMemcpy(triple_array, triple_array_output_d, no_of_out_pairs * sizeof(Triple), cudaMemcpyDeviceToHost);
/*
        cudaFree(triple_array_d);
        cudaFree(y_triple_array_d);
        cudaFree(rmsd_array_d);
*/
 
        
        int m;
        for(m = 0; m < no_of_out_pairs; ++m) {
            Insert(triple_array[m], *heap);
        }
    
 
        y_triple_array_pos += cnt_y;
        
        cudaFree(triple_array_output_d);
        
    }
    
    
    cudaFree(x_d);
    cudaFree(y_d);
    cudaFree(x_cm_d);
    cudaFree(y_cm_d);
    cudaFree(x_triple_array_d);
    cudaFree(y_triple_array_d);
    
    cudaFree(triple_array_d);

    cudaFree(rmsd_array_d);
    
    free(triple_array);
    free_fmatrix(x);
    free_fmatrix(x_cm);
    free_fmatrix(y);
    free_fmatrix(y_cm);
   
    return 0;
} 

#ifdef	__cplusplus
}
#endif
