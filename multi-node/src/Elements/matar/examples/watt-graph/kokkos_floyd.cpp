#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <matar.h>
#include <limits.h>
#include <Kokkos_Random.hpp>
#include <time.h>
#include <chrono>

using namespace mtr; // matar namespace

using gen_t = Kokkos::Random_XorShift64_Pool<DefaultExecSpace>;


// k := connect to k nearest neighbors
// n := size of graph
// p := rewire prob. should be 0 <= p <= 1
// G := an empty matrix to place the edges.
// A watts storgatz graph is one where each node is connected to it's k nearest neighbors but each edge has a small 
// rewiring chance. This is to show that a graph with mostly local connections can have short average shortest distances 
void wattsStorgatzGraph(int k, int n, double p, CArrayKokkos<int> &G){
    
    gen_t rand_pool64(5374857);
    CArrayKokkos<double> coins(n, k);
    CArrayKokkos<int> offsets(n, k);
    FOR_ALL(i , 0, n,
             j, 1, k+1, {
                // Get a random number state from the pool for the active thread
                gen_t::generator_type rand_gen = rand_pool64.get_state();

                // generate random numbers in the range (0,10]
                coins(i, j-1) = ((double) rand_gen.urand64(10000)) / 10000.0;
                offsets(i, j-1) = rand_gen.urand64(n - 2*k) + k;
            });
    Kokkos::fence();
    FOR_ALL(i, 0, n,
           j, 1, k+1,{ 
                // Give the state back, which will allow another thread to acquire it
                int idx_forw = (i + j) % n;
                int idx_back = (i + n - j) % n;
                G(i,idx_back) = 1;
                if(coins(i, j-1) < p ){
                    int random_idx = offsets(i, j-1);
                    random_idx = (i + random_idx) % n;
                    G(i, random_idx) = 1;
                } else {
                    G(i, idx_forw) = 1;
                }
             });
    Kokkos::fence();
    return;
}


void floydW(CArrayKokkos<int> G, CArrayKokkos<float> &res, int n_nodes){
	int k;
    FOR_ALL(i, 0, n_nodes,
            j, 0, n_nodes, {
                if(G(i,j) == 1){
                    res(i,j) = 1;
                }
            });
    for(k = 0; k < n_nodes; k++){
        FOR_ALL(i, 0, n_nodes, 
                j, 0, n_nodes, {
                 if(i != j){
			
				    float dist1 = res(i,k) + res(k,j);
			        float dist2 = res(i,j);	
                    res(i,j) = (dist1 < dist2) ? dist1 : dist2;
                  }   
        });
    }
}

double averageDistance(CArrayKokkos<float> G, int n){
    double total = 0;
    double loc_sum;
    REDUCE_SUM(i, 0, n,
               j, 0, n,
               loc_sum, {
                    loc_sum += ((double) G(i,j)) /n ;
               }, total); 
    return  total /(n  );
}

void printer(CArrayKokkos<int> G, int n){
        int x = 0;
        x++;
        FOR_ALL(i, 0, n,
            j, 0, n,{
                printf("%d, %d) : %d \n", i, j, G(i,j));
        });
}


int main(int argc, char** argv){
    int node_size = 4000;
    double rewire_p = 0.0;
    int k_nearest = 6;
    if((argc > 4) || (argc == 1)){
        printf("Usage is ./test_kokkoks_floyd <number of nodes> <rewire prob.>\n");
        printf("Using default values: [number of nodes: %d] [rewire_prob : %.2f] [k_nearest : %d]\n", node_size, rewire_p, k_nearest); 
    } else {
       node_size = atoi(argv[1]);
       rewire_p = atof(argv[2]);
       k_nearest = atoi(argv[3]); 
    }
    printf("%d, %.5f, %d", node_size, rewire_p, k_nearest);
    Kokkos::initialize(); {
    
    auto start = std::chrono::high_resolution_clock::now(); // start clock
    CArrayKokkos<int> G(node_size, node_size);
    CArrayKokkos<float> results(node_size, node_size);
    FOR_ALL(i, 0, node_size,
            j, 0, node_size,{
            G(i,j) = 0;
            if(i == j){
                results(i,j) = 0;
            } else {
                results(i,j) =  std::numeric_limits<double>::infinity();
            }
        });
    
    wattsStorgatzGraph(k_nearest, node_size, rewire_p, G);
    
    auto lap = std::chrono::high_resolution_clock::now(); // start clock
    
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(lap-start);
    printf(", %.2f,", elapsed.count() * 1e-9);
    floydW(G, results, node_size);
    
    auto lap2 = std::chrono::high_resolution_clock::now(); // start clock
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(lap2-lap);
    printf("%.2f,", elapsed.count() * 1e-9);
    double average_steps = averageDistance(results, node_size);
    
    auto lap3 = std::chrono::high_resolution_clock::now(); // start clock
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(lap3-lap2);
    auto elapsed2 = std::chrono::duration_cast<std::chrono::nanoseconds>(lap3-start);

    printf("%.2f, %.2f, ", elapsed.count() * 1e-9, elapsed2.count() * 1e-9);
    printf("%f\n", average_steps);	
    
    }
    Kokkos::finalize(); 
}
