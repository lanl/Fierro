#include <stdio.h>
#include <math.h>
#include <matar.h>
#include <limits.h>
#include <time.h>
#include <chrono>

using namespace mtr; // matar namespace

//Helper function to prefill a graph, seems generically useful
void graphFiller(int n, int diag, int off_diag, CArray<int> &G){
    int i,j;
    for(i = 0; i < n; i++){
        for(j =0; j < n; j++){                
            if(i == j){
		        G(i,j) = diag;
            } else {
			    G(i,j) = off_diag;
			}
        }
    }
}

// k := connect to k nearest neighbors
// n := size of graph
// p := rewire prob. should be 0 <= p <= 1
// G := an empty matrix to place the edges.
// A watts storgatz graph is one where each node is connected to it's k nearest neighbors but each edge has a small 
// rewiring chance. This is to show that a graph with mostly local connections can have short average shortest distances 
void wattsStorgatzGraph(int k, int n, double p, CArray<int> &G){
	int idx_forw, idx_back, random_edge, random_idx;
	double coin_flip;
	for(int i = 0; i < n; i++){
		for(int j = 1 ; j <= k; j++){
            idx_forw = (i + j) % n;
			idx_back = (i + n - j) % n;
			coin_flip = ((double) std::rand())/RAND_MAX;
            if(coin_flip < p){
	            // ramdom number from (k, n-k) 
                random_edge = rand() % (n - 2*k) + (k);
                // index is i + the above offset
				random_idx = (i + random_edge) % n;
                G(i, random_idx) = 1;
            }else{
				G(i,idx_forw) = 1;
			}
			G(i,idx_back) = 1;
		}
	}
}


void floydW(CArray<int> G, CArray<int> &res, int n_nodes){
	int i,j,k;
	graphFiller(n_nodes, 0, INT_MAX, res);
	for(i = 0; i<n_nodes; i++){
		for(j =0; j<n_nodes; j++){
			if(G(i,j) == 1){
				res(i,j) = 1;
			}
		}
	}
	for(k=0; k < n_nodes; k++){
		for(i = 0; i < n_nodes; i++){
			for(j = 0; j < n_nodes ; j++){
				if(i == j){
					continue;
				}
				int dist1 = res(i,k) + res(k,j);
			    int dist2 = res(i,j);	
				if(dist1 < 0){
                    dist1 = INT_MAX;
                }
                if(dist2 < 0){
                    dist2 = INT_MAX;
                }
                if(dist1 < dist2){
					res(i,j) = dist1;
				}
			}
		}
	}		
}

double averageDistance(CArray<int> G, int n){
    int i,j;
    double total = 0.0;
    for(i = 0; i<n; i++){
        for(j = 0; j<n; j++){
            if(G(i,j) > n){
                printf("Ohh dear, this shouldn't happen\n");
            }
            total += ((double) G(i,j)) / n;
        }
    }
    return total / ( n  );
}

void printer(CArray<int> G, int n){
    int i,j;
    for(i = 0; i< n; i++){
        for(j = 0; j<n; j++){
            printf("%d ",G(i,j));
        }
        printf("\n");
    }
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
    printf("%d, %.3f, %d,", node_size, rewire_p, k_nearest);

    auto start = std::chrono::high_resolution_clock::now(); // start clock
    auto G = CArray<int>(node_size, node_size);
    auto results = CArray<int> (node_size, node_size);
    wattsStorgatzGraph(k_nearest, node_size, rewire_p, G);
    auto lap = std::chrono::high_resolution_clock::now(); // start clock
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(lap-start);
    printf("%.2f,", elapsed.count() * 1e-9);
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
