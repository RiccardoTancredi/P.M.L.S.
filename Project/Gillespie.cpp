#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <vector>

using namespace std;

/* our custom libraries */
#include "save_data.h" // Include the custom header file

double sum_vector(vector<double> &x);

void set_matrix_to_zero(double **matrix, int N, int M){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            matrix[i][j] = 0.;
        }
    }
}


int main(){
    /*
    =======================
    = Exercise: BVP model =
    =======================
    */

    // Set up the rnd generator
    int iseed=0; 
    srand(static_cast<unsigned>(iseed));

    // Initial conditions
    double a = 0.7;
    double b = 0.8;
    double c = 3.;
    double z = -0.4;

    vector<double> Omega = {1.*pow(10, 3)};

    int t_tot = 2000;
    double dt = pow(10, -3);

    double t = 0.;
    double lambda_C = 0.;
    double u = 0.;
    double tau = 0.; int remain_time;
    int new_state = 0;
    vector <double> p_j = {0., 0., 0., 1.};
    vector<double> rates = {0., 0., 0., 0.};
    int M = 2; int N = t_tot/dt;
    // cout << "N = " << N << endl;
    double **population = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        population[i] = (double *)malloc(M * sizeof(double));
    }

    string s = "BVP_res_limit_cycle.txt";
    int k = 0;
    vector<double> C = {0., 0.};
    
    for (double omega : Omega){
        C[0] = 3.2*omega; C[1] = 1.375*omega;
        set_matrix_to_zero(population, N, M);
        k = 0; t = 0;
        population[k][0] = C[0]; population[k][1] = C[1];
        k += 1; t += dt;

        cout << "Omega = " << omega << endl;

        while (t < t_tot){
            rates[0] = c*C[1] + 2*c*C[0]*C[0]/omega; 
            rates[1] = 3*c*C[0] + 4*c*omega/3. - c*z*omega + c*C[0]*C[0]*C[0]/(3.*omega*omega); 
            rates[2] = (a+2+2*b)*omega/c; 
            rates[3] = b*C[1]/c + C[0]/c;
            lambda_C = sum_vector(rates); 
            // cout << "lambda_C = " << lambda_C << endl;
            // Extract tau from the exponential distribuation of mean 1/lambda_C
            u = (double)rand() / (RAND_MAX);
            tau = -log(1-u)/lambda_C; 
            // Remain in configuration C for t = tau seconds
            // cout << "tau = " << tau << endl;
            if (tau > dt){
                remain_time = min(int(tau/dt), N-k);
                // cout << "Remain time = " << remain_time << endl;
                for (int i = k; i < k + remain_time; i++){
                    population[i][0] = C[0]; population[i][1] = C[1];
                }
                k += remain_time;
                t += dt*remain_time;
            }
            
            // Then, pick a state j != C
            // For the BVP there are 5 reachable states
            p_j[0] = rates[0]/lambda_C; 
            p_j[1] = (rates[0] + rates[1])/lambda_C; 
            p_j[2] = (rates[0] + rates[1] + rates[2])/lambda_C;
            
            // u = (double)rand() / (RAND_MAX);
            if (u < p_j[0]){
                new_state = 0;
            }
            else if (u < p_j[1] && u > p_j[0]){
                new_state = 1;
            }
            else if (u < p_j[2] && u > p_j[1]){
                new_state = 2;
            }
            else{
                new_state = 3;
            }
            // Change state
            // cout << "State changed in " << new_state << endl;
            if (new_state == 0){
                C[0] += 1;
            }
            else if (new_state == 1){
                C[0] -= 1;
            }
            else if (new_state == 2){
                C[1] += 1;
            }
            else{
                C[1] -= 1;
            }   
            // cout << "The new population is C' = " << C[0] << ", " << C[1] << endl; 

            if (k < N){
                population[k][0] = C[0]; population[k][1] = C[1];
            }
            
            t += dt; k += 1;
        }
        // cout << "p_j[2] = " << p_j[2] << endl;
        // cout << "lambda_C = " << lambda_C << endl;
        // cout << "Population C = " << C[0] << " " << C[1] << endl;
        // cout << "Rates :" << rates[0] << " " << rates[1] << " " << rates[2] << " " << rates[3] << " " << endl; 
        saveMatrixToFile(population, N, M, s);

    }

    free(population);

    return 0;
}



double sum_vector(vector<double> &x){
    double res = 0.;
    for(auto val : x){
        res += val;
    }
    return res;
}