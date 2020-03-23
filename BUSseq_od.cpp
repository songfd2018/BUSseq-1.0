#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h> // pow, sqrt, lgamma
#include <cmath> //
#include "omprng.h"
# include <chrono> // time
#include <unistd.h> // getopt 
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <sys/stat.h>  //mkdir

using namespace std;

int rand_cate(double* _prop, omprng _rng){
    int res = 0;
    double u = _rng.runif();
    while(u > _prop[res]){
        u = u - _prop[res];
        res++;
    }
    return res;
}

int rand_Ber(double _prob, omprng _rng){
    int res = 1;
    double u = _rng.runif();
    if(u > _prob){
        res = 0;
    }
    return res;
}


int rand_NB(double _r, double _mu, omprng _rng){
    
    double lambda, nu, p, u;
    double STEP = 500;// steps for large lambda in Poisson distribution
    int res = 0;

    // sample lambda from gamma distribution
    lambda = _rng.rgamma(_r, _mu/_r);
        
    // sample x from Poisson distribution
    nu = lambda;
    p = 1.0;

    do{
        res ++;
        u = _rng.runif();
        p = p * u;
        while(p < 1 & nu > 0){
            if(nu > STEP){
                p = p * exp(STEP);
                nu = nu - STEP;
            }else{
                p = p * exp(nu);
                nu = 0;
            }
        }
    }while(p > 1);

    res--;
    return res;
}

void rand_Dir(double *_xi, int _K, omprng _rng, double *_pi){

    double *rn_gam = new double[_K];
    double sum_gam = 0.0;
    for(int k = 0; k < _K; k++){
        rn_gam[k] = _rng.rgamma(_xi[k], 1.0);
        sum_gam += rn_gam[k];
    }
    for(int k = 0; k < _K; k++){
        _pi[k] = rn_gam[k] / sum_gam;
    }

    delete [] rn_gam;
}

void _update_logmu(int _B, int* _nb, 
                int *_W, double _alpha, double* _beta, double* _nu, double* _delta,//parameter
                double* _logmu){

    int cell_index = 0;
    int k;
    //auxiliry
	for(int b = 0; b < _B; b++) {
		for (int i = 0; i < _nb[b]; i++) {
            k = _W[cell_index];
			_logmu[cell_index] = _alpha + _beta[k] + _nu[b] + _delta[cell_index];
            cell_index ++;
		}
	}
}

void _update_zx(int _B, int* _nb,
                double **_gamma, double *_phi, double *_logmu,
                int *_Y, omprng _rng, 
                int *_X, int *_Z){

    //int ind, ind_nu;
    //int ind_n = 0; //index the row of (b,i)	
    int cell_index = 0;
    double log_rat, u, temp_max, acc_rat;
    int temp_x;


    for (int b = 0; b < _B; b ++){
            for(int i = 0; i < _nb[b]; i ++){
			//ind = ind_n + j * _N;
		    //ind_nu = b + j * _B;
            // printf("Sampling for cell %d.\n", cell_index + 1);
            // if(cell_index == 66){
            //     printf("Gamma_b0 = %f.\n", _gamma[b][0]);
            //     printf("mu_bi = %f\n", exp(_logmu[cell_index]));
            //     printf("phi = %f\n", _phi[b]);
            // }
            
			if (_Y[cell_index] == 0) {
				//update z_{big}
				if (_X[cell_index] == 0) {
                    // printf("Sampling for Z.\n");
				    log_rat = _gamma[b][0];
					_Z[cell_index] = rand_Ber(1 / (1 + exp(-log_rat)),_rng);
				}else {
					_Z[cell_index] = 1;
				}// end of if X == 0

				//Rprintf("z %d %d %d = %d\n", b, i, j, Dropout[index]);
				//update x_{big}
				if (_Z[cell_index] == 1) {// No dropout event
                    // printf("Sampling for X.\n");
					//MH sampling
					temp_x = rand_NB(_phi[b], exp(_logmu[cell_index]), _rng);
                    // printf("Sampling for X.\n");
					u = _rng.runif();
                    //make the exponential be compuational
                    temp_max = 0.0;
                    if(temp_max < - _gamma[b][0] - _gamma[b][1] * temp_x){
                        temp_max = - _gamma[b][0] - _gamma[b][1] * temp_x;
                    }
                    if(temp_max < - _gamma[b][0] - _gamma[b][1] * _X[cell_index]){
                        temp_max = - _gamma[b][0] - _gamma[b][1] * _X[cell_index];
                    }

					acc_rat = (exp(-temp_max) + exp( - _gamma[b][0] - _gamma[b][1] * _X[cell_index] - temp_max));
                    acc_rat = acc_rat / (exp(-temp_max) + exp( - _gamma[b][0] - _gamma[b][1] * temp_x - temp_max));
					if (u < acc_rat) {
						_X[cell_index] = temp_x;
					}
				}else {
					_X[cell_index] = 0;
				}// end of if Z == 0
			}// end of if Y == 0
			cell_index ++;
		}// end of i for
	}// end of b for
}

void _update_zx_optional(int _B, int* _nb, bool* drop,
                double **_gamma, double *_phi, double *_logmu,
                int *_Y, omprng _rng, 
                int *_X, int *_Z){

    //int ind, ind_nu;
    //int ind_n = 0; //index the row of (b,i)	
    int cell_index = 0;
    double log_rat, u, temp_max, acc_rat;
    int temp_x;


    for (int b = 0; b < _B; b ++){
        if(drop[b]){
            for(int i = 0; i < _nb[b]; i ++){
			//ind = ind_n + j * _N;
		    //ind_nu = b + j * _B;
            // printf("Sampling for cell %d.\n", cell_index + 1);
            // if(cell_index == 66){
            //     printf("Gamma_b0 = %f.\n", _gamma[b][0]);
            //     printf("mu_bi = %f\n", exp(_logmu[cell_index]));
            //     printf("phi = %f\n", _phi[b]);
            // }
            
			if (_Y[cell_index] == 0) {
				//update z_{big}
				if (_X[cell_index] == 0) {
                    // printf("Sampling for Z.\n");
				    log_rat = _gamma[b][0];
					_Z[cell_index] = rand_Ber(1 / (1 + exp(-log_rat)),_rng);
				}else {
					_Z[cell_index] = 1;
				}// end of if X == 0

				//Rprintf("z %d %d %d = %d\n", b, i, j, Dropout[index]);
				//update x_{big}
				if (_Z[cell_index] == 1) {// No dropout event
                    // printf("Sampling for X.\n");
					//MH sampling
					temp_x = rand_NB(_phi[b], exp(_logmu[cell_index]), _rng);
                    // printf("Sampling for X.\n");
					u = _rng.runif();
                    //make the exponential be compuational
                    temp_max = 0.0;
                    if(temp_max < - _gamma[b][0] - _gamma[b][1] * temp_x){
                        temp_max = - _gamma[b][0] - _gamma[b][1] * temp_x;
                    }
                    if(temp_max < - _gamma[b][0] - _gamma[b][1] * _X[cell_index]){
                        temp_max = - _gamma[b][0] - _gamma[b][1] * _X[cell_index];
                    }

					acc_rat = (exp(-temp_max) + exp( - _gamma[b][0] - _gamma[b][1] * _X[cell_index] - temp_max));
                    acc_rat = acc_rat / (exp(-temp_max) + exp( - _gamma[b][0] - _gamma[b][1] * temp_x - temp_max));
					if (u < acc_rat) {
						_X[cell_index] = temp_x;
					}
				}else {
					_X[cell_index] = 0;
				}// end of if Z == 0
			}// end of if Y == 0
			cell_index ++;
		}// end of i for
        }else{
            cell_index = cell_index + _nb[b];
        }
	}// end of b for
}

double _update_alpha(int _B, int* _nb,//dimension
                    double _mu_a, double _sigma_a,//prior
                    int *_W, double _alpha, double *_beta, double *_nu, double *_delta, double *_phi, //parameter
                    int *_X, omprng _rng){
    //int ind, ind_beta, ind_nu;
    //int ind_n; //index the row of (b,i)

    int cell_index = 0;

    //proposal
    double alpha_iter = _rng.rnorm(_alpha, 0.1);
    double logr = 0.0;
    int k;
    double res;

    //prior
    logr = logr - pow(alpha_iter - _mu_a, 2.0) / 2 / pow(_sigma_a , 2.0) + pow(_alpha - _mu_a, 2.0) / 2 / pow(_sigma_a, 2.0);

    for (int b = 0; b < _B; b++) {
		for (int i = 0; i < _nb[b]; i++) {
            k = _W[cell_index];

            //numerator
			logr = logr + alpha_iter * _X[cell_index] - (_phi[b] + _X[cell_index]) * log(_phi[b] + exp(alpha_iter + _beta[k] + _nu[b] + _delta[cell_index]));
			//denomerator
			logr = logr - _alpha * _X[cell_index] + (_phi[b] + _X[cell_index]) * log(_phi[b] + exp(_alpha + _beta[k] + _nu[b] + _delta[cell_index]));

            cell_index ++;
        }
    }

    if (logr > log(_rng.runif())) {
		res = alpha_iter;
	}else{
        res = _alpha;
    }

    return res;
}

void _update_l(int _K,//dimension
                double _p, double _tau0, double _tau1,//prior
                double *_beta, omprng _rng,//parameter
                int *_L){
    //int ind_beta;
    double log_rat;
	for (int k = 1; k < _K; k++) {
		//ind_beta = j + k * _G;
	    log_rat = 0.0; //the odds ratio of L_k = 1
 		log_rat += log(_p) - log(1 - _p);
		log_rat += - log(_tau1)/2.0 + log(_tau0)/2.0;
		log_rat += - pow(_beta[k], 2) / 2.0 / _tau1;
        log_rat += pow(_beta[k], 2) / 2.0 / _tau0;
		_L[k] = rand_Ber(1.0 / (1.0 + exp(-log_rat)), _rng);
	}
	
}

void _update_beta(int _B, int *_nb, int _K,//dimension 
                double _tau0, double _tau1, int *_L,//prior
                int *_W, double _alpha, double *_nu, double *_delta, double *_phi, //parameter 	
                int *_X, omprng _rng,//latent variable
                double *_beta){
        
    //index the row of (b,i)
    int cell_index, k;
    double *beta_iter = new double[_K];
    double *logr = new double[_K];
 
    for (k = 1; k < _K; k++) {
        //symmetric proposal
        beta_iter[k] = _rng.rnorm(_beta[k], 0.1);
        logr[k] = 0.0;

        //prior
        if (_L[k] == 1) {
			logr[k] = logr[k] - pow(beta_iter[k], 2.0) / 2 / _tau1 + pow(_beta[k], 2.0) / 2 / _tau1;
		}else{
			logr[k] = logr[k] - pow(beta_iter[k], 2.0) / 2 / _tau0 + pow(_beta[k], 2.0) / 2 / _tau0;
		}
    }

    cell_index = 0;
	for (int b = 0 ; b < _B; b++){
		for (int i = 0; i < _nb[b]; i++) {
            k = _W[cell_index];        
			//numerator
			logr[k] += beta_iter[k] * _X[cell_index]; 
            logr[k] += - (_phi[b] + _X[cell_index]) * log(_phi[b] + exp(_alpha + beta_iter[k] + _nu[b] + _delta[cell_index]));
			//denomerator
			logr[k] += - _beta[k] * _X[cell_index];
            logr[k] += (_phi[b] + _X[cell_index]) * log(_phi[b] + exp(_alpha + _beta[k] + _nu[b] + _delta[cell_index]));
					
            cell_index ++;
		}	
	}

    for (k = 1; k < _K; k++) {
        if (logr[k] > log(_rng.runif())) {
			_beta[k] = beta_iter[k];
		}
    }
    delete [] logr;
    delete [] beta_iter;
}

void _update_nu(int _B, int *_nb,
                double *_mu_c, double _sigma_c,//prior
                int *_W, double _alpha, double *_beta, double *_delta, double *_phi,//parameter
                int *_X, omprng _rng, //latent variable
                double *_nu){

    //int ind, ind_beta, ind_nu, ind_n;

	int cell_index = _nb[0];
    double nu_iter;
    double logr;
    int k;

	for(int b = 1; b < _B; b ++){
		//ind_nu = b + j * _B;
				
        //proposal
        nu_iter = _rng.rnorm(_nu[b], 0.1);
        logr = 0.0;

        //prior
        logr = logr - pow(nu_iter - _mu_c[b], 2.0) / 2 / pow(_sigma_c, 2.0) + pow(_nu[b] - _mu_c[b], 2.0) / 2 / pow(_sigma_c, 2.0);

        for(int i = 0; i < _nb[b]; i ++){

            //ind_beta = j + _w[ind_n] * _G;
			//ind = j * _N + ind_n;
            k = _W[cell_index];
            
            //numerator
			logr += nu_iter * _X[cell_index];
			logr += - (_phi[b] + _X[cell_index]) * log(_phi[b] + exp(_alpha + _beta[k] + nu_iter + _delta[cell_index]));

            //denomerator
			logr += - _nu[b] * _X[cell_index];
            logr += (_phi[b] + _X[cell_index]) * log(_phi[b] + exp(_alpha + _beta[k] + _nu[b] + _delta[cell_index]));
            
            cell_index ++;

            //ind_n = ind_n + 1;
        }

        if (logr > log(_rng.runif())) {
			_nu[b] = nu_iter;
		}			
	}
}

void _update_phi(int _B, int *_nb,
                double *_phi_prior,//prior
                double *_logmu,//parameter
                int *_X, omprng _rng, //latent variable
                double *_phi){
    //int ind, ind_nu, ind_n;
    int cell_index = 0;
    double phi_iter, logr;

	for (int b = 0; b < _B; b++) {
		
        phi_iter = _rng.rgamma(_phi[b], 1);
		logr = 0.0;

		for (int i = 0; i < _nb[b]; i++) {

			//numerator
			logr += lgamma(phi_iter + _X[cell_index]) - lgamma(phi_iter);
			logr += phi_iter * log(phi_iter) - (phi_iter + _X[cell_index]) * log(phi_iter + exp(_logmu[cell_index]));
			
            //denomerator
			logr += - lgamma(_phi[b] + _X[cell_index]) + lgamma(_phi[b]);
			logr += - _phi[b] * log(_phi[b]) + (_phi[b] + _X[cell_index]) * log(_phi[b] + exp(_logmu[cell_index]));
			
            cell_index ++;
		}
				
        //prior
		logr += (_phi_prior[0] - 1) * log(phi_iter) - _phi_prior[1] * phi_iter;
		logr += - (_phi_prior[0] - 1) * log(_phi[b]) + _phi_prior[1] * _phi[b];
		
        //proposal
		logr += (phi_iter - 1.0) * log(_phi[b]) + phi_iter - lgamma(phi_iter);
		logr += - (_phi[b] - 1.0) * log(phi_iter) - _phi[b] + lgamma(_phi[b]);

		if (logr > log(_rng.runif())) {
			_phi[b] = phi_iter;
		}
	}
}

int main(int argc, char** argv) {

    //////////////////////////////////////////////////////////
    //  0. Input the project and version from command line  //
    //////////////////////////////////////////////////////////
    string proj, wk_dir, readcount_dir;
    int ver, K, rep;
    int iter_max, iter_out; // the overall iteration number and the number of iterations to print the posterior sampling 
    int seed; // seed for random number generator 
    int nc; // number of codes for parallel
    int opt;

    //while ((opt = getopt(argc, argv, "p:v:N:G:B:K:r:i:s:")) != -1) {
    while ((opt = getopt(argc, argv, ":d::r::p:v:K:i:o:s:c:")) != -1) {
        switch (opt) {
        case 'd':
            if(optarg){
                wk_dir.assign(optarg);
            }else{
                wk_dir.assign("./");
            }
            cout << "The working directory is " << wk_dir << endl;
            break;
	case 'r':
            if(optarg){
                readcount_dir.assign(optarg);
            }else{
                readcount_dir.assign("./");
            }
            cout << "The directory to load read count data is " << readcount_dir << endl;
            break;        
	case 'p':
            proj.assign(optarg);
            break;
        case 'v':
            ver = atoi(optarg);
            break;
        case 'K':
            K = atoi(optarg);
            break;
        case 'i':
            iter_max = atoi(optarg);
            break;
        case 'o':
            iter_out = atoi(optarg);
            break;
        case 's':
            seed = atoi(optarg);
            break;
        case 'c':
            nc = atoi(optarg);
            break;
        case '?':
            cout << "Unknown option " << optarg << endl;
            break;
        case ':':
            cout << "Missing option for " << optarg << endl;
            break;
        case 1:
            cout << "Non-option arg: " << optarg << endl;
            break;
        default: /* '?' */
            cout << "Error Usage: " << argv[0] <<" [-d WorkingDir] [-p ProjectName] [-v Version] [-K CellTypeNumber] [-r ReplicationIndex] [-i IterationNumber] [-o NumberPerOutput] [-s Seed] [-c Cores] name" << endl;
            exit(EXIT_FAILURE);
        }
    }

    // set count data matrix and the number of cells in each batch
	int iter_noupdate = 500;
    string fcount, fdim;
    double tau1 = 50.0;
    fcount = readcount_dir + "count_data_";
    fcount = fcount + proj; 
    fcount = fcount + "_v";
    fcount = fcount + to_string(ver);
    fcount = fcount + ".txt";
    fdim = readcount_dir + "dim_";
    fdim = fdim + proj;
    fdim = fdim + "_v";
    fdim = fdim + to_string(ver);
    fdim = fdim + ".txt";

    // Create directory to store the results
    //string output_dir = wk_dir + "result/";
    //int check = mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    //if (!check) 
    //    cout << "Directory " << output_dir << " is created." << endl; 
    //else { 
    //    cout << "Directory " << output_dir << " already exists." << endl;
    //} 

    auto start_overall = chrono::system_clock::now();

	// Set the number of cores for parallel
    omp_set_num_threads(nc);

    ////////////////////////
    // 1. Load count data //
    ////////////////////////
    // allocate memory for the observation
    int N, G, B;

    ifstream countFile, dimFile;
    dimFile.open(fdim.c_str());
    if(!dimFile){
        cout << "Unable to open the file of dimensions from " << fdim << endl;
        exit(1); // terminate with error
    }

    // the number of cells
    dimFile >> N;
    cout << "There are totally " << N << " cells and ";

    // the number of genes
    dimFile >> G;
    cout << G << " genes in ";

    // the number of batches
    dimFile >> B;
    cout << B << " batches." << endl;
    cout << "Loading the number of cells in each batch from " << fdim << endl;

    // the number of cells in each batch
    int *nb = new int[B];
    for(int b = 0; b < B; b++){
        dimFile >> nb[b];
        cout << "The " << (b+1) << "-th batch has "<< nb[b] <<" cells." << endl;
    }

    // whether there are dropout events in each batch or not
    bool *Drop_ind = new bool[B];
    bool All_Drop = true;
    int di;
    for(int b = 0; b < B; b++){
        dimFile >> di;
        if(di == 1){
            Drop_ind[b] = true;
            cout << "The " << (b+1) << "-th batch has dropout events." << endl;
        }else if(di == 0){
            Drop_ind[b] = false;
            All_Drop = false;
            cout << "The " << (b+1) << "-th batch does not have any dropout event." << endl;
        }else{
            cout << "Unable to figure out whether there are dropout events in the " << (b+1) << "-th batch. Please enter 0 or 1." << endl;
            exit(1); // terminate with error
        }
    }
    

    dimFile.close();

    int **Y = new int*[G];
    for(int g = 0; g < G; g++){
        Y[g] = new int[N];// for parallel G
    }

    cout << "Loading the count data matrix from " << fcount << endl;
    
    countFile.open(fcount.c_str());
    if (!countFile) {
        cout << "Unable to open count data file" << endl;
        exit(1); // terminate with error
    }
    for(int g = 0; g < G; g++){
        for(int i = 0; i < N; i++){
            countFile >> Y[g][i];
        }
    }
    countFile.close();



    ///////////////////////////////////////////////////////
    // 2. Set hyper-parameters and Initialize parameters //
    ///////////////////////////////////////////////////////
    cout << "Set initial values." << endl;
    omprng MCMC_Rng;
	// cout << "Initializing rng." << endl;
    MCMC_Rng.fixedSeed(seed);
    
    // allocate memory for parameters
    // cell type proportion
	// cout << "Initializing pi." << endl;
    double temp_sum = (K + 1) * K / 2;// sum up 1 to K
    double **prop = new double*[B];
    for(int b = 0; b < B; b++){
        prop[b] = new double[K];
        //if (b == 0) cout << "Take a look at pi[0]: ";
        for(int k = 0; k < K; k ++){
            prop[b][k] = (K - k)/temp_sum;
            //prop[b][k] = 1.0 / K;
        }
        //if (b == 0) cout << endl;
        
    }

    // check rand_cate
    //cout << "Checking the rand_cate function: " << endl; 
    //for(int i = 0; i < 30; i ++){
    //    cout << rand_cate(prop[0],MCMC_Rng) << " ";
    //}
    //cout << endl;

    // cell type indicator
    int *W = new int[N];
    int *pointer_w = W;
	// cout << "Initializing W." << endl;
    // cout << "W_init is: " << endl; 
    for(int b = 0; b < B; b ++){
        //if (b == 0) cout << "Take a look at W: ";
        for(int i = 0; i < nb[b]; i++){
            *pointer_w = rand_cate(prop[b],MCMC_Rng);
            //if( b==0 ){
            //    cout << *pointer_w << " ";
            //}
            pointer_w++;
        }
        //if (b == 0) cout << endl;
    }
    // cout << endl;
    
    // cell-specific effect
    // cout << "delta_init is: " << endl;
	// cout << "Initializing delta." << endl;
    // cout << "The empirical mean of delta: " << endl;
    double *delta = new double[N];
	double *mu_d = new double[N];
    double *pointer_del = delta;
    int cell_index = 0;
    for(int b = 0; b < B; b ++){
        double sum_y0;
        for(int i = 0; i < nb[b]; i++){
            
            if(i == 0){
                *pointer_del = 0.0;// the first cell in each batch as reference
                sum_y0 = 0.0;// sum_g log(1+Y_b0g)
                for(int g = 0; g < G; g++){
                    sum_y0 += Y[g][cell_index];
                }
                pointer_del++;
            }else{
                double sum_y;
                sum_y = 0.0;// sum_g log(1+Y_big)
                for(int g = 0; g < G; g++){
                    sum_y += Y[g][cell_index];
                }
                *pointer_del =  log(sum_y) - log(sum_y0);
				mu_d[cell_index] = *pointer_del;
                // cout << "mu_d[" << cell_index << "] = " << mu_d[cell_index] << endl;
                pointer_del++;
				 
            }
            cell_index ++;
        }
    }
     
    // baseline effects
    double *alpha = new double[G];
	// cout << "Setting mu_a" << endl;
	double *mu_a = new double[G];
	// cout << "Setting beta" << endl;
    // cell-type specific effects
    double **beta = new double*[G];
    for(int g = 0; g < G; g++){
        beta[g] = new double[K];
    }
    
    // cout << "Calculate the empirical mean" << endl;
    for(int g = 0; g < G; g++){
         
        double *sum_logy = new double[K];
        int *count_type = new int[K];
        for(int k = 0; k < K; k++){
            sum_logy[k] = 0.0;
            count_type[k] = 0;
        }
        cell_index = 0;
        
        // only in the first batch
        for(int i = 0; i < nb[0]; i++){
            int k = W[cell_index];
            sum_logy[k] += log(1+Y[g][cell_index]/exp(delta[cell_index]));
            count_type[k] ++;
            cell_index ++;
        }

        alpha[g] = sum_logy[0]/(count_type[0] + 1);// + 1 to prevent count_type[0] = 1
		mu_a[g] = alpha[g];
        beta[g][0] = 0.0;
        for(int k = 1; k < K; k++){
            beta[g][k] = sum_logy[k]/(count_type[k] + 1) - alpha[g];
        }

        delete [] sum_logy;
        delete [] count_type;
    }

    // batch effects
    double **nu = new double*[G];
	double *mu_c = new double[B];
    for(int b = 0; b < B; b++){
        mu_c[b] = 0.0;
    }
    // cout << "The empirical mean of nu: " << endl;
    for(int g = 0; g < G; g++){
        nu[g] = new double[B];
        nu[g][0] = 0.0;
        double sum_logy0;
        sum_logy0 = 0.0;
        cell_index = 0;
        for(int i = 0; i < nb[0]; i++){
            sum_logy0 += log(1+Y[g][cell_index]/exp(delta[cell_index]));
            cell_index ++;
        }

        for(int b = 1; b < B; b++){
            double sum_logy = 0.0;
            for(int i = 0; i < nb[b]; i++){
                sum_logy += log(1+Y[g][cell_index]/exp(delta[cell_index]));
                cell_index ++;
            }
            nu[g][b] = sum_logy / nb[b] - sum_logy0 / nb[0];
			mu_c[b] += nu[g][b];
        }
		
    }
    for(int b = 0; b < B; b++){
        mu_c[b] = mu_c[b] / G;
        // cout << "mu_c[" << b << "] = " << mu_c[b] << endl;
    }
    
    // over-dispersion parameters
	// cout << "Initializing phi." << endl;
    double **phi = new double*[G];
    for(int g = 0; g < G; g++){
        phi[g] = new double[B];
        for(int b = 0; b < B; b++){
            phi[g][b] = 5.0;
        }
    }

    // intercept and odds ratio of dropout events
	// cout << "Initializing gamma." << endl;
    double **gamma = new double*[B];
    for(int b = 0; b < B; b++){
        gamma[b] = new double[2];
        gamma[b][0] = 0.0;
        if(Drop_ind[b]){
            gamma[b][1] = -0.1;
        }else{
            gamma[b][1] = 0.0;
        }

    } 

    // underlying true read count
    int **X = new int*[G];
    for(int g = 0; g < G; g++){
        X[g] = new int[N];
        for(int i = 0; i < N; i++){
            X[g][i] = Y[g][i];
        }
    }
    
    // dropout indicator
    int **Z = new int*[G];
    for(int g = 0; g < G; g++){
        Z[g] = new int[N];
        for(int i = 0; i < N; i++){
            Z[g][i] = 0;
        }
    }

    // intrinsic gene proportion and variance of non-intrinsic genes
    double p, tau0; // p and tau0 as well as tau1
    p = MCMC_Rng.runif(0,0.5);
    tau0 = 0.005;
    // for beta spike and slab prior
    

    // intrinsic gene indicators
    int **L = new int*[G];
#pragma omp parallel for
    for(int g = 0; g < G; g++){
        L[g] = new int[K];
        L[g][0] = 0;
        for(int k = 1; k < K; k++){
            double log_rat;  // log odds of gene g in the cell type k being an intrinsic gene
            log_rat = log(p/(1-p));
            log_rat = log_rat - log(tau1)/2 - pow(beta[g][k], 2) / 2 / tau1;
            log_rat = log_rat + log(tau0)/2 + pow(beta[g][k], 2) / 2 / tau0;
            L[g][k] = rand_Ber(1/(1+exp(-log_rat)), MCMC_Rng);
        }
    }

    // hyper-parameters
    // for prop Dir prior
    double xi = 2.0;

    // for alpha Normal prior
    double sigma_a = sqrt(5);

    // for nu Normal prior
    double sigma_c = sqrt(5);

    // for delta Normal prior
    double sigma_d = sqrt(5);

    // for gamma_b0 Normal prior
    double sigma_z = sqrt(3);
    
    // for gamma_b1 Gamma prior
    double gamma_prior[2];
    gamma_prior[0] = 0.001;
    gamma_prior[1] = 0.01;

    // for phi Gamma prior
    double phi_prior[2];
    phi_prior[0] = 1.0;
    phi_prior[1] = 0.1;

    // for p  Beta prior
    double p_prior[2];
    p_prior[0] = 1.0;
    p_prior[1] = 3.0;

    // for tau0 Gamma prior
    double tau0_prior[2];
    tau0_prior[0] = 2;
    tau0_prior[1] = 0.01;

    cout << "Allocate memory for recording the posterior sampling." << endl;
    // Allocate memory to store posterior sampling
    double **alpha_post = new double*[iter_out];
    for(int iter = 0; iter < iter_out; iter ++){
        alpha_post[iter] = new double[G];
    }

    double **beta_post = new double*[iter_out];
    for(int iter = 0; iter < iter_out; iter ++){
        beta_post[iter] = new double[G * K];
    }

    double **nu_post = new double*[iter_out];
    for(int iter = 0; iter < iter_out; iter ++){
        nu_post[iter] = new double[G * B];
    }

    double **delta_post = new double*[iter_out];
    for(int iter = 0; iter < iter_out; iter++){
        delta_post[iter] = new double[N];
    }
    
    double **gamma_post = new double*[iter_out];
    for(int iter = 0; iter < iter_out; iter++){
        gamma_post[iter] = new double[B * 2];
    }

    double **phi_post = new double*[iter_out];
    for(int iter = 0; iter < iter_out; iter++){
        phi_post[iter] = new double[G * B];
    }

    double **pi_post = new double*[iter_out];
    for(int iter = 0; iter < iter_out; iter++){
        pi_post[iter] = new double[B * K];
    }

    int **w_post = new int*[iter_out];
    for(int iter = 0; iter < iter_out; iter++){
        w_post[iter] = new int[N];
    }
    
    double *p_post = new double[iter_out];
    double *tau0_post = new double[iter_out];

    int **l_post = new int*[iter_out];
    for(int iter = 0; iter < iter_out; iter++){
        l_post[iter] = new int[G * K];
    }
    
    //////////////////////
    // 3. MCMC sampling //
    //////////////////////
    // auxiliary variable
    cout << "Start MCMC sampling." << endl;   

    double **logmu = new double*[G];
    for(int g = 0; g < G; g++){
        logmu[g] = new double[N];
    }
    
    // set initial value of logmu_{big} = alpha_g + beta_{gw_{bi}} + nu_{bg} + delta_{bi} 
    // cout << "Generate auxiliary logmu." << endl; 
    // auto start_logmu = chrono::system_clock::now();
#pragma omp parallel for
    for(int g = 0; g < G; g++){
        _update_logmu(B, nb, 
                W, alpha[g], beta[g], nu[g], delta,//parameter
                logmu[g]);
    }
    // auto end_logmu = chrono::system_clock::now(); 
    // chrono::duration<double> elapsed_seconds_logmu = end_logmu-start_logmu;
    // cout << "elapsed time of logmu is: " << elapsed_seconds_logmu.count() << "s" << endl;
    
    // divide the overall iterations into several parts with length iter_out
    // to control the RAM occupation
    int out_times = (iter_max - 1) / iter_out + 1;
    int iter_last = iter_max - iter_out * (out_times - 1);
    cout << "Output the posterior sampling to txt file per " << iter_out << " iterations." << endl;
    int IND_UPDATE_PTAU0 = 0;
	if(iter_noupdate > iter_max / 10 * 3){
		iter_noupdate = iter_max / 10 * 3;
	} 

    // auxiliary variable to update gamma
    double gamma_iter, logr;//temp_pres, temp_iter, logr;
    double sigma_zsq = pow(sigma_z, 2.0);

    // auxiliary variable to update p and tau0
    double *p_postdist = new double[2];
    double *tau0_postdist = new double[2];

    // auxiliary variable  to update delta
    double delta_iter;

    // auxiliary variable to update w
    double *proposal_pi = new double[K];
	for (int k = 0; k < K; k++) {
		proposal_pi[k]=1.0 / K;
	}
    int w_proposal, w_current;
    double log_proposal, log_current;

    // auxiliary variable to update prop
    double *count_w = new double[K];

    // auxiliary variable to store posterior sampling
    int q;
    
    // file variable to output the posterior sampling
    ofstream post_File;

    // result/simulation_v0_r1
    string output_dir = wk_dir + "MCMC_sampling";
    //output_dir = output_dir + "_v";
    //output_dir = output_dir + to_string(ver);
    output_dir = output_dir + "_K";
    output_dir = output_dir + to_string(K);
    //output_dir = output_dir + "_r";
    //output_dir = output_dir + to_string(rep);
    //output_dir = output_dir + "_c";
    //output_dir = output_dir + to_string(nc);
    output_dir = output_dir + "/";
    int check = mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (!check) 
        cout << "Directory " << output_dir << " is created." << endl; 
    else { 
        cout << "Directory " << output_dir << " already exists." << endl;
    } 
    string out_file;


    for(int t = 0; t < out_times; t ++){
        int iter_num = iter_out;
        if(t == out_times - 1){
            iter_num = iter_last;
        }

        if(iter_noupdate > iter_num){
            iter_noupdate = iter_noupdate - iter_num;
        }

        cout << "Starting " << t * iter_out + 1 << "-" << t * iter_out + iter_num << " iterations for the " << t + 1 << "-th output." << endl;

        auto start_MCMC = chrono::system_clock::now();
        for(int iter = 0; iter < iter_num; iter ++){

            // cout << "Start " << iter + 1 << " iterations in the " << t + 1 << "-th output." << endl;
        if(All_Drop){
/////////////////////////////////////
            //  1) update z_{big} and x_{big}  //
            // cout << "Update z and x." << endl;
            // auto start_zx = chrono::system_clock::now();
    #pragma omp parallel for
            for(int g = 0; g < G; g++){
                _update_zx(B, nb, 
                        gamma, phi[g], logmu[g],
                        Y[g], MCMC_Rng,
                        X[g], Z[g]);
            }
            
            // auto end_zx = chrono::system_clock::now(); 
            // chrono::duration<double> elapsed_seconds_zx = end_zx-start_zx;
            //if(iter == 0){
            // cout << "elapsed time of updating Z and X is: " << elapsed_seconds_zx.count() << "s" << endl;
            //}
    
            //////////////////////////////////////////
            //  2) update gamma_{b0} and gamma_{b1} // 
            // cout << "Update gamma." << endl;
            // auto start_gamma = chrono::system_clock::now();
            // gamma_{b0}        
            cell_index = 0;
            for (int b = 0; b < B; b++) {
                gamma_iter = MCMC_Rng.rnorm(gamma[b][0], 0.1);

                //prior
                logr = - pow(gamma_iter, 2.0) / 2 / sigma_zsq + pow(gamma[b][0], 2.0) / 2 / sigma_zsq;
            
                for(int i = 0; i < nb[b]; i++ ){
        #pragma omp parallel
            {
                    // auxiliary variable to update gamma
                    double temp_pres, temp_iter;
                    double logr_thread = 0.0;
                    int* X_thread, *Z_thread;
        #pragma omp for
                    for(int g = 0; g < G; g++){
                        X_thread = X[g];
                        Z_thread = Z[g];
                        //numerator
                        temp_iter = gamma_iter + X_thread[cell_index] * gamma[b][1];//prevent temp_iter is extremely large
                        if(temp_iter > 0){
                            logr_thread += gamma_iter * Z_thread[cell_index] - temp_iter  - log( 1 + exp(-temp_iter) );
                        }else{
                            logr_thread += gamma_iter * Z_thread[cell_index] - log( 1 + exp(temp_iter) );
                        }
					
				        //denomerator
                        temp_pres = gamma[b][0] + X_thread[cell_index] * gamma[b][1];
                        if(temp_pres>0){
                            logr_thread += - gamma[b][0] * Z_thread[cell_index] + temp_pres  + log(1 + exp(-temp_pres));
                        }else{
                            logr_thread += - gamma[b][0] * Z_thread[cell_index]  + log(1 + exp(temp_pres) );
                        }
                    }
        #pragma omp critical
                {
                    logr += logr_thread;
                }
            }// end of omp parallel
                    cell_index ++;
                }// end of i

                // cout << " logr = " << logr << endl;
                if (logr > log(MCMC_Rng.runif())) {
			        gamma[b][0] = gamma_iter;
                }
            }// end of b
            
            // gamma_{b1}
            cell_index = 0;
            for (int b = 0; b < B; b++) {
		        //pro posal
                gamma_iter = - MCMC_Rng.rgamma( - 10 * gamma[b][1], 0.1);
                
                //if(gamma_iter < 0){
                //prior
                logr = (gamma_prior[0] - 1) * (log( gamma_iter / gamma[b][1] )) + gamma_prior[1] * (gamma_iter - gamma[b][1]);
            
                //proposal
                //numerator
                logr = logr - lgamma( - 10 * gamma_iter) + ( - 10 * gamma_iter - 1) * log( - gamma[b][1]) - 10 * gamma_iter * log(10) + 10 * gamma[b][1];
                //denomerator
                logr = logr + lgamma( - 10 * gamma[b][1]) - ( - 10 * gamma[b][1] - 1) * log( - gamma_iter) + 10 * gamma[b][1] * log(10) - 10 * gamma_iter;
                
                for (int i = 0; i < nb[b]; i++) {

        #pragma omp parallel
            {
                    // auxiliary variable to update gamma
                    double temp_pres, temp_iter;
                    double logr_thread = 0.0;
                    int* X_thread, *Z_thread;
        #pragma omp for
			        for (int g = 0; g < G; g++) {
                        X_thread = X[g];
                        Z_thread = Z[g];
                        //numerator
                        temp_iter = gamma[b][0] + X_thread[cell_index] * gamma_iter;//prevent temp_iter is extremely large
                        if(temp_iter > 0){
                            logr_thread += X_thread[cell_index] * gamma_iter * Z_thread[cell_index] - temp_iter  - log(1 + exp(-temp_iter));
                        }else{
                            logr_thread += X_thread[cell_index] * gamma_iter * Z_thread[cell_index] - log(1 + exp(temp_iter));
                        }
					
				        //denomerator
                        temp_pres = gamma[b][0] + X_thread[cell_index] * gamma[b][1];
                        if(temp_pres>0){
                            logr_thread += - X_thread[cell_index] * gamma[b][1] * Z_thread[cell_index] + temp_pres + log(1 + exp(-temp_pres));
                        }else{
                            logr_thread += - X_thread[cell_index] * gamma[b][1] * Z_thread[cell_index] + log(1 + exp(temp_pres));
                        }
                        //cout << " g = " << g << ", logr_thread = " << logr_thread << endl;
                    }
        #pragma omp critical
                {
                    logr += logr_thread;
                    
                }

            }
                    cell_index ++;
                }

                // cout << " logr = " << logr << endl;
                if (logr > log(MCMC_Rng.runif())) {
                    gamma[b][1] = gamma_iter;
		        }
            }
            // auto end_gamma = chrono::system_clock::now(); 
            // chrono::duration<double> elapsed_seconds_gamma = end_gamma-start_gamma;
            
            //if(iter == 0){
            // cout << "elapsed time of updating gamma is: " << elapsed_seconds_gamma.count() << "s" << endl;
            //}
        }else{
            /////////////////////////////////////
            //  1) update z_{big} and x_{big}  //
            // cout << "Update z and x." << endl;
            // auto start_zx = chrono::system_clock::now();
    #pragma omp parallel for
            for(int g = 0; g < G; g++){
                _update_zx_optional(B, nb, Drop_ind, 
                        gamma, phi[g], logmu[g],
                        Y[g], MCMC_Rng,
                        X[g], Z[g]);
            }
            
            // auto end_zx = chrono::system_clock::now(); 
            // chrono::duration<double> elapsed_seconds_zx = end_zx-start_zx;
            //if(iter == 0){
            // cout << "elapsed time of updating Z and X is: " << elapsed_seconds_zx.count() << "s" << endl;
            //}
    
            //////////////////////////////////////////
            //  2) update gamma_{b0} and gamma_{b1} // 
            // cout << "Update gamma." << endl;
            // auto start_gamma = chrono::system_clock::now();
            // gamma_{b0}        
            cell_index = 0;
            for (int b = 0; b < B; b++) {
                if(Drop_ind[b]){
                gamma_iter = MCMC_Rng.rnorm(gamma[b][0], 0.1);

                //prior
                logr = - pow(gamma_iter, 2.0) / 2 / sigma_zsq + pow(gamma[b][0], 2.0) / 2 / sigma_zsq;
            
                for(int i = 0; i < nb[b]; i++ ){
        #pragma omp parallel
            {
                    // auxiliary variable to update gamma
                    double temp_pres, temp_iter;
                    double logr_thread = 0.0;
                    int* X_thread, *Z_thread;
        #pragma omp for
                    for(int g = 0; g < G; g++){
                        X_thread = X[g];
                        Z_thread = Z[g];
                        //numerator
                        temp_iter = gamma_iter + X_thread[cell_index] * gamma[b][1];//prevent temp_iter is extremely large
                        if(temp_iter > 0){
                            logr_thread += gamma_iter * Z_thread[cell_index] - temp_iter  - log( 1 + exp(-temp_iter) );
                        }else{
                            logr_thread += gamma_iter * Z_thread[cell_index] - log( 1 + exp(temp_iter) );
                        }
					
				        //denomerator
                        temp_pres = gamma[b][0] + X_thread[cell_index] * gamma[b][1];
                        if(temp_pres>0){
                            logr_thread += - gamma[b][0] * Z_thread[cell_index] + temp_pres  + log(1 + exp(-temp_pres));
                        }else{
                            logr_thread += - gamma[b][0] * Z_thread[cell_index]  + log(1 + exp(temp_pres) );
                        }
                    }
        #pragma omp critical
                {
                    logr += logr_thread;
                }
            }// end of omp parallel
                    cell_index ++;
                }// end of i

                // cout << " logr = " << logr << endl;
                if (logr > log(MCMC_Rng.runif())) {
			        gamma[b][0] = gamma_iter;
		        }else{
                    cell_index = cell_index + nb[b];
                }
                }
            }// end of b
            
            // gamma_{b1}
            cell_index = 0;
            for (int b = 0; b < B; b++) {
                if(Drop_ind[b]){
		        //pro posal
                gamma_iter = - MCMC_Rng.rgamma( - 10 * gamma[b][1], 0.1);
                
                //if(gamma_iter < 0){
                //prior
                logr = (gamma_prior[0] - 1) * (log( gamma_iter / gamma[b][1] )) + gamma_prior[1] * (gamma_iter - gamma[b][1]);
            
                //proposal
                //numerator
                logr = logr - lgamma( - 10 * gamma_iter) + ( - 10 * gamma_iter - 1) * log( - gamma[b][1]) - 10 * gamma_iter * log(10) + 10 * gamma[b][1];
                //denomerator
                logr = logr + lgamma( - 10 * gamma[b][1]) - ( - 10 * gamma[b][1] - 1) * log( - gamma_iter) + 10 * gamma[b][1] * log(10) - 10 * gamma_iter;
                
                for (int i = 0; i < nb[b]; i++) {

        #pragma omp parallel
            {
                    // auxiliary variable to update gamma
                    double temp_pres, temp_iter;
                    double logr_thread = 0.0;
                    int* X_thread, *Z_thread;
        #pragma omp for
			        for (int g = 0; g < G; g++) {
                        X_thread = X[g];
                        Z_thread = Z[g];
                        //numerator
                        temp_iter = gamma[b][0] + X_thread[cell_index] * gamma_iter;//prevent temp_iter is extremely large
                        if(temp_iter > 0){
                            logr_thread += X_thread[cell_index] * gamma_iter * Z_thread[cell_index] - temp_iter  - log(1 + exp(-temp_iter));
                        }else{
                            logr_thread += X_thread[cell_index] * gamma_iter * Z_thread[cell_index] - log(1 + exp(temp_iter));
                        }
					
				        //denomerator
                        temp_pres = gamma[b][0] + X_thread[cell_index] * gamma[b][1];
                        if(temp_pres>0){
                            logr_thread += - X_thread[cell_index] * gamma[b][1] * Z_thread[cell_index] + temp_pres + log(1 + exp(-temp_pres));
                        }else{
                            logr_thread += - X_thread[cell_index] * gamma[b][1] * Z_thread[cell_index] + log(1 + exp(temp_pres));
                        }
                        //cout << " g = " << g << ", logr_thread = " << logr_thread << endl;
                    }
        #pragma omp critical
                {
                    logr += logr_thread;
                    
                }

            }
                    cell_index ++;
                }

                // cout << " logr = " << logr << endl;
                if (logr > log(MCMC_Rng.runif())) {
                    gamma[b][1] = gamma_iter;
		        }
                //}
                }else{
                    cell_index = cell_index + nb[b];
                }

            }
        }
        
            ///////////////////////////////
            //  3) update alpha_g by MH  //
            // cout << "Update alpha." << endl;
            // auto start_alpha = chrono::system_clock::now();
        #pragma omp parallel for
            for(int g = 0; g < G; g++){
                alpha[g] = _update_alpha(B, nb,//dimension
                                    mu_a[g], sigma_a,//prior
                                    W, alpha[g], beta[g], nu[g], delta, phi[g], //parameter
                                    X[g], MCMC_Rng);
            
            }
            // auto end_alpha = chrono::system_clock::now(); 
            // chrono::duration<double> elapsed_seconds_alpha = end_alpha-start_alpha;
            
            //if(iter == 0){
            // cout << "elapsed time of updating alpha is: " << elapsed_seconds_alpha.count() << "s" << endl;
            //}
  
            //////////////////////
            //  4) update L_gk  //
            // cout << "Update L." << endl;
            // auto start_L = chrono::system_clock::now();
        #pragma omp parallel for
            for(int g = 0; g < G; g++){
                _update_l(K,//dimension
                        p, tau0, tau1,//prior
                        beta[g], MCMC_Rng,//parameter
                        L[g]);
            }
            // auto end_L = chrono::system_clock::now(); 
            // chrono::duration<double> elapsed_seconds_L = end_L-start_L;
            
            //if(iter == 0){
            // cout << "elapsed time of updating beta is: " << elapsed_seconds_L.count() << "s" << endl;
            //}
            //////////////////////////////////////
            //  5) update p and 6) update tau0  //
        
            if (IND_UPDATE_PTAU0 == 1) {
                // cout << "Update p and tau0." << endl;
				// auto start_pt = chrono::system_clock::now();
		        int sum_L = 0;
		        double sum_beta = 0.0;
        #pragma omp parallel
            {   
                int sum_L_thread = 0;
                double sum_beta_thread = 0.0;
                int *L_thread;
                double *beta_thread;
        #pragma omp for
                for(int g = 0; g < G; g++){
                    L_thread = L[g];
                    beta_thread = beta[g];
                    for(int k = 1; k < K; k++){
                        sum_L_thread += L_thread[k];
                        if(L_thread[k] == 0){
                            sum_beta_thread += pow(beta_thread[k] , 2.0);
                        }
                    }
                }
        #pragma omp critical
                {
                    sum_L += sum_L_thread;
                    sum_beta += sum_beta_thread;
                }
            }
                
                p_postdist[0] = p_prior[0] + sum_L;
                p_postdist[1] = p_prior[1] + G * (K - 1) - sum_L;
                p = MCMC_Rng.rbeta(p_postdist[0],p_postdist[1]);
                
                
                tau0_postdist[0] = tau0_prior[0] + (G * (K - 1) - sum_L) / 2.0;
                tau0_postdist[1] = tau0_prior[1] + sum_beta / 2.0;
				tau0 = 1.0 / MCMC_Rng.rgamma(tau0_postdist[0],1.0 / tau0_postdist[1]);
                // auto end_pt = chrono::system_clock::now(); 
                // chrono::duration<double> elapsed_seconds_pt = end_pt-start_pt;
            
                //if(iter == 0){
                // cout << "elapsed time of updating p and tau0 is: " << elapsed_seconds_pt.count() << "s" << endl;
                //}
			}
            
            
            /////////////////////
            // 7) update beta  //
            // auto start_beta = chrono::system_clock::now();
        #pragma omp parallel for
            for(int g = 0; g < G; g++){     
                _update_beta(B, nb, K,//dimension 
                    tau0, tau1, L[g],//prior
                    W, alpha[g], nu[g], delta, phi[g], //parameter 	
                    X[g], MCMC_Rng,//latent variable
                    beta[g]);
            }
            // auto end_beta = chrono::system_clock::now(); 
            // chrono::duration<double> elapsed_seconds_beta = end_beta-start_beta;

            //if(iter == 0){
            // cout << "elapsed time of updating beta is: " << elapsed_seconds_beta.count() << "s" << endl;
            //}

            ///////////////////
            // 8) update nu  //
            // auto start_nu = chrono::system_clock::now();
        #pragma omp parallel for
            for(int g = 0; g < G; g++){
                _update_nu(B, nb,
                    mu_c, sigma_c,//prior
                    W, alpha[g], beta[g], delta, phi[g],//parameter
                    X[g], MCMC_Rng, //latent variable
                    nu[g]);
            }
            // auto end_nu = chrono::system_clock::now(); 
            // chrono::duration<double> elapsed_seconds_nu = end_nu-start_nu;
            
            //if(iter == 0){
            // cout << "elapsed time of updating nu is: " << elapsed_seconds_nu.count() << "s" << endl;
            //}
            /////////////////////
            // 9) update delta //
            // cout << "Update delta." << endl;
            // auto start_delta = chrono::system_clock::now();
            cell_index = 0;
            for(int b = 0; b < B; b++){
			    for(int i = 0; i < nb[b]; i ++){
                    if(i > 0){//let ind_n be the cell index
                        delta_iter = MCMC_Rng.rnorm(delta[cell_index], 0.1);
                        logr = 0.0;

                        //prior
                        logr += - pow(delta_iter - mu_d[cell_index] , 2.0) / 2 / pow(sigma_d, 2.0); 
                        logr += pow(delta[cell_index] - mu_d[cell_index] , 2.0) / 2 / pow(sigma_d, 2.0);
        #pragma omp parallel
            {
                        double logr_thread = 0.0;
                        int *X_thread;
                        double *beta_thread;
                        double *nu_thread;
                        double *phi_thread;
        #pragma omp for
                        for(int g = 0; g < G;g ++){
                            X_thread = X[g];
                            beta_thread = beta[g];
                            nu_thread = nu[g];
                            phi_thread = phi[g];
                            int k = W[cell_index];
                            //numerator
					        logr_thread += delta_iter * X_thread[cell_index];
                            logr_thread += - (phi_thread[b] + X_thread[cell_index]) * log(phi_thread[b] + exp(alpha[g] + beta_thread[k] + nu_thread[b] + delta_iter));
					        //denomerator
					        logr_thread += - delta[cell_index] * X_thread[cell_index];
                            logr_thread += (phi_thread[b] + X_thread[cell_index]) * log(phi_thread[b] + exp(alpha[g] + beta_thread[k] + nu_thread[b] + delta[cell_index]));

                        }
        #pragma omp critical
                {
                        logr += logr_thread;
                }
            }
                        if (logr > log(MCMC_Rng.runif())) {
				            delta[cell_index] = delta_iter;
			            }	
                    }
                    cell_index++;
                }
            }
    
            // auto end_delta = chrono::system_clock::now(); 
            // chrono::duration<double> elapsed_seconds_delta = end_delta-start_delta;
            
            //if(iter == 0){
            // cout << "elapsed time of updating delta is: " << elapsed_seconds_delta.count() << "s" << endl;
            //}
        
        #pragma omp parallel for
            for(int g = 0; g < G; g++){
                _update_logmu(B, nb, 
                    W, alpha[g], beta[g], nu[g], delta,//parameter
                    logmu[g]);
            }
        
            /////////////////////
            // 10) update phi  //
            // cout << "Update phi." << endl;
            // auto start_phi = chrono::system_clock::now();
        #pragma omp parallel for
            for(int g = 0; g < G; g++){
                _update_phi(B, nb,
                    phi_prior,//prior
                    logmu[g],//parameter
                    X[g], MCMC_Rng, //latent variable
                    phi[g]);
            }
        
            // auto end_phi = chrono::system_clock::now(); 
            // chrono::duration<double> elapsed_seconds_phi = end_phi-start_phi;
        
            //if(iter == 0){
            // cout << "elapsed time of updating phi is: " << elapsed_seconds_phi.count() << "s" << endl;
            //}

            ///////////////////
            // 11) update w  //

            // cout << "Update w." << endl;
            // auto start_w = chrono::system_clock::now();
            cell_index = 0;
            for(int b = 0; b < B; b++){
                for(int i = 0; i < nb[b]; i++){
                    //get a proposal
                    // print proposal_pi
                    // cout << "Take a look the proposal_pi: ";
                    // for(int k = 0; k < K; k++){
                    //     cout << proposal_pi[k] << " ";
                    // }
                    // cout << endl;

				    w_proposal = rand_cate(proposal_pi, MCMC_Rng);
				    w_current = W[cell_index];

					if(w_proposal != w_current){
                    

				    log_proposal = log(prop[b][w_proposal]);
				    log_current = log(prop[b][w_current]);
				    //calculate the posterior ratio in log scale
            #pragma omp parallel
            {
                        //double logr_thread = 0.0;
                        double log_proposal_thread, log_current_thread;
                        log_proposal_thread = 0.0;
                        log_current_thread = 0.0;
                        double temp_logmu;
                        int *X_thread;
                        double *beta_thread;
                        double *nu_thread;
                        double *phi_thread;
            #pragma omp for
				    for (int g = 0; g < G; g++) {
					    X_thread = X[g];
                        beta_thread = beta[g];
                        nu_thread = nu[g];
                        phi_thread = phi[g];

					    temp_logmu = alpha[g] + beta_thread[w_proposal] + nu_thread[b] + delta[cell_index];
					    log_proposal_thread += beta_thread[w_proposal] * X_thread[cell_index];
					    log_proposal_thread += - (phi_thread[b] + X_thread[cell_index]) * log(phi_thread[b] + exp(temp_logmu));

					    //ind_beta = j + w_current * _G;
                        temp_logmu = alpha[g] + beta_thread[w_current] + nu_thread[b] + delta[cell_index];
					    log_current_thread += beta_thread[w_current] * X_thread[cell_index];
					    log_current_thread += - (phi_thread[b] + X_thread[cell_index]) * log(phi_thread[b] + exp(temp_logmu));
					
				    }
            #pragma omp critical
                {
                        log_proposal += log_proposal_thread;
                        log_current += log_current_thread;
                }
            }
				    logr = log_proposal - log_current;
				    if (logr > log(MCMC_Rng.runif())) {
					    W[cell_index] = w_proposal;
				    }
					}
                    //if(cell_index == 11){
                    //    cout << "The cell type index of cell " << cell_index << " is " << W[cell_index] << endl;
                    //    cout << "w_proposal = " << w_proposal << endl;
                    //    cout << "w_current = " << w_current << endl;
                    //    cout << "log_proposal = " << log_proposal << endl;
                    //    cout << "log_current = " << log_current << endl;
                    //    
                    //}
                    cell_index ++;
                }
            }
            // auto end_w = chrono::system_clock::now(); 
            // chrono::duration<double> elapsed_seconds_w = end_w-start_w;
            
            //if(iter == 0){
            // cout << "elapsed time of updating w is: " << elapsed_seconds_w.count() << "s" << endl;
            //}
        
        #pragma omp parallel for
            for(int g = 0; g < G; g++){
                _update_logmu(B, nb, 
                    W, alpha[g], beta[g], nu[g], delta,//parameter
                    logmu[g]);
            }

        
            ////////////////////////
            //  12) update pi_bk  //
            // cout << "Update pi." << endl;
            // auto start_pi = chrono::system_clock::now();
            cell_index = 0;
            for(int b = 0; b < B; b++){
			    for(int k = 0; k < K; k++){
                    count_w[k] = xi;
                }
			    for (int i = 0; i < nb[b]; i++) {
				    count_w[W[cell_index]] = count_w[W[cell_index]] + 1.0;
				    cell_index ++;
			    }
			    rand_Dir(count_w, K, MCMC_Rng, prop[b]);
		    }

            if(iter == iter_noupdate - 1){
                IND_UPDATE_PTAU0 = 1;
            }
            // auto end_pi = chrono::system_clock::now(); 
            // chrono::duration<double> elapsed_seconds_pi = end_pi-start_pi;
        
            //if(iter == 0){
            // cout << "elapsed time of updating pi is: " << elapsed_seconds_pi.count() << "s" << endl;
            //}

            /////////////////////////////////////////
            //  13) Record the posterior sampling  //
            // cout << "Record posterior sampling." << endl;
            // auto start_record = chrono::system_clock::now();
        #pragma omp parallel for
            for(int g = 0; g < G; g ++){
                alpha_post[iter][g] = alpha[g];
            }

            
        #pragma omp parallel
            {
                int q_thread;
        #pragma omp for
            for(int g = 0; g < G; g ++){
                q_thread = g * K;
                for(int k = 0; k < K; k++){
                    beta_post[iter][q_thread] = beta[g][k];
                    q_thread++;
                }
            }
            }

        #pragma omp parallel
            {
                int q_thread;
        #pragma omp for
            for(int g = 0; g < G; g ++){
                q_thread = g * B;
                for(int b = 0; b < B; b++){
                    nu_post[iter][q_thread] = nu[g][b];
                    q_thread++;
                }
            }
            }

            for(int i = 0; i < N; i ++){
                delta_post[iter][i] = delta[i];
            }
            
            q = 0;
            for(int b = 0; b < B; b ++){
                gamma_post[iter][q] = gamma[b][0];
                q++;
                gamma_post[iter][q] = gamma[b][1];
                q++;
            }

        #pragma omp parallel
            {
                int q_thread;
        #pragma omp for
            for(int g = 0; g < G; g ++){
                q_thread = g * B;
                for(int b = 0; b < B; b++){
                    phi_post[iter][q_thread] = phi[g][b];
                    q_thread++;
                }
            }
            }

            q = 0;
            for(int b = 0; b < B; b ++){
                for(int k = 0; k < K; k++){
                    pi_post[iter][q] = prop[b][k];
                    q++;
                }
            }

            for(int i = 0; i < N; i ++){
                w_post[iter][i] = W[i];
            }

            p_post[iter] = p;
            tau0_post[iter] = tau0;

        #pragma omp parallel
            {
                int q_thread;
        #pragma omp for
            for(int g = 0; g < G; g ++){
                q_thread = g * K;
                for(int k = 0; k < K; k++){
                    l_post[iter][q_thread] = L[g][k];
                    q_thread++;
                }
            }
            }

            // auto end_record = chrono::system_clock::now(); 
            // chrono::duration<double> elapsed_seconds_record = end_record-start_record;
        
            //if(iter == 0){
            // cout << "elapsed time of recording posterior sampling is: " << elapsed_seconds_record.count() << "s" << endl;
            //}

        }

        auto end_MCMC = chrono::system_clock::now(); 
        chrono::duration<double> elapsed_seconds_MCMC = end_MCMC-start_MCMC;
        cout << "elapsed time of "<< iter_num <<" iterations of MCMC sampling for the "<< t+1 <<"-th output is: " << elapsed_seconds_MCMC.count() << "s" << endl;
    
        ///////////////////////////////////
        // output the posterior sampling //
        cout << "Writing posterior sampling into the directory " << output_dir << endl;

        out_file = output_dir + "alpha_post.txt";
        // cout << "Writing alpha_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(int g = 0; g < G; g++){
                post_File<< alpha_post[iter][g];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "beta_post.txt";
        // cout << "Writing beta_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(q = 0; q < G * K;q ++){
                post_File<< beta_post[iter][q];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();
        
        out_file = output_dir + "nu_post.txt";
        // cout << "Writing nu_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(q = 0; q < G * B; q ++){
                post_File<< nu_post[iter][q];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "delta_post.txt";
        // cout << "Writing delta_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(int i = 0; i < N; i ++){
                post_File<< delta_post[iter][i];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "gamma_post.txt";
        // cout << "Writing gamma_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(q = 0; q < B * 2; q ++){
                post_File<< gamma_post[iter][q];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "phi_post.txt";
        // cout << "Writing phi_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(q = 0; q < G * B; q ++){
                post_File<< phi_post[iter][q];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "pi_post.txt";
        // cout << "Writing pi_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(q = 0; q < B * K; q ++){
                post_File<< pi_post[iter][q];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "w_post.txt";
        // cout << "Writing w_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(int i = 0; i < N; i ++){
                post_File << w_post[iter][i];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "p_post.txt";
        // cout << "Writing p_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            post_File<< p_post[iter];
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "tau0_post.txt";
        // cout << "Writing tau0_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            post_File<< tau0_post[iter];
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "l_post.txt";
        // cout << "Writing w_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(q = 0; q < G * K; q ++){
                post_File<< l_post[iter][q];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();
        
        // output the imputed true read count
        if(t == out_times - 1){

            out_file = output_dir + "x_imputed.txt";
            // cout << "Writing x_imputed into " << out_file << endl;
            post_File.open(out_file.c_str(), ios::out | ios::app);
            for(int g = 0; g < G; g++){
                for(int i = 0; i < N; i ++){
                    post_File<< X[g][i];
                    post_File << " ";
                }
                post_File << endl;
            }
            post_File.close();
        }

    }

    auto end_overall = chrono::system_clock::now(); 
    chrono::duration<double> elapsed_seconds_overall = end_overall-start_overall;
    cout << "elapsed time of the overall algorithm is: " << elapsed_seconds_overall.count() << "s" << endl;
    
    //free the memory
    delete [] count_w;
    delete [] proposal_pi;
 
    delete [] p_postdist;
    delete [] tau0_postdist;

    for(int iter = 0; iter < iter_out; iter ++){
        delete [] alpha_post[iter];
    }
    delete [] alpha_post;

    for(int iter = 0; iter < iter_out; iter ++){
        delete [] beta_post[iter];
    }
    delete [] beta_post;

    for(int iter = 0; iter < iter_out; iter ++){
        delete [] nu_post[iter];
    }
    delete [] nu_post;

    for(int iter = 0; iter < iter_out; iter++){
        delete [] delta_post[iter];
    }
    delete [] delta_post;

    for(int iter = 0; iter < iter_out; iter++){
        delete [] gamma_post[iter];
    }
    delete [] gamma_post;

    for(int iter = 0; iter < iter_out; iter++){
        delete [] phi_post[iter];
    }
    delete [] phi_post; 

    for(int iter = 0; iter < iter_out; iter++){
        delete [] pi_post[iter];
    }
    delete [] pi_post;

    for(int iter = 0; iter < iter_out; iter++){
        delete w_post[iter];
    }
    delete [] w_post;
     
    delete [] p_post;
    delete [] tau0_post;

    for(int iter = 0; iter < iter_out; iter++){
        delete [] l_post[iter];
    }
    delete [] l_post;


    for(int g = 0; g < G; g++){
        delete [] X[g];
    }
    delete [] X;
    
    // dropout indicator
    for(int g = 0; g < G; g++){
        delete [] Z[g];
    }
    delete [] Z;

    // cell type indicator
    delete [] W;

    // intrinsic gene indicators
    for(int g = 0; g < G; g++){
        delete [] L[g];
    }
    delete [] L;
    
    // allocate memory for parameters
    // cell type proportion
    for(int b = 0; b < B; b++){
        delete [] prop[b];
    }
    delete [] prop;

    // baseline effects
    delete [] alpha;
	delete [] mu_a;

    // cell-type specific effects
    for(int g = 0; g < G; g++){
        delete [] beta[g];
    }
    delete [] beta;

    // batch effects
    for(int g = 0; g < G; g++){
        delete [] nu[g];
    }
    delete [] nu;
	delete [] mu_c;

    // cell-specific effect
    delete [] delta;
	delete [] mu_d;

    // over-dispersion parameters
    for(int g = 0; g < G; g++){
        delete [] phi[g];
    }
    delete [] phi;

    for(int b = 0; b < B; b++){
        delete [] gamma[b];
    } 
    delete [] gamma;


    for(int g = 0; g < G; g++){
        delete [] Y[g];
    }
    delete [] Y;
    delete [] nb;
    delete [] Drop_ind;
/*
	////////////////////////////
    // 4. Posterior inference //
    ////////////////////////////

	////////////////////////////////////////////////
	//  1) load posterior sampling of parameters  //
	ifstream load_File;
	string load_name;
	int burnin = 3000;
	
	// calculate alpha_est
	double *alpha_est = new double[G];
	for(int g = 0; g < G; g++){
		alpha_est[g] = 0.0;
	}
	load_name = output_dir + "alpha_post.txt";
	load_File.open(load_name);
	int iter = 0;
	
	string row_dropped;
	while (iter < burnin){
		getline(load_File, row_dropped);
		iter ++;
	}
	double temp;
	for(; iter < iter_max; iter++){
        for(int g = 0; g < G; g++){
            load_File >> temp;
			alpha_est[g] += temp;
        }
    }
	for(int g = 0; g < G; g++){
		alpha_est[g] = alpha_est[g] / (iter_max - burnin);
    }
	load_File.close();


	cout << "Writing posterior sampling into the directory " << output_dir << endl;

        out_file = output_dir + "alpha_post.txt";
        // cout << "Writing alpha_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(int g = 0; g < G; g++){
                post_File<< alpha_post[iter][g];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "beta_post.txt";
        // cout << "Writing beta_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(q = 0; q < G * K;q ++){
                post_File<< beta_post[iter][q];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();
        
        out_file = output_dir + "nu_post.txt";
        // cout << "Writing nu_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(q = 0; q < G * B; q ++){
                post_File<< nu_post[iter][q];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "delta_post.txt";
        // cout << "Writing delta_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(int i = 0; i < N; i ++){
                post_File<< delta_post[iter][i];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "gamma_post.txt";
        // cout << "Writing gamma_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(q = 0; q < B * 2; q ++){
                post_File<< gamma_post[iter][q];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "phi_post.txt";
        // cout << "Writing phi_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(q = 0; q < G * B; q ++){
                post_File<< phi_post[iter][q];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "pi_post.txt";
        // cout << "Writing pi_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(q = 0; q < B * K; q ++){
                post_File<< pi_post[iter][q];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "w_post.txt";
        // cout << "Writing w_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(int i = 0; i < N; i ++){
                post_File << w_post[iter][i];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "p_post.txt";
        // cout << "Writing p_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            post_File<< p_post[iter];
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "tau0_post.txt";
        // cout << "Writing tau0_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            post_File<< tau0_post[iter];
            post_File << endl;
        }
        post_File.close();

        out_file = output_dir + "l_post.txt";
        // cout << "Writing w_post into " << out_file << endl;
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for(int iter = 0; iter < iter_num; iter++){
            for(q = 0; q < G * K; q ++){
                post_File<< l_post[iter][q];
                post_File << " ";
            }
            post_File << endl;
        }
        post_File.close();
        
        // output the imputed true read count
        if(t == out_times - 1){

            out_file = output_dir + "x_imputed.txt";
            // cout << "Writing x_imputed into " << out_file << endl;
            post_File.open(out_file.c_str(), ios::out | ios::app);
            for(int g = 0; g < G; g++){
                for(int i = 0; i < N; i ++){
                    post_File<< X[g][i];
                    post_File << " ";
                }
                post_File << endl;
            }
            post_File.close();
        }
*/
    return 0;
}