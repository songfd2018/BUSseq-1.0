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
#include <algorithm>    // sort
#include <vector>  //

using namespace std;

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

//Max value in a vector
double vec_max(double* vec,int n){
	double res = vec[0];
	for (int i = 1; i < n; i++) {
		if (res < vec[i]) {
		    res = vec[i];
		}
	}
	return res;
}

// FDR calculation
double fdrDEindicator(double **_PPI, double _kappa, int _G, int _K){
    
    double xi, fdr;
    double sum_xi = 0.0;
    int count_intri = 0;
    if(_kappa > 0){
        for(int g = 0; g < _G; g++){
            for(int k = 1; k < _K; k++){
                xi = 1 - _PPI[g][k];
                if(xi <= _kappa){
                    sum_xi += xi;
                    count_intri ++;
                }
            }
        }
        fdr = sum_xi/count_intri;
    }else{
        fdr = 0.0;
    }
        
    return(fdr);
}

bool descreasing(double i,double j) { return (i>j); }

// calculate the DE posterior probability threshold
double postprob_DE_thr_fun(double **_PPI, double _fdr_threshold, int _G, int _K){
    
    double kappa;
    double postprob_DE_thr = 0.5;
    double fdr = 0.0;
    vector<double> vec_PPI;
    for(int g = 0; g < _G; g++){
        copy(_PPI[g]+1, _PPI[g] + _K, back_inserter(vec_PPI));
    }
    // sort all PPIs decreasingly
	// cout << "The length of vec_PPI is " << vec_PPI.size() << endl;
    sort(vec_PPI.begin(),vec_PPI.end(),descreasing);
	 
    // unique values in the PPI vector
	// cout << "The length of vec_PPI is " << vec_PPI.size() << endl;
    vector<double>::iterator it;
    it = unique (vec_PPI.begin(), vec_PPI.end());                                                           //                ^
    vec_PPI.resize(distance(vec_PPI.begin(),it));
    //cout << "The length of vec_PPI is " << vec_PPI.size() << endl;

    int i = 0;
    while(fdr <= _fdr_threshold & i < vec_PPI.size()){
        kappa = 1 - vec_PPI[i];
        fdr = fdrDEindicator(_PPI, kappa, _G, _K);
        // cout << "kappa = " << kappa << ", fdr = " << fdr << "." << endl;
        i++;
    }
	double PPI_thres;
	if(i < vec_PPI.size()){
		i = i - 2;// the index of last i such that fdr <= _fdr_threshold to control FDR
		PPI_thres = vec_PPI[i];
	}else{
		i = 0; // Even the largest
		PPI_thres = postprob_DE_thr;
	}
    
    // cout << "vec_PPI[i] = " << vec_PPI[i] << endl;
    if(vec_PPI[i] > postprob_DE_thr){
        postprob_DE_thr = vec_PPI[i];
    }
    return postprob_DE_thr;
}

int IG_index(double *_PPI, double _postprob_DE_thr, int _K){

    int IG_indicator = 0;
    for(int k = 1; k < _K; k++){
        if(_PPI[k] >= _postprob_DE_thr){
            IG_indicator = 1;
        }
        // cout << "PPI_" << k << " = " << _PPI[k] << ", ";
    }
    // cout << "IG_indicator = " << IG_indicator << "." << endl;
    
    return IG_indicator;
}


int main(int argc, char** argv) {

    //////////////////////////////////////////////////////////
    //  0. Input the project and version from command line  //
    //////////////////////////////////////////////////////////
    string proj, wk_dir, readcount_dir;
    int ver, K, rep;
    int iter_max, iter_burnin; // the overall iteration number and the number of iterations to print the posterior sampling 
    int seed; // seed for random number generator 
    int nc; // number of codes for parallel
    int opt;

    //while ((opt = getopt(argc, argv, "p:v:N:G:B:K:r:i:s:")) != -1) {
    while ((opt = getopt(argc, argv, ":d::r::p:v:K:i:b:c:")) != -1) {
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
        case 'b':
            iter_burnin = atoi(optarg);
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
        default: /* '?' */
            cout << "Error Usage: " << argv[0] <<" [-p ProjectName] [-v Version] [-K CellTypeNumber] [-r ReplicationIndex] [-i IterationNumber] [-b NumberBurnin] [-c Cores] name" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    // set count data matrix and the number of cells in each batch
    
    string fcount, fdim, output_dir, infer_dir;
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
    
    auto start_overall = chrono::system_clock::now();

    // result/Inference_simulation_v0_K3_r1/
    output_dir = wk_dir; // + "result/";
    
    infer_dir = output_dir + "Inference";
    // infer_dir = infer_dir + proj;
    // infer_dir = infer_dir + "_v";
    // infer_dir = infer_dir + to_string(ver);
    infer_dir = infer_dir + "_K";
    infer_dir = infer_dir + to_string(K);
    // infer_dir = infer_dir + "_r";
    // infer_dir = infer_dir + to_string(rep);
    infer_dir = infer_dir + "/";
    // result/simulation_v0_K3_r1_c4/
    output_dir = output_dir + "MCMC_sampling";
    // output_dir = output_dir + "_v";
    // output_dir = output_dir + to_string(ver);
    output_dir = output_dir + "_K";
    output_dir = output_dir + to_string(K);
    // output_dir = output_dir + "_r";
    // output_dir = output_dir + to_string(rep);
    // output_dir = output_dir + "_c";
    // output_dir = output_dir + to_string(nc);
    output_dir = output_dir + "/";


    int check = mkdir(infer_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (!check) 
        cout << "Directory " << infer_dir << " is created." << endl; 
    else { 
        cout << "Directory " << infer_dir << " already exists." << endl;
    } 
    cout << "The output directory of posterior sampling is: " << infer_dir << endl;

    ////////////////////////
    // 1. Load count data //
    ////////////////////////
    int N, G, B;

    ifstream countFile, dimFile;
    dimFile.open(fdim.c_str());
    if(!dimFile){
        cout << "Unable to open the file of dimensions from " << fdim << endl;
        exit(1); // terminate with error
    }
    dimFile >> N;
    cout << "There are totally " << N << " cells and ";
    dimFile >> G;
    cout << G << " genes in ";
    dimFile >> B;
    cout << B << " batches." << endl;
    cout << "Loading the number of cells in each batch from " << fdim << endl;
    int *nb = new int[B];
    for(int b = 0; b < B; b++){
        dimFile >> nb[b];
        cout << "The " << (b+1) << "-th batch has "<< nb[b] <<" cells." << endl;
    }

    // whether there are dropout events in each batch or not
    bool *Drop_ind = new bool[B];
    int di;
    for(int b = 0; b < B; b++){
        dimFile >> di;
        if(di == 1){
            Drop_ind[b] = true;
            cout << "The " << (b+1) << "-th batch has dropout events." << endl;
        }else if(di == 0){
            Drop_ind[b] = false;
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

	////////////////////////////
    // 2. Posterior inference //
    ////////////////////////////

	////////////////////////////////////////////////
	//  1) load posterior sampling of parameters  //
	ifstream load_File;
    ofstream est_File;
	string load_name, est_name;
    double temp;
    int iter;
    string row_dropped;
	// int burnin = 3000;
	
	// calculate alpha_est
	double *alpha_est = new double[G];
	for(int g = 0; g < G; g++){
		alpha_est[g] = 0.0;
	}
	load_name = output_dir + "alpha_post.txt";
    cout << "Load alpha_post from " << load_name << endl;
	load_File.open(load_name);
	iter = 0;
	
	while (iter < iter_burnin){
		getline(load_File, row_dropped);
		iter ++;
	}
	for(; iter < iter_max; iter++){
        for(int g = 0; g < G; g++){
            load_File >> temp;
			alpha_est[g] += temp;
        }
    }
	for(int g = 0; g < G; g++){
		alpha_est[g] = alpha_est[g] / (iter_max - iter_burnin);
    }
	load_File.close();

    // calculate beta_est
	double **beta_est = new double*[G];
	
    for(int g = 0; g < G; g++){
        beta_est[g] = new double[K];
        for(int k = 0; k < K; k++){
            beta_est[g][k] = 0.0;
        }
    }
	load_name = output_dir + "beta_post.txt";
	cout << "Load beta_post from " << load_name << endl;
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin){
		getline(load_File, row_dropped);
		iter ++;
	}
	for(; iter < iter_max; iter++){
        for(int g = 0; g < G; g++){
            for(int k = 0; k < K; k++){
                load_File >> temp;
			    beta_est[g][k] += temp;
            }
        }
    }
	for(int g = 0; g < G; g++){
        for(int k = 0; k < K; k++){
            beta_est[g][k] = beta_est[g][k] / (iter_max - iter_burnin);
        }
    }
	load_File.close();

    // calculate nu_est
	double **nu_est = new double*[G];
    for(int g = 0; g < G; g++){
        nu_est[g] = new double[B];
        for(int b = 0; b < B; b++){
            nu_est[g][b] = 0.0;
        }
    }
	load_name = output_dir + "nu_post.txt";
	cout << "Load nu_post from " << load_name << endl;
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin){
		getline(load_File, row_dropped);
		iter ++;
	}
	for(; iter < iter_max; iter++){
        for(int g = 0; g < G; g++){
            for(int b = 0; b < B; b++){
                load_File >> temp;
			    nu_est[g][b] += temp;
            }
        }
    }
	for(int g = 0; g < G; g++){
        for(int b = 0; b < B; b++){
            nu_est[g][b] = nu_est[g][b] / (iter_max - iter_burnin);
        }
    }
	load_File.close();

    // calculate delta_est
	double *delta_est = new double[N];
    for(int i = 0; i < N; i++){
        delta_est[i] = 0.0;
    }
	load_name = output_dir + "delta_post.txt";
	cout << "Load delta_post from " << load_name << endl;
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin){
		getline(load_File, row_dropped);
		iter ++;
	}
	for(; iter < iter_max; iter++){
        for(int i = 0; i < N; i++){
            load_File >> temp;
			delta_est[i] += temp;
        }
    }
	for(int i = 0; i < N; i++){
        delta_est[i] = delta_est[i] / (iter_max - iter_burnin);
    }
	load_File.close();

    // calculate gamma_est
	double **gamma_est = new double*[B];
    for(int b = 0; b < B; b++){
        gamma_est[b] = new double[2];
        gamma_est[b][0] = 0.0;
        gamma_est[b][1] = 0.0;
    }
	load_name = output_dir + "gamma_post.txt";
	cout << "Load gamma_post from " << load_name << endl;
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin){
		getline(load_File, row_dropped);
		iter ++;
	}
	for(; iter < iter_max; iter++){
        for(int b = 0; b < B; b++){
            load_File >> temp;
			gamma_est[b][0] += temp;
            load_File >> temp;
			gamma_est[b][1] += temp;
        }
    }
	for(int b = 0; b < B; b++){
        gamma_est[b][0] = gamma_est[b][0] / (iter_max - iter_burnin);
        gamma_est[b][1] = gamma_est[b][1] / (iter_max - iter_burnin);
    }
	load_File.close();

    // calculate phi_est
	double **phi_est = new double*[G];
    for(int g = 0; g < G; g++){
        phi_est[g] = new double[B];
        for(int b = 0; b < B; b++){
            phi_est[g][b] = 0.0;
        }
    }
	load_name = output_dir + "phi_post.txt";
	cout << "Load phi_post from " << load_name << endl;
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin){
		getline(load_File, row_dropped);
		iter ++;
	}
	for(; iter < iter_max; iter++){
        for(int g = 0; g < G; g++){
            for(int b = 0; b < B; b++){
                load_File >> temp;
			    phi_est[g][b] += temp;
            }
        }
    }
	for(int g = 0; g < G; g++){
        for(int b = 0; b < B; b++){
            phi_est[g][b] = phi_est[g][b] / (iter_max - iter_burnin);
        }
    }
	load_File.close();

    // calculate pi_est
	double **pi_est = new double*[B];
    for(int b = 0; b < B; b++){
        pi_est[b] = new double[K];
        for(int k = 0; k < K; k++){
            pi_est[b][k] = 0.0;
        }
    }
	load_name = output_dir + "pi_post.txt";
	cout << "Load pi_post from " << load_name << endl;
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin){
		getline(load_File, row_dropped);
		iter ++;
	}
	for(; iter < iter_max; iter++){
        for(int b = 0; b < B; b++){
            for(int k = 0; k < K; k++){
                load_File >> temp;
			    pi_est[b][k] += temp;
            }
        }
    }
	for(int b = 0; b < B; b++){
        for(int k = 0; k < K; k++){
            pi_est[b][k] = pi_est[b][k] / (iter_max - iter_burnin);
        }
    }
	load_File.close();

    // calculate w_est
    int w_temp;
	int *w_est = new int[N];
    for(int i = 0; i < N; i++){
        w_est[i] = 0;
    }
    int **count_w = new int*[N]; // count the frequency of cell types after burnin
    for(int i = 0; i < N; i++){
        count_w[i] = new int[K];
        for(int k = 0; k < K; k ++){
            count_w[i][k] = 0;
        }
    }
	load_name = output_dir + "w_post.txt";
	cout << "Load w_post from " << load_name << endl;
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin){
		getline(load_File, row_dropped);
		iter ++;
	}
	for(; iter < iter_max; iter++){
        for(int i = 0; i < N; i++){
            load_File >> w_temp;
            count_w[i][w_temp] ++;
        }
    }

    for(int i = 0; i < N; i++){
        for(int k = 1; k < K; k++){
            if(count_w[i][k] > count_w[i][w_est[i]]){
                w_est[i] = k; // obtain the mode from burnin to the last iteration
            }
        }
    }
	load_File.close();

    // calculate p_est
	double p_est = 0.0;
	load_name = output_dir + "p_post.txt";
	cout << "Load p_post from " << load_name << endl;
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin){
		getline(load_File, row_dropped);
		iter ++;
	}
	for(; iter < iter_max; iter++){
        load_File >> temp;
		p_est += temp;

    }
    p_est = p_est / (iter_max - iter_burnin);
	load_File.close();

    // calculate tau0_est
	double tau0_est = 0.0;
	load_name = output_dir + "tau0_post.txt";
	cout << "Load tau0_post from " << load_name << endl;
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin){
		getline(load_File, row_dropped);
		iter ++;
	}
	for(; iter < iter_max; iter++){
        load_File >> temp;
		tau0_est += temp;

    }
    tau0_est = tau0_est / (iter_max - iter_burnin);
	load_File.close();

    // calculate PPI_est
	double **PPI_est = new double*[G];
	
    for(int g = 0; g < G; g++){
        PPI_est[g] = new double[K];
        for(int k = 0; k < K; k++){
            PPI_est[g][k] = 0.0;
        }
    }
	load_name = output_dir + "l_post.txt";
	cout << "Load l_post from " << load_name << endl;
	load_File.open(load_name);
	iter = 0;
	int temp_int;
	
	while (iter < iter_burnin){
		getline(load_File, row_dropped);
		iter ++;
	}
	for(; iter < iter_max; iter++){
        for(int g = 0; g < G; g++){
            for(int k = 0; k < K; k++){
                load_File >> temp_int;
			    PPI_est[g][k] += temp_int;
            }
        }
		// cout << "Finish loading " << iter << "-th iteration." << endl;
    }
	for(int g = 0; g < G; g++){
        for(int k = 0; k < K; k++){
            PPI_est[g][k] = PPI_est[g][k] / (iter_max - iter_burnin);
        }
    }
	load_File.close();

    // Extract intrinsic genes
    double fdr_threshold = 0.05, postprob_DE_threshold;
    int count_gene = 0;
    int *D = new int[G];// intrinsic gene indicator
	cout << "The threshold of PPI to identify intrinsic genes is " << postprob_DE_threshold << " to control FDR at level "<< fdr_threshold << "." << endl;
    postprob_DE_threshold = postprob_DE_thr_fun(PPI_est,fdr_threshold, G, K);
    
    for(int g = 0; g < G; g++){
        D[g] = IG_index(PPI_est[g], postprob_DE_threshold, K);
        count_gene += D[g];
    }
    cout << "There are " << count_gene << " identified intrisic genes." << endl;

    // load x_imputed
	int **x_imputed = new int*[G];
    for(int g = 0; g < G; g++){
        x_imputed[g] = new int[N];
    }
	load_name = output_dir + "x_imputed.txt";
	cout << "Load x_imputed from " << load_name << endl;
	load_File.open(load_name);
	for(int g = 0; g < G; g++){
        for(int i = 0; i < N; i++){
            load_File >> x_imputed[g][i];
        }
    }
	load_File.close();

    // Calculate likelihood and BIC
    
    omprng Loglike_Rng;
    Loglike_Rng.fixedSeed(12345);// for Monte carlo estimation of E(exp(gamma_0+gamma_1x)/(1+exp(gamma_0+gamma_1x)))
    int cell_index;
    double loglike_obs;
    // auxiliary variables
    double *lpy = new double[K]; // log(pi_bk prod_{g=1}^G Pr(Y_big = y_big | Theta ))
    double lpy_max, sum_lpy; // lr0_temp, sum_lr0, lpy_max, sum_lpy;
    // Set the number of cores for parallel
    omp_set_num_threads(nc);

    // Calculate BIC 
    auto start_BIC2 = chrono::system_clock::now();
    double BIC;
    loglike_obs = 0.0;
    cell_index = 0;
    for (int b = 0; b < B; b++) {
        if(Drop_ind[b]){
for (int i = 0; i < nb[b]; i++) {
			for (int k = 0; k < K; k++) {
				lpy[k] = log(pi_est[b][k]);
			}
			for (int k = 0; k < K; k++) {
			// auto start_bik = chrono::system_clock::now();
        # pragma omp parallel
            {
                int read, x_max;
                double logmubikg, pbgk, logp, log1mp, lr0_temp, sum_lr0, lpy_thread;
                lpy_thread = 0.0;
        # pragma omp for
				for (int g = 0; g < G; g++) {
						read = Y[g][cell_index];
                        logmubikg = alpha_est[g] + beta_est[g][k] + nu_est[g][b] + delta_est[cell_index];
					    pbgk = exp(logmubikg) / (exp(logmubikg) + phi_est[g][b]);
                        if(pbgk < pow(0.1, 100)){
                            logp = -100;
                            log1mp = log(1-pbgk);
                        }else if(1 - pbgk < pow(0.1, 100)){
                            logp = log(pbgk);
                            log1mp = -100;
                        }else{
                            logp = log(pbgk);
                            log1mp = log(1 - pbgk);
                        }

						if (read > 0) {
							lpy_thread += - log(1 + exp(gamma_est[b][0] + gamma_est[b][1] * read));
							lpy_thread += lgamma(phi_est[g][b] + read) - lgamma(read + 1) - lgamma(phi_est[g][b]);
                            lpy_thread += read * logp + phi_est[g][b] * log1mp;
						}else {
							x_max = (int)3 * exp(logmubikg);
							lr0_temp = phi_est[g][b] * log1mp; //x=0
							sum_lr0 = lr0_temp;

							for (int x = 1; x < x_max; x++) {
								lr0_temp = gamma_est[b][0] + gamma_est[b][1] * x - log(1 + exp(gamma_est[b][0] + gamma_est[b][1] * x));
								lr0_temp += lgamma(phi_est[g][b] + x) - lgamma(x + 1) - lgamma(phi_est[g][b]);
                                lr0_temp += x * logp + phi_est[g][b] * log1mp;
								if (lr0_temp > sum_lr0) {
									sum_lr0 = lr0_temp + log(1 + exp(sum_lr0 - lr0_temp));
								}else {
									sum_lr0 = sum_lr0 + log(1 + exp(lr0_temp - sum_lr0));
								}
							}
							lpy_thread += sum_lr0;
							//Rprintf("y %d %d %d = 0 and sum_lr0 is %f if the cell belongs to %d-th subtype.\n",b, i ,j ,sum_lr0, k);
						}
					}
        # pragma omp critical
                {
                    lpy[k] += lpy_thread;
                }
            }
					//Rprintf("lpy %d = %f for the %d-th cell.\n",k+1,lpy[k],index_n);
					// auto end_bik = chrono::system_clock::now(); 
					// chrono::duration<double> elapsed_seconds_bik = end_bik-start_bik;
					// cout << "elapsed time of calculating the likelihood of Y_{" << b << "," << i << "} regarded in cell type "<< k << " is: " << elapsed_seconds_bik.count() << "s" << endl;
				}// end of k 
				lpy_max = vec_max(lpy, K);
				sum_lpy = 0.0;
				for (int k = 0; k < K; k++) {
					sum_lpy = sum_lpy + exp(lpy[k] - lpy_max);
					//Rprintf("logproby[%d]=%f",k,lpy[k]);
				}
				loglike_obs += lpy_max + log(sum_lpy);
				cell_index ++;
				// printf("Finish the %d-th cell, lpy_max= %f, sum_lpy = %f, loglike = %f.\n", cell_index, lpy_max, sum_lpy, loglike_obs);
			}// end of i
        }else{
            for (int i = 0; i < nb[b]; i++) {
			for (int k = 0; k < K; k++) {
				lpy[k] = log(pi_est[b][k]);
			}
			for (int k = 0; k < K; k++) {
			// auto start_bik = chrono::system_clock::now();
        # pragma omp parallel
            {
                int read, x_max;
                double logmubikg, pbgk, logp, log1mp, lr0_temp, sum_lr0, lpy_thread;
                lpy_thread = 0.0;
        # pragma omp for
				for (int g = 0; g < G; g++) {
						read = Y[g][cell_index];
                        logmubikg = alpha_est[g] + beta_est[g][k] + nu_est[g][b] + delta_est[cell_index];
					    pbgk = exp(logmubikg) / (exp(logmubikg) + phi_est[g][b]);
                        if(pbgk < pow(0.1, 100)){
                            logp = -100;
                            log1mp = log(1-pbgk);
                        }else if(1 - pbgk < pow(0.1, 100)){
                            logp = log(pbgk);
                            log1mp = -100;
                        }else{
                            logp = log(pbgk);
                            log1mp = log(1 - pbgk);
                        }
                        
                        //lpy_thread += - log(1 + exp(gamma_est[b][0] + gamma_est[b][1] * read));
						lpy_thread += lgamma(phi_est[g][b] + read) - lgamma(read + 1) - lgamma(phi_est[g][b]);
                        lpy_thread += read * logp + phi_est[g][b] * log1mp;
					}
        # pragma omp critical
                {
                    lpy[k] += lpy_thread;
                }
            }
					//Rprintf("lpy %d = %f for the %d-th cell.\n",k+1,lpy[k],index_n);
					// auto end_bik = chrono::system_clock::now(); 
					// chrono::duration<double> elapsed_seconds_bik = end_bik-start_bik;
					// cout << "elapsed time of calculating the likelihood of Y_{" << b << "," << i << "} regarded in cell type "<< k << " is: " << elapsed_seconds_bik.count() << "s" << endl;
				}// end of k 
				lpy_max = vec_max(lpy, K);
				sum_lpy = 0.0;
				for (int k = 0; k < K; k++) {
					sum_lpy = sum_lpy + exp(lpy[k] - lpy_max);
					//Rprintf("logproby[%d]=%f",k,lpy[k]);
				}
				loglike_obs += lpy_max + log(sum_lpy);
				cell_index ++;
				// printf("Finish the %d-th cell, lpy_max= %f, sum_lpy = %f, loglike = %f.\n", cell_index, lpy_max, sum_lpy, loglike_obs);
			}// end of i
        }
		}// end of b

    int NumBatchDrop = 0;
    for(int b = 0; b < B; b++){
        NumBatchDrop = NumBatchDrop + Drop_ind[b];
    }
    BIC = - 2.0 * loglike_obs + log(G * N) * ((B + G) * K+ 2 * NumBatchDrop + G * (B * 2 - 1) + N - B);
	// all parameters of interest contain pi_{bk}, gamma_{b0(1)}, alpha_g, beta_{gk}, nu_{bg}, delta_{bi}, phi_{bg} 
    auto end_BIC2 = chrono::system_clock::now(); 
    chrono::duration<double> elapsed_seconds_BIC2 = end_BIC2-start_BIC2;
    cout << "elapsed time of calculating previous BIC is: " << elapsed_seconds_BIC2.count() << "s" << endl;
    

    cout << "Writing posterior mean of alpha into the directory " << infer_dir << endl;
    est_name = infer_dir + "alpha_est.txt";
    est_File.open(est_name.c_str());
    for(int g = 0; g < G; g++){
        est_File << alpha_est[g];
        est_File << " ";
    }
    est_File << endl;
    est_File.close();

    cout << "Writing posterior mean of beta into the directory " << infer_dir << endl;
    est_name = infer_dir + "beta_est.txt";
    est_File.open(est_name.c_str());
    for(int g = 0; g < G; g++){
        for(int k = 0; k < K; k++){
            est_File << beta_est[g][k];
            est_File << " ";
        }
        est_File << endl;
    }
    est_File << endl;
    est_File.close();

    cout << "Writing posterior mean of nu into the directory " << infer_dir << endl;
    est_name = infer_dir + "nu_est.txt";
    est_File.open(est_name.c_str());
    for(int g = 0; g < G; g++){
        for(int b = 0; b < B; b++){
            est_File << nu_est[g][b];
            est_File << " ";
        }
        est_File << endl;
    }
    est_File << endl;
    est_File.close();

    cout << "Writing posterior mean of delta into the directory " << infer_dir << endl;
    est_name = infer_dir + "delta_est.txt";
    est_File.open(est_name.c_str());
    for(int i = 0; i < N; i++){
        est_File << delta_est[i];
        est_File << " ";
    }
    est_File << endl;
    est_File.close();

    cout << "Writing posterior mean of gamma into the directory " << infer_dir << endl;
    est_name = infer_dir + "gamma_est.txt";
    est_File.open(est_name.c_str());
    for(int b = 0; b < B; b++){
        est_File << gamma_est[b][0];
        est_File << " ";
        est_File << gamma_est[b][1];
        est_File << " ";
        est_File << endl;
    }
    est_File << endl;
    est_File.close();

    cout << "Writing posterior mean of phi into the directory " << infer_dir << endl;
    est_name = infer_dir + "phi_est.txt";
    est_File.open(est_name.c_str());
    for(int g = 0; g < G; g++){
        for(int b = 0; b < B; b++){
            est_File << phi_est[g][b];
            est_File << " ";
        }
        est_File << endl;
    }
    est_File << endl;
    est_File.close();

    cout << "Writing posterior mean of pi into the directory " << infer_dir << endl;
    est_name = infer_dir + "pi_est.txt";
    est_File.open(est_name.c_str());
    for(int b = 0; b < B; b++){
        for(int k = 0; k < K; k++){
            est_File << pi_est[b][k];
            est_File << " ";
        }
        est_File << endl;
    }
    est_File << endl;
    est_File.close();

    cout << "Writing posterior mode of w into the directory " << infer_dir << endl;
    est_name = infer_dir + "w_est.txt";
    est_File.open(est_name.c_str());
    for(int i = 0; i < N; i++){
        est_File << w_est[i];
        est_File << " ";
    }
    est_File << endl;
    est_File.close();

    cout << "Writing posterior mean of p into the directory " << infer_dir << endl;
    est_name = infer_dir + "p_est.txt";
    est_File.open(est_name.c_str());
    est_File << p_est;    
    est_File << endl;
    est_File.close();

    cout << "Writing posterior mean of tau0 into the directory " << infer_dir << endl;
    est_name = infer_dir + "tau0_est.txt";
    est_File.open(est_name.c_str());
    est_File << tau0_est;    
    est_File << endl;
    est_File.close();

    cout << "Writing posterior mean of PPI into the directory " << infer_dir << endl;
    est_name = infer_dir + "PPI_est.txt";
    est_File.open(est_name.c_str());
    for(int g = 0; g < G; g++){
        for(int k = 0; k < K; k++){
            est_File << PPI_est[g][k];
            est_File << " ";
        }
        est_File << endl;
    }
    est_File << endl;
    est_File.close();

    cout << "Writing intrinsic gene indicators into the directory " << infer_dir << endl;
    est_name = infer_dir + "IG_est.txt";
    est_File.open(est_name.c_str());
    for(int g = 0; g < G; g++){
        est_File << D[g];
        est_File << endl;
    }
    est_File << endl;
    est_File.close();

    cout << "Writing posterior mean of imputed X into the directory " << infer_dir << endl;
    est_name = infer_dir + "x_imputed.txt";
    est_File.open(est_name.c_str());
    for(int g = 0; g < G; g++){
        for(int i = 0; i < N; i++){
            est_File << x_imputed[g][i];
            est_File << " ";
        }
        est_File << endl;
    }
    est_File << endl;
    est_File.close();

    cout << "Writing BIC into the directory " << infer_dir << endl;
    est_name = infer_dir + "BIC.txt";
    est_File.open(est_name.c_str());
    est_File << BIC;    
    est_File << endl; 
    est_File.close();


    delete [] alpha_est;
    for(int g = 0; g < G; g ++){
        delete [] beta_est[g];
    }
    delete [] beta_est;
    for(int g = 0; g < G; g ++){
        delete [] nu_est[g];
    }
    delete [] nu_est;
    delete [] delta_est;
    for(int b = 0; b < B; b ++){
        delete [] gamma_est[b];
    }
    delete [] gamma_est;
    for(int g = 0; g < G; g ++){
        delete [] phi_est[g];
    }
    delete [] phi_est;
    for(int b = 0; b < B; b ++){
        delete [] pi_est[b];
    }
    delete [] pi_est;
    delete [] w_est;
    for(int i = 0; i < N; i ++){
        delete [] count_w[i];
    }
    delete [] count_w;
    for(int g = 0; g < G; g ++){
        delete [] PPI_est[g];
    }
    delete [] PPI_est;
    delete [] D;
    
    delete [] lpy;

    for(int g = 0; g < G; g++){
        delete [] Y[g];
    }
    delete [] Y;
    delete [] nb;
    delete [] Drop_ind;
 
    return 0;
}