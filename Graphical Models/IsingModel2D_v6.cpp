// IsingModel2D.cpp : This subsystem emules a 2D Ising Model with boundary conditions
// 
// Included Simulated Annealing to Generate Sample Paths (Jan 10th 2021)
//
// Author - Luis Alvaro Correia
// 
// Date - Jan 11th 2021

#include <iostream>
#include<stdio.h>
#include<stdlib.h>
#include <cmath>
#include <ctime>
#include <fstream>
#include <string>
#include <Windows.h>
//#include "matplotlibcpp.h"

//namespace plt = matplotlibcpp;

const int MAX_SIZE_LATTICE = 30;

struct ModelParam
{
    double beta;
    long K;
    int L;
    double J;
} Param;

long    L = 0L,                                         // Default size of the 2D Ising Model
        N = 0L,                                         // Size of the lattice
        K = 0L,                                         // Number of Samples per Gibbs simulation;
        burnin = 0L;                                    // Number of samples to discard  
double  beta = 0.0;                                     // Initial value for Beta
double  lBeta[5] = { 0.01, 0.06, 0.15, 0.25, 0.75 };    // List of possible values of Beta
double  J = 1.0;                                        // Value of Coupling Energy (default = 1.0)
int     NoSim = 4;                                      // Number of Monte Carlo Simulations per beta(= no.sample paths)
int     job = 0;

ModelParam IsingParam { 0.01, 10000L, 5, 1.0 };         // Default Parameters

char Version[] = "_v6-CPP_";

clock_t start_time;

int durationP = 800;     // milliseconds
int freqP = 1000;        // Hz
int durationJob = 1000;  // milliseconds
int freqJob = 1500;      // Hz

//  OPTIMIZED - Returns the nearest neighbors of sigma(i, j).For nearest neighbor interactions,
// these are sigma(i + 1, j), sigma(i - 1, j), sigma(i, j + 1), and sigma(i, j - 1).

int * get_neighbors(int *sigma, int i, int j)
{
    static int nb[4]{};
    
    nb[0] = *(sigma + ((i - 1) < 0 ? (L - 1) : (i - 1)) * L + j); // Up
    nb[1] = *(sigma + i * L + ((j + 1) > (L - 1) ? 0 : (j + 1))); // Right
    nb[2] = *(sigma + ((i + 1) > (L - 1) ? 0 : (i + 1)) * L + j); // Down
    nb[3] = *(sigma + i * L + ((j - 1) < 0 ? (L - 1) : (j - 1))); // Left

    return &nb[0];
}

double Hamiltonian(int *x)
{
    double S = 0.0;
    int* neighbors = NULL;
    
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++) {
            neighbors = get_neighbors(x, i, j);
            S -= J * ((double)*(x + 2 * i + j) * (double)(*neighbors) + 
                (double)*(x + 2 * i + j) * ((double)*(neighbors + 1)) +
                (double)*(x + 2 * i + j) * ((double)*(neighbors + 2)) + 
                (double)*(x + 2 * i + j) * ((double)*(neighbors + 3)));
            }

    return S;

}

//    Returns p = Pr{ sigma(i,j) = +1 } and q = Pr{ sigma(i,j) = -1 }. These
//    are otherwise known as the Boltzmann factors.

double posterior(int * neighbors)
{
    // Hp and Hm are the Hamiltonians for sigma_ij = +/ -1
    double Hp = -1 * J * ((double)(*neighbors) + ((double)*(neighbors + 1)) + 
        ((double)*(neighbors + 2)) + ((double)*(neighbors + 3)));
    double Hm = +1 * J * ((double)(*neighbors) + ((double)*(neighbors + 1)) + 
        ((double)*(neighbors + 2)) + ((double)*(neighbors + 3)));
    // p and q are the probability that sigma_ij is + / -1
    double p = exp(-beta * J * Hp) / (exp(-beta * J * Hp) + exp(-beta * J * Hm));
    // q = np.exp(-beta * J * Hm) / (np.exp(-beta * J * Hp) + np.exp(-beta * J * Hm))
    return p;
}

//    Returns a zero - filled array sigma of dimension(i, j, K) which the Markov chain
//    is stored in.That's not completely true. Before returning the array, the 
//    initial Markov state is defined in sigma[:, : , 0], where each sigma(i, j, k = 0) =
//    +/ -1 with equal probability.

void initialize_sigma(int * sigma_series, double * H_Gibbs)
{
// creates sigma(k = 0), with elements + / -1 with equal probability
//  sigma_series stores sigma for every time - step
    
    int * sigma0 = new int[N];

    for (int k = 0; k < N; k++)
        *(sigma0+k) = (rand() % 2) * 2 - 1;  // initialize sigma with [-1,1] w/ probability p=0.5
    memset(sigma_series, 0, N * sizeof(int)); // np.zeros((L, L, K));
    memset(H_Gibbs, 0, K * sizeof(double));
    memcpy(sigma_series, sigma0, N * sizeof(int));
    *H_Gibbs = Hamiltonian(sigma0);

    // Deallocates Memory
    delete[] sigma0;
}

//   Returns the new state of sigma after updating sigma(i, j) according to the
//   posterior distribution used in Gibbs sampling.
int * Gibbs_update(int * sigma, int i, int j)
{
    int* neighbors = NULL;

    // The nearest neighbors of sigma(i, j) = sigma(i + 1, j), sigma(i - 1, j),
    //  sigma(i, j + 1), and sigma(i, j - 1).
    neighbors = get_neighbors(sigma, i, j);
    //  p and q are the probabilities of being spin + 1 or spin - 1, respectively.
    //  p + q = 1.
    double p = posterior(neighbors);
    //   sets sigma(i, j) according to wp and wm
    if (((double)rand() / (double)(RAND_MAX)) < p)
        *(sigma + i * L + j) *= 1;  // Accept Configuration (maintain sign)
    else
        *(sigma + i * L + j) *= -1; // Reject Configuration (invert sign)
    
    // returns new state for sigma
    return sigma;
}


//    Each iteration, the Gibbs sampler selects one of the L ^ 2 lattice elements
//    randomly, i.e.sigma(i, j).A new value of sigma(i, j) is then drawn from
//    the posterior distribution P[sigma(i, j) | all sigma(I != i, J != j)].
//    The posterior distribution includes only 4 sigma terms because the Ising
//    model assumes nearest neighbor interactions : sigma(I = i + -1, J = j + -1).Note
//    that sigma being updated one(i, j) pair at a time is the characteristic
//    partial resampling feature of Gibbs sampling.
void Gibbs(double* SP, double* H_Gibbs)
{
    //  sigma_series stores each configuration of the Markov chain
    int* prev_sigma_series = new int[N];
    int* sigma_series = new int[N];

    double Mag = 0.0;
    double Cum_Mag = 0.0;
    double G = 0.0;
    int i = 0, j = 0;
    double P = pow(2, N);

    // Process 1st sample separatedly (details on loop)
    initialize_sigma(prev_sigma_series, H_Gibbs);
    *(H_Gibbs) = Hamiltonian(prev_sigma_series);
    Mag = (double)1.0 / exp(-beta * J * (*(H_Gibbs)));
    Cum_Mag = Mag;
    G = Cum_Mag / P;
    *(SP) = log((double)1.0 / G) / N;

    // Selects(i, j) randomly, processes the new Hamiltonian of new configuration, calculating current Gamma 
    for (int k = 1; k < K; k++) {
        
        // Randomly select (i,j) to change configuration
        i = (rand() % L);
        j = (rand() % L);

        // Update the new configuration with updated version
        memcpy(sigma_series, Gibbs_update(prev_sigma_series, i, j), N * sizeof(int));

        // Calculates the Hamiltonian of new configuration
        *(H_Gibbs + k) = Hamiltonian(sigma_series);

        // Processes the calculus of magnetization of current configuration (need to store history???)
        Mag = 1.0 / exp(-beta * J * (*(H_Gibbs + k)));

        // Calculates Cummulative magnetization to calculate Gamma
        Cum_Mag += Mag;

        // Calculates Gamma (G) of current (k+1) samples 
        G = Cum_Mag / (((double)k+(double)1.0) * P);

        // Update the Sample Path of Gamma
        *(SP+k) =  log((double)1.0 / G) / N;

        // Update configuration for next k
        memcpy(prev_sigma_series, sigma_series, N * sizeof(int));
    }

    // Release memory allocated for configurations
    delete[] prev_sigma_series;
    delete[] sigma_series;

}

void simulateFreeEnergy(double * SP_FGibbs, double* H_Gibbs)
{
    std::cout << "\nStarting Simulations...";

    clock_t Gibbs_start_time = clock();

    /// simulates sigma using Gibbs sampler
    std::cout << "Starting Gibbs sampler for Beta=" << beta << " -- " << NoSim << " simulations with " << K << " samples each\n";

    // sigma_Gibbs, SP_FGibbs = Gibbs(SP_FGibbs)           // Generates K Samples using Gibbs Sampler

    Gibbs(SP_FGibbs, H_Gibbs);

    float ExecTime = (float)(clock() - Gibbs_start_time) / CLOCKS_PER_SEC;

    std::cout << "\nSimulations (Exec. time): --- " << ExecTime << " (in seconds) ---\n";

}


int getParametersFromUser(void)
{
    char YN = 'N';

    // Get the value of Beta (beta)
    std::cout << "Please enter the value of beta: ";
    std::cin >> IsingParam.beta;
    std::cout << "You entered "<< IsingParam.beta << ". Is this correct[Y/N]: ";
    std::cin >> YN;

    if (YN != 'Y' && YN != 'y') {
        std::cout << "\nProgram aborted by user!\n";
        return (FALSE);
    }
    
    // Get number of Samples (K)
    std::cout << "\nPlease enter the number of samples: ";
    std::cin >> IsingParam.K;
    std::cout << "You entered " << IsingParam.K << " samples. Is this correct[Y/N]: ";
    std::cin >> YN;

    if (YN != 'Y' && YN != 'y') {
        std::cout << "\nProgram aborted by user!\n";
        return (FALSE);
    }

    // Get size of the lattice (L)
    std::cout << "\nPlease enter the size of the lattice [2-30]: ";
    std::cin >> IsingParam.L;
    std::cout << "You entered " << IsingParam.L << " for lattice size. Is this correct[Y/N]: ";
    std::cin >> YN;

    if (IsingParam.L <= 1 || IsingParam.L > MAX_SIZE_LATTICE) {
        std::cout << "\nERROR: Invalid size of lattice (" << IsingParam.L << "). Program aborted!\n";
        return FALSE;
    }
     
    if (YN != 'Y' && YN != 'y') {
        std::cout << "\nProgram aborted by user!\n";
        return (FALSE);
    }

    // Get size value of J
    std::cout << "\nPlease enter the coupling energy: ";
    std::cin >> IsingParam.J;
    std::cout << "You entered " << IsingParam.J << " for coupling energy. Is this correct[Y/N]: ";
    std::cin >> YN;

    if (YN != 'Y' && YN != 'y') {
        std::cout << "\nProgram aborted by user!\n";
        return (FALSE);
    }

    return(TRUE);
}

long generateJobNumber(void)
{
    return (long)trunc(rand() % 100000);
}

int saveSamplePath(double * SP, int ID, std::string name, unsigned long size, int job)
{
    std::string fullFileName = name + "_ID" + std::to_string(ID) + "_" + std::to_string(job) + ".dat";
    
    std::fstream SPfile = std::fstream(fullFileName, std::ios::out | std::ios::binary);
    
    if (!SPfile) {
        std::cout << "\n>>> ERROR: Cannot open file (file:" << fullFileName << ")\n";
        return -1;
    }

    std::cout << "\nSaving Sample Paths... (file:" << fullFileName << ")\n";

    SPfile.write((char *)SP, (static_cast<std::streamsize>(size)) * sizeof(double));

    SPfile.close();
    
    std::cout << "Concluded!\n";

    return 0;
}

//void plot_sample_paths(double * SP, std::string lb_plot)
//{
//    plt::figure_size(1200, 720);
//    plt::title("Sample Path");
//    
//    // plt::xscale("log");
//    plt::plot(std::vector<double> SP);
    // plt::xlabel('log(Samples)');
//    
//    std::cout << "Saving result to " << lb_plot << std::endl;
//    plt::save(lb_plot);
//    plt::close()
//}


//inline std::string getCurrentDateTime(std::string s) {
//    time_t now = time(0);
//    struct tm  tstruct;
//    // char b[80] = "";
//    char buf[80];
//    char str[26]{};
//    // tstruct = *localtime(&now);
//    asctime_s(str, sizeof str, localtime_s(&tstruct, &now));
//    if (s == "now")
//        strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);
//    else if (s == "date")
//        strftime(buf, sizeof(buf), "%Y-%m-%d", &tstruct);
//    return buf;
//}

inline void Logger(std::string logMsg) {

    //std::string filePath = "log_" + getCurrentDateTime("date") + ".txt";
    // std::string now = getCurrentDateTime("now");
    std::string filePath = "log_2DIsingModel.txt";
    std::ofstream ofs(filePath.c_str(), std::ios_base::out | std::ios_base::app);
    // ofs << now << '\t' << logMsg << '\n';
    ofs << logMsg << '\n';
    ofs.close();
}

int main()
{

    std::cout << "IsingModel2D - Gibbs Sampling generation for Graphical Models - " << Version << "\n";
    std::cout << "-----------------------------------------------------------------------------\n\n";

    if (!getParametersFromUser())
        return -1;  // Program Aborted

    // Populate global variables with parameters received from user
    beta = IsingParam.beta;
    K = IsingParam.K;
    L = IsingParam.L;
    J = IsingParam.J;

    start_time = clock();  // Start counting time

    // Generate Job Number
    srand((unsigned)time(0));  // Initialize Random number generator

    job = generateJobNumber();

    std::cout << "\n--->>> JOB No. ------------ #" << job << " -----\n";
    std::cout << "Job Summary---\n";
    std::cout << "---------------------------------------------------------\n";
    std::cout << "Beta (temperature): " << beta << "\n";
    std::cout << "K (No. of Samples): " << K << "\n";
    std::cout << "L (size of lattice): " << L << "\n";
    std::cout << "J (coupling energy): " << J << "\n";
    std::cout << "---------------------------------------------------------\n\n";

    // Log Processing information

    Logger("---------------------------------------------------------\n"); 
    Logger("\n------------>>> JOB No. ------------ #"+ std::to_string(job));
    Logger("Job Summary---\n");
    Logger("---------------------------------------------------------\n");
    Logger("Beta (temperature): "+ std::to_string(beta));
    Logger("K (No. of Samples): "+ std::to_string(K));
    Logger("L (size of lattice): " + std::to_string(L));
    Logger("J (coupling energy): " + std::to_string(J));
    Logger("---------------------------------------------------------\n");

    // re-Initialize Random number generator for debugging purposes
    srand(1234);  

    double * SP_FGibbs_TMP = new double[K];

    memset(SP_FGibbs_TMP, 0, K*sizeof(double)); // for automatically-allocated arrays

    N = (int)pow(L, 2);

    burnin = (long)(K * 0.1);

    // Include Validation of Beta, K and N before processing

    // unsigned long size = 10000L;  // Size of Sample Path - CALCULATE!!!

    // Transferred from Gibbs function to optimize memory allocation
    double* H_Gibbs = new double[K];

    // Runs 04 Sample - Paths per Beta
    std::cout << "\nBLOCK - Beta = " << beta << " ---\n";
    for (int i = 0; i < NoSim; i++)
    {
        // plot 4 Sample paths per beta
        std::cout << "\n---- Sample Path " << (i + 1) << "- Beta = " << beta << " ---\n";
        simulateFreeEnergy(SP_FGibbs_TMP, H_Gibbs);
        saveSamplePath(SP_FGibbs_TMP, i, "SP_beta" + std::to_string(beta) + Version, K, job);
        Beep(freqP, durationP);
    }

    Beep(freqJob, durationJob);

    clock_t end_time = clock();

    float ExecTime = (float)(end_time - start_time) / CLOCKS_PER_SEC;

    std::cout << "\nGibbs Sampling Execution time: --- %d (in seconds) ---\n" << ExecTime << std::endl;

    // Deallocates memory from Gibbs variables
    delete[] H_Gibbs;

    // Deallocates memory from Sample Paths
    delete[] SP_FGibbs_TMP;

    return 0;
}

