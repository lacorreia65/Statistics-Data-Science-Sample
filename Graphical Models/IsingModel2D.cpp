// IsingModel2D.cpp : This subsystem emules a 2D Ising Model with boundary conditions
// 
// Author - Luis Alvaro Correia
// 
// Date - Jan 7th 2021
//

#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <string>
#include <Windows.h>
// #include "matplotlibcpp.h"

// namespace plt = matplotlibcpp;

const int MAX_SIZE_LATTICE = 15;

struct ModelParam
{
    float beta;
    long K;
} Param;

//struct stGibbs
//{
//    sigma 
//};

long    L = 0L,                                         // Default size of the 2D Ising Model
        N = 0L,                                         // Size of the lattice
        K = 0L,                                         // Number of Samples per Gibbs simulation;
        burnin = 0L;                                    // Number of samples to discard  
float   beta = 0.0f;                                     // Initial value for Beta
float   lBeta[5] = { 0.01f, 0.06f, 0.15f, 0.25f, 0.75f };   // List of possible values of Beta
float   J = 0.01f;                                      // Value of Coupling Energy
int     NoSim = 4;                                      // Number of Monte Carlo Simulations per beta(= no.sample paths)
int     job = 0;

char Version[] = "_v1-CPP_";

clock_t start_time = clock();

int durationP = 800;     // milliseconds
int freqP = 1000;        // Hz
int durationJob = 1000;  // milliseconds
int freqJob = 1500;      // Hz

//int* AdjVert(int* p, int L)
//{
//   static int Vert[4][2]{};
//
//    Vert[0][0] = *p == 0 ? (L - 1) : *p - 1; Vert[0][1] = *(p + 1);
//    Vert[1][0] = *p == (L - 1) ? 0 : *p + 1; Vert[1][1] = *(p + 1);
//    Vert[2][0] = *p; Vert[2][1] = *(p + 1) == (L - 1) ? 0 : *(p + 1) + 1;
//    Vert[3][0] = *p; Vert[3][1] = *(p + 1) == 0 ? (L - 1) : *(p + 1) - 1;
//
//    return &Vert[0][0];
//}

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

float Hamiltonian(int *x)
{
    float S = 0.0f;
    int* neighbors = NULL;
    
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++) {
            neighbors = get_neighbors(x, i, j);
            S -= J * (*(x + 2 * i + j) * (*neighbors) + *(x + 2 * i + j) * (*(neighbors+1)) +
                *(x + 2 * i + j) * (*(neighbors + 2)) + *(x + 2 * i + j) * (*(neighbors + 3)));
            }
    return S;

}

//    Returns p = Pr{ sigma(i,j) = +1 } and q = Pr{ sigma(i,j) = -1 }. These
//    are otherwise known as the Boltzmann factors.

float posterior(int * neighbors, int * sigma)
{
    // Hp and Hm are the Hamiltonians for sigma_ij = +/ -1
    float Hp = -1 * J * ((*neighbors) + (*(neighbors + 1)) + (*(neighbors + 2)) + (*(neighbors + 3)));
    float Hm = +1 * J * ((*neighbors) + (*(neighbors + 1)) + (*(neighbors + 2)) + (*(neighbors + 3)));
    // p and q are the probability that sigma_ij is + / -1
    float p = exp(-beta * J * Hp) / (exp(-beta * J * Hp) + exp(-beta * J * Hm));
    // q = np.exp(-beta * J * Hm) / (np.exp(-beta * J * Hp) + np.exp(-beta * J * Hm))
    return p;
}

//    Returns a zero - filled array sigma of dimension(i, j, K) which the Markov chain
//    is stored in.That's not completely true. Before returning the array, the 
//    initial Markov state is defined in sigma[:, : , 0], where each sigma(i, j, k = 0) =
//    +/ -1 with equal probability.

void initialize_sigma(int * sigma_series, float * H_Gibbs)
{
// creates sigma(k = 0), with elements + / -1 with equal probability
//  sigma_series stores sigma for every time - step
    
    int * sigma0 = new int[N];

    for (int k = 0; k < N; k++)
        *(sigma0+k) = (rand() % 2) * 2 - 1;  // initialize sigma with [-1,1] w/ probability p=0.5
    memset(sigma_series, 0, N * sizeof(int)); // np.zeros((L, L, K));
    memset(H_Gibbs, 0, K * sizeof(float));
    memcpy(sigma_series, sigma0, N * sizeof(int));
    // std::copy(sigma_series = sigma0;
    *H_Gibbs = Hamiltonian(sigma0);
    // return sigma_series, H_Gibbs;
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
    float p = posterior(neighbors, sigma);
    //   sets sigma(i, j) according to wp and wm
    if (((float)rand() / (RAND_MAX)) < p)
        *(sigma + i * L + j) = 1;
    else
        *(sigma + i * L + j) = -1;
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
void Gibbs(float* SP)
{
    //  sigma_series stores sigma for each step of the Markov chain
    float* H_Gibbs = new float[K];
    int* sigma_series = new int[N*K];
    float* Mag = new float[K];
    float Cum_Mag = 0.0f;
    float G = 0.0f;
    int i = 0, j = 0;
    float P = (float)pow(2, N);

    initialize_sigma(sigma_series, H_Gibbs);
    // H_Gibbs = initialize_sigma()
    // Initialize Hamiltonian
    //  selects(i, j) randomly and calls draw_sigma_ij(...) to update sigma(i, j)
    for (int k = 1; k < K; k++) {
        i = (rand() % L) + 1;
        j = (rand() % L) + 1;
        memcpy((sigma_series + k * N), Gibbs_update(sigma_series + (k - 1) * N, i, j), N * sizeof(int));
        *(H_Gibbs + k) = Hamiltonian(sigma_series + k * N);
        *(Mag+k) = (float) 1.0 / exp(-beta * J * (*(H_Gibbs + k)));
        Cum_Mag += *(Mag + k);
        G = Cum_Mag / (k * P);
        *(SP+k) = (float) log(1.0 / G) / N;
    }
    // return sigma_series, SP
}



void simulateFreeEnergy(float * SP_FGibbs)
{
    std::cout << "\nStarting Simulations...";

    clock_t Gibbs_start_time = clock();

    /// simulates sigma using Gibbs sampler
    std::cout << "Starting Gibbs sampler for Beta=" << beta << " -- " << NoSim << " simulations with " << K << " samples each\n";

    // sigma_Gibbs, SP_FGibbs = Gibbs(SP_FGibbs)           // Generates K Samples using Gibbs Sampler

    Gibbs(SP_FGibbs);

    float ExecTime = (float)(clock() - Gibbs_start_time) / CLOCKS_PER_SEC;

    std::cout << "\nSimulations (Exec. time): --- " << ExecTime << " (in seconds) ---\n";
    // return SP_FGibbs;

}


ModelParam getBetaKFromUser(void)
{
    ModelParam P{-1.0, -1L};
    //float value1;
    //long value2;
    char YN1, YN2;

    std::cout << "Please enter the value of beta: ";
    std::cin >> P.beta;
    std::cout << "You entered "<< P.beta << ". Is this correct[Y/N]: ";
    std::cin >> YN1;

    if (YN1 == 'Y' || YN1 == 'y') 
    {
        std::cout << "Please enter the number of samples: ";
        std::cin >> P.K;
        std::cout << "You entered " << P.K << " samples. Is this correct[Y/N]: ";
        std::cin >> YN2;

        if (YN2 == 'Y' or YN2 == 'y')
            return P;
        else
            return { -1.0, -1L };
    }
    else
        return { -1.0, 0L };
}


int getNFromUser(void)
{
    int  value1;
    char YN1;

    std::cout << "Please enter the size of the lattice [1-15]: ";
    std::cin >> value1;
    std::cout << "You entered " << value1 << ". Is this correct[Y/N]: ";
    std::cin >> YN1;
    if (YN1 == 'Y' or YN1 == 'y')
        if (value1 <= 0 || value1 > MAX_SIZE_LATTICE)
            return -1;
        else
            return value1;
    else
        return -1;
}

long generateJobNumber(void)
{
    return (long)trunc(rand() % 100000);
}

int saveSamplePath(float * SP, int ID, std::string name, unsigned long size, int job)
{
    std::string fullFileName = name + "_ID" + std::to_string(ID) + "_" + std::to_string(job) + ".dat";
    
    std::fstream SPfile = std::fstream(fullFileName, std::ios::out | std::ios::binary);
    
    if (!SPfile) {
        std::cout << "\n>>> ERROR: Cannot open file (file:" << fullFileName << ")\n";
        return -1;
    }

    std::cout << "\nSaving Sample Paths... (file:" << fullFileName << ")\n";

    SPfile.write((char *)SP, (static_cast<std::streamsize>(size)) * sizeof(float));

    SPfile.close();
    
    std::cout << "Concluded!\n";

    return 0;
}

void plot_sample_paths(float * SP, std::string lb_plot)
{
    plt::figure(figsize = (12, 10));
    plt::xscale("log");
    plt::plot(SP[0, :]);
    plt::plot(SP[1, :]);
    plt::plot(SP[2, :]);
    plt::plot(SP[3, :]);
    plt::xlabel('log(Samples)');
    plt::savefig(lb_plot);
    plt::close()
}




int main()
{
    //int p[2];
    
    //p[0] = 0;
    //p[1] = 0;

    srand((unsigned)time(0));  // Initialize Random number generator

    // std::cout << p;

    //int myList[5] = { 0, 4, 8, 12, 16 }; //Line 1
    //int yourList[5]; //Line 2

    //for (int index = 0; index < 5; index++)
    //    yourList[index] = myList[index];

    //std::cout myList;

    // int Up[2], Down[2], Right[2], Left[2];

    //int aVert[4][2] = *(AdjVert(p, L));

    //for (int i = 0; i < L; i++)
    //    for (int j = 0; j < 2; j++) {
    //        std::cout << "*(vert + " << i << ";" << j << ") : ";
    //        std::cout << *(aVert + 2*i + j) << ")\n";
    //    }

    //AdjVert(int p[2])


    job = generateJobNumber();

    std::cout << "\n--->>> JOB No. ------------" << job << "\n";


    // memset(myarray, 0, N * sizeof(*myarray)); // for heap-allocated arrays, where N is the number of elements

    Param = getBetaKFromUser();  // Get input for Betaand K from user

    beta = Param.beta;
    K = Param.K;

    float * SP_FGibbs_TMP = new float[K];

    memset(SP_FGibbs_TMP, 0x30, K*sizeof(float)); // for automatically-allocated arrays

    std::cout << "\nBeta from user: " << beta << "\n\n";
    std::cout << "\nK from user: " << K << "\n\n";

    L = getNFromUser();

    N = (int)pow(L, 2);

    burnin = (long)(K * 0.1);

    std::cout << "\nL from User=" << L << "\n\n";

    // Include Validation of Beta, K and N before processing

    // unsigned long size = 10000L;  // Size of Sample Path - CALCULATE!!!

    // Runs 04 Sample - Paths per Beta
    std::cout << "\nBLOCK - Beta = " << beta << " ---\n";
    for (int i = 0; i < NoSim; i++)
    {
        // plot 4 Sample paths per beta
        std::cout << "\n---- Sample Path " << (i + 1) << "- Beta = " << beta << " ---\n";
        simulateFreeEnergy(SP_FGibbs_TMP);
        saveSamplePath(SP_FGibbs_TMP, i, "SP_beta" + std::to_string(beta) + Version, K, job);
        Beep(freqP, durationP);
    }

    Beep(freqJob, durationJob);

    clock_t end_time = clock();

    float ExecTime = (float)(end_time - start_time) / CLOCKS_PER_SEC;

    std::cout << "\nGibbs Sampling Execution time: --- %d (in seconds) ---\n" << ExecTime << std::endl;

    return 0;
}

