/* RESETTING + REFRACTORY PERIOD SIMULATION. OBTAINING THE PDF OF THE  PROPAGATION PHASE AT TIME IN X.*/

#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <initializer_list>
#include <random>
#include <string>
#include <fstream>
#include <algorithm>


std::random_device device;		// Random device engine, usually based on /dev/random on UNIX-like systems
std::mt19937 gen(device());		// Initialize Mersennes' twister using rd to generate the seed

int main() {

	std::uniform_real_distribution<double> dis(0.0, 1.0); // Real Random number between 0 and 1 with uniform distribution
	std::normal_distribution<double> gaussian(0.0, 1.0); // Normal distribution with (mean) mu = 0 and (variance) sigma2 = 1.

	// We are going to simulate the dynamics of a typical one-dimensional brownian particle spreads (propagation phase) but after certain time it resets to a position xr. 
	// After each reset the particle stays motionless (refractory period phase).
	// The duration of each phase follow their probability density fuction, typically a poissonian distribution.


	///////////////////////////////////////////////////	
	// PARAMETERS
	///////////////////////////////////////////////////

	double r1 = 1.0; // Rate of resetting (Rate of the PDF which indicates the duration of the propagation phase)
	double r2 = 2.0; // Rate of refractory period phase (Duration of refractory period phase)
	double D = 1.0; // Diffusion constant
	double xr = 0.0; // Mean position of resetting
	double x0 = 0.0; // Initial position

	double tmax = 1.0e+1; // Simulation ends when the time reachs tmax
	double t = 0.0; // Current time 

	double dt = 1.0e-4; // Step between two consecutive times
	double Dt = 0.10;
	int tcomp = tmax / dt;
	int tout = Dt / dt;
	int tsave = tmax / Dt; 

	int nbins = 101; // Number of bins left/right to the origin (50 + 1 where we set the boundary)
	double a = 10.0; // Position of the boundary

	double dx = a / (nbins - 1); // Width of the bin

	double c1 = sqrt(2 * D * dt); // Amplitude of the white gaussian noise

	int nruns = 10000000; // How many simulations we are carrying out

	///////////////////////////////////////////////////	
	// FILES
	///////////////////////////////////////////////////

	// We want to describe the PDF of the propagation phase for any time. 
	// We are going to obtain the number of counts for each value of x and t and what phase the system is in.

	std::string r1String = "_r1=" + std::to_string(r1).substr(0, 3);
	std::replace(r1String.begin(), r1String.end(), '.', '_');
	std::string r2String = "_r2=" + std::to_string(r2).substr(0, 3);
	std::replace(r2String.begin(), r2String.end(), '.', '_');
	std::string binString = "_bins=" + std::to_string(nbins);

	std::string pDiffFile = "pDiff" + r1String + r2String + binString + "_runs=" + std::to_string(nruns) + ".dat";
	std::string pRefrFile = "pRefr" + r1String + r2String + binString + "_runs=" + std::to_string(nruns) + ".dat";

	std::ofstream output1(pDiffFile.c_str(), std::ofstream::trunc);
	std::ofstream output2(pRefrFile.c_str(), std::ofstream::trunc);

	///////////////////////////////////////////////////	
	// VARIABLES
	///////////////////////////////////////////////////

	double tr = 0.0; 
	double tau = 0.0;
	std::vector<double> tevents;

	double xOld = 0, xNew = 0;

	// Counters
	double tcounter = 0.0; 
	int rcounter = 0; 
	int outcounter = tout;
	int k = 0;
	int p = 0;

	// Output matrix
	std::vector<int> nrefr(tsave + 1, 0);
	std::vector<std::vector<int>> ndiff(tsave + 1, std::vector<int>(2 * nbins, 0));

	///////////////////////////////////////////////////	
	// STOCHASTIC GENERATION OF RESETTING TIMES
	///////////////////////////////////////////////////

	while (p < nruns)
	{
		while (tcounter <= tmax)
		{
			tr = tau - 1 / r1 * log(dis(gen));
			tau = tr -1 / r2 * log(dis(gen));
			tevents.push_back(tr);
			tevents.push_back(tau);

			tcounter = tau;
		}

		// Initial position
		xOld = 0; 
		xNew = 0;

		for (size_t i = 0; i <= tcomp; i++)
		{
			if ( tevents[rcounter] <= t )
			{
				xNew = 0; // Reset to xr ( = 0 )
				rcounter++;
			}
			if ( rcounter % 2 == 0 )
			{
				xNew = xOld + c1 * gaussian(gen);

			}
			if (tout == outcounter)
			{
				if (rcounter % 2 != 0)
				{
					nrefr[k] += 1;
				} 
				else {
					if (xOld <= -a) {
						ndiff[k][0] += 1;
					}
					else if (xOld > a) {
						ndiff[k].back() += 1;
					}
					else {
						ndiff[k][nbins + std::floor(xOld / dx)] += 1;
					}
				}
				outcounter = 0;
				k++;
			}
			xOld = xNew;
			outcounter++;
			t += dt;
		}
		tevents.clear();
		p++;

		tr = 0.0;
		tau = 0.0; 

		tcounter = 0.0;
		rcounter = 0;
		outcounter = tout;
		k = 0;
		t = 0;

		std::cout << " Run = " << p << std::endl;
	}
	for (size_t i = 0; i < ndiff.size(); i++)
	{
		for (size_t j = 0; j < ndiff[i].size(); j++)
		{
			output1 << ndiff[i][j] << " ";
		}
		output1 << std::endl;
	}
	
	output1.close();
	for (size_t i = 0; i < nrefr.size(); i++)
	{
		output2 << nrefr[i] << " ";
	}

	output2.close();
	return 0;

}