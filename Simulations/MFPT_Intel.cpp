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
#include <chrono>


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
	double r2 = 5.0; // Rate of refractory period phase (Duration of refractory period phase)

	std::cout << "r1 = " << r1 << " r2 = " << r2 << std::endl;

	double D = 1.0; // Diffusion constant
	//double xr = 0.0; // Position of resetting
	//double x0 = 0.0; // Initial position

	double tmax = 1.0e+1; // Simulation ends when the time reachs tmax
	double t = 0.0; // Current time 

	double dt = 1.0e-5; // Step between two consecutive times

	double a = 1; // Target

	double c1 = sqrt(2 * D * dt); // Amplitude of the white gaussian noise

	int nruns = 1e6; // How many simulations we are carrying out

	///////////////////////////////////////////////////	
	// OUTPUT FILES
	///////////////////////////////////////////////////

	// We want to describe the PDF of the propagation phase for any time. 
	// We are going to obtain the number of counts for each value of x and t and what phase the system is in.

	std::string r1String = "_r1=" + std::to_string(r1).substr(0, 3);
	std::replace(r1String.begin(), r1String.end(), '.', '_');
	std::string r2String = "_r2=" + std::to_string(r2).substr(0, 3);
	std::replace(r2String.begin(), r2String.end(), '.', '_');

	std::string MFPTFile = "Intel_MFPT" + r1String + r2String + "_runs=" + std::to_string(nruns) + ".dat";

	std::ofstream output1(MFPTFile.c_str(), std::ofstream::trunc);

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

	// Output matrix

	std::vector<double> FPT(nruns);


	///////////////////////////////////////////////////	
	// MEASURING TIME
	///////////////////////////////////////////////////

	// Start measuring time
	auto begin = std::chrono::high_resolution_clock::now();

	///////////////////////////////////////////////////	
	// STOCHASTIC GENERATION OF RESETTING TIMES
	///////////////////////////////////////////////////

	for (size_t p = 0; p < nruns; p++)
	{
		while (xOld < a)
		{
			if (tau <= t)
			{
				tr = tau - 1 / r1 * log(dis(gen));
				tau = tr - 1 / r2 * log(dis(gen));
				tevents.push_back(tr);
				tevents.push_back(tau);

				// tcounter = tau;
			}
			if (tevents[rcounter] <= t)
			{
				xNew = 0; // Reset to xr ( = 0 )
				rcounter++;
			}
			if (rcounter % 2 == 0)
			{
				xNew = xOld + c1 * gaussian(gen);

			}

			xOld = xNew;
			t += dt;
		}
		// std::cout << "xOld = " << xOld << " , tfin = 0" << t << std::endl;

		FPT[p] = t;

		// Reset all the variables

		tevents.clear();
		tr = 0.0;
		tau = 0.0;
		rcounter = 0;
		t = 0;

		// Initial position
		xOld = 0;
		xNew = 0;

		std::cout << p << std::endl;
	}

	for (size_t i = 0; i < FPT.size(); i++)
	{
		output1 << FPT[i] << " ";
	}

	output1.close();

	// Stop measuring time and calculate the elapsed time
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin);

	std::cout << " Time measured: " << elapsed.count() * 1e-9;

	return 0;

}