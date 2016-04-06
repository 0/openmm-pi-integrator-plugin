#ifndef OPI_PIINTEGRATOR_H_
#define OPI_PIINTEGRATOR_H_

#include "openmm/Integrator.h"
#include "openmm/Kernel.h"


// Reduced Planck constant (from 2010 CODATA recommended values from NIST).
#define HBAR (6.35077993e-2) /* kJ ps/mol */


using namespace OpenMM;


namespace PiIntegratorPlugin {

class PiIntegrator : public Integrator {
public:
	PiIntegrator(double stepSize, double beta, int numBeads, double centroidFriction);

	double getBeta() const {
		return beta;
	}

	int getNumBeads() const {
		return numBeads;
	}

	int getCentroidFriction() const {
		return centroidFriction;
	}

	void step(int steps);

protected:
	void initialize(ContextImpl& context);

	void cleanup();

	std::vector<std::string> getKernelNames();

	double computeKineticEnergy();

private:
	Kernel kernel;

	double beta, centroidFriction;
	int numBeads;
};

} // namespace PiIntegratorPlugin


#endif /*OPI_PIINTEGRATOR_H_*/
