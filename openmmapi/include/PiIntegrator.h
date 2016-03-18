#ifndef OPI_PIINTEGRATOR_H_
#define OPI_PIINTEGRATOR_H_

#include "openmm/Integrator.h"
#include "openmm/Kernel.h"


using namespace OpenMM;


namespace PiIntegratorPlugin {

class PiIntegrator : public Integrator {
public:
	PiIntegrator(double stepSize);

	void step(int steps);

protected:
	void initialize(ContextImpl& context);

	void cleanup();

	std::vector<std::string> getKernelNames();

	double computeKineticEnergy();

private:
	Kernel kernel;
};

} // namespace PiIntegratorPlugin


#endif /*OPI_PIINTEGRATOR_H_*/
