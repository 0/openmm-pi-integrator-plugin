#ifndef OPI_PIINTEGRATORKERNELS_H_
#define OPI_PIINTEGRATORKERNELS_H_

#include "openmm/KernelImpl.h"


namespace OpenMM {

class ContextImpl;
class Platform;
class System;

} // namespace OpenMM

using namespace OpenMM;


namespace PiIntegratorPlugin {

class PiIntegrator;

class PiIntegratorStepKernel : public KernelImpl {
public:
	PiIntegratorStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
	}

	static std::string Name() {
		return "PiIntegratorStep";
	}

	virtual void initialize(const System& system, const PiIntegrator& integrator) = 0;

	virtual void execute(ContextImpl& context, const PiIntegrator& integrator) = 0;

	virtual double computeKineticEnergy(ContextImpl& context, const PiIntegrator& integrator) = 0;
};

} // namespace PiIntegratorPlugin


#endif /*OPI_PIINTEGRATORKERNELS_H_*/
