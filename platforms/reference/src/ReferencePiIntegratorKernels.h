#ifndef OPI_REFERENCEPIINTEGRATORKERNELS_H_
#define OPI_REFERENCEPIINTEGRATORKERNELS_H_

#include "PiIntegratorKernels.h"


namespace PiIntegratorPlugin {

class ReferencePiIntegratorStepKernel : public PiIntegratorStepKernel {
public:
	ReferencePiIntegratorStepKernel(std::string name, const Platform& platform) :
			PiIntegratorStepKernel(name, platform) {
	}

	void initialize(const System& system, const PiIntegrator& integrator);

	void execute(ContextImpl& context, const PiIntegrator& integrator);

	double computeKineticEnergy(ContextImpl& context, const PiIntegrator& integrator);
};

} // namespace PiIntegratorPlugin


#endif /*OPI_REFERENCEPIINTEGRATORKERNELS_H_*/
