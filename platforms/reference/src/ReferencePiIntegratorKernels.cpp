#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/SimTKOpenMMRealType.h"

#include "PiIntegrator.h"
#include "ReferencePiIntegratorKernels.h"


using namespace PiIntegratorPlugin;


void ReferencePiIntegratorStepKernel::initialize(const System& system, const PiIntegrator& integrator) {
}

void ReferencePiIntegratorStepKernel::execute(ContextImpl& context, const PiIntegrator& integrator) {
	const RealOpenMM dt = integrator.getStepSize();

	context.setTime(context.getTime() + dt);
}

double ReferencePiIntegratorStepKernel::computeKineticEnergy(ContextImpl& context, const PiIntegrator& integrator) {
	return 0.;
}
