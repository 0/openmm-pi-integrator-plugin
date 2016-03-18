#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"

#include "PiIntegrator.h"
#include "PiIntegratorKernels.h"


using namespace PiIntegratorPlugin;


PiIntegrator::PiIntegrator(double stepSize) {
	setStepSize(stepSize);
}

void PiIntegrator::step(int steps) {
	if (context == NULL) {
		throw OpenMMException("This Integrator is not bound to a context");
	}

	for (int i = 0; i < steps; i++) {
		kernel.getAs<PiIntegratorStepKernel>().execute(*context, *this);
	}
}

void PiIntegrator::initialize(ContextImpl& contextRef) {
	if (owner != NULL && &contextRef.getOwner() != owner) {
		throw OpenMMException("This Integrator is already bound to a context");
	}

	if (contextRef.getSystem().getNumConstraints() > 0) {
		throw OpenMMException("PiIntegrator cannot be used with Systems that include constraints");
	}

	context = &contextRef;
	owner = &contextRef.getOwner();
	kernel = context->getPlatform().createKernel(PiIntegratorStepKernel::Name(), contextRef);

	kernel.getAs<PiIntegratorStepKernel>().initialize(contextRef.getSystem(), *this);
}

void PiIntegrator::cleanup() {
	kernel = Kernel();
}

std::vector<std::string> PiIntegrator::getKernelNames() {
	std::vector<std::string> names;
	names.push_back(PiIntegratorStepKernel::Name());

	return names;
}

double PiIntegrator::computeKineticEnergy() {
	return kernel.getAs<PiIntegratorStepKernel>().computeKineticEnergy(*context, *this);
}
