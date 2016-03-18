#include "openmm/reference/ReferencePlatform.h"
#include "openmm/Platform.h"
#include "openmm/OpenMMException.h"

#include "ReferencePiIntegratorKernelFactory.h"
#include "ReferencePiIntegratorKernels.h"


using namespace OpenMM;
using namespace PiIntegratorPlugin;


extern "C" void registerPlatforms() {
}

extern "C" void registerKernelFactories() {
	for (int i = 0; i < Platform::getNumPlatforms(); i++) {
		Platform& platform = Platform::getPlatform(i);

		if (dynamic_cast<ReferencePlatform*>(&platform) == NULL) {
			continue;
		}

		ReferencePiIntegratorKernelFactory* factory = new ReferencePiIntegratorKernelFactory();
		platform.registerKernelFactory(PiIntegratorStepKernel::Name(), factory);
	}
}


KernelImpl* ReferencePiIntegratorKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
	if (name == PiIntegratorStepKernel::Name()) {
		return new ReferencePiIntegratorStepKernel(name, platform);
	}

	throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '") + name + "'").c_str());
}
