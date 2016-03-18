#ifndef OPI_REFERENCEPIINTEGRATORKERNELFACTORY_H_
#define OPI_REFERENCEPIINTEGRATORKERNELFACTORY_H_

#include "openmm/KernelFactory.h"


namespace OpenMM {

class KernelImpl;

} // namespace OpenMM

using namespace OpenMM;


namespace PiIntegratorPlugin {

class ReferencePiIntegratorKernelFactory : public KernelFactory {
public:
	KernelImpl* createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const;
};

} // namespace PiIntegratorPlugin


#endif /*OPI_REFERENCEPIINTEGRATORKERNELFACTORY_H_*/
