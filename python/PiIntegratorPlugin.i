%module PiIntegratorPlugin

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

%{
#include "openmm/RPMDIntegrator.h"
#include "OpenMM.h"
#include "OpenMMDrude.h"

#include "PiIntegrator.h"
%}


namespace PiIntegratorPlugin {

class PiIntegrator : public OpenMM::Integrator {
public:
    PiIntegrator(double stepSize);
    void step(int steps);
};

} // namespace PiIntegratorPlugin
