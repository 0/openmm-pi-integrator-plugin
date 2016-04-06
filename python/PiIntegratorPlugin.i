%module PiIntegratorPlugin

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

%{
#include "openmm/RPMDIntegrator.h"
#include "OpenMM.h"
#include "OpenMMDrude.h"

#include "PiIntegrator.h"
%}


%pythoncode %{
import simtk.unit as unit
%}

%pythonappend PiIntegratorPlugin::PiIntegrator::getBeta() const %{
    val = unit.Quantity(val, 1/unit.kilojoule_per_mole)
%}

%pythonappend PiIntegratorPlugin::PiIntegrator::getCentroidFriction() const %{
    val = unit.Quantity(val, 1/unit.picosecond)
%}


namespace PiIntegratorPlugin {

class PiIntegrator : public OpenMM::Integrator {
public:
    PiIntegrator(double stepSize, double beta, int numBeads, double centroidFriction);

    double getBeta() const;

    int getNumBeads() const;

    int getCentroidFriction() const;

    void step(int steps);
};

} // namespace PiIntegratorPlugin
