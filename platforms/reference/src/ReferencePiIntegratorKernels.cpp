#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/reference/SimTKOpenMMRealType.h"
#include "openmm/reference/SimTKOpenMMUtilities.h"

#include "PiIntegrator.h"
#include "ReferencePiIntegratorKernels.h"


using namespace std;
using namespace OpenMM;
using namespace PiIntegratorPlugin;


static vector<RealVec>& extractPositions(ContextImpl& context) {
	ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
	return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractVelocities(ContextImpl& context) {
	ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
	return *((vector<RealVec>*) data->velocities);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
	ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
	return *((vector<RealVec>*) data->forces);
}


/**
 * Convert all the Cartesian coordinates to normal mode coordinates.
 *
 * This uses a direct application of the DFT instead of calling an FFT package.
 */
static void cartesianToNormalMode(vector<RealVec>& x, vector<RealVec>& nm, int P) {
	int numPhysicalParticles = x.size() / P;

	// M2 is the index of the first mode to use sin. For an odd number of
	// beads, M1 is ignored; for an even number of beads, it is the index of
	// the strictly real mode in the middle.
	int M1 = (P+1)/2;
	int M2 = (P+2)/2;

	for (int n = 0; n < numPhysicalParticles; n++) {
		// Offset for the start of the nth physical particle.
		int i = n*P;

		// Clear the result.
		for (int k = 0; k < P; k++) {
			nm[i+k] *= 0.;
		}

		for (int j = 0; j < P; j++) {
			nm[i] += x[i+j]*(RealOpenMM)SQRT(1./P);

			for (int k = 1; k < M1; k++) {
				nm[i+k] += x[i+j]*SQRT(2./P)*COS(2.*M_PI*j*k/P);
			}

			if (M1 < M2) {
				nm[i+M1] += x[i+j]*SQRT(1./P)*POW(-1, j);
			}

			for (int k = M2; k < P; k++) {
				nm[i+k] += x[i+j]*SQRT(2./P)*SIN(2.*M_PI*j*k/P);
			}
		}
	}
}

/**
 * Convert all the normal mode coordinates to Cartesian coordinates.
 *
 * This uses a direct application of the DFT instead of calling an FFT package.
 */
static void normalModeToCartesian(vector<RealVec>& x, vector<RealVec>& nm, int P) {
	int numPhysicalParticles = x.size() / P;

	// M2 is the index of the first mode to use sin. For an odd number of
	// beads, M1 is ignored; for an even number of beads, it is the index of
	// the strictly real mode in the middle.
	int M1 = (P+1)/2;
	int M2 = (P+2)/2;

	for (int n = 0; n < numPhysicalParticles; n++) {
		// Offset for the start of the nth physical particle.
		int i = n*P;

		// Clear the result.
		for (int j = 0; j < P; j++) {
			x[i+j] *= 0.;
		}

		for (int j = 0; j < P; j++) {
			x[i+j] += nm[i]*(RealOpenMM)SQRT(1./P);

			for (int k = 1; k < M1; k++) {
				x[i+j] += nm[i+k]*SQRT(2./P)*COS(2.*M_PI*j*k/P);
			}

			if (M1 < M2) {
				x[i+j] += nm[i+M1]*SQRT(1./P)*POW(-1, j);
			}

			for (int k = M2; k < P; k++) {
				x[i+j] += nm[i+k]*SQRT(2./P)*SIN(2.*M_PI*j*k/P);
			}
		}
	}
}

/**
 * Apply the friction and random force of the thermostat.
 */
static void applyThermostat(vector<RealVec>& nm_vel, vector<RealOpenMM>& masses, int P, RealOpenMM tau, RealOpenMM dt, RealOpenMM centroidFriction) {
	int numPhysicalParticles = nm_vel.size() / P;

	for (int n = 0; n < numPhysicalParticles; n++) {
		// Offset for the start of the nth physical particle.
		int i = n*P;

		// Don't move particles with zero mass.
		if (masses[i] == 0.) {
			continue;
		}

		for (int k = 0; k < P; k++) {
			RealOpenMM gamma;

			if (k == 0) {
				gamma = centroidFriction;
			} else {
				RealOpenMM omega_k = 2.*SIN(M_PI*k/P)/(HBAR*tau);
				gamma = 2.*omega_k;
			}

			RealOpenMM c1 = EXP(-0.5*dt*gamma);
			RealOpenMM c2 = SQRT((1-c1*c1)/(masses[i+k]*tau));

			for (int dim = 0; dim < 3; dim++) {
				nm_vel[i+k][dim] *= c1;
				nm_vel[i+k][dim] += c2*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
			}
		}
	}
}

/*
 * Apply the accelerations due to the force field.
 */
static void applyForces(vector<RealVec>& vel, vector<RealVec>& f, vector<RealOpenMM>& masses, RealOpenMM dt) {
	for (int k = 0; k < vel.size(); k++) {
		// Don't move particles with zero mass.
		if (masses[k] == 0.) {
			continue;
		}

		vel[k] += f[k]*dt/(2.*masses[k]);
	}
}

/*
 * Apply the free particle propagation.
 */
static void propagateFreeParticle(vector<RealVec>& nm_pos, vector<RealVec>& nm_vel, vector<RealOpenMM>& masses, int P, RealOpenMM tau, RealOpenMM dt) {
	int numPhysicalParticles = nm_pos.size() / P;

	for (int n = 0; n < numPhysicalParticles; n++) {
		// Offset for the start of the nth physical particle.
		int i = n*P;

		// Don't move particles with zero mass.
		if (masses[i] == 0.) {
			continue;
		}

		// Centroid mode.
		nm_pos[i] += nm_vel[i]*dt;

		// Other modes.
		for (int k = 1; k < P; k++) {
			RealOpenMM omega_k = 2.*SIN(M_PI*k/P)/(HBAR*tau);
			RealOpenMM c = COS(omega_k*dt);
			RealOpenMM s = SIN(omega_k*dt);

			RealVec temp = nm_pos[i+k];

			nm_pos[i+k] *= c;
			nm_pos[i+k] += nm_vel[i+k]*s/omega_k;
			nm_vel[i+k] *= c;
			nm_vel[i+k] -= temp*s*omega_k;
		}
	}
}


void ReferencePiIntegratorStepKernel::initialize(const System& system, const PiIntegrator& integrator) {
	numParticles = system.getNumParticles();
	numBeads = integrator.getNumBeads();
	tau = integrator.getBeta() / numBeads;
	centroidFriction = integrator.getCentroidFriction();
}

void ReferencePiIntegratorStepKernel::execute(ContextImpl& context, const PiIntegrator& integrator) {
	const System& system = context.getSystem();
	const RealOpenMM dt = integrator.getStepSize();

	vector<RealOpenMM> masses(numParticles);

	for (int j = 0; j < numParticles; j++) {
		masses[j] = system.getParticleMass(j);
	}

	vector<RealVec>& pos = extractPositions(context);
	vector<RealVec>& vel = extractVelocities(context);
	vector<RealVec>& f = extractForces(context);
	vector<RealVec> nm_pos(numParticles);
	vector<RealVec> nm_vel(numParticles);

	// Take a single step of the PILE integrator, as described in: Ceriotti et
	// al., J. Chem. Phys. 133, 124104 (2010). Velocities are used instead of
	// momenta to match OpenMM.

	// Half-step 1.
	cartesianToNormalMode(vel, nm_vel, numBeads);
	applyThermostat(nm_vel, masses, numBeads, tau, dt, centroidFriction);
	normalModeToCartesian(vel, nm_vel, numBeads);

	// Half-step 2.
	computeForces(context, integrator);
	applyForces(vel, f, masses, dt);

	// Step 3.
	cartesianToNormalMode(pos, nm_pos, numBeads);
	cartesianToNormalMode(vel, nm_vel, numBeads);
	propagateFreeParticle(nm_pos, nm_vel, masses, numBeads, tau, dt);
	normalModeToCartesian(pos, nm_pos, numBeads);
	normalModeToCartesian(vel, nm_vel, numBeads);

	// Half-step 2.
	computeForces(context, integrator);
	applyForces(vel, f, masses, dt);

	// Half-step 1.
	cartesianToNormalMode(vel, nm_vel, numBeads);
	applyThermostat(nm_vel, masses, numBeads, tau, dt, centroidFriction);
	normalModeToCartesian(vel, nm_vel, numBeads);

	context.setTime(context.getTime() + dt);
}

double ReferencePiIntegratorStepKernel::computeKineticEnergy(ContextImpl& context, const PiIntegrator& integrator) {
	const System& system = context.getSystem();
	vector<RealVec>& vel = extractVelocities(context);
	double energy = 0.;

	// Sum m*v*v for all degrees of freedom.
	for (int i = 0; i < numParticles; i++) {
		double mass = system.getParticleMass(i);
		energy += mass*vel[i].dot(vel[i]);
	}

	return 0.5*energy;
}

void ReferencePiIntegratorStepKernel::computeForces(ContextImpl& context, const PiIntegrator& integrator) {
	context.updateContextState();
	context.calcForcesAndEnergy(true, false);
}
