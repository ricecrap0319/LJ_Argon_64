#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

const int N = 64;                 // Number of Argon molecules
const double dt = 0.0025;       // Time step for the simulation
const double boxLength = 4.0;     // Adjusted to fit a 4x4x4 lattice
const double epsilon = 1.0;
const double sigma = 1.0;
const double mass = 1.0;          // Mass of each Argon atom
const double kb = 1.0;            // Boltzmann constant
const double temp = 1.0;          // Initial temperature
const int numSteps = 2000000;      // Total number of simulation steps
const int rescaleFreq = 100;      // Frequency of velocity rescaling
const double targetTemperature = 1.0; // Target temperature for rescaling

struct Argon {
    std::vector<double> v;  // Velocity components {vx, vy, vz}
    std::vector<double> r;  // Position components {rx, ry, rz}
    std::vector<double> a;  // Acceleration components {ax, ay, az}
};

std::vector<double> calcForce(const std::vector<double>& r1, const std::vector<double>& r2) {
    std::vector<double> force(3);
    std::vector<double> dr(3);
    double r = 0.0;

    for (int i = 0; i < 3; i++) {
        dr[i] = r1[i] - r2[i];
        // Apply periodic boundary conditions
        dr[i] -= boxLength * round(dr[i] / boxLength);
    }

    r = sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
    if (r == 0) return force;  // Avoid division by zero

    double forceMagnitude = 24 * epsilon * (2 * pow(sigma / r, 13) - pow(sigma / r, 7)) ;
    for (int i = 0; i < 3; i++) {
        force[i] = forceMagnitude * dr[i] / r; // Up for debate
    }
    return force;
}

double calcPotentialEnergy(const std::vector<Argon>& molecules) {
    double totalPotentialEnergy = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            std::vector<double> dr(3);
            double r = 0.0;
            for (int k = 0; k < 3; ++k) {
                dr[k] = molecules[i].r[k] - molecules[j].r[k];
                dr[k] -= boxLength * round(dr[k] / boxLength);
                r += dr[k] * dr[k];
            }
            r = sqrt(r);
            if (r != 0) {
                totalPotentialEnergy += 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));
            }
        }
    }
    return totalPotentialEnergy;
}

double calcKineticEnergy(const std::vector<Argon>& molecules) {
    double totalKineticEnergy = 0.0;
    for (const auto& molecule : molecules) {
        double v_squared = 0.0;
        for (int i = 0; i < 3; i++) {
            v_squared += molecule.v[i] * molecule.v[i];
        }
        totalKineticEnergy += 0.5 * mass * v_squared;
    }
    return totalKineticEnergy;
}

double calcCurrentKineticTemperature(const std::vector<Argon>& molecules) {
    double totalKineticEnergy = calcKineticEnergy(molecules);
    return (2.0 / (3.0 * N * kb)) * totalKineticEnergy;
}

void rescaleVelocities(std::vector<Argon>& molecules, double targetTemp, double currentTemp) {
    double scalingFactor = sqrt(targetTemp / currentTemp);
    for (auto& molecule : molecules) {
        for (int i = 0; i < 3; i++) {
            molecule.v[i] *= scalingFactor;
        }
    }
}

int main() {
    std::vector<Argon> molecules(N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0, sqrt(kb * temp / mass));

    // Initialize positions in a 4x4x4 lattice and velocities
    for (int i = 0; i < N; i++) {
        molecules[i].r = {double(i % 4), double((i / 4) % 4), double(i / 16)};
        molecules[i].v = {dist(gen), dist(gen), dist(gen)};
        molecules[i].a = {0, 0, 0};
    }
 

    // Start of velocity verlet 
    for (int step = 0; step < numSteps; step++) {
        // Calculate forces based on current positions
        for (int i = 0; i < N; i++) {
            std::vector<double> totalForce(3, 0.0);
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    std::vector<double> force = calcForce(molecules[i].r, molecules[j].r);
                    for (int k = 0; k < 3; k++) {
                        totalForce[k] += force[k];
                    }
                }
            }
            molecules[i].a = totalForce;
        }
    //First update of velocity 
    for (int i = 0; i < N; i++) {
            for (int k = 0; k < 3; k++) {
                molecules[i].v[k] += 0.5 * (molecules[i].a[k]) * dt;
            }
        }

        // Update positions
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < 3; k++) {
                molecules[i].r[k] += molecules[i].v[k] * dt;
                // Apply periodic boundary conditions
                molecules[i].r[k] = fmod(molecules[i].r[k] + boxLength, boxLength);
            }
        }

        // Calculate new forces based on updated positions
        for (int i = 0; i < N; i++) {
            std::vector<double> totalForce(3, 0.0);
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    std::vector<double> force = calcForce(molecules[i].r, molecules[j].r);
                    for (int k = 0; k < 3; k++) {
                        totalForce[k] += force[k];
                    }
                }
            }
            molecules[i].a = totalForce;
        }

        // Update velocities using the new accelerations
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < 3; k++) {
                molecules[i].v[k] += 0.5 * (molecules[i].a[k]) * dt;
            }
        }

        // Calculate energies and temperature
        double totalKineticEnergy = calcKineticEnergy(molecules);
        double totalPotentialEnergy = calcPotentialEnergy(molecules);
        double totalEnergy = totalKineticEnergy + totalPotentialEnergy;
        double currentTemperature = (2.0 / (3.0 * N * kb)) * totalKineticEnergy;

        // Output energies and temperature every 10000 steps
        if (step % 1000 == 0) {
            std::cout << "Step " << step << ":" << std::endl;
            std::cout << "Total Kinetic Energy: " << totalKineticEnergy << " Joules" << std::endl;
            std::cout << "Total Potential Energy: " << totalPotentialEnergy << " Joules" << std::endl;
            std::cout << "Total Energy: " << totalEnergy << " Joules" << std::endl;
            std::cout << "Current Temperature: " << currentTemperature << " Kelvin" << std::endl;
            std::cout << "----------------------------------------" << std::endl;
        }

        // Rescale velocities if needed
        if ((step + 1) % rescaleFreq == 0 && step < 100) {
            rescaleVelocities(molecules, targetTemperature, currentTemperature);
        }
    }

    return 0;
}
