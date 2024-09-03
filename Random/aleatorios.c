#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// Define a structure to represent a particle
typedef struct {
    double x;
    double y;
    double radius;
    int synapses;
} Particle;

// Function to calculate the distance between two particles
double distance(Particle p1, Particle p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

// Function to check if a particle overlaps with any existing particles
int isOverlapping(Particle *particles, int numParticles, Particle newParticle) {
    for (int i = 0; i < numParticles; ++i) {
        if (distance(particles[i], newParticle) < 2 * newParticle.radius) {
            return 1;
        }
    }
    return 0;
}

//Function to set the random particles without overlaping within the walls
float isInside(float max, float min){
    float scale = ((float)rand() / RAND_MAX); // Generate a random float between 0 and 1
    return min + scale * (max - 2 * min); // Scale and shift the value to the desired range
}

// Function to determine connection probability based on distance
double connectionProbability(double dist, double maxDistance) {
    if (dist >= maxDistance) {
        return 0.0;
    }
    return 1.0 - (dist / maxDistance);
}

int main() {
    int numParticles, N_connections, N_synapses;
    double L_x, L_y, radius;
    char particleFile[100], connectionFile[100];

    // Seed the random number generator
    srand(time(NULL));

    // Read the number of particles, and the dimensions of the box, radius, connections and max number of synapses
    printf("Enter the number of particles: "); scanf("%d", &numParticles);
    printf("Enter the width (L_x) of the box: "); scanf("%lf", &L_x);
    printf("Enter the height (L_y) of the box: "); scanf("%lf", &L_y);
    printf("Enter the radius of the particles: "); scanf("%lf", &radius);
    printf("Enter the number N of connections: "); scanf("%d", &N_connections);
    printf("Enter the number N of synapses per particle: "); scanf("%d", &N_synapses);
    printf("Enter the filename to save particle data (.dat): "); scanf("%s", particleFile);
    printf("Enter the filename to save connection data (connections.dat): "); scanf("%s", connectionFile);

    // Allocate memory for the particles
    Particle *particles = (Particle *)malloc(numParticles * sizeof(Particle));
    if (particles == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        
        return 1;
    }

    // Initialize particle positions randomly within the box without overlap
    int placedParticles = 0;
    while (placedParticles < numParticles) {
        Particle newParticle;
        newParticle.x = isInside(L_x,radius);
        newParticle.y = isInside(L_y,radius);
        newParticle.radius = radius;

        if (!isOverlapping(particles, placedParticles, newParticle)) {
            particles[placedParticles] = newParticle;
            placedParticles++;
        }
    }

    // Open the file for writing
    FILE *particleFilePtr = fopen(particleFile, "w");
    if (particleFilePtr == NULL) {
        fprintf(stderr, "Failed to open file for writing\n");
        free(particles);
        
        return 1;
    }

    // Write the particle positions and radii to the file
    for (int i = 0; i < numParticles; ++i) {
        fprintf(particleFilePtr, "%lf %lf %lf\n", particles[i].x, particles[i].y, particles[i].radius);
    }

    // Close the file
    fclose(particleFilePtr);

    // Open the file for writing
    FILE *connectionFilePtr = fopen(connectionFile, "w");
    if (connectionFilePtr == NULL) {
        fprintf(stderr, "Failed to open connection data file for writing\n");
        free(particles);
        return 1;
    }

    // Initialize current particle
    int currentParticle = rand() % numParticles;
    double maxDistance = sqrt(pow(L_x - 2*radius,2) + pow(L_y - 2*radius,2)); // Maximum possible distance within the box

    for (int i = 0; i < N_connections; ++i) {
        int nextParticle = -1;

        // Attempt to find a valid next particle
        for (int attempt = 0; attempt < numParticles; ++attempt) {
            int candidateParticle = rand() % numParticles;

            if (candidateParticle != currentParticle) {
                double dist = distance(particles[currentParticle], particles[candidateParticle]);
                double prob = connectionProbability(dist, maxDistance);
                // int synapse = particles[candidateParticle].synapses; //Check for the number of connections
                
                // Use probability to decide whether to connect
                if (((double)rand() / RAND_MAX) < prob 
                    && particles[candidateParticle].synapses < N_synapses - 1 //Secures the next loop choice 
                        && particles[currentParticle].synapses < N_synapses) { 
                            nextParticle = candidateParticle;
                            break;
                }
            }
        }

        // If a valid next particle was found, create a connection
        if (nextParticle != -1) {
            fprintf(connectionFilePtr, "%lf %lf %lf %lf\n",
                    particles[currentParticle].x, particles[currentParticle].y,
                    particles[nextParticle].x, particles[nextParticle].y);
             
            particles[nextParticle].synapses += 1; // Add a synapse if the connection is made
            particles[currentParticle].synapses += 1;
            currentParticle = nextParticle; // Move to the next particle
        } else {
            // If no valid next particle is found, stop connecting
            fprintf(stderr, "No valid connection found from particle %d\n", currentParticle);
            break;
        }
    }

    // Close the file
    fclose(connectionFilePtr);

    // Free the allocated memory
    free(particles);

    printf("Particle data saved to %s and connection data saved to %s\n", particleFile, connectionFile);

    return 0;
}
