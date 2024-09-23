#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


// Define a structure to represent a neuron
typedef struct 
{
    double x;
    double y;
    double radius;
    int synapses;
} Particle;

//Function prototypes
double distance(Particle p1, Particle p2);
int isOverlapping(Particle *particles, int numParticles, Particle newParticle);
float isInside(float max, float min);
double connectionProbability(double dist, double maxDistance, double intDistFactor);
double interDistance(Particle p1, Particle p2, Particle p3);


int main() {
    int numParticles, N_connections, N_synapses;
    double L_x, L_y, radius;
    char particleFile[100], connectionFile[100], modo;

    // Seed the random number generator
    srand(time(NULL));

    // Read the number of neurons, and the dimensions of the box, radius, connections and max number of synapses
    printf("Random or Lattice (R/L): "); scanf("%c", &modo);
    printf("Enter the number of neurons: "); scanf("%d", &numParticles);
    printf("Enter the width (L_x) of the box: "); scanf("%lf", &L_x);
    printf("Enter the height (L_y) of the box: "); scanf("%lf", &L_y);
    printf("Enter the radius of the neurons: "); scanf("%lf", &radius);
    printf("Enter the number N of connections: "); scanf("%d", &N_connections);
    printf("Enter the number N of synapses per neuron: "); scanf("%d", &N_synapses);
    printf("Enter the filename to save neuron data (.dat): "); scanf(" %s", particleFile);
    printf("Enter the filename to save connection data (.dat): "); scanf(" %s", connectionFile);

    // Allocate memory for the neurons
    Particle *particles = (Particle *)malloc(numParticles * sizeof(Particle));
    if (particles == NULL) 
    {
        fprintf(stderr, "Memory allocation failed\n");
        free(particles);
        return 1;
    }

    //Random positioning or Lattice
    if (modo == 'R' || modo == 'r')
    {
        // Initialize neuron positions randomly within the box without overlap
        int placedParticles = 0;
        while (placedParticles < numParticles) 
        {
            Particle newParticle;
            newParticle.x = isInside(L_x,radius);
            newParticle.y = isInside(L_y,radius);
            newParticle.radius = radius;

            if (!isOverlapping(particles, placedParticles, newParticle)) {
                particles[placedParticles] = newParticle;
                placedParticles++;
            }
        }
    }
    else if (modo == 'L' || modo == 'l')
    {
        // Initialize particle positions within the box without overlap
        int numRows = (int)sqrt(numParticles);
        int numCols = (numParticles + numRows - 1) / numRows; // This handles cases where numParticles is not a perfect square

        double dx = (L_x - 2 * radius) / (numCols - 1); // Spacing between particles in the x direction
        double dy = (L_y - 2 * radius) / (numRows - 1); // Spacing between particles in the y direction

        int placedParticles = 0;
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                if (placedParticles < numParticles) {
                    particles[placedParticles].x = radius + j * dx;
                    particles[placedParticles].y = radius + i * dy;
                    particles[placedParticles].radius = radius;
                    placedParticles++;
                }
            }
        }

    }
    else
    {
        printf("\n Not a valid mode \n");
        return 1;
    }

    // Open the file for writing
    FILE *particleFilePtr = fopen(particleFile, "w");
    if (particleFilePtr == NULL) 
    {
        fprintf(stderr, "Failed to open file for writing\n");
        free(particles);
        return 1;
    }

    // Write the neuron positions and radii to the file
    for (int i = 0; i < numParticles; ++i)
    {
        fprintf(particleFilePtr, "%lf %lf %lf\n", particles[i].x, particles[i].y, particles[i].radius);
    }

    // Close the file
    fclose(particleFilePtr);

    // Open the file for writing
    FILE *connectionFilePtr = fopen(connectionFile, "w");
    if (connectionFilePtr == NULL) 
    {
        fprintf(stderr, "Failed to open connection data file for writing\n");
        free(particles);
        return 1;
    }

    // Initialize current neuron
    int currentParticle = rand() % numParticles;
    double maxDistance = sqrt(pow(L_x - 2*radius,2) + pow(L_y - 2*radius,2)); // Maximum possible distance within the box

    for (int i = 0; i < N_connections; ++i)
    {
        int nextParticle = -1;

        // Attempt to find a valid next neuron
        for (int attempt = 0; attempt < 10*numParticles; ++attempt) 
        {
            int candidateParticle = rand() % numParticles;

            // Secures the next neuron can look for a new synapse 
            if (candidateParticle != currentParticle 
                    && particles[candidateParticle].synapses < N_synapses - 1
                        && particles[currentParticle].synapses < N_synapses) 
            {
                
                //Reduced probability by intermediate neurons
                double intDistFactor = 1.0;
                for(int j = 0; j < numParticles; j++)
                {
                    if (j != currentParticle && j != candidateParticle)
                    {
                        double interdistance = interDistance(particles[currentParticle], particles[candidateParticle], particles[j]);

                        if(interdistance < radius)
                        {
                            if (interdistance > 0.5 * radius)
                            {
                             intDistFactor *= interdistance / radius;
                            }
                            else if (particles[j].synapses < N_synapses - 1
                                && distance(particles[currentParticle], particles[j]) < distance(particles[currentParticle], particles[candidateParticle])) 
                            {
                                intDistFactor = 1.0;
                                candidateParticle = j;
                            }
                        }
                    }
                }
                
                double dist = distance(particles[currentParticle], particles[candidateParticle]);
                double prob = connectionProbability(dist, maxDistance, intDistFactor);
                
                // Use probability to decide whether to connect
                if (((double)rand() / RAND_MAX) < prob)
                { 
                            nextParticle = candidateParticle;
                            break;
                }
            }
        }

        // If a valid next neuron was found, create a connection
        if (nextParticle != -1)
        {
            fprintf(connectionFilePtr, "%lf %lf %lf %lf\n",
                    particles[currentParticle].x, particles[currentParticle].y,
                    particles[nextParticle].x, particles[nextParticle].y);
             
            // Add a synapse if the connection is made 
            particles[nextParticle].synapses += 1; 
            particles[currentParticle].synapses += 1;
            
            // Move to the next neuron
            currentParticle = nextParticle; 
        } 
        else 
        {
            // If no valid next neuron is found, stop connecting
            fprintf(stderr, "No valid connection found from particle %d\n", currentParticle);
            N_connections += 1;
            //break;
        }
    }

    // Close the file
    fclose(connectionFilePtr);
    
    // Free the allocated memory
    free(particles);
    printf("Neuron data saved to %s and connection data saved to %s\n", particleFile, connectionFile);

    return 0;
}


// Function to calculate the distance between two neurons
double distance(Particle p1, Particle p2)
{
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

// Function to check if a neuron overlaps with any existing neuron
int isOverlapping(Particle *particles, int numParticles, Particle newParticle)
{
    for (int i = 0; i < numParticles; ++i) {
        if (distance(particles[i], newParticle) < 2 * newParticle.radius) {
            return 1;
        }
    }
    return 0;
}

//Function to set the random neurons without overlaping within the walls
float isInside(float max, float min)
{
    float scale = ((float)rand() / RAND_MAX); // Generate a random float between 0 and 1
    return min + scale * (max - 2 * min); // Scale and shift the value to the desired range
}

// Function to determine connection probability based on distance
double connectionProbability(double dist, double maxDistance, double intDistFactor)
{
    double decayFactor = 10.0;
    if (dist >= maxDistance) 
    {
        return 0.0;
    }
    return (1 - dist/maxDistance) * exp(- decayFactor * (dist/maxDistance)) * intDistFactor;
}

//Function to calculate the distance between a neuron and the line connecting a pair of neurons
double interDistance(Particle p1, Particle p2, Particle p3)
{
    double lineLength = sqrt(pow(p2.y - p1.y, 2) + pow(p2.x - p1.x, 2));
    
    if (lineLength < 1e-9) // Avoid floating-point issues
    {
        return distance(p1, p3); // If p1 and p2 are essentially the same point, return distance to p1
    }

    return fabs((p2.y - p1.y) * p3.x - (p2.x - p1.x) * p3.y + p2.x * p1.y - p2.y * p1.x) / lineLength;
}

