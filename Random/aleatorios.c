#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/*  +=================+
    || Define macros ||  
    +=================+  */
#define PI 3.141592654

/*  +===============================================+
    ||  Define a structure to represent a neuron   ||  
    +===============================================+  */
typedef struct 
{
    double x;
    double y;
    double radius;
    int synapses;
} Particle;

/*  +==========================+
    ||  FUNCTION PROTOTYPES   ||  
    +==========================+  */
double packingFraction(int numParticles, double L_x, double L_y, double radius);
double distance(Particle p1, Particle p2);
int isOverlapping(Particle *particles, int numParticles, Particle newParticle);
float setInside(float max, float min);
double connectionProbability(double dist, double maxDistance, double intDistFactor);
double interDistance(Particle p1, Particle p2, Particle p3);


/*  +====================+
    ||  MAIN FUNCTION   ||  
    +====================+  */
int main() {
    int numParticles, N_connections, N_synapses;
    double L_x, L_y, radius, totalDistance, packingfraction;
    char modo, particleFile[100], connectionFile[100];

    // Seed the random number generator
    srand(time(NULL));

    // Read the number of neurons, and the dimensions of the box, radius, connections and max number of synapses
    printf("Random or Lattice mode (R/L): "); scanf("%c", &modo);
    printf("Enter the number of neurons: "); scanf("%d", &numParticles);
    printf("Enter the width (L_x) of the box: "); scanf("%lf", &L_x);
    printf("Enter the height (L_y) of the box: "); scanf("%lf", &L_y);
    printf("Enter the radius of the neurons: "); scanf("%lf", &radius);

    // Check the packing fraction, if too dense, it cant locate all the neurons
    packingfraction = packingFraction(numParticles, L_x, L_y, radius);
    if (packingfraction >= 0.65)
    {
        printf("Packing fraction: %f is too dense. \n", packingfraction);
        return 1;
    }
    else{
        printf("Packing fraction: %.2f \n", packingfraction);
    }

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

    if (modo == 'R' || modo == 'r') // Random positioning for the neurons
    {
        // Initialize neuron positions randomly within the box without overlap
        int placedParticles = 0;
        while (placedParticles < numParticles) 
        {
            Particle newParticle;
            newParticle.x = setInside(L_x,radius);
            newParticle.y = setInside(L_y,radius);
            newParticle.radius = radius;

            if (isOverlapping(particles, placedParticles, newParticle)) 
            {
                particles[placedParticles] = newParticle;
                placedParticles++;
            }
        }
    }
    else if (modo == 'L' || modo == 'l') // Lattice positioning for the neurons
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

    // Open the file for writing neuron positions and making titles
    FILE *particleFilePtr = fopen(particleFile, "w");
    if (particleFilePtr == NULL) 
    {
        fprintf(stderr, "Failed to open file for writing\n");
        free(particles);
        return 1;
    }
    fprintf(particleFilePtr, "xlabel ylabel radii\n");

    // Write the neuron positions and radii to the file
    for (int i = 0; i < numParticles; ++i)
    {
        fprintf(particleFilePtr, "%lf %lf %lf\n", particles[i].x, particles[i].y, particles[i].radius);
    }
    // Close the file
    fclose(particleFilePtr);

    // Open the file for writing neuron connections and making titles
    FILE *connectionFilePtr = fopen(connectionFile, "w");
    if (connectionFilePtr == NULL) 
    {
        fprintf(stderr, "Failed to open connection data file for writing\n");
        free(particles);
        return 1;
    }
    fprintf(connectionFilePtr, "X1 Y1 X2 Y2\n");

    // Initialize current neuron
    int currentParticle = rand() % numParticles;
    double maxDistance = sqrt(pow(L_x - 2*radius,2) + pow(L_y - 2*radius,2)); // Maximum possible distance within the box

    for (int i = 0; i < N_connections; ++i)
    {
        int nextParticle = -1;

        // Attempt to find a valid next neuron
        for (int attempt = 0; attempt < numParticles; ++attempt) 
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
                        // Calculate the distance between a neuron and the line connecting a pair of neurons
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
                                //If the neuron is practically intersecting, this particle is the new candidate
                                intDistFactor = 1.0;
                                candidateParticle = j;
                            }
                        }
                    }
                }
                
                // Use probability to decide whether to connect
                double dist = distance(particles[currentParticle], particles[candidateParticle]);
                double prob = connectionProbability(dist, maxDistance, intDistFactor);
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
            totalDistance += distance(particles[currentParticle], particles[nextParticle]);

            fprintf(connectionFilePtr, "%lf %lf %lf %lf\n",
                        particles[currentParticle].x, particles[currentParticle].y, particles[nextParticle].x, particles[nextParticle].y);
             
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
        }
    }

    // Print the end-to-end distance to the file and close the file
    fprintf(connectionFilePtr, "  \n \n Total Distance: %.2lf \n", totalDistance); fclose(connectionFilePtr);
    
    // Free the allocated memory
    free(particles);
    printf("Neuron data saved to %s and connection data saved to %s\n", particleFile, connectionFile);

    return 0;
}

/*  +================+
    ||  FUNCTIONS   ||  
    +================+  */

/*  +-----------------------------------------------------------------------------------------+  
    |       FUNCTION TO CALCULATE THE PACKING FRACTION                                        |
    |                                                                                         |
    | The packing fraction is a measure of density in a given space. Its the fraction of the  |
    | volume (or area in 2D) of the space that is occupied by the objects.                    |
    +-----------------------------------------------------------------------------------------+  */ 
double packingFraction(int numParticles, double L_x, double L_y, double radius)
{
    return (numParticles * PI * pow(radius, 2)) / (L_x * L_y);
}
/*  +--------------------------------------------------------------------------------------------+  
    |       FUNCTION TO CALCULATE THE DISTANCE BETWEEN TWO POINTS                                |
    |                                                                                            |
    | Takes (x1 , y1) form particle 1 and (x2 , y2) from particle 2 and calculates its distance. |
    +--------------------------------------------------------------------------------------------+  */ 
double distance(Particle p1, Particle p2)
{
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

/*  +-----------------------------------------------------------------------------------------------+  
    |       FUNCTION TO CHECK IF A NEURON OVERLAPS WITH ANY EXISTING NEURON                         |
    |                                                                                               |
    | For a new particle selected, checks for overlaps with any of the already numParticles placed. |
    | If the distance between them is less than the sum of its radii, they overlap (True = 0).      |
    +-----------------------------------------------------------------------------------------------+  */ 
int isOverlapping(Particle *particles, int numParticles, Particle newParticle)
{
    for (int i = 0; i < numParticles; ++i) {
        if (distance(particles[i], newParticle) < 2 * newParticle.radius) {
            return 0;
        }
    }
    return 1;
}

/*  +-----------------------------------------------------------------------------------------------+
    |       FUNCTION TO SET RANDOM NEURONS WITHOUT OVERLAPING WITHIN THE WALLS                      | 
    |                                                                                               |
    | Takes the range values (maximum and minimum) and generates random points within those values. |
    +-----------------------------------------------------------------------------------------------+   */
float setInside(float max, float min)
{
    float scale = ((float)rand() / RAND_MAX); // Generate a random float between 0 and 1
    return min + scale * (max - 2 * min); // Scale and shift the value to the desired range
}

/*  +--------------------------------------------------------------------------------------------------+
    |       FUNCTION TO DETERMINE CONNECTION PROBABILITY BASED ON DISTANCE                             |
    |                                                                                                  |
    | Takes the distance between the current neuron and the candidate neuron, the maximum distance     |
    | possible within the walls and a reduction factor by intermediate intersecting neurons.           |
    | The maximum probability is 1 and the minimum is 0, reduced by distance and intermediate neurons. |
    | The reduced connection probability is modulated by a decayFactor.                                |
    +--------------------------------------------------------------------------------------------------+  */
double connectionProbability(double dist, double maxDistance, double intDistFactor)
{
    double decayFactor = 10.0;
    if (dist >= maxDistance) 
    {
        return 0.0;
    }
    return (1 - dist/maxDistance) * exp(- decayFactor * (dist/maxDistance)) * intDistFactor;
}

/*  +---------------------------------------------------------------------------------------------------------+
    |       FUNCTION TO CALCULATE THE DISTANCE FROM A NEURON TO THE LINE SEGMENT CONNECTING TWO OTHER NEURONS |
    | Takes the curren (p1), candidate (p2) and any other neuron (p3) placed.                                 |
    | Calculates first the distance between p1 and p2 to avoid floating-point issues (division by 0).         |
    | If p1 and p2 are essentially the same point, return distance to p1. Else, do the calculation.           | 
    +---------------------------------------------------------------------------------------------------------+  */
double interDistance(Particle p1, Particle p2, Particle p3)
{
    double lineLength = distance(p1 , p2);
   
    if (lineLength < 1e-9) 
    {
        return distance(p1, p3);
    }
    return fabs((p2.y - p1.y) * p3.x - (p2.x - p1.x) * p3.y + p2.x * p1.y - p2.y * p1.x) / lineLength;
}