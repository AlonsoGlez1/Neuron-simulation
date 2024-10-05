#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

/*  +=========================================================================================================+
    ||                                          DEFINE MACROS                                                ||  
    +=========================================================================================================+  */

#define PI 3.141592654          // Constant for the value of pi
#define DECAYFACTOR 20.0        // Factor for the exponential decay in connection probability
#define THRESHOLD_RADIUS 20.0   // Radius within consider nearby particles
#define PACKING_FRACTION 0.5    // Packing fraction of the system (area occupied / total area)
#define INT_DIST_POWER 2        // Power of the factor of interdistance neurons reduction of probability 
#define TRUE 0
#define FALSE 1


/*  +=========================================================================================================+
    ||                                DEFINE A STRUCTURE TO REPRESENT A NEURON                               ||  
    +=========================================================================================================+  */

typedef struct 
{
    double x;       // x-coordinate of the neuron
    double y;       // y-coordinate of the neuron
    double radius;  // Radius of the neuron
    int synapses;   // Number of synapses (connections) the neuron has
} Particle;


/*  +=========================================================================================================+
    ||                                          FUNCTION PROTOTYPES                                          ||  
    +=========================================================================================================+  */

void packingFraction(int numParticles, double L_x, double L_y, double radius);
void placeRandom(int numParticles, double L_x, double L_y, double radius, Particle *particles);
void placeLattice(int numParticles, double L_x, double L_y, double radius, Particle *particles);
double Distance(Particle p1, Particle p2);
int isOverlapping(Particle *particles, int numParticles, Particle newParticle);
float setInside(float max, float min);
int isCloser(Particle p1, Particle p2, Particle p3);
double connectionProbability(double dist, double maxDistance, double intDistFactor);
double interDistance(Particle p1, Particle p2, Particle p3);
void findNearbyParticles(Particle *particles, int numParticles, int ***nearbyParticles, int **nearbyCounts);
void freeNearbyParticles(int **nearbyParticles, int numParticles, int *nearbyCounts);


/*  +=========================================================================================================+
    ||                                             MAIN FUNCTION                                             ||  
    +=========================================================================================================+  */

int main() {
    int numParticles, timeSteps, N_synapses;          // Number of particles, timeSteps for connections and synapses
    double L_x, L_y, radius, totalDistance;               // Physical dimensions
    char modo, particleFile[100], connectionFile[100];    // Mode (random|lattice) and filenames


    // Seed the random number generator
    srand(time(NULL));

    // Read the number of neurons, the dimensions of the box and neuron radii
    printf("Random or Lattice mode (R/L): ");       scanf("%c", &modo);
    printf("Enter the number of neurons: ");        scanf("%d", &numParticles);
    printf("Enter the width (L_x) of the box: ");   scanf("%lf", &L_x);
    printf("Enter the height (L_y) of the box: ");  scanf("%lf", &L_y);
    printf("Enter the radius of the neurons: ");    scanf("%lf", &radius);

    // Check the packing fraction, if too dense, it cannot locate all neurons
    packingFraction(numParticles, L_x, L_y, radius);


    //Read timeSteps, N_synapses and file names
    printf("Enter the time steps: ");                                 scanf("%d", &timeSteps);
    printf("Enter the number N of synapses per neuron: ");          scanf("%d", &N_synapses);
    printf("Enter the filename to save neuron data (.dat): ");      scanf(" %s", particleFile);
    printf("Enter the filename to save connection data (.dat): ");  scanf(" %s", connectionFile);


    // Allocate memory for the neurons
    Particle *particles = (Particle *)malloc(numParticles * sizeof(Particle));
    if (particles == NULL) 
    {
        fprintf(stderr, "Memory allocation failed\n");
        free(particles);
        return EXIT_FAILURE;
    }

/*  +=========================================================================================================+
    ||                                      INITIALIZE NEURON POSITIONS                                      ||  
    +=========================================================================================================+  */

    if (modo == 'R' || modo == 'r') // Random positioning for the neurons
    {
        placeRandom(numParticles, L_x, L_y, radius, particles);
    }
    else if (modo == 'L' || modo == 'l') // Lattice positioning for the neurons
    {
        placeLattice(numParticles, L_x, L_y, radius, particles);
    }
    else
    {
        fprintf(stderr, "\n Not a valid mode \n");
        free(particles);
        return EXIT_FAILURE;
    }


    // Open the file for writing neuron positions
    FILE *particleFilePtr = fopen(particleFile, "w");
    if (particleFilePtr == NULL) 
    {
        fprintf(stderr, "\n Failed to open file for writing \n");
        free(particles);
        return EXIT_FAILURE;
    }
    fprintf(particleFilePtr, "xlabel ylabel radii \n");

    // Write the neuron positions and radii to the file
    for (int i = 0; i < numParticles; ++i)
    {
        fprintf(particleFilePtr, "%lf %lf %lf\n", particles[i].x, particles[i].y, particles[i].radius);
    }
    // Close the file
    fclose(particleFilePtr);


/*  +=========================================================================================================+
    ||                                     INITIALIZE NEURON CONNECTIONS                                     ||  
    +=========================================================================================================+  */

    // Open the file for writing neuron connections
    FILE *connectionFilePtr = fopen(connectionFile, "w");
    if (connectionFilePtr == NULL) 
    {
        fprintf(stderr, "\nFailed to open connection data file for writing\n");
        free(particles);
        return 1;
    }
    fprintf(connectionFilePtr, "X1 Y1 X2 Y2\n");


    // Create an array of arrays (list of nearby neurons)
    int **nearbyParticles;
    int *nearbyCounts;
    findNearbyParticles(particles, numParticles, &nearbyParticles, &nearbyCounts);


    int currentParticle = rand() % (numParticles + 1);                        // Initialize the first neuron randomly
    double maxDistance = sqrt(pow(L_x - 2*radius,2) + pow(L_y - 2*radius,2)); // Maximum distance within the box
    int failedConnectionAttempt = 0;                                          // Counter of connection tried and failed


    // Start trying connections
    for (int i = 0; i < timeSteps; ++i)
    {
        int nextParticle = -1;

        // Attempt to find a valid next neuron from nearby neurons
        for (int attempt = 0; attempt < nearbyCounts[currentParticle]; ++attempt) 
        {
            int candidateParticle = nearbyParticles[currentParticle][rand() % (nearbyCounts[currentParticle] + 1)];

            // Ensure the candidate has room for new synapses 
            if (particles[candidateParticle].synapses < N_synapses - 1 
                && particles[currentParticle].synapses < N_synapses) 
            {
                double intDistFactor = 1.0; // Reduction factor for intermediate neurons

                // Loop through common nearby neurons for interdistance calculation
                for(int j = 0; j < nearbyCounts[currentParticle]; j++)
                {
                    int intermediateParticle = nearbyParticles[currentParticle][j];

                    // Check if the neuron is also nearby to the candidate
                    for(int k = 0; k < nearbyCounts[candidateParticle]; k++)
                    {
                        if(intermediateParticle == nearbyParticles[candidateParticle][k] 
                           && intermediateParticle != currentParticle)
                        {
                            // Calculate the interdistance
                            double interdistance = interDistance(particles[currentParticle], particles[candidateParticle], particles[j]);

                            // Adjust the connection probability based on interdistance
                            if(interdistance < radius)
                            {
                                if (interdistance <= (0.5 * radius)
                                    && isCloser(particles[currentParticle], particles[intermediateParticle], particles[candidateParticle])
                                    && particles[intermediateParticle].synapses < N_synapses - 1) 
                                {
                                    // The intermediate neuron becomes the new candidate
                                    candidateParticle = intermediateParticle;
                                    intDistFactor = 1.0;
                                    j = 0; 
                                    k = 0; 
                                }
                                else
                                {
                                    intDistFactor *= pow(interdistance / radius , INT_DIST_POWER);
                                }
                            }
                        }
                    }
                }
                

                // Probability check for connection based on distance
                double dist = Distance(particles[currentParticle], particles[candidateParticle]);
                if (connectionProbability(dist, maxDistance, intDistFactor) > (float)rand() / RAND_MAX)
                { 
                    nextParticle = candidateParticle;
                    break;
                }
            }
        }


        // If a valid next neuron was found, create a connection
        if (nextParticle != -1)
        {   
            totalDistance += Distance(particles[currentParticle], particles[nextParticle]);

            fprintf(connectionFilePtr, "%lf %lf %lf %lf\n", 
                    particles[currentParticle].x, particles[currentParticle].y,
                    particles[nextParticle].x, particles[nextParticle].y);
                
            // Add a synapse if the connection is made 
            particles[nextParticle].synapses += 1; 
            particles[currentParticle].synapses += 1;
            
            // Move to the next neuron
            printf("Connection made from: %d\n", currentParticle);
            currentParticle = nextParticle; 
        } 
        else 
        {
            // If no valid next neuron is found, try again
            fprintf(stderr, "No valid connection found from particle %d\n", currentParticle);
            failedConnectionAttempt += 1;
        }
    }
    // Print the end-to-end distance and close the file
    fprintf(connectionFilePtr, "  \n \n Total Distance: %.2lf \n", totalDistance); fclose(connectionFilePtr);

    // Free allocated memory
    free(particles);
    freeNearbyParticles(nearbyParticles, numParticles, nearbyCounts);

    printf("Neuron data saved to %s and connection data saved to %s\n", particleFile, connectionFile);
    
    return EXIT_SUCCESS;
}

/*  +=========================================================================================================+
    ||                                               FUNCTIONS                                               ||  
    +=========================================================================================================+  */


/**************************************************************************************************************
* @brief Computes the packing fraction based on the number of neurons, box dimensions and neuron radius.
*
* @param numParticles The total number of particles.
* @param L_x          The width of the box.
* @param L_y          The height of the box.
* @param radius       The radius of each particle.
*
* @return If the packing fraction is too big, the program ends.
***************************************************************************************************************/
void packingFraction(int numParticles, double L_x, double L_y, double radius)
{
    double packingfraction = (numParticles * PI * pow(radius, 2)) / (L_x * L_y);
    if (packingfraction >= PACKING_FRACTION)
    {
        printf("Packing fraction: %.5f is too dense. \n", packingfraction);
        exit(EXIT_FAILURE);
    }
    printf("Packing fraction: %.5f \n", packingfraction);
}



/**************************************************************************************************************  
* @brief Places neurons inside a box in random positions avoiding overlaping.
*
* @param numParticles The total number of particles.
* @param L_x          The width of the box.
* @param L_y          The height of the box.
* @param radius       The radius of each particle.
* @param particles    Array of structures of particles.   
**************************************************************************************************************/
void placeRandom(int numParticles, double L_x, double L_y, double radius, Particle *particles)
{
    int placedParticles = 0;
        while (placedParticles < numParticles) 
        {
            Particle newParticle;
            newParticle.x = setInside(L_x,radius);
            newParticle.y = setInside(L_y,radius);
            newParticle.radius = radius;
            newParticle.synapses = 0;

            //Check for overlap with existing neurons
            if (isOverlapping(particles, placedParticles, newParticle)) 
            {
                particles[placedParticles] = newParticle;
                placedParticles++;
            }
        }
}




/**************************************************************************************************************  
* @brief Places neurons inside a box in random positions avoiding overlaping.
*
* @param numParticles The total number of particles.
* @param L_x          The width of the box.
* @param L_y          The height of the box.
* @param radius       The radius of each particle.
* @param particles    Array of structures of particles.   
**************************************************************************************************************/
void placeLattice(int numParticles, double L_x, double L_y, double radius, Particle *particles)
{
    int numRows = (int)sqrt(numParticles);
        int numCols = (numParticles + numRows - 1) / numRows; // Handles non-perfect squares

        double dx = (L_x - 2 * radius) / (numCols - 1);       // Spacing in the x direction
        double dy = (L_y - 2 * radius) / (numRows - 1);       // Spacing in the y direction

        int placedParticles = 0;
        for (int i = 0; i < numRows; ++i) 
        {
            for (int j = 0; j < numCols; ++j) 
            {
                if (placedParticles < numParticles) 
                {
                    particles[placedParticles].x = radius + j * dx;
                    particles[placedParticles].y = radius + i * dy;
                    particles[placedParticles].radius = radius;
                    placedParticles++;
                }
            }
        }
}



/************************************************************************************************************** 
* @brief Computes the Euclidean distance between two neurons.
*
* @param p1 The first neuron.
* @param p2 The second neuron.
*
* @return The Euclidean distance between two neurons.
**************************************************************************************************************/
double Distance(Particle p1, Particle p2)
{
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}



/**************************************************************************************************************  
* @brief Computes and compares the distance between two sets of neurons.
*
* @param p1 The current neuron.
* @param p2 A possible intermediate neuron.
* @param p3 The current candidate neuron.
*
* @returns True (0) if p3 is closer to p1 than p3, else False (1).
**************************************************************************************************************/
int isCloser(Particle p1, Particle p2, Particle p3)
{
    if (Distance(p1, p2) < Distance(p1, p3));
    {
        return TRUE;
    }
    return FALSE;
}



/**************************************************************************************************************  
* @brief Checks if a new neuron overlaps with any of the existing neurons.
*
* @param particles    Array of structures of existing particles.
* @param numParticles The number of existing particles.
* @param newParticle  The new particle to check for overlap.
*
* @return True (0) if there is no overlap, false (1) if there is overlap.
**************************************************************************************************************/
int isOverlapping(Particle *particles, int numParticles, Particle newParticle)
{
    for (int i = 0; i < numParticles; ++i) {
        if (Distance(particles[i], newParticle) < 2 * newParticle.radius) {
            return TRUE;
        }
    }
    return FALSE;
}



/**************************************************************************************************************  
* @brief Generates a random position inside the specified boundaries.
*
* @param max The upper limit of the box (either width or height).
* @param min The lower limit for valid placement, typically the radius.
*
* @return A random position within the bounds [min, max - min].
**************************************************************************************************************/
float setInside(float max, float min)
{
    float scale = ((float)rand() / RAND_MAX); // Generate a random float between 0 and 1
    return min + scale * (max - 2 * min); // Scale and shift the value to the desired range
}



/**************************************************************************************************************  
* @brief Calculates the probability of forming a connection between two neurons based on distance.
*
* @param dist          The distance between the two neurons.
* @param maxDistance   The maximum possible distance in the box.
* @param intDistFactor A factor based on intermediate distances between neurons.
*
* @return A value between 0 and 1, where larger distances reduces the connection probability.
**************************************************************************************************************/
double connectionProbability(double dist, double maxDistance, double intDistFactor)
{
    if (dist >= maxDistance) 
    {
        return 0.0;
    }
    return (1 - dist/maxDistance) * exp(- DECAYFACTOR * (dist/maxDistance)) * intDistFactor;
}



/**************************************************************************************************************  
* @brief Computes the distance between a neuron and the line that connects a pair of neurons.
*
* @param p1 The first neuron, it connects to the second one
* @param p2 The second neuron, it connects to the first one.
* @param p3 The third neuron.
*
* @return The calculated interdistance between the three neurons.
**************************************************************************************************************/
double interDistance(Particle p1, Particle p2, Particle p3)
{
    double lineLength = Distance(p1 , p2);
   
    if (lineLength < 1e-9) 
    {
        return Distance(p1, p3);
    }
    return fabs((p2.y - p1.y) * p3.x - (p2.x - p1.x) * p3.y + p2.x * p1.y - p2.y * p1.x) / lineLength;
}



/**************************************************************************************************************
* @brief Generates a list of nearby neurons for each neuron, based on a threshold radius.
*
* @param particles       Array of particles.
* @param numParticles    The total number of particles.
* @param nearbyParticles A pointer to an array of lists, each containing the indices of nearby neurons for each neuron.
* @param nearbyCounts    A pointer to an array storing the number of nearby particles for each particle.
**************************************************************************************************************/
void findNearbyParticles(Particle *particles, int numParticles, int ***nearbyParticles, int **nearbyCounts) {
    // Allocate memory for nearby particles list and counts
    *nearbyParticles = (int **)malloc(numParticles * sizeof(int *));
    *nearbyCounts = (int *)malloc(numParticles * sizeof(int));

    for (int i = 0; i < numParticles; ++i) {
        (*nearbyCounts)[i] = 0;
        (*nearbyParticles)[i] = (int *)malloc(numParticles * sizeof(int)); // allocate max size initially
    }

    // Populate nearby particles list for each particle
    for (int i = 0; i < numParticles; ++i) {
        for (int j = 0; j < numParticles; ++j) {
            if (i != j) {
                double dist = Distance(particles[i], particles[j]);
                if (dist <= THRESHOLD_RADIUS) {
                    (*nearbyParticles)[i][(*nearbyCounts)[i]] = j; // Store index of nearby particle
                    (*nearbyCounts)[i]++;
                }
            }
        }
    }
}



/**************************************************************************************************************
* @brief Frees the memory allocated for the nearby particles lists.
*
* @param nearbyParticles The list of nearby particles for each particle.
* @param numParticles    The total number of particles.
* @param nearbyCounts    The array storing the count of nearby particles for each particle.
**************************************************************************************************************/
void freeNearbyParticles(int **nearbyParticles, int numParticles, int *nearbyCounts) {
    for (int i = 0; i < numParticles; ++i) {
        free(nearbyParticles[i]);
    }
    free(nearbyParticles);
    free(nearbyCounts);
}