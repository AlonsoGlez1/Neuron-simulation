#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define PI 3.141592654       // Constant for the value of pi
#define DECAYFACTOR 15.0     // Factor for the exponential decay in connection probability
#define CUTOFF_RADIUS 25.0   // Radius within consider nearby neurons
#define PACKING_FRACTION 0.5 // Packing fraction of the system (area occupied / total area)
#define INT_DIST_POWER 2     // Power of the interdistance factor reduction
#define BRANCHES 3           // How many branches per initial neuron
#define TRUE 1
#define FALSE 0


// Define a structure to represent a neuron
typedef struct
{
    double x;       // x-coordinate of the neuron
    double y;       // y-coordinate of the neuron
    double radius;  // Radius of the neuron
    int synapses;   // Number of synapses (connections)
} Neuron;

/*
+=========================================================================================================+
||                                          FUNCTION PROTOTYPES                                          ||
+=========================================================================================================+  
*/

void packingFraction(int N_neurons, double L_x, double L_y, double radius);
void placeRandom(int N_neurons, double L_x, double L_y, double radius, Neuron *neurons);
void placeLattice(int N_neurons, double L_x, double L_y, double radius, Neuron *neurons);
double Distance(Neuron p1, Neuron p2);
int isOverlapping(Neuron *neurons, int N_neurons, Neuron newNeuron);
float setInside(float max, float min);
int isBetween(Neuron p1, Neuron p2, Neuron p3);
double connectionProbability(double dist, double maxDistance, double intDistFactor);
double interDistance(Neuron p1, Neuron p2, Neuron p3);
void findNearbyParticles(Neuron *neurons, int N_neurons, int ***nearbyNeurons, int **nearbyCounts);
void freeNearbyParticles(int **nearbyNeurons, int N_neurons, int *nearbyCounts);


int main()
{
   int N_neurons, timeSteps, N_synapses;      // Number of neurons, timeSteps for connections and synapses
   double L_x, L_y, radius;                   // Physical dimensions
   char neuronFile[100], connectionFile[100]; // File names
   char modo;                                 // Mode (random|lattice)


   // Seed the random number generator
   srand(time(NULL));


   // Read the number of neurons, the dimensions of the box and neuron radii
   printf("Random or Lattice mode (R/L): ");      scanf(" %c", &modo);
   printf("Enter the number of neurons: ");       scanf("%d", &N_neurons);
   printf("Enter the width (L_x) of the box: ");  scanf("%lf", &L_x);
   printf("Enter the height (L_y) of the box: "); scanf("%lf", &L_y);
   printf("Enter the radius of the neurons: ");   scanf("%lf", &radius);


   // Check the packing fraction, if it's too dense, it cannot allocate all the neurons
   packingFraction(N_neurons, L_x, L_y, radius);


   //Read timeSteps, N_synapses and file names
   printf("Enter the time steps: ");                              scanf("%d", &timeSteps);
   printf("Enter the max number of synapses per neuron: ");       scanf("%d", &N_synapses);
   printf("Enter the filename to save neuron data (.dat): ");     scanf(" %99s", neuronFile);
   printf("Enter the filename to save connection data (.dat): "); scanf(" %99s", connectionFile);


   // Append "_dat" to the user input
    strcat(neuronFile, "_dat"); 
    strcat(connectionFile, "_dat"); 


   // Allocate memory for the neurons
   Neuron *neurons = (Neuron *)malloc(N_neurons * sizeof(Neuron));
   if(neurons == NULL)
   {
      fprintf(stderr, "Memory allocation failed\n");
      return EXIT_FAILURE;
   }

/*
+=========================================================================================================+
||                                      INITIALIZE NEURON POSITIONS                                      ||
+=========================================================================================================+
*/

   if (modo == 'R' || modo == 'r')
   {
      placeRandom(N_neurons, L_x, L_y, radius, neurons); // Random positioning for the neurons
   }
   else if(modo == 'L' || modo == 'l')
   {
      placeLattice(N_neurons, L_x, L_y, radius, neurons); // Lattice positioning for the neurons
   }
   else
   {
      fprintf(stderr, "\nNot a valid mode\n");
      free(neurons);
      return EXIT_FAILURE;
   }


   // Open the file for writing neuron positions
   FILE *neuronFilePtr = fopen(neuronFile, "w");
   if(neuronFilePtr == NULL)
   {
      fprintf(stderr, "\nFailed to open neuron file for writing\n");
      free(neurons);
      return EXIT_FAILURE;
   }

   // Write the neuron positions and radii to the file
   fprintf(neuronFilePtr, "'Xlabel Ylabel Radius'\n");
   for(int i = 0; i < N_neurons; ++i)
   {
      fprintf(neuronFilePtr, "%lf %lf %lf\n", neurons[i].x, neurons[i].y, neurons[i].radius);
   }
   fclose(neuronFilePtr);


/* 
+=========================================================================================================+
||                                     INITIALIZE NEURON CONNECTIONS                                     ||
+=========================================================================================================+
*/

   // Open the file for writing neuron connections
   FILE *connectionFilePtr = fopen(connectionFile, "w");
   if (connectionFilePtr == NULL)
   {
      fprintf(stderr, "\nFailed to open connection data file for writing\n");
      free(neurons);
      return 1;
   }


   // Create an array of arrays (list of nearby neurons)
   int **nearbyNeurons;
   int *nearbyCounts;
   findNearbyParticles(neurons, N_neurons, &nearbyNeurons, &nearbyCounts);


   int initialNeuron = rand() % (N_neurons + 1); // Initialize randomly the first neuron
   double maxDistance = CUTOFF_RADIUS;           // Maximum distance possible

   for(int branch = 0; branch < BRANCHES; branch++)
   {
      // Restart values at the beginning of each branching
      int currentNeuron = initialNeuron;
      int failedConnectionAttempt = 0;
      double totalDistance = 0.0;
      double end_to_end_Distance = 0.0;
      int N_connections = 0;

      // Print the initial neuron in the connection file
      fprintf(connectionFilePtr, "'X Y Radius'\n");
      fprintf(connectionFilePtr, "%lf %lf %lf\n", neurons[initialNeuron].x, neurons[initialNeuron].y, neurons[initialNeuron].radius);

      // Start connection attempts
      for(int i = 0; i < timeSteps; ++i)
      {
         int nextNeuron = -1;

         // Attempt to find a next valid neuron from nearby neurons list
         for(int attempt = 0; attempt < nearbyCounts[currentNeuron]; ++attempt)
         {
            int candidateNeuron = nearbyNeurons[currentNeuron][rand() % (nearbyCounts[currentNeuron] + 1)];

            // Ensure the candidate has room for more synapses AND TEMPORAL: NOT RECONNECT TO THE INITIAL NEURON
            if(candidateNeuron != initialNeuron
               && neurons[candidateNeuron].synapses < N_synapses - 1
               && neurons[currentNeuron].synapses < N_synapses)
            {
               double intDistFactor = 1.0; // Reduction factor by intermediate neurons

               // Loop through common nearby neurons for interdistance calculation
               for(int j = 0; j < nearbyCounts[currentNeuron]; j++)
               {
                  int interNeuron = nearbyNeurons[currentNeuron][j];

                  // Check if the neuron is also nearby to the candidate
                  for(int k = 0; k < nearbyCounts[candidateNeuron]; k++)
                  {
                     if(interNeuron == nearbyNeurons[candidateNeuron][k]
                        && interNeuron != currentNeuron)
                     {
                        // Calculate the interdistance
                        double interdistance = interDistance(neurons[currentNeuron], neurons[candidateNeuron], neurons[j]);

                        // Adjust the connection probability based on intermediate neurons
                        if(interdistance < radius
                           && isBetween(neurons[currentNeuron], neurons[interNeuron], neurons[candidateNeuron]))
                        {
                           if(interdistance <= (0.5 * radius)
                              && interNeuron != initialNeuron
                              && neurons[interNeuron].synapses < N_synapses - 1)
                           {
                              // The intermediate neuron becomes the new candidate
                              candidateNeuron = interNeuron;
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
               double dist = Distance(neurons[currentNeuron], neurons[candidateNeuron]);
               if(connectionProbability(dist, maxDistance, intDistFactor) > (float)rand() / RAND_MAX)
               {
                  nextNeuron = candidateNeuron;
                  break;
               }
            }
         }


         // If a valid next neuron was found, create a connection
         if(nextNeuron != -1)
         {
            // Calculate the total distance through the connections
            totalDistance += Distance(neurons[currentNeuron], neurons[nextNeuron]);
            
            //Keep track of N_connections
            N_connections += 1;

            // Write the connected neuron
            fprintf(connectionFilePtr, "%lf %lf %lf\n", neurons[nextNeuron].x, neurons[nextNeuron].y, neurons[nextNeuron].radius);

            // Add a synapse if the connection is made
            neurons[nextNeuron].synapses += 1;
            neurons[currentNeuron].synapses += 1;

            // Move to the next neuron
            fprintf(stdout, "Connection made from: %d\n", currentNeuron);
            currentNeuron = nextNeuron;
         }
         else
         {
            // If no valid connection is found, save the number of attempts and try again
            fprintf(stderr, "No valid connection found from neuron %d\n", currentNeuron);
            failedConnectionAttempt += 1;

            // Print space to the file to not make a new connection in the timeStep for the gif animation
            fprintf(connectionFilePtr, "%lf %lf %lf\n", neurons[currentNeuron].x, neurons[currentNeuron].y, neurons[currentNeuron].radius);
         }


         // When reaching the last time step, calculate the end-to-end distance
         if(i + 1 == timeSteps)
         {
            end_to_end_Distance = Distance(neurons[initialNeuron], neurons[currentNeuron]);
         }
      }

      // Print the end-to-end distance and total distance
      double connectedFraction = (double)N_connections / N_neurons;
      fprintf(connectionFilePtr, "\n\nEnd-to-end   Total   connectedFraction\n");
      fprintf(connectionFilePtr, "%lf %lf %lf\n\n", end_to_end_Distance, totalDistance, connectedFraction);
   }

   // Close connection file
   fclose(connectionFilePtr);

   // Free allocated memory
   free(neurons);
   freeNearbyParticles(nearbyNeurons, N_neurons, nearbyCounts);

   printf("Neuron data saved to %s and connection data saved to %s\n", neuronFile, connectionFile);

   return EXIT_SUCCESS;
}


/*
+=========================================================================================================+
||                                               FUNCTIONS                                               ||
+=========================================================================================================+  
*/

/**************************************************************************************************************
* @brief Computes the packing fraction based on the number of neurons, box dimensions and neuron radius.
*
* @param N_neurons The total number of neurons.
* @param L_x       The width of the box.
* @param L_y       The height of the box.
* @param radius    The radius of each neuron.
*
* @return If the packing fraction is too big, the program ends.
***************************************************************************************************************/
void packingFraction(int N_neurons, double L_x, double L_y, double radius)
{
   double packingfraction = (N_neurons * PI * pow(radius, 2)) / (L_x * L_y);
   if(packingfraction >= PACKING_FRACTION)
   {
      fprintf(stderr, "Packing fraction: %.5f is too dense.\n", packingfraction);
      exit(EXIT_FAILURE);
   }
   fprintf(stdout, "Packing fraction: %.5f\n", packingfraction);
}



/**************************************************************************************************************
* @brief Places neurons inside a box in random positions avoiding overlaping.
*
* @param N_neurons The total number of neurons.
* @param L_x       The width of the box.
* @param L_y       The height of the box.
* @param radius    The radius of each neuron.
* @param neurons   Array of structures of neurons.
**************************************************************************************************************/
void placeRandom(int N_neurons, double L_x, double L_y, double radius, Neuron *neurons)
{
   int placedParticles = 0;
      while(placedParticles < N_neurons)
      {
         Neuron newNeuron;
         newNeuron.x = setInside(L_x,radius);
         newNeuron.y = setInside(L_y,radius);
         newNeuron.radius = radius;
         newNeuron.synapses = 0;

         //Check for overlap with existing neurons
         if(!isOverlapping(neurons, placedParticles, newNeuron))
         {
            neurons[placedParticles] = newNeuron;
            placedParticles++;
         }
      }
}



/**************************************************************************************************************
* @brief Places neurons inside a box in random positions avoiding overlaping.
*
* @param N_neurons The total number of neurons.
* @param L_x       The width of the box.
* @param L_y       The height of the box.
* @param radius    The radius of each neuron.
* @param neurons   Array of structures of neurons.
**************************************************************************************************************/
void placeLattice(int N_neurons, double L_x, double L_y, double radius, Neuron *neurons)
{   
   int numRows = (int)sqrt(N_neurons);
   int numCols = (N_neurons + numRows - 1) / numRows; // Handles non-perfect squares

   double dx = (L_x - 2 * radius) / (numCols - 1);    // Spacing in the x direction
   double dy = (L_y - 2 * radius) / (numRows - 1);    // Spacing in the y direction

   double distance = sqrt(dx * dx + dy* dy);          // Handles overlapping
   if(distance < 2 * radius)
   {
      fprintf(stderr, "Error: unavoidable overlapping");
      exit(EXIT_FAILURE);
   }

   int placedParticles = 0;
   for(int i = 0; i < numRows; ++i)
   {
       for(int j = 0; j < numCols; ++j)
       {
         if(placedParticles < N_neurons)
         {
            neurons[placedParticles].x = radius + j * dx;
            neurons[placedParticles].y = radius + i * dy;
            neurons[placedParticles].radius = radius;
            neurons[placedParticles].synapses = 0;
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
double Distance(Neuron p1, Neuron p2)
{
   double distance = sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
   if(distance < 1e-4)
   {
     return 0.0;
   }
   return distance;
}



/**************************************************************************************************************
* @brief Compares the distance between three neurons and checks if there is an intermediate one.
*
* @param p1 The current neuron.
* @param p2 A possible intermediate neuron.
* @param p3 The current candidate neuron.
*
* @returns True (1) if p2 is between p1 and p3, else, False (0).
**************************************************************************************************************/
int isBetween(Neuron p1, Neuron p2, Neuron p3)
{
   double distance = Distance(p1, p3);
   if(distance > Distance(p1, p2) && distance > Distance(p2, p3))
   {
     return TRUE;
   }
   return FALSE;
}



/**************************************************************************************************************
* @brief Checks if a new neuron overlaps with any of the existing neurons.
*
* @param neurons   Array of structures of existing neurons.
* @param N_neurons The number of existing neurons.
* @param newNeuron The new neuron to check for overlap.
*
* @return True (1) if there is no overlap, False (0) if there is overlap.
**************************************************************************************************************/
int isOverlapping(Neuron *neurons, int N_neurons, Neuron newNeuron)
{
   for(int i = 0; i < N_neurons; ++i)
   {
      if(Distance(neurons[i], newNeuron) < 2 * newNeuron.radius) 
      {
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
   return min + scale * (max - 2 * min);     // Scale and shift the value to the desired range
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
   if(dist >= maxDistance)
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
double interDistance(Neuron p1, Neuron p2, Neuron p3)
{
   double lineLength = Distance(p1 , p2);
   if(lineLength < 1e-9)
   {
     return Distance(p1, p3);
   }
   return fabs((p2.y - p1.y) * p3.x - (p2.x - p1.x) * p3.y + p2.x * p1.y - p2.y * p1.x) / lineLength;
}



/**************************************************************************************************************
* @brief Generates a list of nearby neurons for each neuron, based on a cutoff radius.
*
* @param neurons       Array of neurons.
* @param N_neurons     The total number of neurons.
* @param nearbyNeurons A pointer to an array of lists, containing the indices of nearby neurons for each neuron.
* @param nearbyCounts  A pointer to an array storing the number of nearby neurons for each neuron.
**************************************************************************************************************/
void findNearbyParticles(Neuron *neurons, int N_neurons, int ***nearbyNeurons, int **nearbyCounts)
{
   // Allocate memory for nearby neurons list and counts
   *nearbyNeurons = (int **)malloc(N_neurons * sizeof(int *));
   *nearbyCounts = (int *)malloc(N_neurons * sizeof(int));

   for(int i = 0; i < N_neurons; ++i)
   {
      (*nearbyCounts)[i] = 0;
      (*nearbyNeurons)[i] = (int *)malloc(N_neurons * sizeof(int)); // allocate max size initially
   }

   // Populate nearby neurons list for each neuron
   for(int i = 0; i < N_neurons; ++i)
   {
      for(int j = 0; j < N_neurons; ++j)
      {
         if(i != j)
         {
            if(Distance(neurons[i], neurons[j]) <= CUTOFF_RADIUS)
            {
               (*nearbyNeurons)[i][(*nearbyCounts)[i]] = j; // Store index of nearby neuron
               (*nearbyCounts)[i]++;
            }
         }
      }
   }
}



/**************************************************************************************************************
* @brief Frees the memory allocated for the nearby neurons lists.
*
* @param nearbyNeurons The list of nearby neurons for each neuron.
* @param N_neurons     The total number of neurons.
* @param nearbyCounts  The array storing the count of nearby neurons for each neuron.
**************************************************************************************************************/
void freeNearbyParticles(int **nearbyNeurons, int N_neurons, int *nearbyCounts)
{
   for (int i = 0; i < N_neurons; ++i)
   {
      free(nearbyNeurons[i]);
   }
   free(nearbyNeurons);
   free(nearbyCounts);
}
