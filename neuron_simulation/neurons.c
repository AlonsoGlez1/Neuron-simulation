#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>


#define PI 3.141592654            // Constant for the value of pi
#define DECAY_FACTOR 50.0         // Factor for the distance decay in connection probability
#define INT_DIST_FACTOR 2.0       // Factor for the interdistance decay in connection probability
#define SCALE_FACTOR 1.0          // Global factor to scale the probability distribution
#define CUTOFF_RADIUS 20.0        // Maximum radius within which to consider nearby neurons
#define PACKING_FRACTION 0.5      // Maximum packing fraction of the system (area occupied / total area)
#define BRANCHES 2                // Number of branches per initial neuron
#define IMPORT_NEURONS_LIST false // Select false to generate neurons, select true to read data from neurons_dat
#define ONLY_GENERATE_LIST false  // Select false to run the whole program, select true to only generate neurons.


// Define a structure to represent a neuron
typedef struct
{
   double x;       // x-coordinate of the neuron
   double y;       // y-coordinate of the neuron
   double radius;  // Radius of the neuron
   int synapses;   // Number of synapses (connections)
} Neuron;


// Function prototypes
void packingFraction(int N_neurons, double L_x, double L_y, double radius);
void importNeurons(char* filename, Neuron** neurons, int* N_neurons);
int countLinesInFile(FILE *file);
void readNeuronData(FILE *file, Neuron* neurons, int N_neurons);
void placeRandom(int N_neurons, double L_x, double L_y, double radius, Neuron *neurons);
void placeLattice(int N_neurons, double L_x, double L_y, double radius, Neuron *neurons);
void initializeNeurons(Neuron *neurons, int N_neurons, double L_x, double L_y, double radius, char mode);
void writeNeuronData(FILE *file, Neuron *neurons, int N_neurons);
double Distance(Neuron p1, Neuron p2);
bool isOverlapping(Neuron *neurons, int N_neurons, Neuron newNeuron);
float setInside(float max, float min);
bool isBetween(int *currentNeuron, int *interNeuron, int *candidateNeuron, double **distanceMatrix);
double connectionProbability(double *dist, double *intDistFactor, double *maxDistance);
double interDistance(int *currentNeuron, int *candidateNeuron, int *candidateNeuronIndex, int *interNeuron, int *interNeuronIndex, double **distanceMatrix, double ***interDistanceCache);
void findNearbyNeurons(int N_neurons, int ***nearbyNeurons, int **nearbyCounts, double **distanceMatrix);
void freeNearbyNeurons(int **nearbyNeurons, int N_neurons);
double*** allocate3DArray(int x, int y, int z);
void free3DArray(double ***array, int x, int y);
double** allocate2DArray(int rows, int cols);
void free2DArray(double** array, int rows);
int** allocate2DIntArray(int rows, int cols);
void free2DIntArray(int** array, int rows);
double*** initializeInterdistanceCache(int N_neurons, int *nearbyCounts);
void freeInterdistanceCache(double ***array, int N_neurons, int *nearbyCounts);
double** initializeDistanceMatrix(Neuron* neurons, int N_neurons);
double radiusOfGyration(Neuron *neurons, int  *N_connections, int **neuronsList, int currentBranch, int currentTime, double *maxExtension);


int main()
{
   int N_neurons, timeSteps, N_synapses; // Number of neurons, timeSteps for connections and synapses
   double L_x, L_y, radius;              // Physical dimensions
   char mode;                            // Mode (random|lattice)
   Neuron *neurons = NULL;               // List of structures of neurons

   // Seed the random number generator
   srand(time(NULL));

/*
+=========================================================================================================+
||                                     INITIALIZE NEURON POSITIONS                                       ||
+=========================================================================================================+
*/
   if (IMPORT_NEURONS_LIST) 
   {
      importNeurons("neurons_dat", &neurons, &N_neurons);

      printf("Enter the width (L_x) of the box: ");            scanf("%lf", &L_x);
      L_y = L_x;
      radius = neurons[0].radius;
      packingFraction(N_neurons, L_x, L_y, radius);

      printf("Enter the time steps: ");                        scanf("%d", &timeSteps);
      printf("Enter the max number of synapses per neuron: "); scanf("%d", &N_synapses);
   }
   else
   {
      // Read the number of neurons, the dimensions of the box and neuron radii
      printf("Random or Lattice mode (R/L): ");     scanf(" %c", &mode);
      // Check the packing fraction to ensure it isn't too dense
      printf("Enter the number of neurons: ");      scanf("%d", &N_neurons);
      printf("Enter the width (L_x) of the box: "); scanf("%lf", &L_x);
      printf("Enter the radius of the neurons: ");  scanf("%lf", &radius);
      L_y = L_x;
      packingFraction(N_neurons, L_x, L_y, radius);

      // Allocate memory for the neurons
      neurons = (Neuron *)calloc(N_neurons, sizeof(Neuron));
      if(!neurons)
      {
         fprintf(stderr, "Memory allocation failed\n");
         return EXIT_FAILURE;
      }
      
      //Generate and store the positions list
      initializeNeurons(neurons, N_neurons, L_x, L_y, radius, mode);
      
      // Write the neuron positions and radii to a file
      FILE *neuronFilePtr = fopen("neurons_dat", "w");
      writeNeuronData(neuronFilePtr, neurons, N_neurons);
      fclose(neuronFilePtr);

      if(ONLY_GENERATE_LIST)
      {
         printf("Neuron data saved to neurons_dat\n");
         return EXIT_SUCCESS;
      }

      printf("Enter the time steps: ");                        scanf("%d", &timeSteps);
      printf("Enter the max number of synapses per neuron: "); scanf("%d", &N_synapses);
   }  

/*
+=========================================================================================================+
||                                     INITIALIZE NEURON CONNECTIONS                                     ||
+=========================================================================================================+
*/
   double maxDistance = sqrt(2.0) * (L_x - 2*radius); // Diagonal of the box, maximum distance possible

   // Calculate distances between all neurons and store them
   double **distanceMatrix = initializeDistanceMatrix(neurons, N_neurons);

   // Create an array of arrays (list of nearby neurons)
   int **nearbyNeurons = NULL;
   int *nearbyCounts = NULL;
   findNearbyNeurons(N_neurons, &nearbyNeurons, &nearbyCounts, distanceMatrix);

   // Save intedistance calculations for every neuron to avoid redundant calculations
   double ***interDistanceCache = initializeInterdistanceCache(N_neurons, nearbyCounts);

   /* Create and allocate memory for a 3D array to store data for each branch
      Data: 0: N_Connections. 1: Total Distance. 2: Connected Fraction. 4: End-To-EndDistance. 5: Gyration Radius. 6: Maximum Extension.
      Array gets stored as systemData[Data][TimeSteps][Branch] */
   double ***eachBranchData = allocate3DArray(6, timeSteps, BRANCHES);
   int **neuronsList = allocate2DIntArray(BRANCHES, timeSteps + 1);

   // Allocate arrays to keep track of each branch's progress
   int *currentNeurons = (int *)calloc(BRANCHES, sizeof(int));
   int *N_connections = (int *)calloc(BRANCHES, sizeof(int));
   double *totalDistances = (double *)calloc(BRANCHES, sizeof(double));
   double *endToEndDistances = (double *)calloc(BRANCHES, sizeof(double));
   double *connectedFractions = (double *)calloc(BRANCHES, sizeof(double));

   if(!currentNeurons || !N_connections || !totalDistances || !endToEndDistances || !connectedFractions)
   {
      fprintf(stderr, "Memory allocation failed for branch data\n");
      free(neurons);
      freeNearbyNeurons(nearbyNeurons, N_neurons);
      return EXIT_FAILURE;
   }


   // Initialize all branches with the same starting random neuron
   int initialNeuron = rand() % (N_neurons + 1);
   for(int branch = 0; branch < BRANCHES; ++branch)
   {
      currentNeurons[branch] = initialNeuron;
      neuronsList[branch][0] = initialNeuron;
   }
   neurons[initialNeuron].synapses = N_synapses - BRANCHES;

   // Start looking for candidates. Every branch in the same time step.
   for(int time = 0; time < timeSteps; ++time)
   {
      for(int branch = 0; branch < BRANCHES; branch++)
      {
         int currentNeuron = currentNeurons[branch];
         int nextNeuron = -1;

         // Attempt to find the next valid neuron from nearby neurons list
         for(int attempt = 0; attempt < nearbyCounts[currentNeuron]; ++attempt)
         {  
            // Select a candidate with its local index in the nearbyNeurons list
            int candidateNeuronNearbyIndex = rand() % (nearbyCounts[currentNeuron] + 1);
            // Select a candidate with its global index in the neurons list of structures
            int candidateNeuron = nearbyNeurons[currentNeuron][candidateNeuronNearbyIndex];

            // Ensure the candidate is valid for connection
            if(candidateNeuron != initialNeuron
               && neurons[candidateNeuron].synapses < N_synapses - 1
               && neurons[currentNeuron].synapses < N_synapses)
            {
               double intDistFactor = 1.0; // Reduction factor by intermediate neurons

               // Loop through nearby neurons for interdistance calculation
               for(int j = 0; j < nearbyCounts[currentNeuron]; j++)
               {
                  // Select a possible interNeuron with its local index in the nearbyNeurons list
                  int interNeuronNearbyIndex = j;
                  // Select a possible interNeuron with its global index in the neurons list of structures
                  int interNeuron = nearbyNeurons[currentNeuron][j];
                  if(interNeuron == initialNeuron)
                  {
                     continue; // Skip current loop
                  }

                  // Check if the interNeuron is in the space between the current and the candidate neurons
                  if(isBetween(&currentNeuron, &interNeuron, &candidateNeuron, distanceMatrix))
                  {
                     // Calculate the interdistance
                     double interdistance = interDistance(&currentNeuron, &candidateNeuron, &candidateNeuronNearbyIndex, &interNeuron, &interNeuronNearbyIndex, distanceMatrix, interDistanceCache);

                     // Adjust the connection probability based on intermediate neurons
                     if(interdistance < radius)
                     {
                        if(interdistance <= (0.5 * radius) && neurons[interNeuron].synapses < N_synapses - 1)
                        {
                           // The intermediate neuron becomes the new candidate
                           candidateNeuron = interNeuron;
                           intDistFactor = 1.0;
                           j = 0;
                        }
                        else
                        {
                           intDistFactor += interdistance / radius;
                        }
                     }
                  }
               }

               // Probability check for connection based on distance
               double dist = distanceMatrix[currentNeuron][candidateNeuron];
               if(connectionProbability(&dist, &intDistFactor, &maxDistance) > (float)rand() / RAND_MAX)
               {
                  nextNeuron = candidateNeuron;
                  break;
               }
            }
         }

         // If a connection was made, update the branch data
         if(nextNeuron != -1)
         {
            N_connections[branch]++;
            totalDistances[branch] += distanceMatrix[currentNeuron][nextNeuron];
            connectedFractions[branch] = (double)N_connections[branch] / N_neurons;
            neurons[nextNeuron].synapses += 1;
            neurons[currentNeuron].synapses += 1;

            // Update the current neuron for the next iteration
            fprintf(stdout, "Connection made from: %d\n", currentNeuron);
            currentNeurons[branch] = nextNeuron;
            neuronsList[branch][time + 1] = nextNeuron;
         }
         else
         {
            neuronsList[branch][time + 1] = currentNeuron;
            fprintf(stderr, "No valid connection found from neuron %d\n", currentNeuron);
         }
            // Store data in multiple branches on each time step
            eachBranchData[0][time][branch] = N_connections[branch];
            eachBranchData[1][time][branch] = totalDistances[branch];
            eachBranchData[2][time][branch] = connectedFractions[branch];
            eachBranchData[3][time][branch] = distanceMatrix[initialNeuron][currentNeurons[branch]];

            double maximumExtension = 0;
            eachBranchData[4][time][branch] = radiusOfGyration(neurons, N_connections, neuronsList, branch, time + 1, &maximumExtension);
            eachBranchData[5][time][branch] = 2.0 * maximumExtension / (sqrt(2.0) * maxDistance);
      }
   }


   // Open the file for writing neuron connections
   FILE *connectionFilePtr = fopen("connections_dat", "w");
   if (!connectionFilePtr)
   {
      fprintf(stderr, "\nFailed to open connection data file for writing\n");

      // Free allocated memory
      free(neurons);
      freeNearbyNeurons(nearbyNeurons, N_neurons);
      free3DArray(eachBranchData, 6, timeSteps);
      freeInterdistanceCache(interDistanceCache, N_neurons, nearbyCounts);
      free2DIntArray(neuronsList, BRANCHES);
      free2DArray(distanceMatrix, N_neurons);
      free(currentNeurons);
      free(N_connections);
      free(totalDistances);
      free(endToEndDistances);
      free(connectedFractions);

      return EXIT_FAILURE;
   }

   // Connection format string for file output
   const char *formatConnections = "%8.4lf %8.4lf %6.2lf %10.4lf %10.0lf %15.4lf %15.4lf %15.4lf %15.4lf %12d\n";
   for(int branch = 0; branch < BRANCHES; branch++)
   {
      fprintf(connectionFilePtr, "Xlabel   Ylabel   Radius   EndToEnd   N_connections   TotalDist   ConnectedFraction   GyrationRadius   NormalizedMaxExt   TimeStep\n");
      fprintf(connectionFilePtr, formatConnections,
              neurons[initialNeuron].x, neurons[initialNeuron].y, neurons[initialNeuron].radius,
              0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0);

      for(int time = 0; time < timeSteps; time++)
      {
         // Write connected neuron to the file
         fprintf(connectionFilePtr, formatConnections,
                 neurons[neuronsList[branch][time + 1]].x, neurons[neuronsList[branch][time + 1]].y, neurons[neuronsList[branch][time + 1]].radius,
                 eachBranchData[3][time][branch], eachBranchData[0][time][branch], eachBranchData[1][time][branch],
                 eachBranchData[2][time][branch], eachBranchData[4][time][branch], eachBranchData[5][time][branch],  time + 1);
      }
      //Line break to separate data
      fprintf(connectionFilePtr, "\n\n\n");
   }
   // Close connection file
   fclose(connectionFilePtr);


   // Write in a file the results for multiple branches in a single time step
   FILE *resultsFilePtr = fopen("results_dat", "w");
   if(!resultsFilePtr)
   {
      fprintf(stderr, "\nFailed to open results data file for writing\n");

      // Free allocated memory
      free(neurons);
      freeNearbyNeurons(nearbyNeurons, N_neurons);
      free3DArray(eachBranchData, 6, timeSteps);
      freeInterdistanceCache(interDistanceCache, N_neurons, nearbyCounts);
      free2DArray(distanceMatrix, N_neurons);
      free2DIntArray(neuronsList, BRANCHES);
      free(currentNeurons);
      free(N_connections);
      free(totalDistances);
      free(endToEndDistances);
      free(connectedFractions);

      return EXIT_FAILURE;
   }

   const char *formatResults = "%5.0lf %18.4lf %15.4lf %15d\n";
   fprintf(resultsFilePtr, "N_Connections   TotalDist   ConnectedFraction   TimeStep\n");
   fprintf(resultsFilePtr, formatResults, 0, 0.000000, 0.000000, 0);

   // Create and allocate memory for a 2D array that merges the multiple branches data into the same time step
   double** mergedBranchesData = allocate2DArray(3, timeSteps);

   // Merge the data of the 3D array into the same time step
   for(int time = 0; time < timeSteps; time++)
   {
      for(int branch = 0; branch < BRANCHES; branch++)
      {
         mergedBranchesData[0][time] += eachBranchData[0][time][branch];
         mergedBranchesData[1][time] += eachBranchData[1][time][branch];
         mergedBranchesData[2][time] += eachBranchData[2][time][branch];
      }
   }

   // Print data per time step
   for(int time = 0; time < timeSteps; time++)
   {
      fprintf(resultsFilePtr, formatResults, mergedBranchesData[0][time], mergedBranchesData[1][time], mergedBranchesData[2][time], time + 1);
   }

   fclose(resultsFilePtr);

   // Free allocated memory
   free(neurons);
   freeNearbyNeurons(nearbyNeurons, N_neurons);
   freeInterdistanceCache(interDistanceCache, N_neurons, nearbyCounts);
   free3DArray(eachBranchData, 6, timeSteps);
   free2DArray(mergedBranchesData, 3);
   free2DArray(distanceMatrix, N_neurons);
   free2DIntArray(neuronsList, BRANCHES);
   free(currentNeurons);
   free(N_connections);
   free(totalDistances);
   free(endToEndDistances);
   free(connectedFractions);

   printf("Neuron data saved to neurons_dat, connection data saved to connections_dat and results_dat\n");
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

   // Check if the packing fraction exceeds the predefined limit
   if(packingfraction >= PACKING_FRACTION)
   {
      fprintf(stderr, "Packing fraction: %.5f is too dense.\n", packingfraction);
      exit(EXIT_FAILURE);
   }
   fprintf(stdout, "Packing fraction: %.5f\n", packingfraction);
}


/**************************************************************************************************************
* @brief Function to import the file data positions into an array of structures.
*
* @param filename  The filename where the data is stored.
* @param neurons   Array of structures of neurons.
* @param N_neurons The total number of neurons.
***************************************************************************************************************/
void importNeurons(char* filename, Neuron** neurons, int* N_neurons)
{
   FILE *file = fopen(filename, "r");
   if(!file)
   {
      perror("Error opening file for reading.\n");
      exit(EXIT_FAILURE);
   }

   char header[100];
   if(!fgets(header, sizeof(header), file))
   {
      perror("Error reading header");
      fclose(file);
      exit(EXIT_FAILURE);
   }

   *N_neurons = countLinesInFile(file);
   *neurons = (Neuron *)calloc(*N_neurons, sizeof(Neuron));
   if(!*neurons)
   {
      fprintf(stderr, "Memory allocation failed\n");
      exit(EXIT_FAILURE);   
   }  
   readNeuronData(file, *neurons, *N_neurons);

   fclose(file);
}


/**************************************************************************************************************
* @brief Function to count the number of lines in the file, ignoring the headers.
*
* @param file The file name where the data is stored.
*
* @returns The number of data lines in the file (number of neurons).
**************************************************************************************************************/
int countLinesInFile(FILE *file)
{
   int lines = 0;
   char buffer[100];
   while (fgets(buffer, sizeof(buffer), file))
   {
     lines++;
   }
   rewind(file); // Reset file pointer to the beginning
   
   return lines;
}


/**************************************************************************************************************
* @brief Function that reads the data in the file and stores it in an array of structures.
*
* @param file      The file name where the data is stored.
* @param neurons   An array of structures of neurons.
* @param N_neurons Total number of neurons in the file.
**************************************************************************************************************/
void readNeuronData(FILE *file, Neuron* neurons, int N_neurons)
{
   char header[100];
   fgets(header, sizeof(header), file); // Skip header
   for (int i = 0; i < N_neurons; i++)
   {
     if(fscanf(file, "%lf %lf %lf", &neurons[i].x, &neurons[i].y, &neurons[i].radius) != 3)
      {
         fprintf(stderr, "Error reading line %d\n", i + 1);
         free(neurons);
         fclose(file);
         
         exit(EXIT_FAILURE);
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
void placeRandom(int N_neurons, double L_x, double L_y, double radius, Neuron *neurons)
{
   int placedNeurons = 0;
      while(placedNeurons < N_neurons)
      {
         Neuron newNeuron;
         newNeuron.x = setInside(L_x,radius);
         newNeuron.y = setInside(L_y,radius);
         newNeuron.radius = radius;
         newNeuron.synapses = 0;

         //Check for overlap with existing neurons
         if(!isOverlapping(neurons, placedNeurons, newNeuron))
         {
            neurons[placedNeurons] = newNeuron;
            placedNeurons++;
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

   int placedNeurons = 0;
   for(int i = 0; i < numRows; ++i)
   {
       for(int j = 0; j < numCols; ++j)
       {
         if(placedNeurons < N_neurons)
         {
            neurons[placedNeurons].x = radius + j * dx;
            neurons[placedNeurons].y = radius + i * dy;
            neurons[placedNeurons].radius = radius;
            neurons[placedNeurons].synapses = 0;
            placedNeurons++;
         }
       }
   }
}


/**************************************************************************************************************
* @brief Function that prints into a file the positions of generated neurons.
* @param file      The file name.
* @param neurons   Array of structures of neurons.
* @param N_neurons Total amount of neurons.
**************************************************************************************************************/
void writeNeuronData(FILE *file, Neuron *neurons, int N_neurons)
{
   if (!file)
   {
      fprintf(stderr, "Error opening file for writing.\n");
      exit(EXIT_FAILURE);
   }

   // Write header
   fprintf(file, "Xlabel   Ylabel   Radius\n");

   // Write neuron data
   for(int i = 0; i < N_neurons; ++i)
   {
      fprintf(file, "%8.4lf %8.4lf %6.2lf\n", neurons[i].x, neurons[i].y, neurons[i].radius);
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
   double dx = p1.x - p2.x;
   double dy = p1.y - p2.y;

   return sqrt(dx * dx + dy * dy);
}


/**************************************************************************************************************
* @brief Compares the distance between three neurons and checks if there is an intermediate one.
*
* @param currentNeuron The current neuron.
* @param interNeuron A possible intermediate neuron.
* @param candidateNeuron The current candidate neuron.
* @param distanceMatrix  A 2D matrix that stores the distance between every pair of neurons.
*
* @returns True (1) if p2 is between p1 and p3, else, False (0).
**************************************************************************************************************/
bool isBetween(int *currentNeuron, int *interNeuron, int *candidateNeuron, double **distanceMatrix)
{
   if(distanceMatrix[*currentNeuron][*candidateNeuron] > distanceMatrix[*currentNeuron][*interNeuron] 
      && distanceMatrix[*currentNeuron][*candidateNeuron] > distanceMatrix[*interNeuron][*candidateNeuron])
   {
     return true;
   }
   return false;
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
bool isOverlapping(Neuron *neurons, int N_neurons, Neuron newNeuron)
{
   for(int i = 0; i < N_neurons; ++i)
   {
      if(Distance(neurons[i], newNeuron) < 2 * newNeuron.radius)
      {
         return true;
      }
   }
   return false;
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
* @param intDistFactor A factor based on intermediate distances between neurons.
* @param maxDistance   The maximum possible distance within the box.
*
* @return A value between 0 and 1, where larger distances reduces the connection probability.
**************************************************************************************************************/
double connectionProbability(double *dist, double *intDistFactor, double *maxDistance)
{
   return SCALE_FACTOR * (1 - (*dist) / (*maxDistance)) * exp(- DECAY_FACTOR * ((*dist) / (*maxDistance)) - INT_DIST_FACTOR * (*intDistFactor));
}


/**************************************************************************************************************
* @brief Computes the perpendicular distance between a neuron and a pair of neurons.
*
* @param currentNeuron        The first neuron, it connects to the candidate one
* @param candidateNeuron      The second neuron, it connects to the first currentNeuron.
* @param candidateNeuronIndex The index to the second neuron, local to the nearbyNeurons list.
* @param interNeuron          The third neuron, it could be between the first and the second one.
* @param interNeuronIndex     The index to the third neuron, local to the nearbyNeurons list.
* @param distanceMatrix       A 2D matrix that stores the distance between every pair of neurons.
* @param interDistanceCache   A 3D matrix that stores the interdistance between neurons.
*
* @return The calculated interdistance between the three neurons.
**************************************************************************************************************/
double interDistance(int *currentNeuron, int *candidateNeuron, int *candidateNeuronIndex, int *interNeuron, int *interNeuronIndex, double **distanceMatrix, double ***interDistanceCache)
{
   // Check if the interdistance is already calculated
   if(interDistanceCache[*currentNeuron][*candidateNeuronIndex][*interNeuronIndex] > 0)
   {
      return interDistanceCache[*currentNeuron][*candidateNeuronIndex][*interNeuronIndex];
   }

   // Calculate the interdistance
   double d12 = distanceMatrix[*currentNeuron][*candidateNeuron];
   double d13 = distanceMatrix[*currentNeuron][*interNeuron];
   double d23 = distanceMatrix[*candidateNeuron][*interNeuron];

   double s = (d12 + d13 + d23) / 2.0;
   double area = sqrt(s * (s - d12) * (s - d13) * (s - d23));
   double interdistance = (2.0 * area) / d12;

   // Store the result in the cache
   interDistanceCache[*currentNeuron][*candidateNeuronIndex][*interNeuronIndex] = interdistance;

   return interdistance;
}


/**************************************************************************************************************
* @brief Generates a list of nearby neurons for each neuron, based on a cutoff radius.
*
* @param N_neurons      The total number of neurons.
* @param nearbyNeurons  A pointer to an array of lists, containing the indices of nearby neurons for each neuron.
* @param nearbyCounts   A pointer to an array storing the number of nearby neurons for each neuron.
* @param distanceMatrix A 2D matrix that stores the distance between every pair of neurons.
**************************************************************************************************************/
void findNearbyNeurons(int N_neurons, int ***nearbyNeurons, int **nearbyCounts, double **distanceMatrix)
{
   // Allocate memory for nearby neurons list and counts
   *nearbyNeurons = (int **)calloc(N_neurons, sizeof(int *));
   *nearbyCounts = (int *)calloc(N_neurons, sizeof(int));

   if (!*nearbyNeurons || !*nearbyCounts)
   {
      fprintf(stderr, "Failed to allocate memory for nearby neurons or counts\n");
      exit(EXIT_FAILURE);
   }

   // First pass: Count the number of nearby neurons for each neuron
   for(int i = 0; i < N_neurons; ++i)
   {
      (*nearbyCounts)[i] = 0;  // Initialize the count for each neuron
      for(int j = 0; j < N_neurons; ++j)
      {
         if(i != j && distanceMatrix[i][j] <= CUTOFF_RADIUS)
         {
            (*nearbyCounts)[i]++;  // Increment the nearby count for neuron i
         }
      }
   }

   // Allocate memory for the actual number of neighbors
   for(int i = 0; i < N_neurons; ++i)
   {
      (*nearbyNeurons)[i] = (int *)calloc((*nearbyCounts)[i] + 1, sizeof(int));  // Allocate exact memory size

      if(!(*nearbyNeurons)[i])
      {
         fprintf(stderr, "Failed to allocate memory for nearby neurons list\n");
         exit(EXIT_FAILURE);
      }
   }

   // Second pass: Populate the list of nearby neurons
   for(int i = 0; i < N_neurons; ++i)
   {
      int count = 0;  // Index for storing nearby neurons
      for(int j = 0; j < N_neurons; ++j)
      {
         if(i != j && distanceMatrix[i][j] <= CUTOFF_RADIUS)
         {
            (*nearbyNeurons)[i][count] = j; // Store index of nearby neuron
            count++;
         }
      }
   }
}


/**************************************************************************************************************
* @brief Frees the memory allocated for the nearby neurons lists.
*
* @param nearbyNeurons The list of nearby neurons for each neuron.
* @param N_neurons     The total number of neurons.
**************************************************************************************************************/
void freeNearbyNeurons(int **nearbyNeurons, int N_neurons)
{
   for(int i = 0; i < N_neurons; ++i)
   {
      free(nearbyNeurons[i]);
   }
   free(nearbyNeurons);
}


/**************************************************************************************************************
* @brief Allocates memory for a dynamic 2D array with integer elements and initialize them to 0.
*
* @param rows Rows of the array.
* @param cols Columns of the array.
*
* @return The array with memory allocated and all elements initialized to 0.
**************************************************************************************************************/
int** allocate2DIntArray(int rows, int cols)
{
   // Validate parameters
   if(rows <= 0 || cols <= 0)
   {
      fprintf(stderr, "Invalid dimensions: rows and cols must be positive values\n");
      exit(EXIT_FAILURE);
   }

   // Allocate memory for the array of pointers
   int** array = (int **)calloc(rows, sizeof(int*));
   if(!array)
   {
      fprintf(stderr, "Failed to allocate memory for the rows\n");
      exit(EXIT_FAILURE);
   }

   // Allocate memory for each row and initialize to 0
   for(int i = 0; i < rows; i++)
   {
      array[i] = (int *)calloc(cols, sizeof(int)); // Allocate and initialize to 0
      if(!array[i])
      {
         fprintf(stderr, "Failed to allocate memory for row %d\n", i);

         // Free previously allocated memory before exiting
         for(int j = 0; j < i; j++)
         {
            free(array[j]);
         }
         free(array);

         exit(EXIT_FAILURE);
      }
   }

   return array; // Return the allocated 2D array
}


/**************************************************************************************************************
* @brief Frees allocated memory for a dynamic 2D array with integer elements.
*
* @param array The array to be freed (Freedom!).
* @param cols  Columns of the array.
**************************************************************************************************************/
void free2DIntArray(int** array, int rows)
{
   if(!array)
   {
      return; // If the array is NULL, do nothing
   }

   for(int i = 0; i < rows; i++)
   {
      free(array[i]); // Free each row
   }
   free(array); // Free the array of pointers
}


/**************************************************************************************************************
* @brief Allocates memory for a dynamic 2D array with float elements and initialize them to 0.
*
* @param rows Rows of the array.
* @param cols Columns of the array.
*
* @return The array with memory allocated and all elements initialized to 0.
**************************************************************************************************************/
double** allocate2DArray(int rows, int cols)
{
   // Allocate memory for the array of pointers
   double** array = (double **)calloc(rows, sizeof(double*));
   if(!array)
   {
      fprintf(stderr, "Failed to allocate memory for the rows\n");
      exit(EXIT_FAILURE);
   }

   // Allocate memory for each row and initialize to 0
   for(int i = 0; i < rows; i++)
   {
      array[i] = (double *)calloc(cols, sizeof(double)); // Allocate and initialize to 0
      if (!array[i])
      {
         fprintf(stderr, "Failed to allocate memory for row %d\n", i);

         // Free previously allocated memory before exiting
         for(int j = 0; j < i; j++)
         {
            free(array[j]);
         }
         free(array);
         exit(EXIT_FAILURE);
      }
   }

   return array; // Return the allocated 2D array
}


/**************************************************************************************************************
* @brief Frees allocated memory for a dynamic 2D array with float elements.
*
* @param array The array to be freed (Freedom!).
* @param cols  Columns of the array.
**************************************************************************************************************/
void free2DArray(double** array, int rows)
{
   if(!array)
   {
      return; // If the array is NULL, do nothing
   }

   for(int i = 0; i < rows; i++)
   {
      free(array[i]); // Free each row
   }
   free(array); // Free the array of pointers
}


/**************************************************************************************************************
* @brief Allocates memory for a dynamic 3D array and initialize all elements to 0.
*
* @param x First dimension.
* @param y Second direction.
* @param z Third direction.
*
* @return The array with memory allocated and all elements initialized to 0.
**************************************************************************************************************/
double*** allocate3DArray(int x, int y, int z)
{
   // Allocate memory for the array of pointers to 2D arrays
   double ***array = calloc(x, sizeof(double **));
   if (!array)
   {
     fprintf(stderr, "Failed to allocate memory for the first dimension\n");
     exit(EXIT_FAILURE);
   }

   // Allocate memory for each 2D array and initialize to 0
   for(int i = 0; i < x; i++)
   {
      array[i] = calloc(y, sizeof(double *));
      if(!array[i])
      {
         fprintf(stderr, "Failed to allocate memory for the second dimension\n");
         // Free the already allocated memory
         for(int k = 0; k < i; k++)
         {
            free(array[k]);
         }
         free(array);
         exit(EXIT_FAILURE);
      }

      // Allocate memory for each row of the 2D arrays and initialize to 0
      for(int j = 0; j < y; j++)
      {
         array[i][j] = calloc(z, sizeof(double));
         if (!array[i][j])
         {
            fprintf(stderr, "Failed to allocate memory for the third dimension\n");
            // Free already allocated memory
            for(int k = 0; k <= i; k++)
            {
                for(int l = 0; l < y; l++)
                {
                    free(array[k][l]);
                }
                free(array[k]);
            }
            free(array);
            exit(EXIT_FAILURE);
         }
      }
   }

   return array; // Return the allocated 3D array
}


/**************************************************************************************************************
* @brief Frees allocated memory for a dynamic 3D array.
*
* @param array The array to be freed (Freedom!).
* @param x     First dimention.
* @param y     Second dimention.
**************************************************************************************************************/
void free3DArray(double ***array, int x, int y)
{
   for(int i = 0; i < x; i++)
   {
      for(int j = 0; j < y; j++)
      {
         free(array[i][j]); // Free each row of the 2D arrays
      }
      free(array[i]); // Free each 2D array
   }
   free(array); // Free the array of pointers to 2D arrays
}


/**************************************************************************************************************
* @brief Function that allocates a 3D array of interdistances to avoid redundant calculations.
*
* @param N_neurons    The total number of neurons.
* @param nearbyCounts A pointer to an array storing the number of nearby neurons for each neuron.
* 
* @return The allocated 3D array.
**************************************************************************************************************/
double ***initializeInterdistanceCache(int N_neurons, int *nearbyCounts)
{
   // Allocate a 3D array to store interdistances
   double ***interdistanceCache = (double ***)calloc(N_neurons, sizeof(double **));
   if(!interdistanceCache)
   {
      fprintf(stderr, "Memory allocation failed for interdistance cache.\n");
      exit(EXIT_FAILURE);
   }

   for (int i = 0; i < N_neurons; i++)
   {
      interdistanceCache[i] = (double **)calloc(nearbyCounts[i] + 1, sizeof(double *));
      if(!interdistanceCache[i])
      {
         fprintf(stderr, "Memory allocation failed for interdistance cache row %d.\n", i);
         // Free previously allocated memory
         for(int k = 0; k < i; k++)
         {
            free(interdistanceCache[k]);
         }
         free(interdistanceCache);
         exit(EXIT_FAILURE);
      }
      for(int j = 0; j < nearbyCounts[i] + 1; j++)
      {
         interdistanceCache[i][j] = (double *)calloc(nearbyCounts[i] + 1, sizeof(double));
         if(!interdistanceCache[i][j])
         {
            fprintf(stderr, "Memory allocation failed for interdistance cache column.\n");
            // Free previously allocated memory in this row
            for(int k = 0; k < j; k++)
            {
               free(interdistanceCache[i][k]);
            }
            free(interdistanceCache[i]);
            // Free all previously allocated rows
            for(int k = 0; k < i; k++)
            {
               for(int l = 0; l < nearbyCounts[i] + 1; l++)
               {
                  free(interdistanceCache[k][l]);
               }
               free(interdistanceCache[k]);
            }
            free(interdistanceCache);
            exit(EXIT_FAILURE);
         }
      }
   }
   
   return interdistanceCache;
}


/**************************************************************************************************************
* @brief Function that frees the 3D array allocated in the initializeInterdistanceCache function.
*
* @param interdistanceCache The 3D array that stores the calculated interdistances.
* @param N_neurons          The total number of neurons.
* @param nearbyCounts       A pointer to an array storing the number of nearby neurons for each neuron.
**************************************************************************************************************/
void freeInterdistanceCache(double ***interdistanceCache, int N_neurons, int *nearbyCounts)
{
   // Free the 3D array memory
   for(int i = 0; i < N_neurons; i++)
   {
      // Free each pointer in the row
      for(int j = 0; j < nearbyCounts[i] + 1; j++)
      {
         free(interdistanceCache[i][j]); // Free the double pointer
      }
      // Free the row itself
      free(interdistanceCache[i]);
   }
   free(interdistanceCache);

   free(nearbyCounts);
}


/**************************************************************************************************************
* @brief Initialices neuron positions in random or lattice distribution
*
* @param neurons   Array of structures of neurons.
* @param N_neurons Total amount of neurons.
* @param L_x       Width of the box.
* @param L_y       Height of the box.
* @param radius    Radius of the neurons.
* @param mode      Random (R/r) or Lattice (L/l).
**************************************************************************************************************/
void initializeNeurons(Neuron *neurons, int N_neurons, double L_x, double L_y, double radius, char mode) 
{
   if(mode == 'R' || mode == 'r')
   {
      placeRandom(N_neurons, L_x, L_y, radius, neurons);
   }
   else if (mode == 'L' || mode == 'l')
   {
      placeLattice(N_neurons, L_x, L_y, radius, neurons);
   }
   else
   {
      fprintf(stderr, "Invalid mode\n");
      exit(EXIT_FAILURE);
   }
}


/**************************************************************************************************************
* @brief Function that initializes and calculates the distances matrix between every pair or neurons.
*
* @param neurons   Array of structures of neurons.
* @param N_neurons The total number of neurons.
**************************************************************************************************************/
double** initializeDistanceMatrix(Neuron* neurons, int N_neurons)
{
   double** distanceMatrix = allocate2DArray(N_neurons, N_neurons);

   for (int i = 0; i < N_neurons; i++)
   {
      for (int j = i + 1; j < N_neurons; j++)
      {
         double distance = Distance(neurons[i], neurons[j]);
         if(distance < 1e-5)
         {
            distance = 0.0;
         }

         distanceMatrix[i][j] = distance;
         distanceMatrix[j][i] = distance;  // Symmetric entry
      }
      distanceMatrix[i][i] = 0.0; // Distance from a neuron to itself is 0
   }
   return distanceMatrix;
}


/**************************************************************************************************************
* @brief Function that calculates the mean square radius of gyration and the normalized maximum extension.
*
* @param neurons       Array of structures of neurons.
* @param N_neurons     The total number of neurons.
* @param neuronsList   A list of ordered neuron connections per branch per time step.
* @param currentBranch The current branch.
* @param currentTime   The current time step.
* @param maxExtension  Maximum extension of the branch, from the center of mass. Passed by reference.
**************************************************************************************************************/
double radiusOfGyration(Neuron *neurons, int  *N_connections, int **neuronsList, int currentBranch, int currentTime, double *maxExtension)
{
   // It needs at least 3 time steps or 2 connections to connect 3 neurons
   if(currentTime < 2 || N_connections[currentBranch] < 2)
   {
      return 0.0;
   }

   // Store center of mass coordinates
   double x_cm = 0.0;
   double y_cm = 0.0;
   int N = N_connections[currentBranch] + 1;


   // Calculates the center of mass of the chain
   for(int time = 0; time <= currentTime; time++)
   {
      if(time > 0 && neuronsList[currentBranch][time] == neuronsList[currentBranch][time - 1])
      {
         continue; // If the current neuron is the same, skip the calculations.
      }
      x_cm += neurons[neuronsList[currentBranch][time]].x;
      y_cm += neurons[neuronsList[currentBranch][time]].y;
   }
   x_cm /= N;
   y_cm /= N;

   // Calculates square distances
   double sumSquareDistances = 0.0;
   double maxSquareDistance = 0.0; // Track the maximum square distance

   for(int time = 0; time <= currentTime; time++)
   {
      if(time > 0 && neuronsList[currentBranch][time] == neuronsList[currentBranch][time - 1])
      {
         continue; // If the current neuron is the same, skip the calculations.
      }

      double dx = neurons[neuronsList[currentBranch][time]].x - x_cm;
      double dy = neurons[neuronsList[currentBranch][time]].y - y_cm;

      double squareDistance = dx * dx + dy * dy;
      sumSquareDistances += squareDistance;

      // Check for the fartest neuron from the center of mass
      if(squareDistance > maxSquareDistance)
      {
        maxSquareDistance = squareDistance;
      }
   }

   // Normalized maximum extension
   *maxExtension = sqrt(maxSquareDistance);

    // Mean square radius of gyration
    return sqrt(sumSquareDistances / N);
}