//
//  IsingSystem.cpp
//

#include "IsingSystem.h"
#include <fstream>
#include <string>

// colors
namespace colours {
	// note the f's here avoid warnings by telling C++ to read the numbers as floats and not doubles 
	GLfloat blue[] = { 0.1f, 0.3f, 0.9f, 1.0f };   // blue
	GLfloat red[] = { 1.0f, 0.2f, 0.1f, 0.2f };   // red
	GLfloat green[] = { 0.3f, 0.6f, 0.3f, 1.0f };   // green
}


// constructor
IsingSystem::IsingSystem(Window *set_win) {
	cout << "creating system, gridSize " << gridSize << endl;
	win = set_win;
	inverseTemperatureBeta = 0.2;
	slowNotFast = 1;
	isActive = 0;
	endSweeps = 50;
	seed = getSeed();
	numRuns = 1;
	endRuns = 20;
	r_correlation = 1;

	// Allocate memory for the grid, remember to free the memory in destructor
	//   the point here is that each row of the grid is an array
	//   the grid itself is a an array of pointers, one for each row
	// Here we allocate the array of pointers
	grid = new int*[gridSize];
	// Now allocate the indidual rows
	for (int i = 0; i<gridSize; i++) {
		grid[i] = new int[gridSize];
	}

	// this sets the temperatre and initialises the spins grid
	Reset();
}

void IsingSystem::Reset() {

	// double initialTemp = inverseTemperatureBeta;
	
	//resets number of sweeps
	numSweeps = 0;

	// setTemperature(initialTemp);

	// set the grid to -1
	for (int i = 0; i<gridSize; i++) {
		for (int j = 0; j<gridSize; j++) {
			// position is (i,j)
			int pos[2] = { i,j };
			// set this spin to state -1
			setGrid(pos, +1);
		}
	}
}


// destructor
IsingSystem::~IsingSystem() {
	// Close the file (if open)
	if (logfile.is_open())
		logfile.close();

	// Delete the window
	delete win;

	// Delete the grid
	// First we delete the individual rows
	for (int i = 0; i<gridSize; i++)
		delete[] grid[i];
	// Finally delete the array of pointers
	delete[] grid;
}

// this draws the system
void IsingSystem::DrawSquares() {

	double drawScale = 2.0 / (gridSize * 1.1);

	// draw the particles
	double halfSize = 0.5;
	int halfGrid = gridSize / 2;
	for (int x = 0; x<gridSize; x++) {
		for (int y = 0; y<gridSize; y++) {

			double vec[2];
			vec[0] = x - halfGrid;
			vec[1] = y - halfGrid;

			// openGL magic
			glPushMatrix();
			// choose a color
			if (grid[x][y] == -1)
				glColor4fv(colours::green);
			else
				glColor4fv(colours::blue);
			// draw a rectangle for the particle
			glRectd(drawScale*(vec[0] - halfSize),
				drawScale*(vec[1] - halfSize),
				drawScale*(vec[0] + halfSize),
				drawScale*(vec[1] + halfSize));
			// openGL magic
			glPopMatrix();
		}
	}

	// print some information (at top left)
	// this ostringstream is a way to make a string with numbers and words (similar to cout << ... )
	ostringstream str;
	str << "beta " << inverseTemperatureBeta << " size " << gridSize;
	win->displayString(str, -0.9, 0.94, colours::red);

}


// attempt N spin flips, where N is the number of spins
void IsingSystem::MCsweep() {
	for (int i = 0; i<gridSize*gridSize; i++)
		attemptSpinFlip();
}

// here we attempt to flip a spin and accept/reject with Metropolis rule
void IsingSystem::attemptSpinFlip() {
	int pos[2];

	// random site
	pos[0] = rgen.randomInt(gridSize);
	pos[1] = rgen.randomInt(gridSize);

	double hloc = computeLocalField(pos);
	
	double dE = 2.0 * hloc * readGrid(pos);
	if (dE<0)
		flipSpin(pos);
	else if (rgen.random01() < exp(-dE))
		flipSpin(pos);

}

// NOTE: this returns the local field *divided by the temperature* (dimensionless quantity)
double IsingSystem::computeLocalField(int pos[]) {
	double result = 0.0;
	for (int i = 0; i<4; i++) {
		int nborPos[2];
		setPosNeighbour(nborPos, pos, i);
		result += readGrid(nborPos);
	}
	result *= inverseTemperatureBeta;
	return result;
}

// set the value of a grid cell for a particular position
void IsingSystem::setGrid(int pos[], int val) {
	grid[pos[0]][pos[1]] = val;
}

// read the grid cell for a given position
int IsingSystem::readGrid(int pos[]) {
	return grid[pos[0]][pos[1]];
}

// flips spin of grid cell for a given position
void IsingSystem::flipSpin(int pos[]) {
	grid[pos[0]][pos[1]] = -grid[pos[0]][pos[1]];
}

//calculates number of spins 
int IsingSystem::numSpins(int gridSize) {
	return gridSize*gridSize;
}

// calculates magnetisation of whole grid
float IsingSystem::getMagnetisation() {
	float M = 0;
	for (int i = 0; i < gridSize; i++)
	{
		for (int j = 0; j < gridSize; j++)
		{
			M += grid[i][j];
		}
	}
	return (M / (gridSize*gridSize));
}

// gets Energy by summing product of nearest neighbours for each particle
float IsingSystem::getEnergy() {
    float E = 0;
    int neighbour[2], current[2];

    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            current[0] = i;
            current[1] = j;
            for (int k = 0; k < 4; k++) {
                setPosNeighbour(neighbour, current, k);
                E += grid[current[0]][current[1]] * grid[neighbour[0]][neighbour[1]];
            }
        }
    }
    return (-E / (gridSize*gridSize));
}

float IsingSystem::getCorrelation(int r) {
    int neighbour[2], current[2];
    float correlation = 0;

    // Iterate over each point on the grid
    for (int i = 0; i < gridSize; i++) {
        for (int j = 0; j < gridSize; j++) {
            // Set current point
            current[0] = i;
            current[1] = j;

            // Find neighbour in column r distance away
            neighbour[0] = (current[0] + r) % gridSize;
            neighbour[1] = current[1];

            // Calculate the product of the spins of the two sites and add to the total correlation
            correlation += grid[current[0]][current[1]] * grid[neighbour[0]][neighbour[1]];
        }
    }

    // Normalize the correlation by the total number of points on the grid
    correlation /= (gridSize * gridSize);

    return correlation;
}


// send back the position of a neighbour of a given grid cell
// NOTE: we take care of periodic boundary conditions, also positions are integers now not doubles
void IsingSystem::setPosNeighbour(int setpos[], int pos[], int val) {
	switch (val) {
	case 0:
		setpos[0] = (pos[0] + 1) % gridSize;
		setpos[1] = pos[1];
		break;
	case 1:
		setpos[0] = (pos[0] - 1 + gridSize) % gridSize;
		setpos[1] = pos[1];
		break;
	case 2:
		setpos[0] = pos[0];
		setpos[1] = (pos[1] + 1) % gridSize;
		break;
	case 3:
		setpos[0] = pos[0];
		setpos[1] = (pos[1] - 1 + gridSize) % gridSize;
		break;
	}
}

int IsingSystem::getSeed() {
	return seed;
}

// gets file name of csv file
std::string IsingSystem::getFileName(std::string indVar, double depVar, int seed) {
	// sets seed to string data type	
	std::string seedAsString = std::to_string(seed);
	std::string betaAsString = std::to_string(depVar);

	// r_correlation as string
	std::string rAsString = std::to_string(r_correlation);

	//creates filename based on inputs
	std::string filename = "task5_2_data/file_" + betaAsString + "_" + rAsString + ".csv";
	return filename;
}

// creates csv file with dependant variable and seed as name
void IsingSystem::csvHeaders(std::string indVar, double depVar, int seed) {
	//sets filename
	std::string filename = getFileName(indVar, depVar, seed);
	//creates file with filename
	std::ofstream file(filename);
	//labels columns of filename
	file << indVar << "," << "beta" << "," << "magnetisation" << "," << "energy" << "," << "G" << "," << "seed" << ","<< "r_correlation" << endl;
    // Close the file
    file.close();
}

// prints data to csv file
void IsingSystem::printCsv(std::string filename, float indVar, double indVar2, float depVar, float depVar2, float depVar3, int seed, int r) {
	//open csv
	std::ofstream logfile(filename, std::ios_base::app);
	//print data to file
	if (logfile.is_open()) {
		//write to file
		logfile << indVar << "," << indVar2 << "," << depVar << "," << depVar2 << "," << depVar3 << "," << seed << ","<< r << std::endl;
		logfile.close();
	}
	else {
		//couldn't open file for writing
		std::cerr << "Error: Unable to open the file for writing." << std::endl;
	}
}

void IsingSystem::calcVars(std::string filename, int numSweeps) {
	//calculate magnetisation and energy after 10 sweeps
	if ((numSweeps % 1) == 0)
	{
		float M = getMagnetisation();
		float E = getEnergy();
		float G = getCorrelation(r_correlation);
		seed = getSeed();
		printCsv(filename, numSweeps, inverseTemperatureBeta, M, E, G, seed, r_correlation);
	}	
}

// ends automation if endRuns is reached
void IsingSystem::keepGoing() {
	
	if (numSweeps == 0)
	{
		//create new file
		csvHeaders("sweeps", inverseTemperatureBeta , seed);
		fileName = getFileName("sweeps", inverseTemperatureBeta, seed);
		correlation = 0;
		r_correlation = 0;

	}
	if (numSweeps <= endSweeps)
	{
		if (numSweeps > 10)
		{
			r_correlation++;
		}
		fileName = getFileName("sweeps", inverseTemperatureBeta, seed);
		calcVars(fileName, numSweeps);
		MCsweep();
		numSweeps++;
	}
	else if (numRuns < endRuns)
	{
		//increment seed, numRuns counter 
		// seed++;
		
		
		numRuns++;
		Reset();
	}
	else {
		pauseRunning();
		correlation = 0;
		cout << "End number of runs reached" << endl;
	}	
}

// this is the update function which at the moment just does one mc sweep
void IsingSystem::Update() {
	//keeps going if endSweeps not yet reached
	keepGoing();
}

