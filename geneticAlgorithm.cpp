#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <fstream>
#include <time.h>
#include <chrono>
#include <climits>
#include <math.h>
#include <random>
#include <iomanip>

#define TAX_OF_MUTATION 0.1 // In percentage
#define POPULATION_SIZE 10
#define MAX_GENERATIONS 1000
#define GENS_TO_BALANCE 30
#define NUM_OF_GENES 1

using namespace std;

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

//	Returns a random integer between INT_MIN and INT_MAX
int randomInt(int lowerBound = INT_MIN, int upperBound = INT_MAX)
{
	return uniform_int_distribution<int>(lowerBound, upperBound)(rng);
}

//	Returns a random number between 0 and 1
float randomPercentage()
{
	return ((double) randomInt() + INT_MAX) / ((double) INT_MAX - INT_MIN); 
}

//	Returns a random sign
int randomSign()
{
	return (randomPercentage() < 0.5) ? -1 : 1;
}

class Individual
{

public : 
	
	//	gene[i] belongs to the interval genesLowerBound[i] to genesUpperBound[i]
	vector<float> genes;
	vector<float> genesLowerBound;
	vector<float> genesUpperBound;
	int fitness;

	Individual(	vector<float> _genes = vector<float>(0),
				vector<float> _genesLowerBound = vector<float>(0),
				vector<float> _genesUpperBound = vector<float>(0)) 
	{
		
		genes = _genes;
		genesLowerBound = _genesLowerBound;
		genesUpperBound = _genesUpperBound;
		fitness = -1;
	}

	int getFitnessScore()
	{
		if(fitness != -1)
			return fitness;

		string commandLine = "./jsonFileWriter";
		for(float gene : genes)
			commandLine += " " + to_string(gene);
			
		//	Running jsonFileWriter to set the genes parameters in the json file
		system(commandLine.c_str());

		FILE *fpipe;
		char answer[256];

		//	Running Hardware-Accelerated-Ant-Colony-Based simulator
		char command[] = "./main";
		if(!(fpipe = (FILE*) popen(command, "r")))
		{
			cerr << "Failed to read pipe";
			exit(1);
		}

		//	Reading the return of the simulation 
		fgets(answer, sizeof(answer), fpipe);
		pclose(fpipe);

		fitness = stof(answer);
		return fitness;
	}

	//	Returns the indexes of the genes that need to mutate
	vector<int> getGenesIndexesToMutate()
	{
		int numberOfIndexes = 0;
		float randPercentage = randomPercentage();

		if (randPercentage < 0.7)
			numberOfIndexes = 1;
		else if (randPercentage < 0.95)
			numberOfIndexes = 2;
		else
			numberOfIndexes = 3;

		numberOfIndexes = min(numberOfIndexes, (int) genes.size());

		vector<int> indexes((int) genes.size());
		for(int i = 0; i < (int) genes.size(); i++)
			indexes[i] = i;
		random_shuffle(indexes.begin(), indexes.end());

		while((int) indexes.size() > numberOfIndexes)
			indexes.pop_back();

		return indexes;
	}

	//	Method that mutates genes of the individual
	void mutateGenes()
	{
		vector<int> indexes = getGenesIndexesToMutate();

		for(int i : indexes)
		{
			float lw = genesLowerBound[i];
			float up = genesUpperBound[i];
			genes[i] += (double) randomSign() * randomPercentage() * TAX_OF_MUTATION * (up - lw);
			genes[i] = min(genes[i], up);
			genes[i] = max(genes[i], lw);
		}
	}

	//	Operator that checks if the individual is less apt than other
	bool operator < (Individual other)
	{	return (*this).getFitnessScore() < other.getFitnessScore();	}

	//	Operator that checks if the individual is more apt than other
	bool operator > (Individual other)
	{	return (*this).getFitnessScore() > other.getFitnessScore();	}

	//	Operator that checks if the individual has the same aptness than other
	bool operator == (Individual other)
	{	return (*this).getFitnessScore() == other.getFitnessScore();	}

	//	Operator that simulates the crossover between two individuals
	Individual operator + (Individual other)
	{
		Individual child = (*this);
		child.fitness = -1;

		for(int i = 0; i < (int) child.genes.size(); i++)
			child.genes[i] = (((*this).genes[i] + other.genes[i]) / 2.0);

		return child;
	}
};

class GeneticAlgorithm
{

public : 

	vector<Individual> individuals;

	GeneticAlgorithm(vector<Individual> _individuals = vector<Individual>(0))
	{
		individuals = _individuals;
	}

	void initializePopulation(int populationSize, Individual base)
	{
		individuals.clear();
		individuals.resize(populationSize);
		for(int i = 0; i < populationSize; i++)
		{
			/*
			for(int j = 0; j < (int) base.genes.size(); j++)
			{
				float lw = base.genesLowerBound[j];
				float up = base.genesUpperBound[j];
				double gene = (double) randomPercentage() * (up - lw);
				base.genes[j] = gene;
				base.genes[j] = min(base.genes[j], up);
				base.genes[j] = max(base.genes[j], lw);
			}
			*/

			individuals[i] = base;
			individuals[i].mutateGenes();
		}
	}

	//	Return the current index of the best individual
	int getBestIndividualIndex()
	{
		int bestIndividualIndex = 0;
		
		for(int i = 1; i < (int) individuals.size(); i++)
			if(individuals[bestIndividualIndex] < individuals[i])
				bestIndividualIndex = i;

		return bestIndividualIndex;
	}

	//	Best individual becomes the first individual of the vector
	void realocateBestIndividual()
	{
		int index = getBestIndividualIndex();
		swap(individuals[0], individuals[index]);
	}

	void elitism()
	{
		realocateBestIndividual();

		for(int i = 1; i < (int) individuals.size(); i++)
		{
			individuals[i] = (individuals[0] + individuals[i]);
			individuals[i].mutateGenes();
		}
	}
};

int main ()
{
	//	Compiling Hardware-Accelerated-Ant-Colony-Based simulator
	system("make all");

	//	Compiling jsonFileWriter
	system("g++ jsonFileWriter.cpp -o jsonFileWriter");

	
	// ----------- Declaring base individual --------- //
	Individual base;

	base.genes.resize(NUM_OF_GENES);
	base.genesLowerBound.resize(NUM_OF_GENES);
	base.genesUpperBound.resize(NUM_OF_GENES);

	//	Velocity Gene
	base.genes[0] = 0.0004;
	base.genesLowerBound[0] = 0.0002;
	base.genesUpperBound[0] = 0.002;

	// //	PlacePheromoneIntensity Gene
	// base.genes[1] = 60;
	// base.genesLowerBound[1] = 0;
	// base.genesUpperBound[1] = 255;

	base.fitness = -1;
	// ----------------------------------------------- //


	GeneticAlgorithm population;
	population.initializePopulation(POPULATION_SIZE, base);

	ofstream myfile;
	myfile.open("evolution.txt");

	vector<float> bestGenes(1, 0);
	int flag = 0;

	for(int i = 0; i < MAX_GENERATIONS; i++)
	{
		for(int j = 0; j < (int) population.individuals.size(); j++)
			population.individuals[j].getFitnessScore();

		population.realocateBestIndividual();

		for(int j = 0; j < (int) population.individuals.size(); j++)
		{
			myfile << i << " ";
			cout << i << " ";
			for(float gene : population.individuals[j].genes)
			{
				myfile << setprecision(6) << gene << " ";
				cout << setprecision(6) << gene << " ";
			}

			myfile << population.individuals[j].getFitnessScore() << "\n";
			cout << population.individuals[j].getFitnessScore() << "\n";
		}

		if(bestGenes == population.individuals[0].genes)
			flag++;
		else
		{
			flag = 0;
			bestGenes = population.individuals[0].genes;
		}

		if(flag >= GENS_TO_BALANCE)
			break;

		population.elitism();
	}

	if (flag == GENS_TO_BALANCE)
	{
		myfile << "Ended by estabilization.\n";
		cout << "Ended by estabilization.\n";
	} 
	else 
	{
		myfile << "Ended by max generations.\n";
		cout << "Ended by max generations.\n";
	}

	myfile.close();

	cout << "The best individual has the following genes:" << endl;
	for(float gene : bestGenes)
		cout << gene << " ";
	cout << "\n";

	return 0;
}