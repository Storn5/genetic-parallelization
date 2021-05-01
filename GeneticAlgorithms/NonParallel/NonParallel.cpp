#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <random>
#include <map>

using namespace std;

#define MAX_GENS 100
#define POP_SIZE 100
#define ALPH_SIZE 26 + 1
#define ASCII_a 97
#define MAX_TEXT_LEN 100
#define MUTATION_CHANCE 50 // Out of 100 individuals, MUTATION_CHANCE will mutate
#define ALPHA 0.65 // Weight of unigrams
#define BETA 0.35 // Weight of bigrams
#define GAMMA 0.0 // Weight of trigrams

const char ALPHABET[ALPH_SIZE+1] = "abcdefghijklmnopqrstuvwxyz ";

/// <summary>
/// An Individual contains a chromosome (25 genes) and fitness score.
/// <para> Each gene describes the position of a letter in the alphabet (default order is alphabetical). </para>
/// <para> Fitness describes the Individual's fitness (higher is better, default 0.0). </para>
/// </summary>
struct Individual
{
    string genes;
    long double fitness = 0;

    // Initialize Individual's genes to be in alphabetical order by default
    Individual()
    {
        genes = ALPHABET;
    }
    // Initialize Individual's genes explicitly
    Individual(string newGenes)
    {
        genes = newGenes;
    }
};

/// <summary>
/// Decipher text using key.
/// </summary>
/// <param name="ciphertext">Text to be deciphered</param>
/// <param name="key">Key (genes) to decipher the text</param>
/// <returns>Deciphered text</returns>
string decipher(const string& ciphertext, const string& key)
{
    string plaintext = string();
    string::size_type i;
    for (i = 0; i < ciphertext.length(); i++)
    {
        if (isalpha(ciphertext[i]))
            plaintext += key[ciphertext[i] - ASCII_a];
        else
            plaintext += ciphertext[i];
    }
    return plaintext;
}

/// <summary>
/// Count occurences of substring in text.
/// </summary>
/// <param name="substring">Substring to count</param>
/// <param name="text">Text to search</param>
/// <returns>Number of substring occurences in text</returns>
long countSubstring(const string& substring, const string& text)
{
    long num = 0;
    for (size_t i = 0; (i = text.find(substring, i)) != string::npos; num++, i++);
    return num;
}

/// <summary>
/// Count the number of all 1, 2 and 3-letter combinations in the text.
/// </summary>
/// <param name="text">Text to search</param>
/// <returns>Map of char combinations and their frequencies</returns>
const map<string, long double> countNgrams(const string& text)
{
    map<string, long double> ngrams;
    for (int i = 0; i < ALPH_SIZE; i++)
    {
        string key = string() + ALPHABET[i];
        long double count = countSubstring(key, text) / (double) text.length();
        ngrams.emplace(key, count);
        for (int j = 0; j < ALPH_SIZE; j++)
        {
            string key = string() + ALPHABET[i] + ALPHABET[j];
            count = countSubstring(key, text) / (double)text.length();
            ngrams.emplace(key, count);
            /*
            for (int k = 0; k < ALPH_SIZE; k++)
            {
                string key = string() + ALPHABET[i] + ALPHABET[j] + ALPHABET[k];
                count = countSubstring(key, text) / (double) text.length();
                ngrams.emplace(key, count);
            }
            */
        }
    }
    return ngrams;
}

/// <summary>
/// Randomly shuffle each Individual's genes.
/// </summary>
/// <param name="pop">Array of Individuals</param>
void initPopulation(Individual population[])
{
    for (int i = 0; i < POP_SIZE; i ++)
        random_shuffle(begin(population[i].genes), end(population[i].genes) - 1);
}

/// <summary>
/// Calculate fitness of every Individual in population.
/// </summary>
/// <param name="population">Array of Individuals</param>
/// <param name="fitnessSum">Sum of fitness values of population</param>
/// <param name="ciphertext">Text to be deciphered</param>
/// <param name="ngrams">Map of char combinations and their frequencies</param>
void calcFitness(Individual population[], long double* fitnessSum, const string& ciphertext, const map<string, long double> ngrams)
{
    for (int idx = 0; idx < POP_SIZE; idx++)
    {
        map<string, long double> cipherNgrams = countNgrams(decipher(ciphertext, population[idx].genes));
        long double fitness = 0.0;

        for (int i = 0; i < ALPH_SIZE; i++)
        {
            string key = string() + ALPHABET[i];
            fitness += ALPHA * (abs(cipherNgrams.at(key) - ngrams.at(key))); // TODO: fix the fitness function (why sometimes 0???)
            for (int j = 0; j < ALPH_SIZE; j++)
            {
                key = string() + ALPHABET[i] + ALPHABET[j];
                fitness += BETA * (abs(cipherNgrams.at(key) - ngrams.at(key)));
                /*
                for (int k = 0; k < ALPH_SIZE; k++)
                {
                    key = string() + ALPHABET[i] + ALPHABET[j] + ALPHABET[k];
                    fitness += GAMMA * (abs(cipherNgrams.at(key) - ngrams.at(key)));
                }
                */
            }
        }
        population[idx].fitness = (long double)1.0 / fitness;
        *fitnessSum += population[idx].fitness;
    }
}

/// <summary>
/// Weighted sampling of individuals into new population.
/// </summary>
/// <param name="population">Array of Individuals</param>
/// <param name="fitnessSum">Sum of fitness values</param>
void selection(Individual population[], long double* fitnessSum)
{
    Individual newPopulation[POP_SIZE];
    for (int i = 0; i < POP_SIZE; i++)
    {
        unsigned long long randomVal = rand() % (unsigned long long)(*fitnessSum * 1000); // TODO: fix the weighted sampling algorithm
        int index = 0;
        while (randomVal > 0)
        {
            randomVal -= (unsigned long long)(population[index].fitness * 1000);
            index ++;
        }
        newPopulation[i] = population[index - 1];
    }
    population = newPopulation;
}

/// <summary>
/// Mating using 2-point crossover between each 2 parents from the population.
/// </summary>
/// <param name="population">Array of Individuals</param>
void crossover(Individual population[])
{
    for (int i = 0; i < POP_SIZE; i += 2)
    {
        string child1(ALPH_SIZE, ' ');
        string child2(ALPH_SIZE, ' ');
        unsigned short gene1 = rand() % ALPH_SIZE;
        unsigned short gene2 = rand() % ALPH_SIZE;
        if (gene1 > gene2)
        {
            unsigned short temp = gene1;
            gene1 = gene2;
            gene2 = temp;
        }
        for (int j = gene1; j < gene2; j++)
        {
            child1[j] = population[i].genes[j];
            child2[j] = population[i + 1].genes[j];
        }
        unsigned short k1 = 0, k2 = 0;
        for (int j = 0; j < ALPH_SIZE; j++)
        {
            if (k1 < ALPH_SIZE && child1.find(population[i + 1].genes[j]) == std::string::npos)
            {
                while (child1[k1] != ' ')
                    k1++;
                child1[k1] = population[i + 1].genes[j];
            }
            if (k2 < ALPH_SIZE && child2.find(population[i].genes[j]) == std::string::npos)
            {
                while (child2[k2] != ' ')
                    k2++;
                child2[k2] = population[i].genes[j];
            }
        }
        population[i].genes = child1;
        population[i+1].genes = child2;
    }
}

/// <summary>
/// Swap 2 chars randomly.
/// </summary>
/// <param name="genes">Text to swap chars from</param>
void randomSwap(string genes)
{
    unsigned short gene1 = rand() % ALPH_SIZE;
    unsigned short gene2 = rand() % ALPH_SIZE;
    char gene = genes[gene1];
    genes[gene1] = genes[gene2];
    genes[gene2] = gene;
}

/// <summary>
/// Randomly choose individuals to be mutated (randomly swap 2 chars in each mutated individual).
/// </summary>
/// <param name="population">Array of Individuals</param>
void mutation(Individual population[])
{
    for (int i = 0; i < POP_SIZE; i++)
        if (rand() % 100 < MUTATION_CHANCE)
            randomSwap(population[i].genes);
}

/// <summary>
/// Print statistics about a generation.
/// </summary>
/// <param name="population">Array of Individuals</param>
/// <param name="fitnessSum">Sum of fitness values of population</param>
/// <param name="ciphertext">Text to be deciphered</param>
void printStats(Individual population[], long double* fitnessSum, const string& ciphertext)
{
    cout << "Avg fitness:\t\t" << *fitnessSum / POP_SIZE << "\n";
    long double maxFitness = 0;
    string bestGenes = ALPHABET;
    for (int i = 0; i < POP_SIZE; i++)
    {
        if (population[i].fitness > maxFitness)
        {
            maxFitness = population[i].fitness;
            bestGenes = population[i].genes;
        }
    }
    cout << "Max fitness:\t\t" << maxFitness << "\n";
    cout << "Best guess:\t\t" << decipher(ciphertext, bestGenes) << "\n";
    cout << "Best genes:\t\t" << bestGenes << "\n";
}

int main()
{
    // Initialization
    int gen = 0;
    Individual population[POP_SIZE];
    initPopulation(population);

    // Reading input from files
    ifstream ciphertextFile("ciphertext.txt"); // Text to be deciphered
    stringstream ciphertextBuffer;
    ciphertextBuffer << ciphertextFile.rdbuf();
    string ciphertext = ciphertextBuffer.str();
    ciphertextFile.close();

    ifstream sampletextFile("sampletext.txt"); // Text to calculate letter frequencies
    stringstream sampletextBuffer;
    sampletextBuffer << sampletextFile.rdbuf();
    string sampletext = sampletextBuffer.str();
    transform(sampletext.begin(), sampletext.end(), sampletext.begin(), ::tolower);
    sampletext.erase(std::remove_if(sampletext.begin(), sampletext.end(), [](const unsigned& c) { return !isspace(c) && !isalpha(c); }), sampletext.end());
    sampletextFile.close();

    // Calculating n-gram frequencies
    const map<string, long double> ngrams = countNgrams(sampletext);

    // Training loop
    while (gen < MAX_GENS)
    {
        gen += 1;
        cout << "Generation " << gen << "\n";

        long double fitnessSum = 0;
        calcFitness(population, &fitnessSum, ciphertext, ngrams);
        selection(population, &fitnessSum);
        crossover(population);
        mutation(population);

        printStats(population, &fitnessSum, ciphertext);
    }

    return 0;
}
