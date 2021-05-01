#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <random>
#include <map>

using namespace std;

#define MAX_GENS 25
#define POP_SIZE 20
#define ALPH_SIZE 26
#define ASCII_a 97
#define MAX_TEXT_LEN 100
#define MUTATION_CHANCE 50 // Out of 100 individuals, MUTATION_CHANCE will mutate

const char ALPHABET[ALPH_SIZE+1] = "abcdefghijklmnopqrstuvwxyz";

/// <summary>
/// An Individual contains a chromosome (25 genes) and fitness score.
/// <para> Each gene describes the position of a letter in the alphabet (default order is alphabetical). </para>
/// <para> Fitness describes the Individual's fitness (higher is better, default 0.0). </para>
/// </summary>
struct Individual
{
    char genes[ALPH_SIZE+1];
    unsigned long long fitness = 0;

    // Initialize every Individual's genes to be in alphabetical order by default
    Individual()
    {
        snprintf(genes, ALPH_SIZE+1, ALPHABET);
    }
};

/// <summary>
/// Decipher text using key
/// </summary>
/// <param name="ciphertext">Text to be deciphered</param>
/// <param name="key">Key (genes) to decipher the text</param>
/// <returns>Deciphered text</returns>
string decipher(const string& ciphertext, char* key)
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
/// Count occurences of substring in text
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
/// Count the number of all 2- and 3- letter combinations in the text
/// </summary>
/// <param name="text">Text to search</param>
/// <returns>Map of char combinations and their frequencies</returns>
const map<string, unsigned long> countNgrams(const string& text)
{
    map<string, unsigned long> ngrams;
    for (int i = 0; i < ALPH_SIZE; i++)
    {
        for (int j = 0; j < ALPH_SIZE; j++)
        {
            string key = string() + ALPHABET[i] + ALPHABET[j];
            long count = countSubstring(key, text);
            ngrams.emplace(key, count);
            for (int k = 0; k < ALPH_SIZE; k++)
            {
                string key = string() + ALPHABET[i] + ALPHABET[j] + ALPHABET[k];
                count = countSubstring(key, text);
                ngrams.emplace(key, count);
            }
        }
    }
    return ngrams;
}

/// <summary>
/// Randomly shuffle each Individual's genes
/// </summary>
/// <param name="pop">Array of Individuals</param>
void initPopulation(Individual population[])
{
    for (int i = 0; i < POP_SIZE; i ++)
        random_shuffle(begin(population[i].genes), end(population[i].genes) - 1);
}

/// <summary>
/// Calculate fitness of every Individual in population
/// </summary>
/// <param name="population">Array of Individuals</param>
/// <param name="fitnessSum">Sum of fitness values of population</param>
/// <param name="ciphertext">Text to be deciphered</param>
/// <param name="ngrams">Map of char combinations and their frequencies</param>
void calcFitness(Individual population[], unsigned long long* fitnessSum, const string& ciphertext, const map<string, unsigned long> ngrams)
{
    for (int i = 0; i < POP_SIZE; i++)
    {
        map<string, unsigned long> cipherNgrams =
            countNgrams(decipher(ciphertext, population[i].genes));

        unsigned long long fitness = 0;
        for (map<string, unsigned long>::iterator iter = cipherNgrams.begin(); iter != cipherNgrams.end(); ++iter)
        {
            fitness += iter->second * ngrams.at(iter->first);
        }
        population[i].fitness = fitness;
        *fitnessSum += fitness;
    }
}

void selection(Individual population[], unsigned long long* fitnessSum)
{
    Individual newPopulation[POP_SIZE];
    for (int i = 0; i < POP_SIZE; i++)
    {
        long long randomVal = rand() % *fitnessSum;
        int index = 0;
        while (randomVal > 0)
        {
            randomVal -= population[index].fitness;
            index ++;
        }
        newPopulation[i] = population[index - 1];
    }
    population = newPopulation;
}

void crossover(Individual population[])
{
    Individual newPopulation[POP_SIZE];

    population = newPopulation;
}

void swap(char* genes)
{
    unsigned short gene1 = rand() % ALPH_SIZE;
    unsigned short gene2 = rand() % ALPH_SIZE;
    char gene = genes[gene1];
    genes[gene1] = genes[gene2];
    genes[gene2] = gene;
}

void mutation(Individual population[])
{
    for (int i = 0; i < POP_SIZE; i++)
        if (rand() % 100 < MUTATION_CHANCE)
            swap(population[i].genes);
}

/// <summary>
/// Print statistics about a generation
/// </summary>
/// <param name="population">Array of Individuals</param>
/// <param name="fitnessSum">Sum of fitness values of population</param>
/// <param name="ciphertext">Text to be deciphered</param>
void printStats(Individual population[], unsigned long long* fitnessSum, const string& ciphertext)
{
    cout << "Avg fitness:\t\t" << *fitnessSum / POP_SIZE << "\n";
    unsigned long long maxFitness = 0;
    char* bestGenes = (char*)ALPHABET;
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
    sampletextFile.close();

    // Calculating n-gram frequencies
    const map<string, unsigned long> ngrams = countNgrams(sampletext);

    // Training loop
    while (gen < MAX_GENS)
    {
        gen += 1;
        cout << "Generation " << gen << "\n";

        unsigned long long fitnessSum = 0;
        calcFitness(population, &fitnessSum, ciphertext, ngrams);
        selection(population, &fitnessSum);
        crossover(population);
        mutation(population);

        printStats(population, &fitnessSum, ciphertext);
    }

    return 0;
}
