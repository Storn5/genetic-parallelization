#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <random>
#include <time.h>
#include <map>

using namespace std;

#define ASCII_a 97 // DON'T CHANGE (ASCII character a = 97)

#define MAX_GENS 50 // Number of generations
#define POP_SIZE 250 // Number of individuals in population
#define CROSSOVER_CHANCE 50 // Out of 100 parents, CROSSOVER_CHANCE will crossover (the rest will just be copied)
#define MUTATION_CHANCE 50 // Out of 100 children, MUTATION_CHANCE will mutate (the rest will just be copied)
#define ALPHA 0.5 // Weight of unigrams
#define BETA 0.25 // Weight of bigrams
#define GAMMA 0.25 // Weight of trigrams
#define PROPORTIONAL_SELECTION true // Selection proportional to weights (if false, then just kill bottom 50% of population)
#define CUSTOM_FREQ false // Whether frequencies should be calculated from sample text file, or from predefined numbers (false is MUCH faster)
#define FULL_BIGRAMS false // Whether to use the larger full_bigrams.txt file or the smaller bigrams.txt (false is faster)

#define ALPH_SIZE 10 // Length of the key/alphabet to guess
const char ALPHABET[ALPH_SIZE + 1] = "abcdefghij"; // alphabet/key to guess (original state, not correct)

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
        if (isalpha(ciphertext[i]) && ciphertext[i] - ASCII_a < ALPH_SIZE)
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
long double countSubstring(const string& substring, const string& text)
{
    long double num = 0;
    for (size_t i = 0; (i = text.find(substring, i)) != string::npos; num++, i++);
    return num;
}

/// <summary>
/// Count number of n-grams in text.
/// </summary>
/// <param name="text">Text to search</param>
/// <param name="n">n for the n-gram</param>
/// <returns>Number of n-grams in text</returns>
long double numNgrams(const string& text, unsigned short n)
{
    n--;
    long double num = 0;
    for (long i = 0; i < text.length() - n; i++)
    {
        num++;
        for (long j = i; j < i + n + 1; j++)
        {
            if (isspace(text[j]))
                num--;
        }
    }
    return num;
}

/// <summary>
/// Count the number of all 1, 2 and 3-letter combinations in the text.
/// </summary>
/// <param name="text">Text to search</param>
/// <param name="sampleNgrams">Map of char combinations and their frequencies (from sample text)</param>
/// <returns>Map of char combinations and their frequencies</returns>
const map<string, long double> countNgrams(const string& text, map<string, long double> sampleNgrams, long double num_unigrams, long double num_bigrams, long double num_trigrams)
{
    map<string, long double> ngrams;
    for (map<string, long double>::iterator it = sampleNgrams.begin(); it != sampleNgrams.end(); it++)
    {
        if (it->first.length() == 1)
            ngrams.emplace(it->first, (long double)(countSubstring(it->first, text) / num_unigrams));
        else if (it->first.length() == 2)
            ngrams.emplace(it->first, (long double)(countSubstring(it->first, text) / num_bigrams));
        else
            ngrams.emplace(it->first, (long double)(countSubstring(it->first, text) / num_trigrams));
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
void calcFitness(Individual population[], long double* fitnessSum, const string& ciphertext, map<string, long double> ngrams, long double num_unigrams, long double num_bigrams, long double num_trigrams)
{
    for (int idx = 0; idx < POP_SIZE; idx++)
    {
        map<string, long double> cipherNgrams = countNgrams(decipher(ciphertext, population[idx].genes), ngrams, num_unigrams, num_bigrams, num_trigrams);
        long double fitness = 0.0;
        for (map<string, long double>::iterator it = ngrams.begin(); it != ngrams.end(); it++)
        {
            if (it->first.length() == 1)
                fitness += ALPHA * abs(cipherNgrams.at(it->first) - ngrams.at(it->first));
            else if (it->first.length() == 2)
                fitness += BETA * abs(cipherNgrams.at(it->first) - ngrams.at(it->first));
            else
                fitness += GAMMA * abs(cipherNgrams.at(it->first) - ngrams.at(it->first));
        }

        population[idx].fitness = (long double)pow((long double)pow(1 - (fitness / 4), 8), 3);
        *fitnessSum += population[idx].fitness;
    }
}

/// <summary>
/// Weighted sampling of individuals into new population.
/// </summary>
/// <param name="population">Array of Individuals</param>
/// <param name="fitnessSum">Sum of fitness values</param>
void selection(Individual population[POP_SIZE], long double* fitnessSum)
{
    Individual newPopulation[POP_SIZE];
    if (PROPORTIONAL_SELECTION)
    {
        for (int i = 0; i < POP_SIZE; i++)
        {
            int randomVal = rand() % (int)(*fitnessSum * 10000);
            int index = 0;
            while (randomVal > 0 && index < POP_SIZE)
            {
                randomVal -= (int)(population[index].fitness * 10000);
                index++;
            }
            if (index == 0)
                index++;
            newPopulation[i] = population[index - 1];
        }
    }
    else
    {
        long double avgFitness = *fitnessSum / (long double)POP_SIZE;
        int i = 0, j = 0;
        while (i < POP_SIZE)
        {
            for (j = 0; j < POP_SIZE; j++)
            {
                if (population[j].fitness > avgFitness)
                {
                    newPopulation[i] = population[j];
                    i++;
                    if (i >= POP_SIZE)
                        break;
                }
                if (i >= POP_SIZE)
                    break;
            }
        }
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
        if (rand() % 100 >= CROSSOVER_CHANCE)
            continue;
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
void printStats(Individual population[], long double* fitnessSum, const string& ciphertext, const string& plaintext)
{
    //cout << "Avg fitness:\t\t" << *fitnessSum / POP_SIZE << "\n";
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
    const string deciphered = decipher(ciphertext, bestGenes);
    //cout << "Max fitness:\t\t" << maxFitness << "\n";
    //cout << "Best guess:\t\t" << deciphered << "\n";
    cout << "Best genes:\t\t" << bestGenes << "\n";

    long double bestAccuracy = 0.0;
    string::size_type i;
    for (i = 0; i < deciphered.length(); i++)
        if (deciphered[i] == plaintext[i])
            bestAccuracy ++;
    bestAccuracy /= deciphered.length();
    cout << "Best accuracy:\t\t" << bestAccuracy << "\n";
}

int main()
{
    // Initialization
    int gen = 0;
    Individual population[POP_SIZE];
    initPopulation(population);
    map<string, long double> ngrams;
    long double num_unigrams, num_bigrams, num_trigrams;

    // Reading input from files
    ifstream plaintextFile("../../plaintext.txt"); // Target text
    stringstream plaintextBuffer;
    plaintextBuffer << plaintextFile.rdbuf();
    string plaintext = plaintextBuffer.str();
    plaintext.erase(std::remove_if(plaintext.begin(), plaintext.end(), [](const unsigned& c) { return !isspace(c) && !isalpha(c); }), plaintext.end());
    plaintextFile.close();
    ifstream ciphertextFile("../../ciphertext.txt"); // Text to be deciphered
    stringstream ciphertextBuffer;
    ciphertextBuffer << ciphertextFile.rdbuf();
    string ciphertext = ciphertextBuffer.str();
    ciphertext.erase(std::remove_if(ciphertext.begin(), ciphertext.end(), [](const unsigned& c) { return !isspace(c) && !isalpha(c); }), ciphertext.end());
    ciphertextFile.close();
    num_unigrams = numNgrams(ciphertext, 1);
    num_bigrams = numNgrams(ciphertext, 2);
    num_trigrams = numNgrams(ciphertext, 3);

    if (CUSTOM_FREQ) // Calculating frequencies from sample text file
    {
        ifstream sampletextFile("../../sampletext.txt"); // Text to calculate letter frequencies
        stringstream sampletextBuffer;
        sampletextBuffer << sampletextFile.rdbuf();
        string sampletext = sampletextBuffer.str();
        transform(sampletext.begin(), sampletext.end(), sampletext.begin(), ::tolower);
        sampletext.erase(std::remove_if(sampletext.begin(), sampletext.end(), [](const unsigned& c) { return !isspace(c) && !isalpha(c); }), sampletext.end());
        sampletextFile.close();

        // Calculating n-gram frequencies
        for (int i = 0; i < ALPH_SIZE; i++)
        {
            string key = string() + ALPHABET[i];
            ngrams.at(key) = 0.0;
            for (int j = 0; j < ALPH_SIZE; j++)
            {
                key = string() + ALPHABET[i] + ALPHABET[j];
                ngrams.at(key) = 0.0;
                for (int k = 0; k < ALPH_SIZE; k++)
                {
                    key = string() + ALPHABET[i] + ALPHABET[j] + ALPHABET[k];
                    ngrams.at(key) = 0.0;
                }
            }
        }
        ngrams = countNgrams(sampletext, ngrams, num_unigrams, num_bigrams, num_trigrams);
    }
    else // Reading frequencies from predefined numbers
    {
        string filenames[] = { "../../unigrams.txt", "../../bigrams.txt", "../../trigrams.txt" };
        if (FULL_BIGRAMS)
            filenames[1] = "../../full_bigrams.txt";
        for (string filename : filenames)
        {
            ifstream frequencyFile(filename);
            string key;
            long double number, frequency;
            if (FULL_BIGRAMS && filename == "../../full_bigrams.txt")
                while (frequencyFile >> key >> frequency)
                    ngrams.emplace(key, frequency / 100.0);
            else
                while (frequencyFile >> key >> number >> frequency)
                    ngrams.emplace(key, frequency / 100.0);
            frequencyFile.close();
        }
    }

    // Training loop
    clock_t time;
    time = clock();
    while (gen < MAX_GENS)
    {
        gen += 1;
        cout << "Generation " << gen << "\n";

        long double fitnessSum = 0;
        calcFitness(population, &fitnessSum, ciphertext, ngrams, num_unigrams, num_bigrams, num_trigrams);
        selection(population, &fitnessSum);
        crossover(population);
        mutation(population);

        printStats(population, &fitnessSum, ciphertext, plaintext);
    }
    time = clock() - time;
    cout << "Time for " << MAX_GENS << " generations and " << POP_SIZE << " population: " << (float)time / CLOCKS_PER_SEC << " sec.\n";

    return 0;
}
