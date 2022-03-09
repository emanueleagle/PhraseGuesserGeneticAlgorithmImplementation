
#add the ability to input a word?
#do some graphing to determine how many generations it takes to guess the word

import random
import os
import time 

class GeneticAlgo:

    def __init__(self, userInput, input=[]):
        try:
            self.input = list(open(input, 'r').readlines())
        except TypeError:
            self.input = input
        self.userInput = userInput
        self.inputLength = len(input)
        self.parentChromosomes = []
        self.bestChromosome = [] 
        self.bestChromosomeFitness = []
        self.averageFitness = []
        self.childChromosomes = []
        self.chromosomeSize = len(self.userInput)
        self.numOffSpringProduced = 4
        self.transposableRate = 1
        self.transposonSize = 0
        self.populationSize = self.numOffSpringProduced*200
        self.crossoverProb = 1
        self.mutationProb = 0.0001
        self.chromosomeFitnessScores = [0]*self.populationSize
        self.chromosomeFitnessRatios = [0]*self.populationSize
    
    def get_parents(self):
        """get parent chromosomes."""
        return self.parentChromosomes

    def fitness_function(self):
        """fitness function used for measuring performance of individual chromosome."""
        binary = self.userInput
        for i in range(len(self.chromosomeFitnessScores)):
            score = 0
            for b in range(len(binary)):
                if(binary[b] == self.parentChromosomes[i][b]):
                    score += 2
                elif(binary[b] != self.parentChromosomes[i][b]):
                    score -= 1
            self.chromosomeFitnessScores[i] = score

    def set_parents(self):
        """Creates a random population of chromosomes if no input is given."""
        if(self.inputLength <= 1):
            for i in range(self.populationSize):
                chromosome = []
                for i in range(self.chromosomeSize):
                    if(random.random() >= 0.5):
                        chromosome.append(1)
                    else:
                        chromosome.append(0)
                self.parentChromosomes.append(chromosome)
        else:
            for chromosome in self.input:
                self.parentChromosomes.append(chromosome)

    def fitness_ratios(self):
        """Gets the list of fitness scores and calculates the fitness ratios for each chromosome."""
        totalFitnessScore = sum(self.chromosomeFitnessScores)
        for i in range(len(self.chromosomeFitnessScores)):
            ratio = round((self.chromosomeFitnessScores[i]/totalFitnessScore)*100,2)
            self.chromosomeFitnessRatios[i] = ratio

    def mating(self):
        """Mates two chromsomes and generates their offspring."""
        self.childChromsomes = []
        fitnessRatios = self.chromosomeFitnessRatios #get the fitness ratios
        offSpring = [] #list stores the offspring for this generation
        mates = 2 #number of mates that will be reproducing
        rouletteWheel = [] #list to store the roulette wheel ratios
        for i in range(len(fitnessRatios)): #for each fitness ratio in the list of ratios
            if(i == 0): #if the first index
                portion = fitnessRatios[i] #the portion is just equal to the first fitness ratio
                rouletteWheel.append(round(portion,3)) #round and append to the roulette wheel list
            else:
                portion = rouletteWheel[i-1] + fitnessRatios[i] #if not the first index, take the prior index and add the current ratio to create rnages
                rouletteWheel.append(round(portion,3)) #round and append to the roulette wheel
        numOffSpring = len(offSpring)
        while(numOffSpring+2 <= self.populationSize):
            pairOfChromosomes = [] #list to store the pair of chromosomes
            for i in range(mates): #for each mating chromosome
                selector = random.randrange(0,100) #set a random value between 0 and 100
                for ratio in rouletteWheel: #for each ratio in the roulette wheel
                    if selector < ratio: #if the random value is less than the ratio but greater than the prior ratio
                        chromosomeNum = rouletteWheel.index(ratio) #get the index of that ratio, which is the chromosome's number
                        pairOfChromosomes.append(chromosomeNum) #append chromsome number to the pair of chromosomes
                        break
            c1 = self.parentChromosomes[pairOfChromosomes[0]] #define c1 as the first chromosome, index into parent chromosomes 
            c2 = self.parentChromosomes[pairOfChromosomes[1]] #define c2 as the second chromsome, index into the parent chromosomes
            c1 = self.mutation(c1) #run mutation operator
            c2 = self.mutation(c2) #run mutation operator
            crossoverInds = self.crossover(c1, c2) #run the crossover function, return the two new chromosomes
            newOffSpring = self.incompleteDominance(crossoverInds[0],crossoverInds[-1])
            for i in range(len(newOffSpring)):
                offSpring.append(newOffSpring[i])
            #offSpring.append(newOffSpring[0])
            #offSpring.append(newOffSpring[1])
            numOffSpring = len(offSpring)
        self.childChromosomes = offSpring 

    def incompleteDominance(self, c1, c2):
        offSpring = [0]*self.numOffSpringProduced
        for i in range(self.numOffSpringProduced):
            if(i == 0):
                offSpring[i] = c1
            elif(i == (self.numOffSpringProduced-1)):
                offSpring[i] = c2 
            else:
                blendedIndividual = [0]*len(c1)
                for a in range(len(blendedIndividual)):
                    if(c1[a] == c2[a]):
                        blendedIndividual[a] = c1[a]
                    else:
                        blendedIndividual[a] = random.randrange(0,2)
                offSpring[i] = blendedIndividual
        return offSpring

    def crossover(self, c1, c2):
        """Applies crossover operator to two chromosomes: c1 = chromosome1, c2 = chromosome2."""
        randomValue = random.random()
        if randomValue <= self.crossoverProb:
            #crossover
            breakPoint = random.randrange(0, len(c1)) #find a random point at which to break chromomsomes
            while breakPoint == 0:
                breakPoint = random.randrange(0,len(c1)) #if the breakpoint is zero, cloning occurs, so continue to generate values while value is zero
            c1CrossoverRegion = c1[breakPoint:] #get the crossover region for the first chromosome
            c1NonCrossoverRegion = c1[:breakPoint] #get the region of the first chromosome that won't be crossing over 
            c2CrossoverRegion = c2[breakPoint:] #get the crossover region for the second chromosome 
            c2NonCrossoverRegion = c2[:breakPoint] #get the region of the second chromosome that wont be crossing over
            offSpring1 = c1NonCrossoverRegion+c2CrossoverRegion #concatenate the noncrossing over region of chromosome 1 with the crossing over region of chromosome 2
            offSpring2 = c2NonCrossoverRegion+c1CrossoverRegion #concatenate the noncrossing over region of chromosome 2 with the crossing over region of chromosome 1 
            return offSpring1, offSpring2
        else:
            #cloning 
            offSpring1 = c1
            offSpring2 = c2
            return offSpring1, offSpring2

    def mutation(self, chromosome):
        """Applies the mutation operator to a chromsome."""
        chance = random.random()
        chance = 0.01
        if chance < self.mutationProb:
            geneToChange = random.randrange(0, len(chromosome))
            if(chromosome[geneToChange] == 1):
                chromosome[geneToChange] = 0
            elif(chromosome[geneToChange] == 0):
                chromosome[geneToChange] = 1
        return chromosome

    def new_generation(self):
        """Method sets the child chromosomes list as the parent chromosomes list."""
        for i in range(len(self.parentChromosomes)):
            self.parentChromosomes[i] = self.childChromosomes[i]
    
    def best_chromosome(self):
        """Appends the best chromosome in each generation to the best chromsomes list."""
        chromosomeScores = self.chromosomeFitnessScores
        bestChromosome = chromosomeScores.index(max(chromosomeScores))
        chromsomeBestScore = chromosomeScores[bestChromosome]
        avgScore = sum(chromosomeScores)/len(chromosomeScores)
        self.bestChromosome.append(self.parentChromosomes[bestChromosome])
        self.bestChromosomeFitness.append(bestChromosome)
        self.averageFitness.append(avgScore)
        return avgScore, chromsomeBestScore

    def get_bestChromsomes(self):
        """Get the list of best chromsomes for each generation"""
        return self.bestChromsomes

    def get_final_results(self):
        chromosomes = self.parentChromosomes
        fitness = self.chromosomeFitnessRatios
        bestChromosome = chromosomes[fitness.index(max(fitness))]
        averageFitness = sum(fitness)/len(fitness)
        bestChromosomeFitness = fitness[fitness.index(max(fitness))]
        return round(averageFitness,3), bestChromosomeFitness, bestChromosome

def toString(binary):
    prior = 0
    letters = []
    word = ""
    binarySize = 8
    for i in range(0, len(binary), binarySize):
        letters.append(binary[i:i + binarySize])
    for letter in letters:
        letterInBinary = ""
        for i in range(len(letter)):
            letterInBinary += str(letter[i])
        asInt = int(letterInBinary,2)
        ascii = chr(asInt)
        word+=ascii
    return word

def toBinary(word):
    binaryResult = [] #list to hold the binary result
    finalResult = []
    for letter in word: #for letter in word
        ascii = ord(letter) #get the ascii value of the letter
        binaryResult.append(bin(ascii)[2:].zfill(8)) #append the binary of the ascii value to the binary list
    for bString in binaryResult:
        bString = str(bString)
        for binary in bString:
            finalResult.append(int(binary))
    return finalResult

def main():
    while True:
        generations = 500
        uInput = input("Input a word, I'll try and guess what the word is: ")
        if(uInput.lower() == "quit"):
            break
        binaryInput = toBinary(uInput)
        GA = GeneticAlgo(userInput=binaryInput)
        GA.set_parents()
        os.system("echo I'm thinking about it...")
        for i in range(generations):
            GA.fitness_function()
            GA.fitness_ratios()
            GA.mating()
            GA.new_generation()
        os.system("echo Alright, I think I got it...")
        time.sleep(0.5)
        results = GA.get_final_results()
        word = toString(results[2]).lower()
        os.system(f"echo I'm guessing the word is.. {word}")


    

main()