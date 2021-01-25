#!/usr/bin/env python3
import runcmd
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 12:13:44 2019

@author: Bruno A. C. Horta, IQ-UFRJ, Brazil
"""

import copy
import multiprocessing
from random import uniform
from numpy import min,max

class VBGA:
    def __init__(self, individualClass, grid, maxSamples, fitnessMethod,selectionMethod,crossoverMethod,crossoverRate,mutationMethod, mutationRate, populationSize=100):
        self.populationSize = populationSize
        self.individualClass = individualClass
        self.fitnessMethod = fitnessMethod
        self.selectionMethod = selectionMethod
        self.crossoverMethod = crossoverMethod
        self.crossoverRate = crossoverRate
        self.mutationRate = mutationRate
        self.mutationMethod = mutationMethod
        self.population = []
        # extra
        self.grid = grid
        self.maxSamples = maxSamples
       
    def createRandomIndividual(self):
        ind = self.individualClass(self.grid,self.maxSamples)
        ind.randomize()
        return ind
        
    def createInitialPopulation(self):
        self.population = []
        for i in range(0,self.populationSize):
            ind = self.createRandomIndividual()
            self.evaluateIndividual(ind)
            self.population.append(ind)
    
    def evaluateIndividual(self, individual):
        individual.fitValue = self.fitnessMethod(individual)
        #print(individual.fitValue)
    
    def evaluatePopulation(self):
        for individual in self.population:
            self.evaluateIndividual(individual)
        
    def applySelectionMethod(self):
        selectedPopulation = self.selectionMethod(self.population)
        return selectedPopulation

    def fillPopulation(self, selectedPopulation):
        nextPopulation = copy.deepcopy(selectedPopulation)
        for i in range(len(selectedPopulation),self.populationSize):
            ind = self.createRandomIndividual()
            # yMHG: why now?
            self.evaluateIndividual(ind)
            nextPopulation.append(ind)
        self.population = nextPopulation

    def averageFitness(self):
        max = 100000
        for individual in self.population:
            #print(individual.fitValue)
            if individual.fitValue < max:
                max = individual.fitValue
        return max

    def applyCrossoverMethod(self, selected):
        populationCross = copy.deepcopy(selected)
        for i in range(0,len(selected)):
            for j in range(i+1, len(selected)):
                dice = uniform(0,100)
                if dice < self.crossoverRate:
                    father = selected[i]
                    mother = selected[j]
                    childOne = self.createRandomIndividual()
                    childTwo = self.createRandomIndividual()
                    children = self.crossoverMethod(father, mother, childOne, childTwo)
                    #Apply mutation:
                    children[0] = self.applyMutationMethod(children[0])
                    children[1] = self.applyMutationMethod(children[1])
                    self.evaluateIndividual(children[0])
                    self.evaluateIndividual(children[1])
                    populationCross = populationCross + children
        return populationCross

    def applyMutationMethod(self, individual):
        dice = uniform(0, 100)
        if dice < self.mutationRate:
            individual = self.mutationMethod(individual)
        return individual

    def getBest(self):
        max = 100000
        best = self.createRandomIndividual()
        self.evaluateIndividual(best)
        for individual in self.population:
            # print(individual.fitValue)
            if individual.fitValue < max:
                best = individual
                max = individual.fitValue
        return best

    def run(self, maxIterations=100):
        self.createInitialPopulation()
        savePreviousFitness = 100000000000000000000

        for i in range(0, maxIterations):
            print ("Generation: %d/%d" % (i+1,maxIterations))
            self.evaluatePopulation()
            selected = self.applySelectionMethod()
            populationCross = self.applyCrossoverMethod(selected)
            self.fillPopulation(populationCross)
            avFitness = self.averageFitness()
            if(avFitness<savePreviousFitness):
                print("New best fitness: " + str(avFitness))
                savePreviousFitness = avFitness
