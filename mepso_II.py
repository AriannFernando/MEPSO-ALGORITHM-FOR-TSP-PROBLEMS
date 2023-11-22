# -*- coding: utf-8 -*-


########################################################################################################################
# A Particle Swarm Optimization algorithm for solving the traveling salesman problem.
# 
#  
# Authors: Rafael Batres, Ariann Fernando Arriaga Alcántara, Juan Carlos Espinoza García, José Carlos Bañales López.
# Institution: Tecnologico de Monterrey
# Date: November 17, 2023
########################################################################################################################


from operator import attrgetter
import random, sys, time, copy
import random
import copy
import csv
import math
import statistics
from datetime import datetime
from datetime import timedelta
from time import process_time
import tsplib95

import gc


# class that represents a graph
class Graph:

  def __init__(self, amount_vertices, cost_table):

    self.costTable = cost_table
    self.graphSize = amount_vertices



  # Calculate the objective function

  def evaluateCost(self, route):
    routeSize = len(route)
    
    # Determine the cost of traversing the cities
    first = -1
    second = -1
    last = -1
    cost = 0.0
    #print(routeSize)
    #print(route)
    i = 0
    while routeSize > 1 and i < routeSize-1:
      first = route[i]
      second = route[i+1]
      cost = cost + self.costTable[first][second]
      i = i + 1
    # Complete Hamiltonian circuit
    last = route[routeSize-1]
    first = route[0]
    depot = 0
    cost = cost + self.costTable[depot][first] + self.costTable[last][depot]
    return cost



  # gets random unique paths - returns a list of lists of paths
 
  def getRandomPaths(self, max_size):
    random_paths = []
    
    for i in range(max_size):
      
      list_temp = self.getRandomSolution(self.graphSize)
      #print("Initial random sol: ", list_temp)

      if list_temp not in random_paths:
        random_paths.append(list_temp)
    return random_paths
    
  # Generate a random sequence and stores it
  # as a Route
  def getRandomSolution(self, graph_size):
    route = Chromosome()
    visited = [0 for j in range(graph_size)]
    city = -1
    cityCount = 0
    while cityCount < len(visited)-1:
      city = random.randint(1, graph_size-1)
      while visited[city] == 1:
        city = random.randint(1, graph_size-1)        
      route.append(city)
      
      visited[city] = 1
      cityCount = cityCount + 1
    #print("route", route, len(route))
    return route

  def getRandomSolution__(self, graph_size):
    route = Chromosome()
    visited = [0 for j in range(graph_size)]
    city = -1
    cityCount = 0
    while cityCount < len(visited):
      city = random.randint(0, graph_size-1)
      while visited[city] == 1:
        city = random.randint(0, graph_size-1)
      route.append(city)
      visited[city] = 1
      cityCount = cityCount + 1
    return route


# class that represents a particle
class Particle:
  __slots__ = ['solution', 'pbest',  'newSolutionCost', 'pbestCost', 'velocity', 'history']
  def __init__(self, solution, cost):

    # current solution
    self.solution = solution

    # best solution it has achieved so far by this particle
    self.pbest = solution

    # set costs
    self.newSolutionCost = cost
    self.pbestCost = cost

    # velocity of a particle is a sequence of 3-tuple
    # (1, 2, 0.8) means SO(1,2), beta = 0.8
    self.velocity = []

    # past positions
    self.history = []

  # set pbest
  def setPBest(self, new_pbest):
    self.pbest = new_pbest

  # returns the pbest
  def getPBest(self):
    return self.pbest

  # set the new velocity (sequence of swap operators)
  def setVelocity(self, new_velocity):
    self.velocity = new_velocity

  # returns the velocity (sequence of swap operators)
  def getVelocity(self):
    return self.velocity

  # set solution
  def setCurrentSolution(self, solution):
    self.solution = solution

  # gets solution
  def getCurrentSolution(self):
    return self.solution

  # set cost pbest solution
  def setCostPBest(self, cost):
    self.pbestCost = cost

  # gets cost pbest solution
  def getCostPBest(self):
    return self.pbestCost

  # set cost current solution
  def setCostCurrentSolution(self, cost):
    self.newSolutionCost = cost

  # gets cost current solution
  def getCurrentSolutionCost(self):
    return self.newSolutionCost

  # removes all elements of the list velocity
  def clearVelocity(self):
    del self.velocity[:]


# PSO algorithm
class Solver:


  def __init__(self, epoch_stop, graph, iterations, maxEpochs, size_population, beta=1, alfa=1):
    self.epoch_stop = epoch_stop
    self.graph = graph # the graph
    self.iterations = iterations # max of iterations
    self.maxEpochs = maxEpochs
    self.populationSize = size_population # size population
    self.particles = Chromosome() # list of particles
    self.beta = beta # the probability that all swap operators in swap sequence (gbest - x(t-1))
    self.alfa = alfa # the probability that all swap operators in swap sequence (pbest - x(t-1))
    self.last_epoch = 0
    self.costConvergence = None
    
    """
    #graph_size = 5
    
    # initialized with a group of random particles (solutions)
    solutions = self.graph.getRandomPaths(self.populationSize)
    print("One initial solution: ", solutions[0])
    
    # Generates a solution with the Savings method
    
    #amountOfGoodSolutions = int(round(size_population*0.01))
    amountOfGoodSolutions = int(round(size_population*0.0))
    for i in range(amountOfGoodSolutions):
      solutions[i] = self.getSavingsSolution()
    
    #print("Solutin 0:", solutions[0])
    #print("Good solution:", goodSolution)
    #solutions[0] = goodSolution

    # checks if exists any solution
    if not solutions:
      print('Initial population empty! Try run the algorithm again...')
      sys.exit(1)
    
    # Select the best random population among 5 populations
    bestSolutions = list(solutions)
        
    mostDiverse = self.evaluateSolutionsDiversity(solutions)
    for i in range(5):
      solutions = self.graph.getRandomPaths(self.populationSize)
      sim = self.evaluateSolutionsDiversity(solutions)
      print("Diversity of the population: ", sim)
      #cost = self.evaluateSolutionsAverageCost(solutions)
      if sim > mostDiverse:
        mostDiverse = sim
        bestSolutions = list(solutions)
      del solutions[:]

    
    ###bestSolutions[0] = goodSolution
    # creates the particles and initialization of swap sequences in all the particles
    for solution in bestSolutions:
      # creates a new particle
      particle = Particle(solution=solution, cost=graph.evaluateCost(solution))
      # add the particle
      self.particles.append(particle)
    
    # updates "populationSize"
    self.populationSize = len(self.particles)
    """

  def initPopulationFirstEpoch(self, population_size):
    self.particles = Chromosome() # list of particles
    solutions = self.graph.getRandomPaths(self.populationSize)
    print("One initial solution: ", solutions[0])
    
    # checks if exists any solution
    if not solutions:
      print('Initial population empty! Try run the algorithm again...')
      sys.exit(1)
    
    # Select the best random population among 5 populations
    bestSolutions = list(solutions)
        
    mostDiverse = self.evaluateSolutionsDiversity(solutions)
    for i in range(5):
      solutions = self.graph.getRandomPaths(self.populationSize)
      sim = self.evaluateSolutionsDiversity(solutions)
      print("Diversity of the population: ", sim)
      #cost = self.evaluateSolutionsAverageCost(solutions)
      if sim > mostDiverse:
        mostDiverse = sim
        bestSolutions = list(solutions)
      del solutions[:]

    
    ###bestSolutions[0] = goodSolution
    # creates the particles and initialization of swap sequences in all the particles
    for solution in bestSolutions:
      # creates a new particle
      particle = Particle(solution=solution, cost=graph.evaluateCost(solution))
      # add the particle
      self.particles.append(particle)
    
    # updates "populationSize"
    self.populationSize = len(self.particles)

  def initPopulation(self, population_size):
    self.particles = Chromosome() # list of particles
    solutions = self.graph.getRandomPaths(population_size)
    self.populationSize = population_size

    
    # Generates a solution with the Savings method
    
    #amountOfGoodSolutions = int(round(size_population*0.01))
    #amountOfGoodSolutions = int(round(population_size*0.0))
    amountOfGoodSolutions = 0
    for i in range(amountOfGoodSolutions):
      solutions[i] = self.mutateGoodSolution(self.getSavingsSolution())
    
    #print("Solutin 0:", solutions[0])
    #print("Good solution:", goodSolution)
    #solutions[0] = goodSolution

    # checks if exists any solution
    if not solutions:
      print('Initial population empty! Try run the algorithm again...')
      sys.exit(1)

    # Select the best random population among 5 populations
    bestSolutions = list(solutions)
    bestCost = self.evaluateSolutionsAverageCost(solutions)
    
    
    for i in range(5):
      solutions = self.graph.getRandomPaths(self.populationSize)
      cost = self.evaluateSolutionsAverageCost(solutions)
      if cost < bestCost:
        bestCost = cost
        bestSolutions = list(solutions)
      del solutions[:]
    
    
    # creates the particles and initialization of swap sequences in all the particles
    for solution in bestSolutions:
      # creates a new particle
      particle = Particle(solution=solution, cost=graph.evaluateCost(solution))
      # add the particle
      self.particles.append(particle)
      
  def hamming(self, route1, route2): 
    distance = 0
    # Loop over the indices of the string
    for i in range(len(route1)):
        # Add 1 to the distance if these two numbers are not equal
        if route1[i] != route2[i]:
            distance += 1
    # Return the final count of differences
    return distance
    
  def hamming_(self, route1, route2): 
    distance = 0
    # Loop over the indices of the string
    for i in range(len(route1)):
        distance += (route1[i] - route2[i])**2
    distance = math.sqrt(distance)
    # Return the final count of differences
    return distance
  
  def evaluateSolutionsDiversity(self, solutions):
    simSum = 0
    count = 0
    for solution1 in solutions:
      for solution2 in solutions:
        if not (solution1 == solution2):
          count += 1
          #sim = self.cosine_similarity(solution1, solution2)
          sim = self.hamming(solution1, solution2)
          simSum += sim
    return simSum / count
    
  def evaluateSolutionsAverageCost(self, solutions):
  
    totalCost = 0.0
    i = 0
    for solution in solutions:
      cost = self.evaluateCost(solution)
      totalCost += cost
      i+=1
    averageCost = totalCost / float(i)

    return averageCost


  # set gbest (best particle of the population)
  def setGBest(self, new_gbest):
    self.gbest = new_gbest

  # returns gbest (best particle of the population)
  def getGBest(self):
    return self.gbest  
  
  # returns the elapsed milliseconds since the start of the program
  def elapsedTime(self, start_time):
    dt = datetime.now() - start_time
    ms = (dt.days * 24 * 60 * 60 + dt.seconds) * 1000 + dt.microseconds / 1000.0
    return ms
  
  def setEpoch(self, last_epoch):
    self.last_epoch = last_epoch
  
  def getEpoch(self):
    return self.last_epoch

  
  def mutation_probability(self, initial_probability: float=0.1, current_epoch: int=0, max_epochs: int=1000) -> float:
    alpha = max_epochs / math.log(initial_probability/9E-01)
    return initial_probability * math.exp(-current_epoch / alpha)

  def run(self):
    # variables for convergence data
    convergenceData = []
    iterationArray = []
    bestCostArray = []
    epochArray = []
    epochBestCostArray = []
    bestCostSampling = []
    costCon = []

    batchSize = 100 # save data every n iterations
    batchCounter = 0

    HISTORY_SIZE = 100
    MAX_NEIGHBORS = 1

    startTime = datetime.now()
    
    # updates gbest (best particle of the population)
    #self.gbest = min(self.particles, key=attrgetter('pbestCost'))

    
    
    epoch = 0
    delta = 0.0
    new_beta = self.beta
    while epoch < self.maxEpochs:
    #for epoch in range(self.maxEpochs):
      population_size = self.populationSize
      if epoch == 0:
        self.initPopulationFirstEpoch(population_size)
        self.gbest = self.particles[0]
        eliteSolution = []
      else:
        self.initPopulation(population_size)
        
        
      #population_size = len(self.particles)
      print("Epoch: ", epoch, "with ", population_size, " particles")
      print("Alfa = ", self.alfa, "Beta = ", self.beta)
      convergencePerEpoch = []


      eliteCost = sys.maxsize

      ##solution = self.graph.getRandomPaths(2)
      ##particle = Particle(solution=copy.deepcopy(solution[0]), cost=graph.evaluateCost(solution[0]))
      ##print("Particle: ", particle.getCurrentSolution())
      ##self.particles.append(particle)
      
      ##population_size += 1


        
      print("Particles: ", len(self.particles))
      if epoch > 0:
        # Insert the best individual into the new population (1% of the population)
        #for i in range(population_size):
        #print("Elite:", eliteSolution)

        if random.uniform(0,1) < 1.0:
          mutated_elite = self.mutateGoodSolution(eliteSolution)
          mutated_elite_cost = self.graph.evaluateCost(mutated_elite)
          self.particles[0]  = Particle(mutated_elite, mutated_elite_cost)
          #self.particles[0]  = Particle(mutated_elite, eliteCost)
          print("Inserted elite solution!")
          
        
    
      # for each time step (iteration)
      if epoch == 0:
          maxIter = 3 * self.iterations
      else:
          maxIter = self.iterations
      for t in range(maxIter):
        convergencePerIteration = []
        batchCounter = batchCounter + 1
        
        # updates gbest (best particle of the population)
        #self.gbest = min(self.particles, key=attrgetter('pbestCost'))
        #costVariance = (self.particles, key=attrgetter('pbestCost'))
        averageCost = statistics.mean(particle.pbestCost for particle in self.particles)
        costStd = statistics.pstdev(particle.pbestCost for particle in self.particles)

        # for each particle in the swarm
        for particle in self.particles:
          
          previousCost = particle.getCurrentSolutionCost()
          
          #particle.clearVelocity() # cleans the speed of the particle
          #velocity = []
          #velocity_alfa = []
          # what is the difference btn. copy.copy and list?
          gbest = list(self.gbest.getPBest()) # gets solution of the gbest solution
          #pbest = particle.getPBest()[:] # copy of the pbest solution
          ##newSolution = particle.getCurrentSolution()[:] # gets copy of the current solution of the particle
          
          # the previous solution
          previousSolution = particle.getCurrentSolution()[:]
          if particle != self.gbest and previousSolution == gbest:
            ####new_particle_sol = self.graph.getRandomSolution(len(previousSolution)+1)
            ####new_particle_sol_cost = self.graph.evaluateCost(new_particle_sol)
            ####particle.setCurrentSolution(new_particle_sol)
            ####particle.setCostCurrentSolution(new_particle_sol_cost)
            mutated_particle_solution = self.mutateGoodSolution(previousSolution)
            mutated_particle_cost = self.graph.evaluateCost(mutated_particle_solution)
            particle.setCurrentSolution(mutated_particle_solution)
            particle.setCostCurrentSolution(mutated_particle_cost)
          
          if len(particle.history) == HISTORY_SIZE:
            particle.history.pop(0)
          
          bestNeighbor = self.mutateGoodSolution(particle.getCurrentSolution())
          #bestNeighbor = particle.getCurrentSolution()[:]
          bestNeighborCost = self.graph.evaluateCost(bestNeighbor)
          
          """
          if previousCost < bestNeighborCost:
            bestNeighbor = previousSolution[:]
            bestNeighborCost = previousCost
          """

            
          newSolution = particle.getCurrentSolution()[:]
          newSolutionCost = particle.getCurrentSolutionCost()

          #random_particle = random.choice(self.particles)
          
          """
          if random.random() <= (1.0 - self.alfa - self.beta):
            newSolution = self.crossover(list(newSolution), previousSolution)
          elif random.random() <= self.beta:
            newSolution = self.crossover(list(newSolution), self.gbest.getPBest())
          elif random.random() <= self.alfa:
            newSolution = self.crossover(list(gbest), random_particle.getPBest())
          """  
          if random.random() <= new_beta:
            new_son_solution, new_daughter_solution = self.crossover(list(newSolution), self.gbest.getPBest())
            new_son_solution_cost = self.graph.evaluateCost(new_son_solution)
            new_daughter_solution_cost = self.graph.evaluateCost(new_daughter_solution)
            if new_son_solution_cost < new_daughter_solution_cost:
                newSolution = new_son_solution
                # gets cost of the current solution
                newSolutionCost = new_son_solution_cost
            else:
                newSolution = new_daughter_solution
                # gets cost of the current solution
                newSolutionCost = new_daughter_solution_cost
          elif random.random() <= delta:
            new_son_solution, new_daughter_solution = self.crossover(self.graph.getRandomSolution(len(previousSolution)+1), self.gbest.getPBest())
            new_son_solution_cost = self.graph.evaluateCost(new_son_solution)
            new_daughter_solution_cost = self.graph.evaluateCost(new_daughter_solution)
            if new_son_solution_cost < new_daughter_solution_cost:
                newSolution = new_son_solution
                # gets cost of the current solution
                newSolutionCost = new_son_solution_cost
            else:
                newSolution = new_daughter_solution
                # gets cost of the current solution
                newSolutionCost = new_daughter_solution_cost          
          elif random.random() <= self.alfa:
            largest_dist = 0
            #less_common = len(gbest)
            for neighbor_particle in self.particles:
              sol = neighbor_particle.getPBest()
              dist = self.hamming(gbest, sol)
              if dist > largest_dist:
                largest_dist = dist
                #less_common = common
                dissimilar_particle = neighbor_particle

            new_son_solution, new_daughter_solution = self.crossover(list(newSolution), dissimilar_particle.getPBest())

            #newSolution = self.crossover(list(gbest), random_particle.getPBest())
            new_son_solution_cost = self.graph.evaluateCost(new_son_solution)
            new_daughter_solution_cost = self.graph.evaluateCost(new_daughter_solution)
            if new_son_solution_cost < new_daughter_solution_cost:
                newSolution = new_son_solution
                # gets cost of the current solution
                newSolutionCost = new_son_solution_cost
            else:
                newSolution = new_daughter_solution
                # gets cost of the current solution
                newSolutionCost = new_daughter_solution_cost

            # gets cost of the current solution
          #newSolutionCost = self.graph.evaluateCost(newSolution)

          if newSolutionCost < bestNeighborCost:
            bestNeighbor = newSolution[:]
            bestNeighborCost = newSolutionCost

          if bestNeighborCost < previousCost and bestNeighbor not in particle.history:
              # updates the current solution  
            particle.setCurrentSolution(bestNeighbor)
              # updates the cost of the current solution
            particle.setCostCurrentSolution(bestNeighborCost)
            particle.history.append(bestNeighbor)


          #"""

          # checks if new solution is pbest solution
          pbCost =  particle.getCostPBest()
          
          """
          if particle.getCurrentSolutionCost() < pbCost:
            particle.setPBest(particle.getCurrentSolution())
            particle.setCostPBest(particle.getCurrentSolutionCost())
          """  

          if bestNeighborCost < pbCost:
              acceptance_prob = 1
          else:
              # If the new solution is worse, calculate an acceptance probability
              acceptance_prob = self.mutation_probability(0.05, epoch, self.maxEpochs)
          
          
          if random.random() < acceptance_prob:
            particle.setPBest(bestNeighbor)
            particle.setCostPBest(bestNeighborCost)

            #temperature = temperature + (pbCost - newSolutionCost) / (math.log(rnd))
          
          gbestCost = self.gbest.getCostPBest()
                   
          
            
            # Using the value of the Boltzmann constant
          
          # check if new solution is gbest solution
         
          if particle.getCurrentSolutionCost() < gbestCost:
            self.gbest = copy.deepcopy(particle)

        eliteSolution = self.getGBest().getPBest()[:]
        eliteCost = self.gbest.getCostPBest()
        if batchCounter > batchSize:
          #print("Sum of acceptance probabilities:", sumAcceptanceProbabilities)
          print(t, "Gbest cost = ", self.gbest.getCostPBest())
          convergencePerIteration.append(t)
          convergencePerIteration.append(self.gbest.getCostPBest())
          convergencePerIteration.append(averageCost)
          convergencePerIteration.append(costStd)
          convergenceData.append(convergencePerIteration)
          iterationArray.append(t)
          bestCostArray.append(self.gbest.getCostPBest())
          batchCounter = 0

        if self.maxEpochs > 1:
          convergencePerEpoch.append(epoch)
          convergencePerEpoch.append(self.gbest.getCostPBest())
          convergenceData.append(convergencePerEpoch)
          epochArray.append(epoch)
          epochBestCostArray.append(self.gbest.getCostPBest())
    
      epoch = epoch + 1
      self.setEpoch(epoch)
      bestCostSampling.append(self.gbest.getCostPBest())
      if epoch > (200 - self.epoch_stop) :
        costCon.append(self.gbest.getCostPBest())
          
      if epoch >= self.epoch_stop :
        std = statistics.pstdev(bestCostSampling[-self.epoch_stop:])
        print("standard deviation: ", std)
      else:
        std = 1000
      
      if std == 0:
        break
      elif std < 100 and std !=0:
        delta = 0.0
        new_beta = self.beta
      elif std > 200:
        delta = self.beta / 2.0
        new_beta = self.beta / 20.0      
    
    print("What's going on?")
    self.costConvergence=costCon


    """
    csvFile = open('pso-convergence.csv', 'w', newline='')  
    with csvFile: 
         writer = csv.writer(csvFile)
         writer.writerows(convergenceData)
    csvFile.close()
    print("Elapsed time: ", self.elapsedTime(startTime))
    """ 
    
  def acceptanceProbability(self, previous_solution, new_solution, temperature, boltzmann):
    # If the new solution is better, accept it
    if new_solution < previous_solution:
      return 1.0
    elif temperature == 0:
      return 0.0
    # if there is no change then accept the bad solution
    #elif new_solution == previous_solution:
    #  return 0.0
    # If the new solution is worse, calculate an acceptance probability
    else:
      try:
        boltzmannProb = math.exp((previous_solution - new_solution) / (temperature * boltzmann))
      except OverflowError:
        boltzmannProb = float('inf')
      
      #print("Boltzmann probability", boltzmannProb)
      #print("Worse solution accepted ...")
      return boltzmannProb

  # Reverse mutation
  def mutate(self, parent, point1, point2):
    chromosome = parent[:]

    # Inverse the genes between point1 and point2    
    while True:
      chromosome[point1],chromosome[point2] = chromosome[point2],chromosome[point1]
      point1 = point1 + 1
      point2 = point2 - 1
      if point1 > point2:
        break
      
    return chromosome

# Use reverse mutation for elite and savings solutions
  def mutateGoodSolution(self, elite_solution):
    chromosome = elite_solution[:]

    point1 = -1
    point2 = -1
    while True:
      point1 = random.randint(0, len(chromosome)-1)
      point2 = random.randint(0, len(chromosome)-1)
      if point1 < point2:
        break
    # Swap the contents of the two points
    #chromosome[point1],chromosome[point2] = chromosome[point2],chromosome[point1]
    
    if random.random() < 0.16:
      temp = elite_solution[point2]
      for i in range(point1 + 1, point2):
        chromosome[i] = elite_solution[i - 1]  
      i = point2 + 1
      while i < len(elite_solution):
        chromosome[i] = elite_solution[i]
        i = i + 1
      chromosome[point2] =  elite_solution[point2 - 1]
      chromosome[point1] = temp
    #print("point1: ", point1)
    #print("point2: ", point2)
    #print("Original: ", elite_solution)
    #print("Inserted: ", chromosome)
    elif random.random() < 0.16:
      chromosome[point1],chromosome[point2] = chromosome[point2],chromosome[point1]
    # Inverse the genes between point1 and point2
    else:
      while True:
        chromosome[point1],chromosome[point2] = chromosome[point2],chromosome[point1]
        point1 = point1 + 1
        point2 = point2 - 1
        if point1 > point2:
          break
    return chromosome

# Crossover operator with mutation
  # This is an ordered crossover in which the center part of the dad chromosome 
  # is passed to the son and the left and right parts come from the mom chromosome.
  # The right part is filled first using the right part of the mom and then
  # the left part of the mom, while omitting the genes already inserted.
  # Similarly, the left part is filled next, using the left part of the mom 
  # and then the center part of the mom, while omitting the genes already inserted.
  # Subsequently, the offspring is reverse-mutated.
  def crossover(self, dadChromosome, momChromosome):
    #print("")
    
    while True:
      point1 = random.randint(0, len(dadChromosome)-1)
      point2 = random.randint(0, len(dadChromosome)-1)
      if point1 < point2 and point1 !=0 and point2 != len(dadChromosome)-1:
        break
    #print("Point 1: ", point1)
    #print("Point 2: ", point2)
    #print("Dad's chromosome: ", dadChromosome)
    #print("Mom's chromosome: ", momChromosome)
    
    sonChromosome = [-1 for j in range(len(dadChromosome))]

    son_center_segment = list()
    son_left_segment = list()
    son_right_segment = list()
  
    daughterChromosome = [-1 for j in range(len(momChromosome))]

    daughter_center_segment = list()
    daughter_left_segment = list()
    daughter_right_segment = list()
    
    # Copy the middle section
    for i in range(point1, point2+1):
      son_center_segment.append(dadChromosome[i])

    for i in range(point1, point2+1):
      daughter_center_segment.append(momChromosome[i])

    # Fill in the right section of the son
    i = point2 + 1
    j = point2 + 1
    gene = -1
    while i < len(sonChromosome):
      if j == len(momChromosome):
        j = 0
      while j < len(momChromosome):
        gene = momChromosome[j]
        j=j+1
        if gene not in son_center_segment:
          son_right_segment.append(gene)
          i = i + 1
          break
    
    # Fill in the right section of the daughter
    i = point2 + 1
    j = point2 + 1
    gene = -1
    while i < len(daughterChromosome):
      if j == len(dadChromosome):
        j = 0
      while j < len(dadChromosome):
        gene = dadChromosome[j]
        j=j+1
        if gene not in daughter_center_segment:
          daughter_right_segment.append(gene)
          i = i + 1
          break
      
    # Fill in the left section of the son
    i=0
    j=0
    gene = -1
    while i < point1 and j < len(momChromosome):
      while j < len(momChromosome):
        gene = momChromosome[j]
        j=j+1
        if gene not in son_center_segment and gene not in son_right_segment:
          son_left_segment.append(gene)
          break
      i = i + 1
    
    # Fill in the left section of the daughter
    i=0
    j=0
    gene = -1
    while i < point1 and j < len(dadChromosome):
      while j < len(dadChromosome):
        gene = dadChromosome[j]
        j=j+1
        if gene not in daughter_center_segment and gene not in daughter_right_segment:
          daughter_left_segment.append(gene)
          break
      i = i + 1

    sonChromosome = son_left_segment + son_center_segment + son_right_segment
    daughterChromosome = daughter_left_segment + daughter_center_segment + daughter_right_segment
    # Swap the contents of the two points
    ##sonChromosome[point1],sonChromosome[point2] = sonChromosome[point2],sonChromosome[point1]
    #print("Son: ", sonChromosome)

    # After the crossover do a reverse mutation
    # Select two random points in the chromosome
    
    point1 = -1
    point2 = -1
    while True:
      point1 = random.randint(0, len(sonChromosome)-1)
      point2 = random.randint(0, len(sonChromosome)-1)
      if point1 != point2:
        break
    # Inverse the genes between point1 and point2    
    while True:
      sonChromosome[point1],sonChromosome[point2] = sonChromosome[point2],sonChromosome[point1]
      point1 = point1 + 1
      point2 = point2 - 1
      if point1 > point2:
        break
    
    point1 = -1
    point2 = -1
    while True:
      point1 = random.randint(0, len(daughterChromosome)-1)
      point2 = random.randint(0, len(daughterChromosome)-1)
      if point1 != point2:
        break
    # Inverse the genes between point1 and point2    
    while True:
      daughterChromosome[point1],daughterChromosome[point2] = daughterChromosome[point2],daughterChromosome[point1]
      point1 = point1 + 1
      point2 = point2 - 1
      if point1 > point2:
        break

    #print("Son: ", sonChromosome)
    #print("*******")
    return sonChromosome, daughterChromosome



  # The ordered crossover operator
  def oxcrossover(self, dadRoute, momRoute):
    #print("")
    
    visitedDad = [False for i in range(len(dadRoute))]
    visitedMom = [False for i in range(len(momRoute))]
    
  
    while True:
      point1 = random.randint(0, len(dadRoute)-1)
      point2 = random.randint(0, len(dadRoute)-1)
      if point1 < point2 and point1 !=0 and point2 != len(dadRoute)-1:
        break
        
    #print("Point 1: ", point1)
    #print("Point 2: ", point2)
    #print("Dad's route: ", dadRoute)
    #print("Mom's route: ", momRoute)
    
    sonRoute = Route()
    
    sonRoute = [-1 for j in range(len(dadRoute))]
    
    # Copy the middle section
    for i in range(point1, point2+1):
      sonRoute[i] = dadRoute[i]
      visitedDad[dadRoute[i]] = True
      visitedMom[momRoute[i]] = True

    # Fill in the right section of the son
    i = point2 + 1
    j = point2 + 1
    gene = -1
    while i < len(sonRoute):
      if j == len(momRoute):
        j = 0
      while j < len(momRoute):
        gene = momRoute[j]
        j=j+1
        if not visitedDad[gene]:
          sonRoute[i] = gene
          visitedDad[gene] = True
          i = i + 1
          break
      
      
    # Fill in the left section of the son
    i=0
    j=0
    gene = -1
    while i < point1 and j < len(momRoute):
      while j < len(momRoute):
        gene = momRoute[j]
        j=j+1
        if not visitedDad[gene]:
          sonRoute[i] = gene
          visitedDad[gene] = True
          break
      i = i + 1

      


    #"""  
    
    # After the crossover do a reverse mutation
    # Select two random points in the route of the parent
    point1 = -1
    point2 = -1
    while True:
      point1 = random.randint(0, len(sonRoute)-1)
      point2 = random.randint(0, len(sonRoute)-1)
      if point1 != point2:
        break
    # Inverse the genes between point1 and point2    
    while True:
      sonRoute[point1],sonRoute[point2] = sonRoute[point2],sonRoute[point1]
      point1 = point1 + 1
      point2 = point2 - 1
      if point1 > point2:
        break
    #"""    
    return sonRoute

  # Heuristic crossover
  def hcrossover(self, dadRoute, momRoute):

    #print("Mom:", momRoute)
    #print("Dad:", dadRoute)
    #print("Cost of mom:", momCost)
    #print("Cost of dad:", dadCost)
    
    sonRoute = [-1 for j in range(len(dadRoute))]
    visited = [False for k in range(len(dadRoute))]
    
    firstGene = sonRoute[0] = dadRoute[0]
    visited[firstGene] = True
    
    i = 1
    while i < len(dadRoute)+2:  
      fromGene = sonRoute[i-1]
      #print("From Gene:", fromGene)
      connectedGeneInMom = self.getConnectedGene(fromGene, momRoute)
      connectedGeneInDad = self.getConnectedGene(fromGene, dadRoute)
      costMomEdge = GRAPH[fromGene][connectedGeneInMom]
      costDadEdge = GRAPH[fromGene][connectedGeneInDad]
      
      sampleSize = 3
      if costDadEdge < costMomEdge:
        best = connectedGeneInDad
        worse = connectedGeneInMom
      else:
        best = connectedGeneInMom
        worse = connectedGeneInDad
      
      if visited[best]:
        if visited[worse]:
          winner = -1
          winnerCost = -1
          for l in range(0, sampleSize-1):
            candidate = self.getUnvisitedGene(dadRoute, visited)
            candidateCost = GRAPH[fromGene][candidate]
            if winner == -1 or candidateCost > winnerCost:
              winner = candidate
              winnerCost = candidateCost
        else:
          winner = worse
      else:
        winner = best
      #print(i, "winner", winner)
      sonRoute[i] = winner
      visited[winner] = True
      i += 1
      if i == len(dadRoute):
        break
      
    #sonCost = self.evaluateCost(sonRoute)
    #print("Cost of heuristic son:", sonCost)
      
    """
    missingGenes = False
    for j in range(0,  len(dadRoute)-1):
      if j not in sonRoute:
        print(j, " is missing")
        missingGenes = True
    if missingGenes:
      sys.exit("Gene is missing in son:")
    """
    
    #print(" ")
    
    #if sonCost < dadCost:
      #print("Heuristic son", sonRoute)
      #print("Son is better:", sonCost)
    return sonRoute
        

  # The alternative edge crossover
  def aexcrossover(self, dadRoute, momRoute):
    
    sonRoute = [-1 for j in range(len(dadRoute))]
    
    visited = [False for i in range(len(dadRoute))]
    
    fromGeneInDad = dadRoute[0]
    sonRoute[0] = fromGeneInDad
    connectedGeneInDad = self.getConnectedGene(fromGeneInDad, dadRoute)
    sonRoute[1] = connectedGeneInDad
    visited[fromGeneInDad] = True
    visited[connectedGeneInDad] = True
    i = 2
    
    while i < len(dadRoute):
      fromGeneInMom = connectedGeneInDad
      connectedGeneInMom = self.getConnectedGene(fromGeneInMom, momRoute)
      if visited[connectedGeneInMom]:
        #print("Gene in Mom already visited:", connectedGeneInMom)
        connectedGeneInMom = self.getUnvisitedGene(momRoute, visited)
        #print("Randomly selected momÂ´s gene:", connectedGeneInMom)
      sonRoute[i] = connectedGeneInMom
      #print("Connected gene in Mom:", connectedGeneInMom)
      visited[connectedGeneInMom] = True
      i += 1
      fromGeneInDad = connectedGeneInMom
      connectedGeneInDad = self.getConnectedGene(fromGeneInDad, dadRoute)
      if visited[connectedGeneInDad]:
        #print("Gene in Dad already visited:", connectedGeneInDad)
        connectedGeneInDad = self.getUnvisitedGene(dadRoute, visited)
        #print("Randomly selected dad's gene:", connectedGeneInDad)
      sonRoute[i] = connectedGeneInDad
      #print("Connected gene in Dad:", connectedGeneInDad)
      visited[connectedGeneInDad] = True
      i += 1

    #print("Son:", sonRoute)
    
    missingGenes = False
    for j in range(0,  len(dadRoute)-1):
      if j not in sonRoute:
        #print(j, " is missing")
        missingGenes = True
    if missingGenes:
      sys.exit("Gene is missing in son:")
      
    # Select two random points in the route of the son
    point1 = -1
    point2 = -1
    while True:
      point1 = random.randint(0, len(sonRoute)-1)
      point2 = random.randint(0, len(sonRoute)-1)
      if point1 != point2:
        break
        
    # Swap the contents of the two points
    #sonRoute[point1],sonRoute[point2] = sonRoute[point2],sonRoute[point1]
    return sonRoute

  def getConnectedGene(self, from_gene, route):
    i = 0
    to_gene = -1
    while i < len(route):
      if route[i] == from_gene and i == len(route)-1:
        to_gene = route[0]
        break
      if route[i] == from_gene:
        to_gene = route[i+1]
        break
      i += 1
    if to_gene == -1:
      print("Error!!!!!!")
      print("Route:", route)
      print("From gene:", from_gene)
      print("Connected gene:", to_gene)
      sys.exit("Error message")
    return to_gene
      
  def getUnvisitedGene(self, route, visited):
    while True:
      point = random.randint(0, len(route)-1)
      gene = route[point]
      if not visited[gene]:
        break
    #print("Visiting status of gene:", visited[gene])
    return gene

  # Generate a good solution based on the Clarke-Wright savings method
  # Then scrambles it a bit by swapping two cities
  def getSavingsSolution(self):
    best = Route()
    best.append(0)
    savingsMatrix = self.generateSavingsMatrix()
    first = savingsMatrix[0]
    best.append(first.get_from_node())
    best.append(first.get_to_node())
    last = best[-1]

    for i in range(GRAPH_SIZE - 3):
      toNode = self.getNodeFromSavings(last, savingsMatrix, best)
      best.append(toNode)
      last = best[-1]
    """  
    # Select two random points in the savings result
    point1 = -1
    point2 = -1
    while True:
      point1 = random.randint(0, len(best)-1)
      point2 = random.randint(0, len(best)-1)
      if point1 != point2:
        break
    # Swap the contents of the two points
    best[point1],best[point2] = best[point2],best[point1]
    """
    
    # Do a reverse mutation
    # Select two random points in the route of the parent
    """
    point1 = -1
    point2 = -1
    while True:
      point1 = random.randint(0, len(best)-1)
      point2 = random.randint(0, len(best)-1)
      if point1 != point2:
        break
    # Inverse the genes between point1 and point2    
    while True:
      best[point1],best[point2] = best[point2],best[point1]
      point1 = point1 + 1
      point2 = point2 - 1
      if point1 > point2:
        break
    """
    return best

  # Calculate the objective function
  def evaluateCost(self, route):
    routeSize = len(route)
    
    # Determine the cost of traversing the cities
    first = -1
    second = -1
    last = -1
    cost = 0.0
    #print(routeSize)
    #print(route)
    i = 0
    while routeSize > 1 and i < routeSize-1:
      first = route[i]
      second = route[i+1]
      cost = cost + GRAPH[first][second]
      i = i + 1
    # Complete Hamiltonian circuit
    last = route[routeSize-1]
    first = route[0]
    cost = cost + GRAPH[last][first]
    return cost  



  def generateSavingsMatrix(self):
    savingsMatrix = []
    for i in range(1, GRAPH_SIZE):
      for j in range(1, GRAPH_SIZE):
        savings = GRAPH[i][0] + GRAPH[0][j] - GRAPH[i][j]
        # saving = CWSA_mtx[i,-1] + CWSA_mtx[j,-1] - CWSA_mtx[i,j]
        savingsMatrix.append(Savings(i, j, savings))
    savingsMatrix = sorted(savingsMatrix, reverse=True)
    return savingsMatrix
    # sorted(student_objects, key=lambda student: student.age)

  def getNodeFromSavings(self, last, savings_matrix  = [], best = []):
    for i in range(len(savings_matrix)):
      node = savings_matrix[i]
      fromNode = node.get_from_node()
      toNode = node.get_to_node()
      #print("Last: ", last)
      #print("fromNode:", fromNode)
      #print("toNode:", toNode)
      if fromNode == last and toNode not in best:
        return toNode

      
  

# Define Route as a subclass of list
class Route(list):
  def __init__(self):
    self.elements = []
  #def size(self):
  #  return len(self.elements)

# An Individual stores its route along with
# its cost and fitness.
class Individual:
  def __init__(self, route):
    self.route = route
    self.cost = 0.0
    self.fitness = 0.0
  
  def getRoute(self):
    return self.route

  def setRoute(self, route):
    self.route = route
    
  def getFitness(self):
    return self.fitness

  def setFitness(self, fitness):
    self.fitness = fitness
    
  def getCost(self):
    return self.cost
    
  def setCost(self, cost):
    self.cost = cost
    
# Test code begins

class Savings(object):
  def __init__(self, from_node, to_node, savings_value):
    self.from_node = from_node
    self.to_node = to_node
    self.savings_value = savings_value
  
  def get_from_node(self):
    return self.from_node

  def get_to_node(self):
    return self.to_node
        
  def __repr__(self):
    return repr((self.from_node, self.to_node, self.savings_value))
    
  def __lt__(self, other):
    return self.savings_value < other.savings_value
  
  def __eq__(self, other):
    return self.savings_value == other.savings_value


class City:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    def set_x(self, x):
        self.x = x
    
    def get_x(self):
        return self.x
    
    def set_y(self, y):
        self.y = y
    
    def get_y(self):
        return self.y

# Define Chromosome as a subclass of list
class Chromosome(list):
  def __init__(self):
    self.elements = []



# Nest nearest integer
def nint(number):
  return int(number + (0.5 if number > 0 else -0.5))
  
# Calculate Euclidean distance for one location
def calculate_distance(x1, y1, x2, y2):
    return math.sqrt(math.pow((x1 - x2), 2) + math.pow((y1 - y2), 2))

def generateDistanceMatrix(file_name):
  # Open input file
  GRAPH = []
  newDist = []
  x = []
  y = []  
  tsp = tsplib95.load(file_name)
  if (tsp.is_explicit() and tsp.edge_weight_format=="LOWER_DIAG_ROW"):
    graph=tsp.edge_weights
    listgraph = [element for row in graph for element in row]
    M = tsp.dimension
    full_matrix = [[0] * M for _ in range(M)]
    k = 0
    for i in range(M):
        for j in range(i + 1):
            full_matrix[i][j] = listgraph[k]
            full_matrix[j][i] = listgraph[k] 
            k += 1
    return full_matrix, x, y
  elif(tsp.is_explicit()):
    graph=tsp.edge_weights
    return graph, x, y
  else:    
    infile = open(file_name, 'r') 

    # Read instance header
    name = infile.readline().strip().split()[1] # NAME
    fileType = infile.readline().strip().split()[1] # TYPE
    comment = infile.readline().strip().split()[1] # COMMENT
    dimension = infile.readline().strip().split()[1] # DIMENSION
    edgeWeightType = infile.readline().strip().split()[1] # EDGE_WEIGHT_TYPE
    infile.readline()

    # Read node list
    nodelist = []
    #N = int(dimension)
    print("Dimension: ", dimension)
    print("Edge weight type: ", edgeWeightType)
    for i in range(0, int(dimension)):
      x,y = infile.readline().strip().split()[1:]    
      newCity = City(float(x), float(y))
      nodelist.append(newCity)
        #nodelist.append([float(x), float(y)])

    # Close input file
    infile.close()
    GRAPH = []
    newDist = []
    x = []
    y = []  

    print("Dimension: ", dimension)
    N = int(dimension)
    for i  in range(0, N):
      x_i = nodelist[i].get_x()
      y_i = nodelist[i].get_y()
      x.append(x_i)
      y.append(y_i)
      newDist = []
      
      for j in range(0, len(nodelist)):
        x_j = nodelist[j].get_x()
        y_j = nodelist[j].get_y()
      
        
        #  edge weight type GEO means geographic coordinates
        if edgeWeightType == "GEO":
          deg = int(x_i)
          min = x_i - deg
          latitudeI = math.pi * (deg + 5.0 * min / 3.0) / 180.0
          
          deg = int(y_i)
          min = y_i - deg
          longitudeI = math.pi * (deg + 5.0 * min / 3.0) / 180.0
          
          deg = int(x_j)
          min = x_j - deg
          latitudeJ = math.pi * (deg + 5.0 * min / 3.0) / 180.0
          
          deg = int(y_j)
          min = y_j - deg
          longitudeJ = math.pi * (deg + 5.0 * min / 3.0) / 180.0
          
          R = 6378.388 # Earth radius
          q1 = math.cos( longitudeI - longitudeJ )
          q2 = math.cos( latitudeI - latitudeJ )
          q3 = math.cos( latitudeI + latitudeJ )
          dist = int((R * math.acos(0.5*((1.0+q1)*q2-(1.0-q1)*q3))+1.0))
          
          if i == j:
            dist = 0
        else:
          dist = nint(calculate_distance(x_i, y_i, x_j, y_j))
        if dist == 0:
          dist = 10000000
        #newDist = []
        #print("Distance 1: ", dist)
        newDist.append(dist)
      GRAPH.append(newDist)
    print("x: ", x)
    print("y: ", y)
      #print(GRAPH)
    return GRAPH, x, y

# coefficients for MT19937
(w, n, m, r) = (32, 624, 397, 31)
a = 0x9908B0DF
(u, d) = (11, 0xFFFFFFFF)
(s, b) = (7, 0x9D2C5680)
(t, c) = (15, 0xEFC60000)
l = 18
f = 1812433253


# make a arry to store the state of the generator
MT = [0 for i in range(n)]
index = n+1
lower_mask = 0x7FFFFFFF #(1 << r) - 1 // That is, the binary number of r 1's
upper_mask = 0x80000000 #lowest w bits of (not lower_mask)


# initialize the generator from a seed
def mt_seed(seed):
    # global index
    # index = n
    MT[0] = seed
    for i in range(1, n):
        temp = f * (MT[i-1] ^ (MT[i-1] >> (w-2))) + i
        MT[i] = temp & 0xffffffff


# Extract a tempered value based on MT[index]
# calling twist() every n numbers
def extract_number(seed):
    mt_seed(seed)
    global index
    if index >= n:
        twist()
        index = 0

    y = MT[index]
    y = y ^ ((y >> u) & d)
    y = y ^ ((y << s) & b)
    y = y ^ ((y << t) & c)
    y = y ^ (y >> l)

    index += 1
    return float(y) / 0xffffffff

# Generate the next n values from the series x_i
def twist():
    for i in range(0, n):
        x = (MT[i] & upper_mask) + (MT[(i+1) % n] & lower_mask)
        xA = x >> 1
        if (x % 2) != 0:
            xA = xA ^ a
        MT[i] = MT[(i + m) % n] ^ xA

if __name__ == "__main__":
  # Read a TSP file and convert x,y coordinates to distance matrix
  dimension = None
  epoch_stop = 0
  tsp_filename = 'ulysses16.tsp'
  # Open the TSP file in read mode
  with open(tsp_filename, "r") as tsp_file:
      # Iterate through the lines of the file
      for line in tsp_file:
          # Look for the line that starts with "DIMENSION:"
          if line.startswith("DIMENSION:"):
              # Split the line into words and extract the numeric value
              dimension = int(line.split(":")[1].strip())
              break  # End the search once the dimension is found

  #    print("La dimensión del problema TSP es:", dimension)
  epoch_stop=int(round(.2*dimension)) 
  
  
  #GRAPH, x, y = generateDistanceMatrix('att48.tsp')
  #GRAPH, x, y = generateDistanceMatrix('berlin52.tsp')
  #GRAPH, x, y  = generateDistanceMatrix('ulysses16.tsp')
  #GRAPH, x, y  = generateDistanceMatrix('ulysses22.tsp')
  GRAPH, x, y  = generateDistanceMatrix(tsp_filename)
  #GRAPH, x, y  = generateDistanceMatrix('a280.tsp')
  #GRAPH, x, y  = generateDistanceMatrix('capp1.tsp')
  #GRAPH, x, y  = generateDistanceMatrix('capp2.tsp')
  #GRAPH_SIZE = len(GRAPH)
  #GRAPH, x, y  = generateDistanceMatrix('clientes.tsp')
  #GRAPH, x, y  = generateDistanceMatrix('fl417.tsp')

  
  """
  # The following data comes from https://developers.google.com/optimization/routing/tsp/tsp
  # They report the following solution:
  # 0 - 7 - 2 - 3 - 4 - 12 - 6 - 8 - 1 - 11 - 5 - 10 - 9 - 0
  # with a cost of 7293 (Total distance in miles)
  GRAPH = [[10000, 2451, 713, 1018, 1631, 1374, 2408, 213, 2571, 875, 1420, 2145, 1972], # New York
  [2451, 10000, 1745, 1524, 831, 1240, 959, 2596, 403, 1589, 1374, 357, 579], # Los Angeles
  [ 713, 1745, 10000, 355, 920, 803, 1737, 851, 1858, 262, 940, 1453, 1260], # Chicago
  [1018, 1524, 355, 10000, 700, 862, 1395, 1123, 1584, 466, 1056, 1280, 987], # Minneapolis
  [1631, 831, 920, 700, 10000, 663, 1021, 1769, 949, 796, 879, 586, 371], # Denver
  [1374, 1240, 803, 862, 663, 10000, 1681, 1551, 1765, 547, 225, 887, 999], # Dallas
  [2408, 959, 1737, 1395, 1021, 1681, 10000, 2493, 678, 1724, 1891, 1114, 701], # Seattle
  [ 213, 2596, 851, 1123, 1769, 1551, 2493, 10000, 2699, 1038, 1605, 2300, 2099], # Boston
  [2571, 403, 1858, 1584, 949, 1765, 678, 2699, 10000, 1744, 1645, 653, 600], # San Francisco
  [ 875, 1589, 262, 466, 796, 547, 1724, 1038, 1744, 10000, 679, 1272, 1162], # St. Louis
  [1420, 1374, 940, 1056, 879, 225, 1891, 1605, 1645, 679, 10000, 1017, 1200], # Houston
  [2145, 357, 1453, 1280, 586, 887, 1114, 2300, 653, 1272, 1017, 10000, 504], # Phoenix
  [1972, 579, 1260, 987, 371, 999, 701, 2099, 600, 1162, 1200, 504, 10000]] # Salt Lake City
  """
  
  
  
  

  GRAPH_SIZE = len(GRAPH)
  
  
  # creates the Graph instance
  graph = Graph(GRAPH_SIZE, GRAPH)
  
  #print("Graph: ", GRAPH)


  # creates a PSO instance
  # beta is the probability for a global best movement
  #pso = PSO(graph, iterations=5000, maxEpochs=1, size_population=150, beta=0.1, alfa=0.9)
  
  """
  pso = PSO(graph, iterations=1000, maxEpochs=50, size_population=15, beta=0.05, alfa=0.9)
  pso.run()
  """

  
  results = ["Solution", "Cost", "Comp. time", "last epoch", "epoch convergence","Convergence COSTS"]
  fileoutput = []
  fileoutput.append(results)
  # def run_solver_with_params(params):
  #   iters = params
  for i  in range(100):
    results = []
    #pso = Solver(epoch_stop, graph, iterations=1000, maxEpochs=200, size_population=pop_size, beta=0.1, alfa=0.25)
    pso = Solver(epoch_stop, graph, iterations=int(round((dimension*600)/7)), maxEpochs=300, size_population=7, beta=0.29, alfa=0.12)
    #Nodos(dimension)*Factor de 600 = iterations
    #pso = Solver(epoch_stop, graph, iterations= 1428, maxEpochs=200, size_population=8, beta=0.29, alfa=0.12)
    #start_time = datetime.now()
    start_process_time = process_time()
    pso.run() # runs the PSO algorithm
    results.append(pso.getGBest().getPBest())
    results.append(pso.getGBest().getCostPBest())
    epoch = pso.getEpoch()
    epoch_convergence = epoch-epoch_stop
    #dt = datetime.now() - start_time
    #ms = (dt.days * 24 * 60 * 60 + dt.seconds) * 1000 + dt.microseconds / 1000.0
    #results.append(ms)
    processing_time = (process_time() - start_process_time) * 1000.0
    results.append(processing_time)
    results.append(epoch)
    results.append(epoch_convergence)
    results.append(pso.costConvergence)
    fileoutput.append(results)
    #del pso
    gc.collect()    

  #parameter_sets = [457,914,1371,1829,2286]
  #parameter_sets = [1486,2971,4457,5943,7429]
  # parameter_sets = [2886,5771,8657,11543,14429]
  # iteration = 0
  # for params in parameter_sets:
  #   iteration += 1
  #   run_solver_with_params(params)
  #   print(f"Iteration {iteration}, Parameters: {params}")


  """
  results = ["alpha", "beta", "mean"]
  fileoutput = []
  fileoutput.append(results)
  alfa_value = 0.1
  while alfa_value < 1:
    beta_value = 0.1
    while beta_value < 1:
      results = []
      cost_sum = 0.0
      runs = 0
      for i in range(10):
        runs += 1
        pso = Solver(graph, iterations=1000, maxEpochs=50, size_population=15, beta=beta_value, alfa=alfa_value)
        pso.run()
        cost_sum = cost_sum + pso.getGBest().getCostPBest()
      i = 0
      average_cost = cost_sum / runs
      results.append(alfa_value)
      results.append(beta_value)
      results.append(average_cost)
      fileoutput.append(results)
      beta_value += 0.1
    alfa_value += 0.1
  """
  
  """
  # the alpha and beta values are obtained with LHS in R"
  results = ["alpha", "beta", "cost", "comp. time", "epoch"]
  fileoutput = []
  fileoutput.append(results)
  #alfa_values = [0.33083770, 0.81286553, 0.28994727, 0.62146038, 0.52186326, 0.17224262, 0.02020398, 0.78311170, 0.90367650, 0.49015838]
  #beta_values = [0.6197808, 0.89753027, 0.11981711, 0.72693961, 0.40612212, 0.94802543, 0.57793579, 0.23444987, 0.39597191, 0.09485114]
  alfa_values = [0.52966316, 0.4145423, 0.19379634, 0.87264983, 0.97362117, 0.06562221, 0.65973607, 0.2151281, 0.461971, 0.38040526, 0.77591784, 0.12027035, 0.82232354, 0.32580516, 0.58453551, 0.73574494, 0.92142149, 0.26890176, 0.61401872, 0.0111655]
  beta_values = [0.58285168, 0.35864099, 0.85454185, 0.13890171, 0.49436501, 0.23125492, 0.72400869, 0.0334176, 0.98079264, 0.82039207, 0.52108801, 0.40230925, 0.33304859, 0.15509154, 0.05296365, 0.92091665, 0.78477693, 0.62063391, 0.25896396, 0.65795874]
  for i in range(20):
    alfa_value = alfa_values[i]
    beta_value = beta_values[i]
    cost_sum = 0.0
    runs = 0
    print("alfa: ", alfa_value, "beta: ", beta_value)
    for i in range(20):
      results = []
      runs += 1
      pso = Solver(graph, iterations=1000, maxEpochs=50, size_population=15, beta=beta_value, alfa=alfa_value)
      start_time = datetime.now()
      pso.run()
      #cost_sum = cost_sum + pso.getGBest().getCostPBest()
      results.append(alfa_value)
      results.append(beta_value)
      results.append(pso.getGBest().getCostPBest())
      epoch = pso.getEpoch()
      dt = datetime.now() - start_time
      ms = (dt.days * 24 * 60 * 60 + dt.seconds) * 1000 + dt.microseconds / 1000.0
      results.append(ms)
      results.append(epoch)
      fileoutput.append(results)
    #average_cost = cost_sum / runs
    #results.append(alfa_value)
    #results.append(beta_value)
    #results.append(average_cost)
  """  
   

  # pso-results.csv
  
  
  csvFile = open('mepso_II.csv', 'w', newline='')  
  with csvFile: 
    writer = csv.writer(csvFile)
    writer.writerows(fileoutput)
  csvFile.close()
  

  # shows the global best particle
  ##bestPath = pso.getGBest().getPBest()
  ##print('gbest: %s | cost: %d\n' % (bestPath, pso.getGBest().getCostPBest()))
  ##plotPoints(x, y)
  ##plotTSP(bestPath, x, y)
  
  '''
  # random graph
  print('Random graph...')
  random_graph = CompleteGraph(amount_vertices=20)
  random_graph.generates()
  pso_random_graph = PSO(random_graph, iterations=10000, size_population=10, beta=1, alfa=1)
  pso_random_graph.run()
  print('gbest: %s | cost: %d\n' % (pso_random_graph.getGBest().getPBest(), 
          pso_random_graph.getGBest().getCostPBest()))
  '''
