# encoding:utf-8

########################################################################################################################
# A Particle Swarm Optimization algorithm for solving the traveling salesman problem.
# The program reuses part of the code of Marco Castro (https://github.com/marcoscastro/tsp_pso)
#  
# Authors: Rafael Batres, Ariann Fernando Arriaga Alcántara, Juan Carlos Espinoza García, José Carlos Bañales López.
# Institution: Tecnologico de Monterrey
# 
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
# import matplotlib.pyplot as plt
# import pandas as pd
from time import process_time
import tsplib95
#random.seed(0)

# PSO algorithm
class PSO:


  def __init__(self, iteration_stop, graph, iterations, size_population, alfa=1, beta=1):
    self.iteration_stop = iteration_stop
    self.deviation_iteration = int(round(iteration_stop/100))
    self.graph = graph # the graph
    self.iterations = iterations # max of iterations
    self.populationSize = size_population # size population
    self.particles = Chromosome() # list of particles
    self.beta = beta # the probability that all swap operators in swap sequence (gbest - x(t-1))
    self.alfa = alfa # the probability that all swap operators in swap sequence (pbest - x(t-1))
    #graph_size = 5
    
    # initialized with a group of random particles (solutions)
    solutions = self.graph.getRandomPaths(self.populationSize)
       
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
  
  def setIter(self, last_iter):
    self.last_iter = last_iter
  
  def getIter(self):
    return self.last_iter

  def run(self):
    # variables for convergence data
    convergenceData = []
    iterationArray = []
    bestCostArray = []
    bestCostSampling = []


    batchSize = 100 # save data every n iterations
    batchCounter = 0
    iteration = 0
    
    startTime = datetime.now()
    
    # updates gbest (best particle of the population)
    self.gbest = self.particles[0]
        
      # for each time step (iteration)
    iterations = self.iterations
    t = 0
    while t < iterations:
    #for t in range(self.iterations):
      convergencePerIteration = []
      batchCounter = batchCounter + 1
      
      # for each particle in the swarm
      for particle in self.particles:
        previousCost = particle.getCurrentSolutionCost()
        particle.clearVelocity() # cleans the speed of the particle
        velocity = []

        gbest = list(self.gbest.getPBest()) # gets solution of the gbest solution
        pbest = particle.getPBest()[:] # copy of the pbest solution
        newSolution = particle.getCurrentSolution()[:] # gets copy of the current solution of the particle
        

        # generates all swap operators to calculate (pbest - x(t-1))
        for i in range(self.graph.graphSize):
          if newSolution[i] != pbest[i]:
            # generates swap operator
            swapOperation = (i, pbest.index(newSolution[i]), self.alfa)
            # append swap operator in the list of velocity
            velocity.append(swapOperation)
            # makes the swap           
            aux = pbest[swapOperation[0]]
            pbest[swapOperation[0]] = pbest[swapOperation[1]]
            pbest[swapOperation[1]] = aux
            
            #pbest[swapOperation[0]], pbest[swapOperation[1]] = pbest[swapOperation[1]], pbest[swapOperation[0]] 


        # generates all swap operators to calculate (gbest - x(t-1))
        for i in range(self.graph.graphSize):
          if newSolution[i] != gbest[i]:
            # generates swap operator
            swapOperation = (i, gbest.index(newSolution[i]), self.beta)
            # append swap operator in the list of velocity
            velocity.append(swapOperation)
            # makes the swap         
            aux = gbest[swapOperation[0]]
            gbest[swapOperation[0]] = gbest[swapOperation[1]]
            gbest[swapOperation[1]] = aux
            
            #gbest[swapOperation[0]], gbest[swapOperation[1]] = gbest[swapOperation[1]], gbest[swapOperation[0]] 
        
        
        # updates velocity
        particle.setVelocity(velocity)
                
        
        # Inverts the cities between the selected points        
        for swapOperation in velocity:
          if random.random() <= swapOperation[2]:
            point1 = swapOperation[0]
            point2 = swapOperation[1]
            while True:
              newSolution[point1],newSolution[point2] = newSolution[point2],newSolution[point1]
              point1 = point1 + 1
              point2 = point2 - 1
              if point1 > point2:
                break
              
        
        # gets cost of the current solution
        newSolutionCost = self.graph.evaluateCost(newSolution)
        
        if newSolutionCost < previousCost:
          # updates the current solution  
          particle.setCurrentSolution(newSolution)
          # updates the cost of the current solution
          particle.setCostCurrentSolution(newSolutionCost)
                

        # checks if new solution is pbest solution
        pbCost =  particle.getCostPBest()
        if newSolutionCost < pbCost:
          particle.setPBest(newSolution)
          particle.setCostPBest(newSolutionCost)
       
        
        gbestCost = self.gbest.getCostPBest()       
        # check if new solution is gbest solution
        if newSolutionCost < gbestCost:
          self.gbest = particle

        convergencePerIteration.append(t)
        convergencePerIteration.append(self.gbest.getCostPBest())
        convergenceData.append(convergencePerIteration)
        iterationArray.append(t)
        bestCostArray.append(self.gbest.getCostPBest())
      

        if batchCounter > batchSize:
          print(t, "Gbest cost = ", self.gbest.getCostPBest())
          #print("Standard deviation: ", std)
          batchCounter = 0
          bestCostSampling.append(self.gbest.getCostPBest()) 
          if t >= self.iteration_stop :
             bestCostSampling = bestCostSampling[1:] 
      t = t + 1
      if t >= self.iteration_stop :
        std = statistics.pstdev(bestCostSampling[-self.deviation_iteration:])
      else:
        std = 1000
      
      if t==self.iterations:
        self.setIter(t)

      if std == 0:
        self.setIter(t)
        break

    
    # df = pd.DataFrame()
    # df['Iteration'] = pd.Series(iterationArray)
    # df['Best cost'] = pd.Series(bestCostArray)
    # plt.xlabel("Iteration No.")
    # plt.ylabel("Best cost")
    # plt.plot(df['Iteration'], df['Best cost'])

    """
    csvFile = open('pso-convergence.csv', 'w', newline='')  
    with csvFile: 
         writer = csv.writer(csvFile)
         writer.writerows(convergenceData)
    csvFile.close()
    
    """ 
    print("Elapsed time: ", self.elapsedTime(startTime))
  
  

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
    cost = cost + self.costTable[last][first]
    return cost



  # gets random unique paths - returns a list of lists of paths
  def getRandomPaths(self, max_size):
    random_paths = []
    
    for i in range(max_size):
      
      list_temp = self.getRandomSolution(self.graphSize)

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
            full_matrix[j][i] = listgraph[k]  # Reflejar los valores en la parte superior de la diagonal
            k += 1
    return full_matrix
  elif(tsp.is_explicit()):
    graph=tsp.edge_weights
    return graph
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
    print("Dimension: ", dimension)
    N = int(dimension)
    for i  in range(0, N):
      x_i = nodelist[i].get_x()
      y_i = nodelist[i].get_y()
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
      #print(GRAPH)
    return GRAPH

if __name__ == "__main__":
  # Read a TSP file and convert x,y coordinates to distance matrix
  # Uncomment the line that corresponds to the TSP file that you want to use
  dimension = None
  iteration_stop = 0
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
  iteration_stop=int(round(((dimension**2)*.2*600)/150)) 
  
  
  #GRAPH = generateDistanceMatrix('clientes2.tsp')
  GRAPH = generateDistanceMatrix(tsp_filename)
  #GRAPH = generateDistanceMatrix('ulysses22.tsp')
  #GRAPH = generateDistanceMatrix('berlin52.tsp')
  #GRAPH = generateDistanceMatrix('eil101.tsp')
  #GRAPH = generateDistanceMatrix('a280.tsp')
  #GRAPH = generateDistanceMatrix('capp1.tsp')
  #GRAPH = generateDistanceMatrix('capp2.tsp')
  #GRAPH = generateDistanceMatrix('fl417.tsp')
  #GRAPH_SIZE = len(GRAPH)
  
  GRAPH_SIZE = len(GRAPH)
  
  
  # creates the Graph instance
  graph = Graph(GRAPH_SIZE, GRAPH)


  # creates a PSO instance
  # alfa is the probabiliy for a movement based on local best
  # beta is the probability for a movement based on the global best
  
 
  
  results = ["Solution", "Cost", "Comp. time", "Max iter","iter convergence"]
  fileoutput = []
  fileoutput.append(results)
  # def run_solver_with_params(params):
  #       alfa_prob, beta_prob = params
  for i  in range(100):
    results = []
    #pso = Solver(graph, iterations=1000, maxEpochs=80, size_population=15, beta=0.80, alfa=0.2)
    #pso = Solver(graph, iterations=1000, maxEpochs=50, size_population=15, beta=0.51, alfa=0.11)
    pso = PSO(iteration_stop, graph, iterations=int(round(600*dimension*((7*200)/150))) , size_population=150, alfa=0.47442915, beta=0.03153496)
    #Nodos(dimension)*Factor de 600 = iterations
    #start_time = datetime.now()
    start_process_time = process_time()
    pso.run() # runs the PSO algorithm
    results.append(pso.getGBest().getPBest())
    results.append(pso.getGBest().getCostPBest())
    iteration = pso.getIter()
    iteration_convergence = iteration-iteration_stop
    # dt = datetime.now() - start_time
    # ms = (dt.days * 24 * 60 * 60 + dt.seconds) * 1000 + dt.microseconds / 1000.0
    # results.append(ms)
    processing_time = (process_time() - start_process_time) * 1000.0
    results.append(processing_time)
    results.append(iteration)
    results.append(iteration_convergence)
    fileoutput.append(results)
    # shows the global best particle
    print("")
    print("alfa: ", pso.alfa, "beta: ", pso.beta)
    print("Cost of gbest: ", pso.getGBest().getCostPBest()) 
    print("gbest: ", pso.getGBest().getPBest())
    #del pso
    
  # parameter_sets = [
  #   (0.01408918, 0.36408546),
  #   (0.83727808, 0.7652775),
  #   (0.04757019, 0.05348702),
  #   (0.68471084, 0.69839838),
  #   (0.52034017, 0.10770506),
  #   (0.07613794, 0.79369761),
  #   (0.2761332, 0.37480791),
  #   (0.21332816, 0.8332606),
  #   (0.66103729, 0.18269137),
  #   (0.72807695, 0.14628294),
  #   (0.93434318, 0.31469016),
  #   (0.98797699, 0.65905174),
  #   (0.93257113, 0.55087066),
  #   (0.62747079, 0.09556384),
  #   (0.47442915, 0.03153496),
  #   (0.15301107, 0.2279023),
  #   (0.45607367, 0.40835087),
  #   (0.19590808, 0.50628222),
  #   (0.88164369, 0.99902693),
  #   (0.59651399, 0.86278655),
  #   (0.73751353, 0.61679432),
  #   (0.24738666, 0.92635521),
  #   (0.32357573, 0.93966287),
  #   (0.38873597, 0.71028227),
  #   (0.54604574, 0.29745743),
  #   (0.10424685, 0.57967287),
  #   (0.83216319, 0.87321059),
  #   (0.42194966, 0.46582887),
  #   (0.35435905, 0.46811074),
  #   (0.77875317, 0.24357927)
  # ]
  
  # iteration = 0
  # for params in parameter_sets:
  #   iteration += 1
  #   run_solver_with_params(params)
  #   print(f"Iteration {iteration}, Parameters: {params}")
  
  csvFile = open('pso-simple-stdev.csv', 'w', newline='')  
  with csvFile: 
    writer = csv.writer(csvFile)
    writer.writerows(fileoutput)
  csvFile.close()


  # shows the global best particle
  #print("")
  #print("alfa: ", pso.alfa, "beta: ", pso.beta)
  #print("Cost of gbest: ", pso.getGBest().getCostPBest()) 
  #print("gbest: ", pso.getGBest().getPBest())

