#imports
from os import mkdir
import math
from statistics import median, stdev
from matplotlib import pyplot as plt
from time import gmtime, strftime, time
from random import uniform, choice, randint
import uuid


class PSO:

    def __init__(self):
        self.pop = [] #population's positions        
        self.m_nmdf = 0.00 #diversity variable
        self.diversity = []
        self.fbest_list = []

    def generateGraphs(self, fbest_list, diversity_list, max_iterations, uid, run):
        plt.plot(range(0, max_iterations), fbest_list, 'r--')
        plt.savefig(str(uid) + '/graphs/run' + str(run) + '_' + 'convergence.png')
        plt.clf()                                                   
        plt.plot(range(0, max_iterations), diversity_list, 'b--')
        plt.savefig(str(uid) + '/graphs/run' + str(run) + '_' + 'diversity.png')
        plt.clf()                                                   
        
    def updateDiversity(self):
        diversity = 0
        aux_1 = 0
        aux2 = 0
        a = 0
        b = 0
        d = 0
        
       
        for a in range(0, len(self.pop)):
            b = a+1
            for i in range(b, len(self.pop)):
                aux_1 = 0
    
                ind_a = self.pop[a]
                ind_b = self.pop[b]
    
                for d in range(0, len(self.pop[0])):
                    aux_1 = aux_1 + (pow(ind_a[d] - ind_b[d], 2).real)
                aux_1 = (math.sqrt(aux_1).real)
                aux_1 = (aux_1 / len(self.pop[0]))
    
                if b == i or aux_2 > aux_1:
                    aux_2 = aux_1
            diversity = (diversity) + (math.log((1.0) + aux_2).real)
            
        if self.m_nmdf < diversity:
            self.m_nmdf = diversity

        return (diversity/self.m_nmdf).real

    #fitness_function
    def fitness(self, individual):
        'to override'
        'rastrigin' 
        result = 0.00
        for dim in individual:
            result += (dim - 1)**2 - 10 * math.cos(2 * math.pi * (dim - 1))
        return (10*len(individual) + result)
    

    #pso_function========================================================================================================================
    def particleSwarmOptimization(self, pop_size, dim, bounds, max_iterations, runs, vmax=0.1, c1=2, c2=2, inertia=0.7, maximize=True):
        #generete execution identifier
        uid = uuid.uuid4()
        mkdir(str(uid))
        mkdir(str(uid) + '/graphs')
        #to record the results
        results = open(str(uid) + '/results.txt', 'a')
        records = open(str(uid) + '/records.txt', 'a')
        results.write('ID: %s\tDate: %s\tRuns: %s\n' % (str(uid ), strftime("%Y-%m-%d %H:%M:%S", gmtime()), str(runs)))
        results.write('=================================================================================================================\n')
        records.write('ID: %s\tDate: %s\tRuns: %s\n' % (str(uid ), strftime("%Y-%m-%d %H:%M:%S", gmtime()), str(runs)))
        records.write('=================================================================================================================\n')
        avr_fbest_r = []
        avr_diversity_r = []
        fbest_r = []
        best_r = []
        elapTime_r = []
        #runs
        for r in range(runs):
            elapTime = []
            start = time()
            records.write('Run: %i\n' % r)
            records.write('Iter\tGbest\tAvrFit\tDiver\tETime\t\n')
            vpop = [] #population's velocities
            best = [] #global best positions
            fbest = 0.00
                    
            #global best fitness
            if maximize == True:
                fbest = 0.00
            else:
                fbest = math.inf

            #initial_generations
            for ind in range(pop_size):
                lp = []
                lv = []
                for d in range(dim):
                    lp.append(uniform(bounds[d][0],bounds[d][1]))
                    lv.append(uniform(bounds[d][0]*vmax,bounds[d][1]*vmax))
                self.pop.append(lp)
                vpop.append(lv)
            
            fpop = []
            for ind in self.pop:
                fpop.append(self.fitness(ind))
            
            ppop = [values for values in self.pop]
            fppop = [values for values in fpop]
            
            for ind in range(pop_size):
                if maximize == True:
                    if fpop[ind] >= fbest:
                        fbest = float(fpop[ind])
                        best = [values for values in self.pop[ind]]
                else:     
                    if fpop[ind] <= fbest:
                        fbest = float(fpop[ind])
                        best = [values for values in self.pop[ind]]
            
        
            #evolution_step
            for iteration in range(max_iterations):
                avrFit = 0.00
                #update_solutions
                for ind in range(pop_size):
                    for d in range(dim):
                        vpop[ind][d] = inertia*vpop[ind][d] + (c1*uniform(0,1)*(ppop[ind][d]-self.pop[ind][d])) + (c2*uniform(0,1)*(best[d]-self.pop[ind][d]))
                        self.pop[ind][d] = self.pop[ind][d] + vpop[ind][d]
                        if self.pop[ind][d] < bounds[d][0]:
                            self.pop[ind][d] = float(bounds[d][0])
                        if self.pop[ind][d] > bounds[d][1]:
                            self.pop[ind][d] = float(bounds[d][1])
                fpop = []
                for ind in self.pop:
                    aux = self.fitness(ind)
                    fpop.append(aux)
                    avrFit += aux
                avrFit = avrFit/pop_size
                self.diversity.append(self.updateDiversity())
                                            
                for ind in range(pop_size):
                    if maximize == True:
                        if fpop[ind] > fbest:
                            fbest = float(fpop[ind])
                            best = [values for values in self.pop[ind]]

                        if fpop[ind] > fppop[ind]:
                            fppop[ind] = float(fpop[ind])
                            ppop[ind] = [values for values in self.pop[ind]]
                    else:
                        if fpop[ind] < fbest:
                            fbest = float(fpop[ind])
                            best = [values for values in self.pop[ind]]

                        if fpop[ind] < fppop[ind]:
                            fppop[ind] = float(fpop[ind])
                            ppop[ind] = [values for values in self.pop[ind]]
                
                self.fbest_list.append(fbest)
                elapTime.append((time() - start)*1000.0)
                records.write('%i\t%.4f\t%.4f\t%.4f\t%.4f\n' % (iteration, round(fbest,4), round(avrFit,4), round(self.diversity[iteration],4), elapTime[iteration]))
            records.write('Pos: %s\n\n' % str(best))
            fbest_r.append(fbest)
            best_r.append(best)
            elapTime_r.append(elapTime[max_iterations-1])
            self.generateGraphs(self.fbest_list, self.diversity, max_iterations, uid, r)
            avr_fbest_r.append(self.fbest_list)
            avr_diversity_r.append(self.diversity)
            
            self.pop = []
            self.m_nmdf = 0.00 
            self.diversity = []
            self.fbest_list = []
        
        
        fbestAux = [sum(x)/len(x) for x in zip(*avr_fbest_r)]
        diversityAux = [sum(x)/len(x) for x in zip(*avr_diversity_r)]
        self.generateGraphs(fbestAux, diversityAux, max_iterations, uid, 'Overall')
        records.write('=================================================================================================================')
        if maximize==False:
            results.write('Gbest Overall: %.4f\n' % (min(fbest_r)))
            results.write('Positions: %s\n\n' % str(best_r[fbest_r.index(min(fbest_r))]))
        else:
            results.write('Gbest Overall: %.4f\n' % (max(fbest_r)))
            results.write('Positions: %s\n\n' % str(best_r[fbest_r.index(max(fbest_r))]))

        results.write('Gbest Average: %.4f\n' % (sum(fbest_r)/len(fbest_r)))
        results.write('Gbest Median: %.4f #probably should use median to represent due probably non-normal distribution (see Shapiro-Wilk normality test)\n' % (median(fbest_r)))
        if runs > 1:
            results.write('Gbest Standard Deviation: %.4f\n\n' % (stdev(fbest_r)))
        results.write('Elappsed Time Average: %.4f\n' % (sum(elapTime_r)/len(elapTime_r)))
        if runs > 1:
            results.write('Elappsed Time Standard Deviation: %.4f\n' % (stdev(elapTime_r)))
        results.write('=================================================================================================================\n')


if __name__ == '__main__': 
    from pso import PSO 
     
    max_iterations = 200 
    pop_size = 20
    dim = 2 
    runs = 10
    bounds = ((-5.12,5.12), (-5.12,5.12))
    p = PSO()
    p.particleSwarmOptimization(pop_size, dim, bounds, max_iterations, runs, maximize=False)

