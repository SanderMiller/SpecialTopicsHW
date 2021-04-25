import picos as pic
from picos import RealVariable
from copy import deepcopy
from heapq import *
import heapq as hq
import numpy as np
import itertools
import math
counter = itertools.count() 

class BBTreeNode():
    def __init__(self, vars = [], constraints = [], objective='', prob=None):
        self.vars = vars
        self.constraints = constraints
        self.objective = objective
        self.prob = prob

    def __deepcopy__(self, memo):
        '''
        Deepcopies the picos problem
        This overrides the system's deepcopy method bc it doesn't work on classes by itself
        '''
        newprob = pic.Problem.clone(self.prob)
        return BBTreeNode(self.vars, newprob.constraints, self.objective, newprob)
    
    def buildProblem(self):
        '''
        Bulids the initial Picos problem
        '''
        prob=pic.Problem()
   
        prob.add_list_of_constraints(self.constraints)    
        
        prob.set_objective('max', self.objective)
        self.prob = prob
        return self.prob

    def is_integral(self):
        '''
        Checks if all variables (excluding the one we're maxing) are integers
        '''
        for v in self.vars[:-1]:
            if v.value == None or abs(round(v.value) - float(v.value)) > 1e-4 :
                return False
        return True

    def branch_floor(self, branch_var):
        '''
        Makes a child where xi <= floor(xi)
        '''
        n1 = deepcopy(self)
        n1.prob.add_constraint( branch_var <= math.floor(branch_var.value) ) # add in the new binary constraint

        return n1

    def branch_ceil(self, branch_var):
        '''
        Makes a child where xi >= ceiling(xi)
        '''
        n2 = deepcopy(self)
        n2.prob.add_constraint( branch_var >= math.ceil(branch_var.value) ) # add in the new binary constraint
        return n2

    def bbsolve(self):
        '''
        Use the branch and bound method to solve an integer program
        This function should return:
            return bestres, bestnode_vars

        where bestres = value of the maximized objective function
              bestnode_vars = the list of variables that create bestres
        '''

        # these lines build up the initial problem and adds it to a heap
        root = self
        res = root.buildProblem().solve(solver='cvxopt')
        heap = [(res, next(counter), root)]
        bestres = -1e20 # a small arbitrary initial best objective value
        bestnode_vars = root.vars # initialize bestnode_vars to the root vars

        #While the heap is non-empty
        while(len(heap)>0):

            #Pop off current problem
            currRes, currCount, currProb = hq.heappop(heap)

            #Try solving current problem
            try:
                res = currProb.prob.solve(solver='cvxopt')

            except:
                pass

            #Base case, if all are integers return
            if currProb.is_integral():
                return bestres, bestnode_vars

            #Otherwise recursively build heap
            else:
                
                for var in currProb.vars:
                    #Check if a given variable is an integer
                    if var.value == None or abs(round(var.value) - float(var.value)) > 1e-4 :

                        #If not an integer generate floor and ceiling case children
                        child1 = currProb.branch_ceil(var)
                        child2 = currProb.branch_floor(var)

                        #Push new cases and corresponding variables to the heap
                        hq.heappush(heap,(child1.vars, next(counter),child1))
                        hq.heappush(heap,(child2.vars, next(counter),child2))
                        break
        

