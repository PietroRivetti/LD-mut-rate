"""

This class implements an experiment that aims to estimate the spontaneous mutation rate 
of cell lines, made by DNA sequencing.

It implements the generation of the population tree, assignements of mutations across the tree
ad the estimation of the mutation rate. This is performed by also accounting for some error models:

- Mutations prior to cloning
- Ploidy
- Mean sequencing error 

"""


# import
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.stats import binom
from anytree import Node, RenderTree


class LD:
    
    def __init__(self, gen = 5, min_size = 0, bases = 100E+6,  
                 death_prob = 0., mu = 3E-9, N0gen = 0,  seq_err = 0.001, 
                 cov = 128, quiet = True, accept_extinct = False):
        """ 
        Cells are represented as nodes in a binary tree structure.
        Population tree is built based on:
        - minimum extant cell number, if a 'min_size' positive argument is passed,
        - fixed number of generations, otherwise.
        
        If 'accept_extinct' is False only non extinct trees are accepted.
        
        """
        
        # Options for the run
        self._quiet = quiet # if 'True' hide some prints
        self._min_size = min_size # if positive is the minimum size to reach
        self._gen = gen # number of generations if no 'min_size' is passed
        self._accept_extinct = accept_extinct
        
        #-----------------------
        # SIMULATION PARAMETERS
        #-----------------------
        
        # Cell's parameters
        self._bases = bases # cell's bases number (realizations of the mutation process) 
        self._death_prob = death_prob # death probability per generation d/(d+b)
        self._mu = mu # real mutation rate per generation
        
        # Parameters for error models
        self._N0gen = N0gen # generations from reference to clone
        self._seq_err = seq_err # sequencing error rate per base
        self._cov = cov # coverage
        self._threshold = 0. # threshold for frequencies
        self._mean_ploidy = 0 # mean ploidy of cells
        
        #-----------------------
        # RESULTS
        #-----------------------
        
        # Mutation rate estimates for different methods
        self._mu_est = 0. # estimated mutation rate (clean)
        self._mu_est_N0gen = 0. # estimated mutation rate considering also mutation prior to cloning
        self._mu_est_N0_threshold = 0 # mutation prior to cloning and threshold
        self._mu_est_seq_err = 0. # estimated mutation rate considering also sequencing errors
        self._mu_est_seq_err_threshold = 0. #sequencing errors and threshold
        self._mu_est_ploidy = 0. # estimated mutation rate considering also ploidy
        self._mu_est_ploidy_with_correction = 0. # ploidy and threshold
        self._mu_est_no_dead = 0. # estimated mutation rate counting only edges that have alive progeny
        
        #-----------------------
        # CELLS' DATA STRUCTURES 
        #-----------------------
        
        # We create a list of fixed initial size to contain population
        # If needed, more elements will be added letely 
        self._static_list_size = self.max_imax(self._gen) + 1
        self._population = [None for _ in np.arange(self._static_list_size)] # List of population nodes
        
        self._extant =  0 # number of nodes in the last generation
        self._num_layers = 0 # number of layers (generations) in the tree 
        self._extinct = True # bool to check if tree is extinct before reaching the right size
        self._imax = 0 # max cell index = number of nodes in the whole tree
        
        self._min_extant_id = 0 # min ID between extant cells
                                # because ID are ordered, every cell with ID in 
                                # [self.min_extant_id, self.imax] is an extant cell
        
        self._count_with_alive_progeny = 0 # number of cells that have alive progeny (descendants)
                                           # in the whole tree (extant cells are counted)
        #---------------------------
        # MUTATIONS DATA STRUCTURES 
        #---------------------------     
        self._mutation_frequencies = [] # mutation frequency vector
        self._count_mutation = 0 # counter for mutiations that are not extict in the last layer
        self._n_mut_prior = 0 # mutation occurring in N0gen from reference to clone zero (root)
        
        
        # Fill population data structure (tree generation)
        
        # Check if extinct trees are accepted
        if self._accept_extinct:
            # Extinct trees are accepted
            self.generate_tree()
        else:
            #Extinct trees are discarded          
            while(self._extinct):
                #generate new tree
                self.generate_tree()
        
        # we can not estimate mut rate for an extinct tree
        if not self._extinct:
            # Generate mutation
            self.generate_mutations()
            # Test "standard" estimator (clean from errors)
            self.test_LD_estimator()
        else:
            if not self._quiet:
                print("Tree is extinct, can not estimate mutation rate.")
       
    def generate_tree(self):
        """
        Tree is built based on:
        - minimum extant cell number, if a 'min_size' positive argument is passed,
        - fixed number of generation, otherwise.
        """
        
        #set (or reset in case of extinction) cells' data structures 
        self._static_list_size = self.max_imax(self._gen) + 1
        self._population = [None for _ in np.arange(self._static_list_size)] # List of population nodes
        
        self._extant =  0 # number of nodes in the last generation
        self._num_layers = 0 # number of layers (generations) in the tree 
        self._imax = 0 # max cell index = number of nodes in the whole tree
        self._min_extant_id = 0
        #set extinction flag to False
        # will be set to 'True' if extintion happens during tree's generation
        self._extinct = False 
        
        # if minimum size parameter is passed
        if self._min_size > 0:
            #build tree until minimum population size
            self.generate_tree_min_size()
        else:
            # there is no minimum size
            # build for fixed number of generation
            self.generate_tree_fixed_gen()
        
    def generate_tree_min_size(self):
        """
        Generate population until fixed minimum size.
        """            
        # print 'start'
        if not self._quiet:
            print("\n\nGenerating new tree...")
            
        # Every cell is represented as a node in binary tree
        # Create cell zero (root)
        self._population[0] = Node(0)
        self._imax = 0
        self._extant = 1
        self._min_extant_id = 0
        
        # current number of layers in the tree is zero
        self._num_layers = 0
        
        #while number of extant cells is less than minimum size
        while(self._extant < self._min_size):
            # if tree is not extinct
            if not self._extinct:
                #add generation (layer)
                self.add_generation() 
            else:
                #tree is estinct
                #stop generations
                break         
    
    def generate_tree_fixed_gen(self):
        """
        Generate population for fixed number of generations.
        """
        # print 'start'
        if not self._quiet:
            print("\n\nGenerating new tree...")
        
        # Every cell is represented as a node in binary tree
        # Create cell zero (root)
        self._population[0] = Node(0)
        self._imax = 0
        self._extant = 1
        self._min_extant_id = 0
        
        # current number of layers in the tree is zero
        self._num_layers = 0
        
        # for each generation from 0 (now) to (fixed gen number - 1)
        for current_gen in range(self._gen):
            # if tree is not extinct
            if not self._extinct:
                # add generation (layer)
                self.add_generation()
            else:
                # tree is extinct
                # stop generations
                break
                
    def add_generation(self):
        """
        Add one generation to the tree.
        """
        # save current max cell's index
        old_imax = self._imax
        
        # for every cell (node) in the last generation (layer)
        for parent in self.get_last_layer():
            # attempt to generate two daughters
            for attempt in np.arange(2):
                # if death event not verified
                if (np.random.random() > self._death_prob): 
                    # make new cell
                    self._imax += 1
                    # add cell to population
                    # check if "static" use of the list is possible
                    # if imax is too large it means that the list is full
                    if self._imax < self._static_list_size:
                        # insert it in the static list 
                        self._population[self._imax] = Node(self._imax, parent = parent)
                    else:
                        # "static" list is full
                        # it happens when building tree for min_size
                        self._population.append(Node(self._imax, parent = parent))
                        
        
        # update layers' counter            
        self._num_layers += 1
        # update extant number
        self._extant = self._imax - old_imax
        self._min_extant_id = old_imax + 1
        
        # Check for extinction:
        # if no cells were added this generation
        if self._extant == 0:
            # tree is extinct
            self._extinct = True
            # zero extant 
            self._extant = 0
            # minimum extant id is not defined
            self._min_extant_id = None
            # print
            if not self._quiet:  
                    print("\n\nEmpty layer! Tree is extinct.")
                    self.print_tree()
        else:
            #tree is not extinct
            self._extinct = False
              
    def generate_mutations(self):
        """
        Generate mutations and assign them to specific cells (mutants).
        For every mutant check how many descendants it has
        in extant cells. 
        Compute mutation frequency as the number of mutants in extant
        cells, for each mutation.
        """
        # We have 'self._extant' extant nodes. 
        # The number of mutation attempts is
        #           imax * bases
        # compute total number of mutations
        attempts = int(self._bases * self._imax)
        n_mutations = np.random.binomial(attempts, self._mu)
        
        if not self._quiet:
            print("\n\nMutation attempts: " + str(attempts))
            print("Mutations number: " + str(n_mutations))
        
        # assign the mutations in the tree
        for _ in np.arange(n_mutations):
            
            # overall index of mutant is sampled uniformly in [1, imax]
            rand_int = 1 + np.random.randint(attempts)
            
            if not self._quiet:
                print("Rand_int: " + str(rand_int) + ", imax: " + str(self._imax) + 
                      ", randint/imax: " + str(rand_int/self._imax) )
                
            # compute cellID of mutant in tree
            cellID = rand_int % self._imax
            
            # check how many mutants there are in progeny
            n_mutants = self.count_alive_progeny(cellID)
            
            # compute frequency = n_mutants/n_extant
            frequency = n_mutants / self._extant 
                
            # check if mutation is not extinct in the tree
            if n_mutants != 0 :         
                # save mutation frequency
                self._mutation_frequencies.append(frequency)
                # increment counter for mutations that are not extict in the last layer
                self._count_mutation += 1
                
            if not self._quiet:
                print("Mutant cell ID: " + str(cellID) + ", Mutants in progeny: " + str(n_mutants) +
                          ", Extant: " + str(self._extant) + ", Frequency: " + str(frequency) )
        if not self._quiet:
            print("\n")
            
    def test_LD_estimator(self):
        """
        Estimate mutation rate with "standard" formula:
        
        P = (# mutations not extinct in last layer)/(# bases)
        
        mu = - log( 1 - P ) / imax
        
        Equation 3, page 13. 
        """
        
        # clean mut rate estimate
        self._mu_est  = - math.log(1 - self._count_mutation/self._bases)
        self._mu_est /= self._imax
        
        if not self._quiet:
            print("Clean mutation count: ", self._count_mutation)
            print("Clean mutation rate estimated: "+str(self._mu_est))
            print("\n")
    
    #----------------------------
    # ERROR MODELS
    #----------------------------
    
    # role of N0_gen
    def test_LD_estimator_N0gen(self, N0gen = 0):
        
        """ 
        Test the estimator in presence of errors due to mutations occurring in 
        'N0gen' generations from reference to clone. 
        Compute attempts as:
                            mutation attempts = N0gen * bases
        
        Mutation realizations are sampled from binomial distribution with 
        probability of success equals to mutation rate.
        This mutation have frequency 1 because they are already present in root cell.
        Total mutations are computed as the sum of prior mutations number and mutation counter.
        Mutation rate mu is than estimated using the "standard" formula:
        
        P = (# total mutations )/(# bases)
        
        mu = - log( 1 - P ) / imax
        
        Equation 3, page 13. 
        
        We also compute mut rate applying a threshold on max mutation frequency.
        As described in section 3.3.1, page 14.
        """
        
        self._N0gen = N0gen
    
        # Compute mutations prior to cloning
        attempts = self._N0gen * self._bases # num attempt as in the original code
        self._n_mut_prior = np.random.binomial(attempts, self._mu)
        
        # load mutation vector with prior mutations (frequency = 1)
        for _ in np.arange(self._n_mut_prior) :
            self._mutation_frequencies.append(1.)
        
        # test LD estimator:
        # total mutations are computed as the sum of prior mutations number and mutation counter
        self._mu_est_N0gen = - math.log(1 - (self._count_mutation + self._n_mut_prior)/self._bases)
        self._mu_est_N0gen /= self._imax
        
        if not self._quiet:
            print("Estimate with prior mutations: "+str(self._mu_est_N0gen))
            
        # test estimator with threshold
        count = 0 # count mutations that are fine with threshold
        
        # simulated thresholding: frequency < 0.95
        for fobs in self._mutation_frequencies:
            if fobs <= 0.95 :
                count += 1
                
        self._mu_est_N0gen_threshold = - math.log(1 - count/self._bases)
        self._mu_est_N0gen_threshold /= self._imax
        
    # role of sequencing errors
    def test_LD_estimator_seq_err(self, seq_err = 0., threshold = 1/32.):
        
        """ 
        Test the estimator in presence of errors due to sequencing.
        As described in section 3.3.2, page 15. 
        """
        
        self._seq_err = seq_err
        self._threshold = threshold
        
        # sampling due to coverage
        # sequencing error epsilon ~ 0.1% per base
        # proposition: could both be modeled as Binomials
        
        # Errors on mutations
        count = 0 # count mutations that are fine with thresholds
        count_no_threshold = 0 # count mutations
        reads = list([]) # estimated reads
        mutobs = list([]) # empirical mutation frequency vector
        
        for frequency in self._mutation_frequencies:
            
            # sampling of coverage sequences
            R1 = np.random.binomial(self._cov, frequency)
            
            # seq error leading to false negatives
            R2 = np.random.binomial(R1, self._seq_err)
            
            #  false positive
            R3 = np.random.binomial(self._cov - R1, self._seq_err)
            
            # estimated reads
            reads.append(R1-R2+R3)
            
            # empirical mutation frequency vector
            fobs = float((R1-R2+R3)/self._cov)
            mutobs.append(fobs)
                
            # simulated thresholding
            if( (fobs >= (self._threshold)) & (fobs<=0.95) ):
                count += 1
            
            # no threshold comparison
            if ( fobs > 0. ) :
                count_no_threshold += 1
        
        if not self._quiet:
            print("Mutation count: " + str(count)) #debug
            print("Mutation count (no threshold): " + str(count_no_threshold)) #debug
            
        # errors on non-mutated bases (frequency is zero)
        
        # these errors give me a false positve if the fraction of reads
        # with a mutation is higher than the threshold (1/extant)

        R1 = self._cov * self._threshold

        # CDF Binomial of SciPy
        err_bias = 1 - binom.cdf(R1+1, self._cov, self._seq_err)
        
        
        # now I assume that it is binomial
        # we need to aclude the mutate bases
        count += np.random.binomial( self._bases - count, err_bias )
        count_no_threshold += np.random.binomial( self._bases - count, err_bias )
        
        # Test LD estimator
        
        self._mu_est_seq_err_threshold = - math.log(1 - (count)/self._bases)
        self._mu_est_seq_err_threshold /= self._imax
        
        self._mu_est_seq_err = - math.log(1 - (count_no_threshold)/self._bases)
        self._mu_est_seq_err /= self._imax
        
        if not self._quiet:
            print("Estimate with sequencing errors and threshold: "+str(self._mu_est_seq_err_threshold))
            print("Estimate with sequencing errors and no threshold: "+str(self._mu_est_seq_err))
            
    # estimate the role of ploidy
    def test_LD_estimator_ploidy(self, mean_ploidy, threshold = 1./32):
        """ 
        Test the estimator in presence of ploidy.
        As described in section 3.3.3, page 16. 
        """
        
        # mean ploidy for the simulation
        self._mean_ploidy = mean_ploidy
        
        self._threshold = threshold
        
        # ploidy is assumed to be distributed according to
        # a Poisson distribution
        # usually ploidy is in [1,5] so we model
        # ploidy ~ 1 + Poisson (mean_ploidy - 1)
        
        lam = mean_ploidy - 1 
        
        mutobs = list([]) # empirical mutation frequency vector
        count = 0 # count mutations 
        count_with_correction = 0 # count mutations considering mean ploidy
        
        for frequency in self._mutation_frequencies:
            
            #extract random ploidy
            ploidy = 1 + np.random.poisson(lam=lam)
            
            # sampling of coverage sequences
            R1 = np.random.binomial(self._cov/ploidy, frequency)
            
            #here we ignore sequencing errors
            
            # empirical mutation frequency vector
            fobs = float(R1/self._cov)
            mutobs.append(fobs)
            
            # simulated thresholding
            if( (fobs>=(self._threshold)) & (fobs<=0.95) ):
                count += 1
                    
            # simulated thresholding with correction
            # fobs = R1/(coverage/ploidy) = R1*ploidy/coverage = fobs*ploidy
            fobs *= ploidy
            
            if( fobs >= self._threshold ):
                count_with_correction += 1

                    
        # Test LD estimator
        self._mu_est_ploidy = - math.log(1 - (count)/self._bases)
        self._mu_est_ploidy /= self._imax
        
        self._mu_est_ploidy_with_correction = - math.log(1 - (count_with_correction)/self._bases)
        self._mu_est_ploidy_with_correction /= self._imax
        
        if not self._quiet:
            print("\n\n")
            print("Mutation count with ploidy: " + str(count)) #debug
            print("Mutation count with ploidy and correction: " + str(count_with_correction)) #debug
            print("Estimate with ploidy: " + str(self._mu_est_ploidy))
            print("Estimate with ploidy and correction: "+str(self._mu_est_ploidy_with_correction))
            
    # estimate mut rate counting only edges that have alive progeny
    def test_LD_estimator_no_dead(self):
        """
        Corrects the "stadard" mutation rate estimator taking into account 
        cell's death. 
        
        In the "standard" estimator the number of attempts (Poisson Process)
        are considered to be equal the total number of cells in the whole tree.
        (max cell's index, imax)
        
        Here we have 'count_with_alive_progeny':
        
                         P = (# total mutations )/(# bases)
        
                       mu = - log(1-P)/(count_with_alive_progeny)
        
        'count_with_alive_progeny' is the number of cells (node) in the whole tree that have 
        alive progeny (i.e. the number of cells in the whole tree that have descendants in 
        extant cells).
        
        Note that every cells is considered to be part of his own progeny, so extant cells 
        are counted in 'count_with_alive_progeny'.
        
        See section 3.2.4.
        """
        
        #Call the function that update 'self._count_with_alive_progeny'
        self.count_with_alive_progeny()
        
        # here we ignore error models
        
        # Test "corrected" LD estimator for death prob
        self._mu_est_no_dead = - math.log(1 - self._count_mutation/self._bases)
        self._mu_est_no_dead /= self._count_with_alive_progeny
        
        if not self._quiet:
            print("Clean mutation rate estimated (no dead progeny): "+str(self._mu_est_no_dead))  
            

    def test_LD_estimator_no_dead_estimated(self):
        """
        Test the analytic correction given by equation 16 in which the 
        denominator is given by equation 15 and 12.
        
        See section 3.5.3 of the continuous time model.
        """
        
        # in the original code are compared the clean estimate
        # and the one counting only edges that have alive progeny
        # so here we ignore error models
        
        # "age ove the tree" = number of generations
        t = self._num_layers 
        
        # we count the cumulative mean number of cell
        # with at least one alive daughter
        # from generation zero to the second-last one
        
        # Equation 12, and N_0 = 1
        count_with_alive_progeny = (1-math.pow(self._death_prob,2))*1. * (math.pow(2*(1.-self._death_prob),t-1) -1) / (math.log(2*(1-self._death_prob)))
        
        # we add also add extant cells number
        count_with_alive_progeny += self._extant
        
        # Test analytic LD estimator corrected for death prob
        self._mu_est_no_dead_estimated = - math.log(1 - self._count_mutation/self._bases)
        self._mu_est_no_dead_estimated /= count_with_alive_progeny
        
        if not self._quiet:
            print("Clean mutation rate estimated (no dead progeny estimate): "+str(self._mu_est_no_dead_estimated))  
    

    def test_LD_estimator_all_errors(self, N0gen, seq_err, mean_ploidy, threshold):
        """ 
        Estimates mut rate accounting for all error models.
        """
        
        self._mean_ploidy = mean_ploidy
        self._N0gen = N0gen
        self._threshold = threshold
        
        # mutation count
        count = 0
        
        # Mutations prior to cloning
    
        # Compute the number of mutations prior to cloning
        attempts = self._N0gen * self._bases # num attempt as in the original code
        self._n_mut_prior = np.random.binomial(attempts, self._mu)
        
        # load mutation vector with prior mutations (frequency = 1)
        # subclonal mutations are loaded when tree is generated
        for _ in np.arange(self._n_mut_prior) :
            self._mutation_frequencies.append(1.)
            
        
        #apply seq error model and mean ploidy
        
        lam = mean_ploidy - 1 #mean for poisson distribution of ploidy
        
        mutobs = [] # observed mut frequencies
        
        for frequency in self._mutation_frequencies:
            
            # extract random ploidy
            ploidy = 1 + np.random.poisson(lam=lam)
            
            # sampling of coverage sequences
            R1 = np.random.binomial(self._cov/ploidy, frequency)
            
            # seq error leading to false negatives
            R2 = np.random.binomial(R1, self._seq_err)
            #  false positive
            R3 = np.random.binomial(self._cov/ploidy - R1, self._seq_err)
            
            reads = R1 - R2 + R3 
            fobs = reads*ploidy/self._cov
            mutobs.append(fobs)
            
            # check with tresholds
            
            # double threshold:
            # fobs < 0.95 for mutation prior to cloning
            # fobs > self._threshold for seq errors
            
            if (fobs<= 0.95) & (fobs >= self._threshold) :
                count += 1
            
        # errors on non mutated bases
        # these errors give me a false positve if the fraction of reads
        # with a mutation is higher than the threshold (1/extant)

        R1 = self._cov*self._threshold
            
        err_bias = 1 - binom.cdf(R1+1, self._cov, self._seq_err)
        
        # now I assume that it is binomial
        # we need to aclude the mutate bases
        count += np.random.binomial( self._bases - count, err_bias )
            
        # Test LD estimator
        self._mu_est_all_errors = - math.log(1 - (count)/self._bases)
        # update the counter with the number of cells with alive progeny
        self.count_with_alive_progeny()
        
        self._mu_est_all_errors /= self._count_with_alive_progeny
            
    #----------------------------
    # UTILTIES
    #----------------------------
    
    def max_imax(self, gen):
        """
        Return max possible cell index given the number of generations. 
        """
        imax = 0
        for n in range(gen):
            imax += 2*2**(n)
        return imax
    
    def get_last_layer(self):
        """
        Return list of nodes (cell) in the last layer
        """
        # last layer is the list of extant nodes
        
        # return list of elements in population list 
        # population = [ Node(0), Node(1), ..., Node(min_extant_id), ..., Node(imax), None, ..., None]
        # population = [ root cell, ..., ..., ..., first extant, ..., ..., last extant, None, ..., None]
        
        #return 
        #[first extant, ..., ..., last extant]
        return self._population[self._min_extant_id: self._imax+1]
        
    def print_tree(self):
        """
        Print "human readable" generational tree.
        """
        for pre, _, node in RenderTree(self._population[0]):
            print("%s%s" % (pre, node.name))
            
    def count_alive_progeny(self, cellID):
        """
        Check how many alive mutants there are in the progeny of cell with specific ID.
        """

        # There are two cases:
        # 1. mutant cell is an extant cell -> no progeny other than itself
        # 2. mutant cell is not an extant cell -> check for alive descendants
        
        if cellID >= self._min_extant_id:
            # Case 1. mutant cell is an extant cell
            # than the only alive mutant is itself
            return 1
        else:
            # Case 2. mutant cell is not an extant cell 
            # than check for alive descendants:
            # get the list of ID of all descendats of cellID
            descendant_IDs = [node.name for node in self._population[cellID].descendants]
            #initialize mutant counter
            count_alive_mutants = 0
            #for each of them
            for descendantID in descendant_IDs:
                #if descendantID is bigger than or equal to 'min_extant_ID' it's an extant cell
                if descendantID >= self._min_extant_id:
                    # count it as an alive mutant
                    count_alive_mutants += 1
            #return counter
            return count_alive_mutants
        
    def count_with_alive_progeny(self, node_list = None):
        """
        Function to update 'count_with_alive_progeny':
        
        'count_with_alive_progeny' is the number of cells (node) in the whole tree that have 
        alive progeny (i.e. the number of cells in the whole tree that have descendants in 
        extant cells).  
        
        Note that every cells is considered to be part of its own progeny, so extant cells 
        are counted in 'count_with_alive_progeny'.
        
        This function operates recursively (one call for each layer) starting from last layer.
        """
        
        # If no list is passed
        if node_list == None :
            # This is the first call to this function 
            self._count_with_alive_progeny = 0
            # Operates on last layer of the tree
            
            # Take the list of nodes (cell) in the last layer
            node_list = self.get_last_layer()
        
        # If node_list is only the cell zero we have finished
        # Don't count cell zero
        # ("node.name" to compare int type)
        if node_list[0].name == 0:
            # End
            return 
            
        ancestors_list = []
        
        for i in np.arange(np.shape(node_list)[0]):
            #for every element of the node_list
            #i is the index of the element in the node_list
            
            # Count this node
            self._count_with_alive_progeny += 1
            
            # Fill ancestors list
            
            # If it is the first element of the layer
            if i == 0:
                # Add parent to ancestors
                ancestors_list.append(node_list[i].parent)
            else:
                # For other elements (i > 0):
                
                # Check if the parent cell is the same 
                # as the previous node
                
                # ('node.parent.name' to compare int type)
                if node_list[i].parent.name == node_list[i-1].parent.name :
                    # if it's the same don't do anything
                    # because parent cell is already in ancestors list
                    pass
                
                else:
                    # different parent from preavious cell
                    # Add parent to ancestors
                    ancestors_list.append(node_list[i].parent)
                    
        # Continue (recursively) for other layers until reach cell zero (root)
        self.count_with_alive_progeny(node_list = ancestors_list)
    
