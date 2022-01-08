import numpy as np

class Material (object): 
    def __init__(self,number_elements, epsilon, E, limit_array , dt, vis ): 
        self.epsilon = epsilon
        self.E = E
        self.limit_array = limit_array
        self.dt = dt
        self.vis = vis
        self.h =  self.dt / 50
        self.iter = 100
        self.Ct = 0
        self.number_elements = number_elements
        self.sigma = np.zeros([self.number_elements,1]) 
        self.eps_pl_k_prev = np.zeros([self.number_elements,1]) 
    
    # sigma =YM*strain
    def sigma_fun (self,strain_check):  
        self.epsilon = strain_check
        self.sigma = self.E * self.epsilon

        return self.sigma 
    
    #it is used to check if any element has reached plastic stage
    # if the stess of element is less than Yield stress, in elastic zone
    # as soon as the element stress exceeeds the Yield stress, we calculatlate the plastic strain 
    # etrue-eyield
    
    def material_condition(self, eps_pl_k, strain_check, f):
        self.epsilon = strain_check

        if abs(self.sigma[0]) < self.limit_array[0] and abs(self.sigma[1]) < self.limit_array[1]: 
            self.sigma = self.E * self.epsilon

            return self.sigma, self.epsilon, self.E, np.zeros([self.number_elements,1]) , "Elastic Condition in Both Rods"


        elif abs(self.sigma[0]) > self.limit_array[0] and abs(self.sigma[1]) < self.limit_array[1]: 

            # why is this used for?
            #ans: 
            for _ in range(int(10000)):
                self.sigma = self.E * ( self.epsilon - eps_pl_k)
                self.sigma[1] = self.E[1] * self.epsilon[1]
                eps_pl_k_next = eps_pl_k + self.h * self.vis * np.sign(self.sigma) * (np.abs(self.sigma/self.limit_array)-1)

                while ( abs(self.eps_pl_k_prev[0]- eps_pl_k_next[0]) >1e-5):
                    self.eps_pl_k_prev =  eps_pl_k_next
                    self.sigma =  self.E * ( self.epsilon - eps_pl_k_next)
                    self.sigma[1] = self.E[1] * self.epsilon[1]
                    eps_pl_k_next = eps_pl_k + self.h * self.vis * np.sign(self.sigma) * (np.abs(self.sigma/self.limit_array)-1)

                self.sigma =  self.E * ( self.epsilon - eps_pl_k_next) 
                self.sigma[1] = self.E[1] * self.epsilon[1]
                eps_pl_k =  eps_pl_k_next

            self.Ct = self.E / (1 + ( self.E * self.dt * self.vis) / self.limit_array)
            self.Ct[1] = self.E[1]
            eps_pl_k[1]  = eps_pl_k_next[1]= 0.0

            return self.sigma, self.epsilon, self.Ct , eps_pl_k_next ,"Rod 1: Plastic, Rod 2: Elastic"

        elif abs(self.sigma[0]) > self.limit_array[0] and abs(self.sigma[1]) > self.limit_array[1]:   

            for _ in range(int(10000)):

                self.sigma = self.E * ( self.epsilon - eps_pl_k)
                eps_pl_k_next = eps_pl_k + self.h * self.vis * np.sign(self.sigma) * (np.abs(self.sigma/self.limit_array)-1)

                while (abs(self.eps_pl_k_prev[0]- eps_pl_k_next[0]) >1e-5 and abs(self.eps_pl_k_prev[1]- eps_pl_k_next[1]) >1e-5):

                    self.eps_pl_k_prev =  eps_pl_k_next
                    self.sigma =  self.E * ( self.epsilon - eps_pl_k_next)
                    eps_pl_k_next = eps_pl_k + self.h * self.vis * np.sign(self.sigma) * (np.abs(self.sigma/self.limit_array)-1)

                self.sigma = self.E * ( self.epsilon - eps_pl_k_next)
                eps_pl_k =  eps_pl_k_next
                self.Ct = self.E / (1 + ( self.E * self.dt * self.vis) / self.limit_array)
                
            return self.sigma, self.epsilon, self.Ct , eps_pl_k_next ,"Plastic Condition in Both Rods"