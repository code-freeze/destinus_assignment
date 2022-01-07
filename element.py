import numpy as np
import sys
class Element():   
    def __init__(self, number_elements, ar_1, len_1, ar_2, len_2, weight):
        self.number_elements = number_elements
        self.ar_1 = ar_1
        self.len_1 = len_1
        self.len_2 = len_2
        self.ar_2 = ar_2
        self.weight = weight
        self.focus_node = 0

    def element_check(self):
        
        if self.number_elements == 1:
            print ("Minimum requirement for elements not met, elements considered =  2")
            return 1
        else : 
            return 1
    
    def parameters(self):

        self.length = []
        self.area = []
        self.focus_node = int(self.number_elements/2) + 1
        if self.number_elements <= 1:
            sys.exit ()
        if self.number_elements == 2:
            self.focus_node = self.element_check()
        for i in range(self.number_elements):
            if i<self.focus_node:   
                self.length.append(self.len_1/self.focus_node)
                self.area.append(self.ar_1)
            else:
                self.length.append(self.len_2/(self.number_elements-self.focus_node))
                self.area.append(self.ar_2)
        self.length = np.array(self.length)
        self.area = np.array(self.area)

        return self.length , self.area

    
    def B_matrix(self):

        self.B = np.array([(-1/self.length),(1/self.length)])
  
    def J_matrix(self):

        self.J = self.length/2
    

    def k_ele(self, Ct):
        
        k=Ct*(self.area/self.length)* np.array([[1,-1],[-1,1]]) 
       
        return k
 
    def A_matrix(self,ele_num): 
       
        a_element=np.zeros([2,self.number_elements+1])
        a_element[0][ele_num-1]=1
        a_element[1][ele_num]=1
        
        return a_element
    
    def internal_forceernal(self, sigma): 
       
        int_force = []
        for ele in range(self.number_elements):
           
            int_force.append((self.weight*self.B[:,ele].reshape(1,self.number_elements)*self.area[ele]*self.J[ele]*sigma[ele]).reshape(self.number_elements,1))

        return int_force

    def external_force(self, f):
       
        ext_force=np.zeros(self.number_elements+1)
        focus_node_check = self.focus_node
        ext_force[focus_node_check]= f
      
        return ext_force.reshape(self.number_elements+1,1)

    def focus_node_return(self):
       
        return self.focus_node

    def strain(self, u_elements):
       
        strain_combined = []
       
        for ele in range(self.number_elements):
          
            strain_combined.append(np.dot(self.B[:,ele],u_elements[ele,:,:]))
       
        return np.array(strain_combined)

