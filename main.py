from material import Material
from element import Element
import numpy as np
import os
import matplotlib.pyplot as plt

# Declaration of All Parameters

#  ???????????????????
weight = 2  # Since gauss point = 1, the weight = 2

f_total = 5100  # Magnitude of Maximum Force (Newton)
vis=1  # Viscosity of the material model (per second)
max_time=1  # Total time for which force is applied (seconds)
number_elements = 2  # Number of Elements the user wants the problem to be discretized 
limit=240  # Yield Stress in MPa (Elastic Limit)
limit_array = limit*np.ones([number_elements,1]) # Array in order to compute for 'number_elements' number of elements
E=Ct=120000*np.ones([number_elements,1])  # Young's Modulus and Tangent Stiffness (GPa), which are initially equal for Plastic Condtion 
steps = 6000 # No. of Steps user wants to focus_nodeide the time and force

del_t = 1/steps
force_increment = np.linspace(0,f_total,steps)
# This is the constructor for Element class in the Element Routine
# The arguements passed are the lengths, cross sectional areas of two elements and the weight
# this is used to get /update the elemental assignement matrix


# NON-LINEAR FINITE ELEMENT ANALYSIS oF SOLIDS AND STRUCTURES
# René de Borst, Mike A. Crisfield

our_element = Element(number_elements, 10,40, 20, 80, weight)
epsilon = np.zeros(number_elements)  #Initialization of strain values for all elements
our_material = Material(number_elements,epsilon, E, limit_array , del_t, vis)  # This is the constructor for Material class in the material routine
length, area = our_element.parameters()  # Calling Element Routine
our_element.B_matrix()
our_element.J_matrix()
# To append values of Stress & Strain
stress_list =[]
save_strain =[]
# Initializing stress arrays
sigma_prev = np.zeros([number_elements,1]) 
sigma_next = np.zeros([number_elements,1]) 
initial_eps_pl_k = np.zeros([number_elements,1])  # Initializing initial plastic strain to '0'
initial_u_sol = 0.0  # Initializing initial displacement of the applied node
# Initializing internal and external force at global level
internal_force = np.zeros([number_elements+1,1])
external_force = np.zeros([number_elements+1,1])
u_glob =  np.zeros([number_elements+1,1])  # Initialzing Global Displacement Matrix
u_elements = np.zeros([number_elements , 2, 1])  # Initializing Element Displacement Matrix (as a 3D Matrix)

u_prev = np.zeros([number_elements,1])
eps_next = 0
total_disp= []
t_disp = "disp.txt"
displacement_data= [t_disp]
try:
    os.remove("disp.txt")
except:
    pass

for f in force_increment: 
    print("Value of Force : ", f)
    # Calling External Force from Element Routine
    external_force = our_element.external_force(f)
    # Initializing Residual
    G=np.zeros([number_elements+1,1])
    
    #used to calculate the final dispalcement occured due to external loading
    #in this case due to incremental loading, we use this in every iteration to minmiza the delta(F)
    
    for Newton_Raphson_Iteration in range(6): #Since Newton-Raphosn will converge within this Criterion
        K_glob = np.zeros([number_elements+1,number_elements+1])
        if f == 0: 
            sigma = 0
        K_ele = np.array(our_element.k_ele(Ct))
        internal_force_e = np.array(our_element.internal_forceernal(sigma_next-sigma_prev))
        for ele in range(number_elements):
            A = our_element.A_matrix(ele+1)
            internal_force = internal_force + np.matmul(A.T,internal_force_e[:,ele])
            K_glob = K_glob + np.matmul( np.matmul(A.T,K_ele),A) 
            
        focus_node = our_element.focus_node_return()
        K_glob_red = K_glob[1:focus_node+1,1:focus_node+1]
        K_glob_red_inv = np.linalg.inv(K_glob_red)
        G = external_force -  internal_force
        G_red = G[1:focus_node+1]
        # del_u = G_red/K_glob_red 
        del_u  = np.matmul(K_glob_red_inv, G_red)
        u_next = u_prev + del_u
        u_glob[1:focus_node+1] = del_u
        
        # print(u_glob)
        for ele in range(number_elements):
            A = our_element.A_matrix(ele+1)
            u_elements[ele,:,:,] = u_elements[ele,:,:,] + np.matmul(A,u_glob)
        strain_check = our_element.strain(u_elements)
        sigma, epsilon, Ct, eps_next , status = our_material.material_condition(initial_eps_pl_k,strain_check, f)
        sigma_prev = sigma_next
        sigma_next = sigma
        error = del_u
        u_prev = u_next
    initial_eps_pl_k =eps_next
    stress_list.append(sigma[0])
    save_strain.append(epsilon[0])
    total_disp.append(u_next)
    plot_steps = np.linspace(0,1,steps)
    print(Ct)
    
    with open(t_disp, 'a') as fd:
        x_write = u_next[-1]
        fd.write("\n t_disp %s %s " %(x_write,status))
        fd.close()

stress_list=np.array(stress_list)
save_strain=np.array(save_strain)

total_disp = np.array(total_disp)
total_disp = total_disp[:,1:2,:].flatten()
# print(total_disp.flatten())

fig, axs = plt.subplots(1,2,figsize=(25,25))
axs[0].set(xlabel = 'Strain')
axs[0].set(ylabel = 'Stress')
axs[1].set(xlabel = 'Time Step')
axs[1].set(ylabel = 'Displacement of Node')
axs[0].plot(save_strain,stress_list)
axs[1].plot(plot_steps ,total_disp)
plt.show()
