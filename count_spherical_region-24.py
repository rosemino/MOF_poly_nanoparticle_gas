# -*- coding: utf-8 -*-
"""

@author: Cecilia
"""

import numpy as np

# This is the number of atoms in the simulation domain (mof, polymer and gas) in
# the given ZIF-8/PVDF/(single-gas) system.
NA = 69769 

# Number of configuratios
number_of_configurations = 20

# This part of the code is responsible for reading the configurations of a LAMMPS
# trajectory, which was collected during the dynamics. The file contain the follow-
# ing attributes: ATOMS id type x y z vx vy vz
trj = np.zeros((1,5))
box_size = np.zeros((1,2))
ofi = open("dump.pos", 'r')
for it_1 in range (0, number_of_configurations):
    add=[]
    tmp = []
    ofi.readline()
    ofi.readline()
    ofi.readline()
    ofi.readline()
    ofi.readline()
    for it_3 in range (0,3):
        line = ofi.readline()
        for f in zip(line.split(' ')):
            # e = list(e)*
            f = np.array(f)
            tmp = np.append(tmp,float(f))
    tmp = np.array(tmp,float)
    tmp = tmp.reshape(3,2)
    box_size = np.concatenate((box_size, tmp), axis = 0)
    ofi.readline()
    for it_2 in range(0, (NA)):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e in zip(dump.split(' ')):
            # e = list(e)*
            e = np.array(e)
            add = np.append(add,float(e))
    add = np.array(add,float)
    add = add.reshape(NA,8)
    # -------------------------------
    # I dont need the velocities, despite having saved them, and so I will delete
    # it.
    add = np.delete(add, 7, 1)
    add = np.delete(add, 6, 1)
    add = np.delete(add, 5, 1)
    # -------------------------------
    trj = np.concatenate((trj, add), axis = 0)

# Lets delete the first line, which was defined as a bunch of columns containing
# the "0" element.
trj = np.delete(trj,0, axis=0)
box_size = np.delete(box_size, 0, 0)

# -------------------------------------------------------------------
# This is the number of MOF chemical bonds. You need to alter accordingly depending
# on which ZIF-8 nanoparticle (CB, RD or SRD) you have embedded in your given system.
number_of_bonds = 6900 

bonds_mof = []
ofi = open("bonds_mof.dat", 'r')
for it_1 in range(0, number_of_bonds):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e,it_2 in zip(dump.split('\t'), range(4)):
            bonds_mof.append(float(e))
bonds_mof = np.array(bonds_mof,float)
bonds_mof = bonds_mof.reshape(number_of_bonds,4)

# --------------------------------------------------------
gas_molecules = 0
zif8_beads1 = 0
zif8_beads2 = 0

# These are the masses of the ZIF-8 beads types 1 and 2
M1 = 65.38
M2 = 81.0

number_of_configurations_considered = 20

for loop_0 in range (0, number_of_configurations_considered):
    
    output = np.zeros((NA,5))
    output[:,:] = trj[int((loop_0)*NA):int((loop_0 + 1)*NA),:]
    
    # This defines the vectors i, j and k for the given configuration using
    # the simulation domain.
    i = np.zeros((1,3))
    i[0,0] = box_size[int(loop_0*3),1] - box_size[int(loop_0*3),0]
    mod_i = ((i[0,0])**2 + (i[0,1])**2 + (i[0,2])**2)**(1/2)
    j = np.zeros((1,3))
    j[0,1] = box_size[int((loop_0*3) + 1),1] - box_size[int((loop_0*3) + 1),0]
    mod_j = ((j[0,0])**2 + (j[0,1])**2 + (j[0,2])**2)**(1/2)
    k = np.zeros((1,3))
    k[0,2] = box_size[int((loop_0*3) + 2),1] - box_size[int((loop_0*3) + 2),0]
    mod_k = ((k[0,0])**2 + (k[0,1])**2 + (k[0,2])**2)**(1/2)
    
    # Since the values of xlo, ylo and zlo are not 0, I will need this in order
    # to have the coordinates mapped into the simulation domain.
    output[:,2] = (output[:,2]) + box_size[int(loop_0*3),0]
    output[:,3] = (output[:,3]) + box_size[int(loop_0*3) + 1,0]
    output[:,4] = (output[:,4]) + box_size[int(loop_0*3) + 2,0]
    
    x_min = box_size[int(loop_0*3),0]
    y_min = box_size[int(loop_0*3) + 1,0]
    z_min = box_size[int(loop_0*3) + 2,0]
    
    # I think the first thing I will do is to shift the (xlo, ylo, zlo) of the simulation box
    # to the origin to facilitate further reasoning and ensure, within my mindstate at the
    # moment that I am absolutely doing things correctly!
    translator = np.array([[x_min, y_min, z_min]])
    output[:,2] = output[:,2] - translator[0,0]
    output[:,3] = output[:,3] - translator[0,1]
    output[:,4] = output[:,4] - translator[0,2]
    
    # I will also to redefine the values of x_min, y_min and z_min defined previously.
    x_min = 0
    y_min = 0
    z_min = 0
    
    # -----------------------------------------------------------------------------------
    # Now I will try to find out which part of the simulation domain (lx, ly and lz are
    # divided into 2) contains most of my nanoparticle. This is because I did not ensure
    # it is always in the center of the simulation domain (easy to do within LAMMPS).
    counter_xi = 0
    counter_yi = 0
    counter_zi = 0
    counter_xf = 0
    counter_yf = 0
    counter_zf = 0
    
    for it_1 in range (0, len(output)):
        # I will only run all of this stuff if I found a superatom of the MOF phase
        if (output[it_1,1] == 1) or (output[it_1,1] == 2):
            position = np.array([[output[it_1,2], output[it_1,3], output[it_1,4]]])
            # ---------------------------------------
            if position[0,0] > i[0,0]/2:
                counter_xf = counter_xf + 1
            if position[0,0] < i[0,0]/2:
                counter_xi = counter_xi + 1
            # ---------------------------------------
            if position[0,1] > j[0,1]/2:
                counter_yf = counter_yf + 1
            if position[0,1] < j[0,1]/2:
                counter_yi = counter_yi + 1
            # ---------------------------------------
            if position[0,2] > k[0,2]/2:
                counter_zf = counter_zf + 1
            if position[0,2] < k[0,2]/2:
                counter_zi = counter_zi + 1
        
    # By the end of this, I will be able to identify, for each of the three directions,
    # in which part of the simulation box most of the atoms are. 
    
    # I will then have to shift the atoms that sit at the portion containing the least a-
    # mount of atoms AND which has been disrupted by periodic boundaries in order to have
    # the nanoparticle completely unwrapped.
    flag = 1
    while flag == 1:
        flag = 0
        # flag will only remain as 0 if, by the end of the it_1 for loop below, I never happen
        # to need to shift atoms due to crossing of the periodic boundary. Note that I am al-
        # ways shifting the superatoms in the direction of the more populated region, so that
        # I spend less time here.
        for it_1 in range (0, len(bonds_mof)):
            given_CG1 = bonds_mof[it_1,2]
            given_CG2 = bonds_mof[it_1,3]
            counter = 0
            # -----------------------------------------------------------------
            # Finding the positions. The variable counter is meant to speed up things
            # in the sense that I am creating a condition to break this for loop so  
            # that I dont need to scan it up to the end to find the position of the
            # two given superatoms.
            for it_2 in range (0, len(output)):
                if output[it_2,0] == given_CG1:
                    given_CG1_pos = np.array([[output[it_2,2], output[it_2,3], output[it_2,4]]])
                    line_CG1 = it_2
                    counter = counter + 1
                if output[it_2,0] == given_CG2:
                    given_CG2_pos = np.array([[output[it_2,2], output[it_2,3], output[it_2,4]]])
                    line_CG2 = it_2
                    counter = counter + 1
                if counter == 2:
                    break
            # -----------------------------------------------------------------
            dist_x = given_CG1_pos[0,0] - given_CG2_pos[0,0]
            dist_y = given_CG1_pos[0,1] - given_CG2_pos[0,1]
            dist_z = given_CG1_pos[0,2] - given_CG2_pos[0,2]
            
            # if the condition below is met, this means that I should be worrying about
            # translating superatoms downwards since the xi side is more populated. The
            # equal condition is included here but could be included in the other also.
            if int(counter_xi/counter_xf) >= 1:
                if (abs(dist_x) > mod_i/2) & (dist_x > 0):
                    output[line_CG1,2] = output[line_CG1,2] - i[0,0]
                    flag = 1
                if (abs(dist_x) > mod_i/2) & (dist_x < 0):
                    output[line_CG2,2] = output[line_CG2,2] - i[0,0]
                    flag = 1
            if int(counter_xi/counter_xf) < 1:
                if (abs(dist_x) > mod_i/2) & (dist_x > 0):
                    output[line_CG2,2] = output[line_CG2,2] + i[0,0]
                    flag = 1
                if (abs(dist_x) > mod_i/2) & (dist_x < 0):
                    output[line_CG1,2] = output[line_CG1,2] + i[0,0]
                    flag = 1
            # ---------------------------------------------------
            # Similarly, to the y direction:
            if int(counter_yi/counter_yf) >= 1:
                if (abs(dist_y) > mod_j/2) & (dist_y > 0):
                    output[line_CG1,3] = output[line_CG1,3] - j[0,1]
                    flag = 1
                if (abs(dist_y) > mod_j/2) & (dist_y < 0):
                    output[line_CG2,3] = output[line_CG2,3] - j[0,1]
                    flag = 1
            if int(counter_yi/counter_yf) < 1:
                if (abs(dist_y) > mod_j/2) & (dist_y > 0):
                    output[line_CG2,3] = output[line_CG2,3] + j[0,1]
                    flag = 1
                if (abs(dist_y) > mod_j/2) & (dist_y < 0):
                    output[line_CG1,3] = output[line_CG1,3] + j[0,1]
                    flag = 1
            # ---------------------------------------------------
            # Similarly, to the z direction:
            if int(counter_zi/counter_zf) >= 1:
                if (abs(dist_z) > mod_k/2) & (dist_z > 0):
                    output[line_CG1,4] = output[line_CG1,4] - k[0,2]
                    flag = 1
                if (abs(dist_z) > mod_k/2) & (dist_z < 0):
                    output[line_CG2,4] = output[line_CG2,4] - k[0,2]
                    flag = 1
            if int(counter_zi/counter_zf) < 1:
                if (abs(dist_z) > mod_k/2) & (dist_z > 0):
                    output[line_CG2,4] = output[line_CG2,4] + k[0,2]
                    flag = 1
                if (abs(dist_z) > mod_k/2) & (dist_z < 0):
                    output[line_CG1,4] = output[line_CG1,4] + k[0,2]
                    flag = 1
                    
    # -----------------------------------------------------------------------------------
    # Once I finish the part below I can finally calculate the com of the mof without worrying
    # that it is broken across boundaries.
    numerator = np.array([[0.0, 0.0, 0.0]])
    denominator = 0
    for it_1 in range (0, len(output)):
        if (output[it_1,1] == 1):
            position = np.array([[output[it_1,2], output[it_1,3], output[it_1,4]]])
            numerator = numerator + M1*position
            denominator = denominator + M1
        if (output[it_1,1] == 2):
            position = np.array([[output[it_1,2], output[it_1,3], output[it_1,4]]])
            numerator = numerator + M2*position
            denominator = denominator + M2
    
    com_mof = numerator/denominator
    
    # -------------------------------------------------------------------------------
    # Now I will also already re-position the polymer superatoms so that they are
    # in the simulation box the closest possible to the center of mass of the MOF,
    # which could be a ghost domain or the original one.
    for it_1 in range (0, NA):
        if output[it_1,1] > 2:
            atom_position = np.array([[output[it_1,2], output[it_1,3], output[it_1,4]]])
            dist_x = atom_position[0,0] - com_mof[0,0]
            dist_y = atom_position[0,1] - com_mof[0,1]
            dist_z = atom_position[0,2] - com_mof[0,2]
            # ------------------------------------------------------------------
            if (abs(dist_x) > mod_i/2) & (dist_x > 0):
                atom_position[0,0] = atom_position[0,0] - i[0,0]
                output[it_1,2] = output[it_1,2] - i[0,0]
            if (abs(dist_x) > mod_i/2) & (dist_x < 0):
                atom_position[0,0] = atom_position[0,0] + i[0,0]
                output[it_1,2] = output[it_1,2] + i[0,0]
            if (abs(dist_y) > mod_j/2) & (dist_y > 0):
                atom_position[0,1] = atom_position[0,1] - j[0,1]
                output[it_1,3] = output[it_1,3] - j[0,1]
            if (abs(dist_y) > mod_j/2) & (dist_y < 0):
                atom_position[0,1] = atom_position[0,1] + j[0,1]
                output[it_1,3] = output[it_1,3] + j[0,1]
            if (abs(dist_z) > mod_k/2) & (dist_z > 0):
                atom_position[0,2] = atom_position[0,2] - k[0,2]
                output[it_1,4] = output[it_1,4] - k[0,2]
            if (abs(dist_z) > mod_k/2) & (dist_z < 0):
                atom_position[0,2] = atom_position[0,2] + k[0,2]
                output[it_1,4] = output[it_1,4] + k[0,2]
    
    # ---------------------------------------------------------------------------------
    # Once I get here, I should have the nanoparticle unwrapped, the coordinates of its
    # com dully calculated at stored in the com_position array *and* I will have moved
    # all the polymer chains to the simulation domain that sits the closest to the COM
    # of the nanoparticle.
    # Now I need to simply count all the gas beads as well as MOF beads that lie within a
    # sphere of 24 angs of the COM of the nanoparticle.
    # Note that you can set whatever value as radius if you wish.
    for it_1 in range (0, NA):
        if output[it_1,1] == 7:
            atom_position = np.array([[output[it_1,2], output[it_1,3], output[it_1,4]]])
            dist_x = com_mof[0,0] - atom_position[0,0]
            dist_y = com_mof[0,1] - atom_position[0,1]
            dist_z = com_mof[0,2] - atom_position[0,2]
            dist_total = (dist_x**2 + dist_y**2 + dist_z**2)**(1/2)
            if dist_total <= 24:
                gas_molecules = gas_molecules + 1
        # ------------------------------------------------------------
        if output[it_1,1] == 1:
            atom_position = np.array([[output[it_1,2], output[it_1,3], output[it_1,4]]])
            dist_x = com_mof[0,0] - atom_position[0,0]
            dist_y = com_mof[0,1] - atom_position[0,1]
            dist_z = com_mof[0,2] - atom_position[0,2]
            dist_total = (dist_x**2 + dist_y**2 + dist_z**2)**(1/2)
            if dist_total <= 24:
                zif8_beads1 = zif8_beads1 + 1
        # ------------------------------------------------------------
        if output[it_1,1] == 2:
            atom_position = np.array([[output[it_1,2], output[it_1,3], output[it_1,4]]])
            dist_x = com_mof[0,0] - atom_position[0,0]
            dist_y = com_mof[0,1] - atom_position[0,1]
            dist_z = com_mof[0,2] - atom_position[0,2]
            dist_total = (dist_x**2 + dist_y**2 + dist_z**2)**(1/2)
            if dist_total <= 24:
                zif8_beads2 = zif8_beads2 + 1
            
gas_molecules = gas_molecules/number_of_configurations_considered
zif8_beads1 = zif8_beads1/number_of_configurations_considered
zif8_beads2 = zif8_beads2/number_of_configurations_considered
           
vector = np.array([[zif8_beads1, zif8_beads2, gas_molecules]])
ending = np.savetxt('sphere-24.txt', (vector))

# --------------------------------------------------------------------------------------
# Provided that you do this for different sets of 20 configurations in the trajectory
# as we did in our work, you will ultimately have multiple output files as you run the
# code using different dump.pos files as input. You can use them to do some statistics.