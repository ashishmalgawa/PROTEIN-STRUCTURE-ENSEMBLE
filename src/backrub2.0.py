# Usage: Command line argument: path to pdb file,segment size, Number of iterations 
from random import randint
import random
import math
import os
import matplotlib.pyplot as plt
import subprocess
import sys
import time
# Global Variables
SCWRL_OUTPUT_DIR= "scwrloutput"
BACKRUB_OUTPUT_DIR= "backrub_output"
ENSEMBLE_OUTPUT_DIR= "ensembles"
SEGMENT_SIZE = int(sys.argv[2])
CURRENT_WORKING_DIRECTORY = os.getcwd()
TEMPERATURE = 302
def output_dir_exists():
    if not os.path.exists(sys.argv[1]):
        print "Invalid path to pdb file...Exiting now!"
        exit(0)
    if not os.path.exists("output/"):
        print "creating output directory"
        os.mkdir("output")
    if not os.path.exists(SCWRL_OUTPUT_DIR+"/"):
        print "creating scwrl output directory"
        os.mkdir(SCWRL_OUTPUT_DIR)
    if not os.path.exists(ENSEMBLE_OUTPUT_DIR+"/"):
        print "creating ensemble output directory"
        os.mkdir(ENSEMBLE_OUTPUT_DIR)
    

def scwrl(pdb_path):
    output_path=SCWRL_OUTPUT_DIR+'/scwrl_'+pdb_path.split("/")[-1]
    command='Scwrl4 -i "'+pdb_path+'" -o '+output_path
    subprocess.check_output(command, shell=True)
    return output_path

def foldx(pdb_path):
    output_dir="/".join(pdb_path.split("/")[:-1])
    os.chdir(CURRENT_WORKING_DIRECTORY+"/"+output_dir)
    command='foldx --command=Stability --pdb='+pdb_path.split("/")[-1]+' | grep "^Total          ="'
    output=subprocess.check_output(command, shell=True)
    total_energy=float(output.split("=")[1].strip())
    os.chdir(CURRENT_WORKING_DIRECTORY)
    return total_energy

def generateNumberString(num):
    string=""
    for i in range(7-len(str(num))):
        string += " "
    return string+str(num)+" "

def updatepdbline(old,new):
    newline=old[:31]
    for co_ordinate in new:
        newline +=generateNumberString(co_ordinate)
    newline +=old[55:]
    return newline

def backrub(filepath,n):
    f = open(filepath, "r")
    temp_set=set()
    cordinate_set_list=list()
    for line in f:
        temp = line.split(" ")
        if line[:3]== "TER":
            cordinate_set_list.append(temp_set)
            temp_set=set()
        if temp[0] == "ATOM":
            curr_int=int(line[23:26])
            temp_set.add(curr_int)
    no_chains=len(cordinate_set_list)
    f.close()
    #select a segment randomly whose consecutive n ATOMs' co-co_ordinate information is available
    while True:
        flag=1
        random.seed(time.time())
        apply_move_on_chain=randint(0,no_chains-1)
        cordinate_list=list(cordinate_set_list[apply_move_on_chain])
        apply_on_segment=random.choice(cordinate_list)
        for i in range(apply_on_segment,apply_on_segment+n):
            if i not in cordinate_list:
                flag=0
                break
        if flag==1:
            print"apply from segment",apply_on_segment
            break
    
    f = open(filepath, "r")
    ter_count=0
    apply_backrub_on_atom_list=list()
    for line in f:
        if ter_count != apply_move_on_chain:
            
            temp = line.split("  ")
            if line[:3] == "TER":
                ter_count+=1
        else:
            temp = line.split("  ")
            if temp[0] == "TER":
                ter_count+=1
                break
            if temp[0] == "ATOM":
                temp_int=int(line[23:26])
                if temp_int in range(apply_on_segment,apply_on_segment+n):
                    apply_backrub_on_atom_list.append(line)
                    if line[13:15]=="CA" and temp_int==apply_on_segment:
                        apply_backrub_on_atom_list=[line] #emptying list at first Ca atom
                        start_ca=line
                    if line[13:15]=="CA" and temp_int==apply_on_segment+n-1:
                        end_ca=line
                        break
    f.close()
    axis_temp=float(end_ca[30:39])
    axis_endca=list()
    axis_endca.append(axis_temp)
    axis_temp=float(end_ca[40:47])
    axis_endca.append(axis_temp)
    axis_temp=float(end_ca[47:54])
    axis_endca.append(axis_temp)

    axis_temp = float(start_ca[30:39])
    axis_startca = list()
    axis_startca.append(axis_temp)
    axis_temp = float(start_ca[40:47])
    axis_startca.append(axis_temp)
    axis_temp = float(start_ca[47:54])
    axis_startca.append(axis_temp)

    axis=list()
    for i in range(len(axis_endca)):
        axis.append( axis_endca[i]-axis_startca[i])

    degree=round(random.uniform(11,40),2)
    print "rotate by degree ",degree

    print "Applying move on chain",apply_move_on_chain

    f = open(filepath, "r")
    out=open("output/updated_"+os.path.basename(filepath),"w")
    ter_count = 0
    for line in f:
        flag=1
        updated_line=line
        if ter_count != apply_move_on_chain:
            temp = line.split("  ")
            if temp[0] == "TER":
                ter_count += 1
        else:
            temp = line.split("  ")
            if temp[0] == "TER":
                ter_count += 1
            if temp[0] == "ATOM":
                temp_int = int(line[23:26])
                if temp_int in range(apply_on_segment, apply_on_segment + n):
                    if line[13:15]=="N " or line[13:15] == "CA" or line[13:15]=="C " or line[13:15]=="O ":
                        axis_temp=float(line[30:39])
                        v=list()
                        v.append(axis_temp)
                        axis_temp=float(line[40:47])
                        v.append(axis_temp)
                        axis_temp = float(line[47:54])
                        v.append(axis_temp)
                        new_co=rotate(v[0],v[1],v[2],axis_startca[0],axis_startca[1],axis_startca[2],axis[0],axis[1],axis[2],degree)
                        new_co = [round(l, 3) for l in new_co]
                        updated_line=updatepdbline(line,new_co)
                    else:
                        flag=0
        if flag!=0:
            out.write(updated_line)
    f.close()
    return "output/updated_"+os.path.basename(filepath),apply_on_segment

def rotate(x,y,z,a,b,c,u,v,w,theta):
    theta = math.radians(theta)
    L=math.pow(u,2)+math.pow(v,2)+math.pow(w,2)
    new_x=(((a*(math.pow(v,2)+math.pow(w,2))-u*(b*v+c*w-u*x-v*y-w*z))*(1-math.cos(theta)))+L*x*math.cos(theta)+math.pow(L,0.5)*(-c*v+b*w-w*y+v*z)*math.sin(theta))/L
    new_y=(((b*(math.pow(u,2)+math.pow(w,2))-v*(a*u+c*w-u*x-v*y-w*z))*(1-math.cos(theta)))+L*y*math.cos(theta)+math.pow(L,0.5)*(c*u-a*w+w*x-u*z)*math.sin(theta))/L
    new_z=(((c*(math.pow(u,2)+math.pow(v,2))-w*(a*u+b*v-u*x-v*y-w*z))*(1-math.cos(theta)))+L*z*math.cos(theta)+math.pow(L,0.5)*(-b*u+a*v-v*x+u*y)*math.sin(theta))/L
    new_out=list()
    new_out.append(new_x)
    new_out.append(new_y)
    new_out.append(new_z)
    return new_out
def addToEnsemble(pdbfile,ensemble_number):
    print "adding to ensemble"
    command = "cat "+pdbfile+" > "+ENSEMBLE_OUTPUT_DIR+"/"+str(ensemble_number)+".pdb"
    subprocess.check_output(command, shell=True)
    return ENSEMBLE_OUTPUT_DIR+"/"+str(ensemble_number)+".pdb"

def main():
    INPUT_PDB = sys.argv[1]
    output_dir_exists()
    scwrl_output = scwrl(INPUT_PDB) #fixing side chains
    print "Chain fixed pdb",scwrl_output
    fittedOriginalEnergy = foldx(scwrl_output)
    graph_energies=[fittedOriginalEnergy]
    print fittedOriginalEnergy,"fitted Original"
    ensemble_number = 0
    energies=["\t".join(["Ensemble Number".rjust(15),"Segment Number".rjust(15),"New Energy".rjust(15)])]
    print "Running for",sys.argv[3]
    for i in range(int(sys.argv[3])): #finding 10 ensembles
        print "Iteration",i,"Number of ensembles found",ensemble_number
        newpdb,apply_on_segment = backrub(INPUT_PDB,SEGMENT_SIZE)
        newpdb=scwrl(newpdb)
        newEnergy = foldx(newpdb)
        print "newEnergy",newEnergy
        # print "original Diff",newEnergy-originalEnergy
        changeInEnergy = newEnergy-fittedOriginalEnergy
        random.seed(time.time())
        random_number = random.uniform(0,1)
        print "random number",random_number
        
        metropolis_criterion=math.exp(-changeInEnergy / float(TEMPERATURE))
        print "metropolis_criterion", metropolis_criterion
        if changeInEnergy <= 0 or random_number > metropolis_criterion:
            ensemble_number += 1
            graph_energies.append(newEnergy)
            fittedOriginalEnergy = newEnergy
            energies.append("\t".join([str(ensemble_number).rjust(15),str(apply_on_segment).rjust(15),str(newEnergy).rjust(15)]))
            INPUT_PDB =addToEnsemble(newpdb,ensemble_number)
        print "---------------------------------------------------------------------------------------------------"
    print "graph energies",graph_energies
    plt.plot(graph_energies)
    plt.ylabel('Energy')
    plt.xlabel('Ensembles')
    plt.axis([0, len(graph_energies),min(graph_energies),max(graph_energies)])
    plt.savefig("graph.png")
    
    
    
    outputfile=open(sys.argv[1].split("/")[-1]+".backrub.log","w")
    outputfile.write("Selected Protein path:"+sys.argv[1]+"\n")
    outputfile.write("Segment Size:"+sys.argv[2]+"\n")
    outputfile.write("Number of iterations:"+sys.argv[3]+"\n")
    outputfile.write("------------------------------------------------------------------"+"\n")
    outputfile.write("\n".join(energies))
if __name__=="__main__":
    main()
