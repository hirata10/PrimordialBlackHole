import sys 
import os

M = sys.argv[1]
l = int(sys.argv[2])
direc = 'M'+M+'-l'+str(l)

os.mkdir(direc)
os.mkdir(direc + '/MasterCopies')
os.mkdir(direc + '/MasterCopies/MinusPlus')
os.mkdir(direc + '/MasterCopies/PlusPlus')

os.system("cp M1e21-l1/*.py " +direc)


os.system("cp M1e21-l1/MasterCopies/MinusPlus/* "+direc+'/MasterCopies/MinusPlus')
os.system("cp M1e21-l1/MasterCopies/PlusPlus/* "+direc+'/MasterCopies/PlusPlus')

import time

#time.sleep(20)


with open('M1e21-l1/MasterCopies/MinusPlus/GeneratorInUpIn.pl', 'r') as file:
    data = file.readlines()
data[4] = '$l = '+str(l)+ '; \n'
with open(direc+'/MasterCopies/MinusPlus/GeneratorInUpIn.pl', 'w') as file:
    file.writelines(data)

with open('M1e21-l1/MasterCopies/MinusPlus/GeneratorInUpUp.pl', 'r') as file:
    data = file.readlines()
data[4] = '$l = '+str(l)+ '; \n'
with open(direc+'/MasterCopies/MinusPlus/GeneratorInUpUp.pl', 'w') as file:
    file.writelines(data)

with open('M1e21-l1/MasterCopies/MinusPlus/GeneratorUpInIn.pl', 'r') as file:
    data = file.readlines()
data[4] = '$l = '+str(l)+ '; \n'
with open(direc+'/MasterCopies/MinusPlus/GeneratorUpInIn.pl', 'w') as file:
    file.writelines(data)
    
with open('M1e21-l1/MasterCopies/MinusPlus/GeneratorUpInUp.pl', 'r') as file:
    data = file.readlines()
data[4] = '$l = '+str(l)+ '; \n'
with open(direc+'/MasterCopies/MinusPlus/GeneratorUpInUp.pl', 'w') as file:
    file.writelines(data)
    
with open('M1e21-l1/MasterCopies/MinusPlus/GeneratorUpUpIn.pl', 'r') as file:
    data = file.readlines()
data[4] = '$l = '+str(l)+ '; \n'
with open(direc+'/MasterCopies/MinusPlus/GeneratorUpUpIn.pl', 'w') as file:
    file.writelines(data)
    
with open('M1e21-l1/MasterCopies/MinusPlus/GeneratorUpUpUp.pl', 'r') as file:
    data = file.readlines()
data[4] = '$l = '+str(l)+ '; \n'
with open(direc+'/MasterCopies/MinusPlus/GeneratorUpUpUp.pl', 'w') as file:
    file.writelines(data)
    
    
with open('M1e21-l1/MasterCopies/PlusPlus/GeneratorInInIn.pl', 'r') as file:
    data = file.readlines()
data[4] = '$l = '+str(l)+ '; \n'
with open(direc+'/MasterCopies/PlusPlus/GeneratorInInIn.pl', 'w') as file:
    file.writelines(data)
    
with open('M1e21-l1/MasterCopies/PlusPlus/GeneratorInInUp.pl', 'r') as file:
    data = file.readlines()
data[4] = '$l = '+str(l)+ '; \n'
with open(direc+'/MasterCopies/PlusPlus/GeneratorInInUp.pl', 'w') as file:
    file.writelines(data)
    
with open('M1e21-l1/MasterCopies/PlusPlus/GeneratorUpInIn.pl', 'r') as file:
    data = file.readlines()
data[4] = '$l = '+str(l)+ '; \n'
with open(direc+'/MasterCopies/PlusPlus/GeneratorUpInIn.pl', 'w') as file:
    file.writelines(data)
    
with open('M1e21-l1/MasterCopies/PlusPlus/GeneratorUpInUp.pl', 'r') as file:
    data = file.readlines()
data[4] = '$l = '+str(l)+ '; \n'
with open(direc+'/MasterCopies/PlusPlus/GeneratorUpInUp.pl', 'w') as file:
    file.writelines(data)
    
with open('M1e21-l1/MasterCopies/PlusPlus/GeneratorUpUpIn.pl', 'r') as file:
    data = file.readlines()
data[4] = '$l = '+str(l)+ '; \n'
with open(direc+'/MasterCopies/PlusPlus/GeneratorUpUpIn.pl', 'w') as file:
    file.writelines(data)
    
with open('M1e21-l1/MasterCopies/PlusPlus/GeneratorUpUpUp.pl', 'r') as file:
    data = file.readlines()
data[4] = '$l = '+str(l)+ '; \n'
with open(direc+'/MasterCopies/PlusPlus/GeneratorUpUpUp.pl', 'w') as file:
    file.writelines(data)
    
