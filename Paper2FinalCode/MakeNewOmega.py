import sys 
import os

omega = float(sys.argv[1])
nplus = int(sys.argv[2])
nminus = int(sys.argv[3])

direc = 'Omega'+str(omega)+'T'

os.mkdir(direc)
os.mkdir(direc+'/MinusPlus')
os.mkdir(direc+'/MinusPlus/InUpIn')
os.mkdir(direc+'/MinusPlus/InUpUp')
os.mkdir(direc+'/MinusPlus/UpInIn')
os.mkdir(direc+'/MinusPlus/UpInUp')
os.mkdir(direc+'/MinusPlus/UpUpIn')
os.mkdir(direc+'/MinusPlus/UpUpUp')
os.mkdir(direc+'/PlusPlus')
os.mkdir(direc+'/PlusPlus/InInIn')
os.mkdir(direc+'/PlusPlus/InInUp')
os.mkdir(direc+'/PlusPlus/UpInIn')
os.mkdir(direc+'/PlusPlus/UpInUp')
os.mkdir(direc+'/PlusPlus/UpUpIn')
os.mkdir(direc+'/PlusPlus/UpUpUp')

os.system("cp MasterCopies/MinusPlus/*.py " +direc+'/MinusPlus')
os.system("cp MasterCopies/MinusPlus/*.pbs "+direc+'/MinusPlus')

os.system("cp MasterCopies/PlusPlus/*.py " +direc+'/PlusPlus')
os.system("cp MasterCopies/PlusPlus/*.pbs "+direc+'/PlusPlus')


with open('MasterCopies/MinusPlus/GeneratorInUpIn.pl', 'r') as file:
    data = file.readlines()
data[0] = '$N = ' + str(nminus)+ '; \n'
data[1] = '$omega = ' + str(omega) + '; \n'
with open(direc+'/MinusPlus/GeneratorInUpIn.pl', 'w') as file:
    file.writelines(data)

with open('MasterCopies/MinusPlus/GeneratorInUpUp.pl', 'r') as file:
    data = file.readlines()
data[0] = '$N = ' + str(nminus)+ '; \n'
data[1] = '$omega = ' + str(omega)+ '; \n'
with open(direc+'/MinusPlus/GeneratorInUpUp.pl', 'w') as file:
    file.writelines(data)

with open('MasterCopies/MinusPlus/GeneratorUpInIn.pl', 'r') as file:
    data = file.readlines()
data[0] = '$N = ' + str(nminus)+ '; \n'
data[1] = '$omega = ' + str(omega) + '; \n'
with open(direc+'/MinusPlus/GeneratorUpInIn.pl', 'w') as file:
    file.writelines(data)
    
with open('MasterCopies/MinusPlus/GeneratorUpInUp.pl', 'r') as file:
    data = file.readlines()
data[0] = '$N = ' + str(nminus)+ '; \n'
data[1] = '$omega = ' + str(omega) + '; \n'
with open(direc+'/MinusPlus/GeneratorUpInUp.pl', 'w') as file:
    file.writelines(data)
    
with open('MasterCopies/MinusPlus/GeneratorUpUpIn.pl', 'r') as file:
    data = file.readlines()
data[0] = '$N = ' + str(nminus)+ '; \n'
data[1] = '$omega = ' + str(omega) + '; \n'
with open(direc+'/MinusPlus/GeneratorUpUpIn.pl', 'w') as file:
    file.writelines(data)
    
with open('MasterCopies/MinusPlus/GeneratorUpUpUp.pl', 'r') as file:
    data = file.readlines()
data[0] = '$N = ' + str(nminus)+ '; \n'
data[1] = '$omega = ' + str(omega) + '; \n'
with open(direc+'/MinusPlus/GeneratorUpUpUp.pl', 'w') as file:
    file.writelines(data)
    
    
with open('MasterCopies/PlusPlus/GeneratorInInIn.pl', 'r') as file:
    data = file.readlines()
data[0] = '$N = ' + str(nplus)+ '; \n'
data[1] = '$omega = ' + str(omega) +'; \n'
with open(direc+'/PlusPlus/GeneratorInInIn.pl', 'w') as file:
    file.writelines(data)
    
with open('MasterCopies/PlusPlus/GeneratorInInUp.pl', 'r') as file:
    data = file.readlines()
data[0] = '$N = ' + str(nplus)+ '; \n'
data[1] = '$omega = ' + str(omega) +'; \n'
with open(direc+'/PlusPlus/GeneratorInInUp.pl', 'w') as file:
    file.writelines(data)
    
with open('MasterCopies/PlusPlus/GeneratorUpInIn.pl', 'r') as file:
    data = file.readlines()
data[0] = '$N = ' + str(nplus)+ '; \n'
data[1] = '$omega = ' + str(omega) +'; \n'
with open(direc+'/PlusPlus/GeneratorUpInIn.pl', 'w') as file:
    file.writelines(data)
    
with open('MasterCopies/PlusPlus/GeneratorUpInUp.pl', 'r') as file:
    data = file.readlines()
data[0] = '$N = ' + str(nplus)+ '; \n'
data[1] = '$omega = ' + str(omega) +'; \n'
with open(direc+'/PlusPlus/GeneratorUpInUp.pl', 'w') as file:
    file.writelines(data)
    
with open('MasterCopies/PlusPlus/GeneratorUpUpIn.pl', 'r') as file:
    data = file.readlines()
data[0] = '$N = ' + str(nplus)+ '; \n'
data[1] = '$omega = ' + str(omega) +'; \n'
with open(direc+'/PlusPlus/GeneratorUpUpIn.pl', 'w') as file:
    file.writelines(data)
    
with open('MasterCopies/PlusPlus/GeneratorUpUpUp.pl', 'r') as file:
    data = file.readlines()
data[0] = '$N = ' + str(nplus)+ '; \n'
data[1] = '$omega = ' + str(omega) +'; \n'
with open(direc+'/PlusPlus/GeneratorUpUpUp.pl', 'w') as file:
    file.writelines(data)
    
