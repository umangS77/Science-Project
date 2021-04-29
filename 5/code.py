import math

charges = {
    'H1':0.417000,
    'H2':0.417000,
    'O':-0.834000
}

Kc = 332.1
A = 582000
B = 595.0

Lx=23.623
Ly=22.406
Lz=27.1759



def take_input():
    global data
    data = open('./starting_config_300k.pdb').readlines()

coordinates = []

def func():
    for i in range(0,3*(len(data)//3),3):
        data[i] = data[i].split()
        O_coord = {'x':float(data[i][5]),'y':float(data[i][6]),'z':float(data[i][7])}
        data[i+1] = data[i+1].split()
        H1_coord = {'x':float(data[i+1][5]),'y':float(data[i+1][6]),'z':float(data[i+1][7])}
        data[i+2] = data[i+2].split()
        H2_coord = {'x':float(data[i+2][5]),'y':float(data[i+2][6]),'z':float(data[i+2][7])}
        coordinates.append({'O':O_coord,'H1':H1_coord,'H2':H2_coord})

def find_potential(i,j):
    p = 0.0
    
    for mol_A,pos_A in coordinates[i].items():
        for mol_B,pos_B in coordinates[j].items():
            Rsquare = (abs(pos_A['x']-pos_B['x']) - Lx*(round((abs(pos_A['x']-pos_B['x']))/Lx)))**2
            Rsquare+= (abs(pos_A['y']-pos_B['y']) - Ly*(round((abs(pos_A['y']-pos_B['y']))/Ly)))**2
            Rsquare+= (abs(pos_A['z']-pos_B['z']) - Lz*(round((abs(pos_A['z']-pos_B['z']))/Lz)))**2
            d = math.sqrt(Rsquare)
            if mol_A=='O' and mol_B=='O':
                Roo = d
            p+=((Kc*charges[mol_A]*charges[mol_B])/d)
    v = ((A/(Roo**12))-(B/(Roo**6)))
    return v,p


charges = {
    'H1':0.417000,
    'H2':0.417000,
    'O':-0.834000
}

def main():
    take_input()
    func()
    num_of_molecules = len(coordinates)
    PE = 0.0 # total potential energy
    VE = 0.0 # total vanderwaal energy

    for i in range(num_of_molecules):
        for j in range(i):
            if i!=j:
                v, p = find_potential(i,j)
                VE+=v
                PE+=p
    print("Electrostatic energy: ",PE)
    print("Vanderwaal energy: ",VE)
    print("Total energy: = Electrostatic energy + Vanderwaal energy = ",VE+PE)

main()