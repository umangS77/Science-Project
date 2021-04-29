import numpy as np
from itertools import combinations
import autograd
import sys
import os
import re

import functools
from autograd import grad
import autograd.numpy as anp
import autograd
from multiprocessing import Pool, Manager
from matplotlib import pyplot as plt

hartree2J = 4.359744e-18
amu2kg = 1.660538782e-27
bohr2m = 0.52917720859e-10
c = 2.99792458E8

def sep(p1, p2):
    return p2 - p1


def pbc_sep(p1, p2):
    arrow = sep(p1, p2)
    rem = np.mod(arrow, 18)  # The image in the first cube
    mic_separation_vector = np.mod(rem+18/2, 18)-18/2

    return np.array(mic_separation_vector)

class Molecule:

    def __init__(self, xyz_file, units="Angstrom"):

        self.units = units
        self.read(xyz_file)
        # Mass of argon.. Sir mentioned its argon..
        self.masses = [float(39.96238312251) for i in self.atoms]

    def read(self, xyz_file):
        geom_str = xyz_file.read()
        self.atoms = []
        geom = []
        for line in geom_str.split('\n')[2:]:
            if line.strip() == '':
                continue
            atom, x, y, z = line.split()[:4]
            self.atoms.append(atom)
            geom.append([float(x), float(y), float(z)])
        self.geom = np.array(geom)

    def __len__(self):
        return len(self.geom)

    def __str__(self):
        out = "{:d}\n{:s}\n".format(len(self), self.units)
        for atom, xyz in zip(self.atoms, self.geom):
            out += "{:2s} {: >15.10f} {: >15.10f} {: >15.10f}\n".format(
                atom, *xyz)
        return out

    def bohr(self):
        if self.units == "Angstrom":
            self.geom *= 1.889725989
            self.units = "Bohr"
        return self.geom

    def angs(self):
        if self.units == "Bohr":
            self.geom /= 1.889725989
            self.units = "Angstrom"
        return self.geom

    def copygeom(self):
        return np.array(self.geom)


class Frequencies:

    def __init__(self, mol, hessString):

        self.mol = mol
        self.hess = hessString
        self.N = mol.__len__()

        m = []
        for i in range(self.N):
            m += [1/(mol.masses[i])**0.5]*3
        self.MM = np.diag(m)
        self.m = m

    def get_MWhessian(self):

        H0 = np.matrix([i.split() for i in self.hess.splitlines()], float)
        mwH = np.dot(self.MM, np.dot(H0, self.MM))
        return mwH

    def get_frequencies(self):

        self.e, self.l = np.linalg.eigh(self.get_MWhessian())
        self.Q = np.matrix(self.MM)*np.matrix(self.l)
        freq = []
        conv = np.sqrt(hartree2J/(amu2kg*bohr2m**2)
                       ) / (c*2*np.pi)  # dimensional analysis
        # print(conv)
        for i in self.e:
            if i < 0:
                freq.append((-i)**0.5*conv)
            else:
                freq.append(i**0.5*conv)

        return freq

    def frequency_output(self, output):

        mol = self.mol
        freq = self.get_frequencies()

        t = open(output, "w")
        for i in range(3*self.N):
            t.write("%d\n%s cm^{-1}\n" % (self.N, str(freq[i])))
            for j in range(self.N):
                atom = mol.atoms[j]
                x, y, z = mol.geom[j, 0], mol.geom[j, 1], mol.geom[j, 2]
                dx, dy, dz = self.Q[3*j, i], self.Q[3*j+1, i], self.Q[3*j+2, i]
                t.write("{:s}{:12.7f}{:12.7f}{:12.7f}\n".format(atom, x, y, z))
            t.write("\n")
        t.close()
        a = np.array(freq)
        plt.hist(a, bins=100)
        plt.title("Frequency Histogram")
        plt.savefig("hist.png")
        # plt.show()

        return None

def a_sep(p1, p2):
    return p2 - p1


def a_pbc_sep(p1, p2):

    arrow = a_sep(p1, p2)
    rem = anp.mod(arrow, 18)  # The image in the first cube
    mic_separation_vector = anp.mod(rem+18/2, 18)-18/2

    return anp.array(mic_separation_vector)


def a_get_energy(geom):
    eps = 0.238
    sigma = 3.4
    te = 0
    tgeom = anp.array(geom)
    pairs = [(a, b) for idx, a in enumerate(tgeom)
             for b in tgeom[idx + 1:]]
    for pair in pairs:
        rij = anp.linalg.norm(a_pbc_sep(pair[0], pair[1]))
        if rij == 0:
            continue
        te += 4*eps*((sigma/rij)**12-(sigma/rij)**6)

    return te


def gradient_descent(alpha, max_its, w):
    # compute gradient module using autograd
    gradient = grad(a_get_energy)
    # print("A")
    # run the gradient descent loop
    weight_history = [w]           # container for weight history
    # container for corresponding cost function history
    cost_history = [a_get_energy(w)]
    # print("B")
    print("In progress: ")
    for k in range(max_its):
        print('\r' + str(k+1) + "%", end = '')
        # print(k)
        # evaluate the gradient, store current weights and cost function value
        #print(k, w)
        # print("It:", k)
        grad_eval = gradient(w)

        # take gradient descent step
        w = w - alpha*grad_eval
        # print(a_get_energy(w))
        # record weight and cost
        weight_history.append(w)
        cost_history.append(a_get_energy(w))
    print()
    return weight_history, cost_history



def get_energy(geom):
    eps = 0.238
    sigma = 3.4
    te = 0
    tgeom = np.array(geom)
    pairs = [(a, b) for idx, a in enumerate(tgeom)
             for b in tgeom[idx + 1:]]
    for pair in pairs:
        rij = np.linalg.norm(pbc_sep(pair[0], pair[1]))
        if rij == 0:
            continue
        te += 4*eps*((sigma/rij)**12-(sigma/rij)**6)

    return te

def get_config():
    P = []
    count = 0
    with open("molecule.xyz", "w") as f:
        print(108, file=f)
        print("\n", end="", file=f)
        print("In progress:")
        while True:
            p = np.random.rand(3,) * 18
            skip = False
            for p_in_box in P:
                if np.linalg.norm(pbc_sep(p, p_in_box)) <= 3.4:
                    skip = True
                    # print("skipping...")
            print(len(P))
            if not skip:
                P.append(p)
                count = 0
                if len(P) == 108:
                    break
            else:
                count += 1
                if count >= 10000:
                    count = 0
                    tempP = []
                    tempP.append(p)
                    for j in range(len(P)):
                        if np.linalg.norm(pbc_sep(P[j], p)) >= 3.4:
                            tempP.append(P[j])
                    P = tempP
            # print('\r' + str(min(len(P)+1,100)) + "%", end = '')

        print()
        for point in P:
            print(f"C {point[0]} {point[1]} {point[2]}", file=f)


class Hessian(object):

    def __init__(self, mol, disp_size=0.005):

        self.mol = mol
        self.N = len(self.mol)
        self.h = disp_size
        self.energy = {}

    def find_E(self, i, j, hi, hj):
        """
        :params i,j: indices of atoms 1,2
        :params hi,hj: displacements of atoms 1,2 (-1, 0, or 1, corresponds to -h, 0, or h)
        """
        key = "X%dX%d_%d%d" % (i, j, hi, hj)
        return self.energy[key]

    def set_energy(self, key, geom, d=None):
        e = get_energy(np.array(geom))
        if d != None:
            d[key] = e
            return
        self.energy[key] = e

    def process(self, i, d):
        h, N, atoms, geom = self.h, self.N, self.mol.atoms, self.mol.geom

        # print(i)
        for j in range(i):
            forward = "X%dX%d_11" % (i, j)
            reverse = "X%dX%d_-1-1" % (i, j)
            geom_copy2 = self.mol.copygeom()

            geom_copy2[i//3, i % 3] = geom_copy2[i//3, i % 3] + h
            geom_copy2[j//3, j % 3] = geom_copy2[j//3, j % 3] + h

            self.set_energy(forward, geom_copy2, d)

            geom_copy2[i//3, i % 3] = geom_copy2[i//3, i % 3] - 2*h
            geom_copy2[j//3, j % 3] = geom_copy2[j//3, j % 3] - 2*h

            self.set_energy(reverse, geom_copy2, d)

    def run_disps(self):

        h, N, atoms, geom = self.h, self.N, self.mol.atoms, self.mol.geom
        self.set_energy("X0X0_00", geom)

        ####   Run single displacements   ####
        for i in range(3*N):
            # print(i)
            forward = "X%dX0_10" % i
            reverse = "X%dX0_-10" % i
            geom_copy = self.mol.copygeom()
            geom_copy[i//3, i % 3] = geom_copy[i//3, i % 3]+h
            self.set_energy(forward, geom_copy)

            geom_copy[i//3, i % 3] = geom_copy[i//3, i % 3]-2*h
            self.set_energy(reverse, geom_copy)
        ####   Run double displacements    ######
        mylist = [*range(3*N)]
        pool = Pool()
        D = Manager().dict()                     # Create a multiprocessing Pool
        pool.map(functools.partial(self.process, d=D), mylist)
        pool.close()
        pool.join()
        self.energy.update(D)

    def make_Hessian(self):

        self.run_disps()

        h, N = self.h, self.N
        E0 = self.find_E(0, 0, 0, 0)
        self.H = np.zeros((3*self.N, 3*self.N))
        for i in range(3*N):
            # print(i)
            for i in range(3*N):
                self.H[i, i] = (self.find_E(i, 0, 1, 0) +
                                self.find_E(i, 0, -1, 0)-2*E0)/(h**2)
                for j in range(0, i):
                    self.H[i, j] = (self.find_E(i, j, 1, 1)+self.find_E(i, j, -1, -1)-self.find_E(
                        i, 0, 1, 0)-self.find_E(j, 0, 1, 0)-self.find_E(j, 0, -1, 0)-self.find_E(i, 0, -1, 0)+2*E0)
                    self.H[i, j] /= 2*h**2
                    self.H[j, i] = self.H[i, j]

    def make_eigh(self):
        w, v = np.linalg.eigh(self.H)
        np.savetxt("eigen_vectors.dat", v, "%15.7f", " ", "\n")
        np.savetxt("eigen_values.dat", w, "%15.7f", " ", "\n")

    def write_Hessian(self):
        self.make_Hessian()
        self.make_eigh()
        np.savetxt("hessian.dat", self.H, "%15.7f", " ", "\n")

def main():
    print("Initiating the program...\nTakes some time to set Configuration. Please wait...")
    # get_config()  # Q1
    init_mol = None
    with open("molecule.xyz") as f:
        print("Initial Configuration Generated for Arg LJ Model")
        print("Initial Configuration for Q1 is saved in molecule.xyz")
        init_mol = Molecule(f, "Angstrom")
        print("Doing Q2 and Q3...")
        e = get_energy(np.array(init_mol.geom))  # Q2
        print("Starting Steepest Descent Algoithm to minimise Energy...")
        g = gradient_descent(0.135, 100, init_mol.geom)
        print("Steepest Descent Algoithm Ended.. Results below\n")
        print("Original Energy: ", e)
        print("New Energy:", np.min(np.array(g[1])), "\n")
        i = np.argmin(np.array(g[1]))
        print("New Configuration will be saved in new_molecule.xyz")
        with open("new_molecule.xyz", "w") as f2:
            print(len(g[0][i]), file=f2)
            print(file=f2)
            for p in g[0][i]:
                print(
                    f"C {p[0]} {p[1]} {p[2]}", file=f2)
    with open("new_molecule.xyz") as f2:
        print("Generating Hessian Matrix... This will take some time.(15-20 mins)")
        mol = Molecule(f2, "Angstrom")
        mol.bohr()
        hessian = Hessian(mol, 0.00001)
        hessian.write_Hessian()
        print("Q4 has been outputted as eigen_vectors.dat, eigen_values.dat, hessian.dat")
    with open("new_molecule.xyz", "r") as f2:
        print("Q5:")
        mol = Molecule(f2)
        hessian = open("hessian.dat", "r").read()
        freq = Frequencies(mol, hessian)
        freq.frequency_output("modes.xyz")
        print("Modes have been written to modes.xyz")
        print("A histogram for the frequencies have been saved as hist.png")
        print("Kindly look at Report for the output format of these files and what they contain..")
        print("End")


main()