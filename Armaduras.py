#coding=utf8
import numpy as np
class Elemento():
    def __init__(self, elemento, Area, ModElasticidad, Coord_xi, Coord_yi, Coord_xf, Coord_yf, vec_coord):
        self.elemento = elemento
        self.A = Area
        self.E = ModElasticidad
        self.xi = Coord_xi
        self.yi = Coord_yi
        self.xf = Coord_xf
        self.yf = Coord_yf
        self.vc = vec_coord
        #Propiedades del elemento
        self.l = self.Longitud()
        self.kloc = self.rig_local()
        self.Lx = self.Lambda_x()
        self.Ly = self.Lambda_y()
        self.Trans_loc_glob = self.Trans()
        self.Rigidez_global = self.Rig_global()

    def __str__(self):
        print(f'Area = {self.A}')
        print(f'Módulo de eslasticidad = {self.E}')
        print(f'Punto inicial= {(self.xi, self.yi)}')
        print(f'Puntos finales = {(self.xf, self.yf)}')
        print(f'Vector de coordenadas = {self.vc}')
        print(f'Longitud del elemento = {self.l}')
        print(f'Matriz de rigidez local = {self.kloc}')
        print(f'Lambda_x = {self.Lx}')
        print(f'Lambda_y = {self.Ly}')
        print(f'Matriz de transformación = {self.Trans_loc_glob}')
        print(f'Matriz de rigidez global: {self.Rigidez_global}')
        return ''

    #Definición de propiedades del elemento
    def Longitud(self):
        l = np.sqrt((self.xf-self.xi)**2+(self.yf-self.yi)**2)
        return l

    def rig_local(self):
        k_loc = (self.A*self.E/self.l)*np.array([[1, -1], [-1, 1]])
        return k_loc

    def Lambda_x(self):
        return (self.xf-self.xi)/self.l
    def Lambda_y(self):
        return (self.yf-self.yi)/self.l

    def Trans(self):
        return np.array([[self.Lx, self.Ly, 0, 0], [0, 0, self.Lx, self.Ly]])

    def Rig_global(self):
        Trt = np.transpose(self.Trans_loc_glob)
        Tt_local = np.matmul(Trt, self.kloc)
        return np.matmul(Tt_local, self.Trans_loc_glob)

class AnalisisMatricial():
    def __init__(self, tbl_Elem, tbl_Nodo, tbl_Frza, tbl_Desp):
        self.tE = tbl_Elem
        self.tN = tbl_Nodo
        self.tF = tbl_Frza
        self.tD = tbl_Desp

        self.nE = len(self.tE)
        self.nN = len(self.tN)
        self.nG1 = self.nN*2
        [self.nr, self.nkf, self.N] = self.VectorCoordenadasGlobales()
        self.pi = np.array(self.Matriz_Pi())
        self.armad = self.Armadura()
        self.Kglob = self.MatrizRigidezGlobal()
        [self.k11, self.k12, self.k21, self.k22] = self.MAtrizRigidezGlobParcial()

    def __str__(self):
        print(f'# elementos: {self.nE}')
        print(f'# nodos: {self.nN}')
        print(f'# grados de libertad: {self.nG1}')
        print(f'Vector coordenadas globales: {self.nr, self.nkf, self.N} ')
        print(f'Matriz pi: {self.pi}')
        print('Elemento 1:', self.armad[0])
        print('KGlobal: ', self.Kglob)
        print('K11: ', self.k11)
        print('K12: ', self.k12)
        print('K21: ', self.k21)
        print('K22: ', self.k22)
        return ''

    def VectorCoordenadasGlobales(self):
        dicnodos = {}
        cont = 1
        g_l= self.nG1
        numReacciones = 0
        for i in self.tN:
            nodo_key = i[0]
            c_x = i[1]
            c_y = i[2]
            tipo = i[3]
            if tipo == 'Libre':
                numReacciones += 0
                dicnodos.setdefault(nodo_key, [c_x, c_y, cont, cont+1])
                cont+=2
            elif tipo == 'Fijo':
                numReacciones += 2
                dicnodos.setdefault(nodo_key, [c_x, c_y, g_l-1, g_l])
                g_l-=2
            elif tipo == 'Dx':
                numReacciones += 1
                dicnodos.setdefault(nodo_key, [c_x, c_y])
            elif tipo == 'Dy':
                numReacciones += 1
                dicnodos.setdefault(nodo_key, [c_x, c_y])
        knownforces = self.nG1 - numReacciones
        return numReacciones, knownforces, dicnodos

    def Matriz_Pi(self):
        pi = []
        for e in self.tE:
            ni, nf = e[3], e[4]
            pi.append([self.N[ni][2], self.N[ni][3], self.N[nf][2], self.N[nf][3]])
        return pi

    def Armadura(self):
        element = []
        for e in self.tE:
            name = e[0]
            area = e[1]
            ME = e[2]
            Ni = e[3]
            Nf = e[4]
            xi = self.N[Ni][0]
            yi = self.N[Ni][1]
            xf = self.N[Nf][0]
            yf = self.N[Nf][1]
            clix = self.N[Ni][2]
            cliy = self.N[Ni][3]
            clfx = self.N[Nf][2]
            clfy = self.N[Nf][3]
            element.append(Elemento(name, area, ME, xi, yi, xf, yf, [clix, cliy, clfx, clfy]))
        return element

    def MatrizRigidezGlobal(self):
        k = np.zeros((self.nG1, self.nG1))
        for i in range(self.nE):
            ke_global = self.armad[i].Rigidez_global
            for e in range(4):
                for j in range(4):
                    a = self.pi[i, e]-1
                    b = self.pi[i, j]-1
                    k[a, b] = ke_global[e,j] + k[a, b]
        return k

    def MAtrizRigidezGlobParcial(self):
        k11 = self.Kglob[0:self.nkf, 0:self.nkf]
        k12 = self.Kglob[0:self.nkf, self.nkf:self.nG1]
        k21 = self.Kglob[self.nkf:self.nG1, 0:self.nkf]
        k22 = self.Kglob[self.nkf:self.nG1, self.nkf:self.nG1]
        return k11, k12, k21, k22


E1 = Elemento('Elemento 1', 1, 1, 0, 0, 4,3, [5,6,1,2])
#print(E1)

tbl_Elementos = [['E1', 1, 1, 'N5', 'N1'],
                 ['E2', 1, 1, 'N2', 'N1'],
                 ['E3', 1, 1, 'N3', 'N2'],
                 ['E4', 1, 1, 'N2', 'N5'],
                 ['E5', 1, 1, 'N3', 'N5'],
                 ['E6', 1, 1, 'N4', 'N5'],
                 ['E7', 1, 1, 'N4', 'N3']]

tbl_Nodos = [['N1', 14, 7, 'Libre'],
             ['N2', 7, 7, 'Libre'],
             ['N3', 0, 7, 'Fijo'],
             ['N4', 0, 0, 'Fijo'],
             ['N5', 7, 0, 'Libre']]

tbl_Frzas = [[-3, 'N2', 'Dy'],
             [-2, 'N1', 'Dy']]

tbl_Desp = [[0, 'N3', 'Dx'],
            [0, 'N3', 'Dy'],
            [0, 'N4', 'Dx'],
            [0, 'N4', 'Dy']]

print(AnalisisMatricial(tbl_Elementos, tbl_Nodos, tbl_Frzas, tbl_Desp))
