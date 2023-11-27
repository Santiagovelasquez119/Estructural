#coding=utf8
import numpy as np
import pandas as pd
class Elemento():
    def __init__(self, elemento, Area, ModElasticidad, Inercia, Coord_xi, Coord_yi, Coord_xf, Coord_yf, vec_coord):
        self.elemento = elemento
        self.A = Area
        self.E = ModElasticidad
        self.I = Inercia
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
        Glob = np.array([[(self.E * self.A) / self.l, 0, 0, -self.E * self.A / self.l, 0, 0],
                      [0, 12 * self.E * self.I / (self.l ** 3), 6 * self.E * self.I / (self.l ** 2), 0, -12 * self.E * self.I / (self.l ** 3), 6 * self.E * self.I / (self.l ** 2)],
                      [0, 6 * self.E * self.I / (self.l ** 2), 4 * self.E * self.I / self.l, 0, -6 * self.E * self.I / (self.l ** 2), 2 * self.E * self.I / self.l],
                      [-self.E * self.A / self.l, 0, 0, self.E * self.A / self.l, 0, 0],
                      [0, -12 * self.E * self.I / (self.l ** 3), -6 * self.E * self.I / (self.l ** 2), 0, 12 * self.E * self.I / (self.l ** 3), -6 * self.E * self.I / (self.l ** 2)],
                      [0, 6 * self.E * self.I / (self.l ** 2), 2 * self.E * self.I / self.l, 0, -6 * self.E * self.I / (self.l ** 2), 4 * self.E * self.I / self.l]])
        data = pd.DataFrame(Glob)
        #data.to_excel('KglobEl.xlsx')
        return Glob

class AnalisisMatricial():
    def __init__(self, tbl_Elem, tbl_Nodo, tbl_Frza, tbl_Desp):
        self.tE = tbl_Elem
        self.tN = tbl_Nodo
        self.tF = tbl_Frza
        self.tD = tbl_Desp
        self.nE = len(self.tE)
        self.nN = len(self.tN)
        self.nG1 = self.nN*3
        [self.nr, self.nkf, self.N] = self.VectorCoordenadasGlobales()
        self.pi = np.array(self.Matriz_Pi())
        self.armad = self.Armadura()
        self.Kglob = self.MatrizRigidezGlobal()

    # def __str__(self):
    #     print(f'# elementos: {self.nE}')
    #     print(f'# nodos: {self.nN}')
    #     print(f'# grados de libertad: {self.nG1}')
    #     print(f'Vector coordenadas globales: {self.nr, self.nkf, self.N} ')
    #     print(f'Matriz pi: {self.pi}')
    #     print('Elemento 3:', self.armad[2])
    #     print('KGlobal: ', self.Kglob)
    #     return ''

    def VectorCoordenadasGlobales(self):
        dicnodos = {}
        cont = 1
        g_l= self.nG1
        numReacciones = 1
        for i in self.tN:
            nodo_key = i[0]
            c_x = i[1]
            c_y = i[2]
            tipo = i[3]
            if tipo == 'Libre':
                numReacciones += 0
                dicnodos.setdefault(nodo_key, [c_x, c_y, cont, cont+1, cont+2])
                cont+=3
            elif tipo == 'Fijo':
                numReacciones += 2
                dicnodos.setdefault(nodo_key, [c_x, c_y, g_l-2, g_l-1, g_l])
                g_l-=3
            elif tipo == 'Dx':
                numReacciones += 1
                dicnodos.setdefault(nodo_key, [c_x, c_y])
            elif tipo == 'Dy':
                numReacciones += 1
                dicnodos.setdefault(nodo_key, [c_x, c_y])
        dicnodos['N3'][4] = 13
        dicnodos['N5'][2] = 18
        knownforces = self.nG1 - numReacciones
        return numReacciones, knownforces, dicnodos

    def Matriz_Pi(self):
        pi = []
        for e in self.tE:
            ni, nf = e[4], e[5]
            pi.append([self.N[ni][2], self.N[ni][3], self.N[ni][4], self.N[nf][2], self.N[nf][3], self.N[nf][4]])
        return pi

    def Armadura(self):
        element = []
        for e in self.tE:
            name = e[0]
            area = e[1]
            ME = e[2]
            I = e[3]
            Ni = e[4]
            Nf = e[5]
            xi = self.N[Ni][0]
            yi = self.N[Ni][1]
            xf = self.N[Nf][0]
            yf = self.N[Nf][1]
            clix = self.N[Ni][2]
            cliy = self.N[Ni][3]
            clmi = self.N[Ni][4]
            clfx = self.N[Nf][2]
            clfy = self.N[Nf][3]
            clmf = self.N[Nf][4]
            element.append(Elemento(name, area, ME, I, xi, yi, xf, yf, [clix, cliy, clmi, clfx, clfy, clmf]))
        return element

    def MatrizRigidezGlobal(self):
        k = np.zeros((self.nG1, self.nG1))
        for e in range(self.nE):
            ke_global = self.armad[e].Rigidez_global
            for i in range(6):
                for j in range(6):
                    a = self.pi[e, i]-1
                    b = self.pi[e, j]-1
                    c = self.pi[e,]
                    print(a,b)
                    k[a, b] = ke_global[i,j] + k[a, b]
        data = pd.DataFrame(k)
        data.to_excel('KGlob.xlsx')
        return k

I = (0.5**3)*0.3/12
A= 0.5*0.3
E = 4700*np.sqrt(28)/1000


tbl_Elementos = [['E1', A, E, I, 'N1', 'N6'],
                 ['E2', A, E, I, 'N6', 'N4'],
                 ['E3', A, E, I, 'N4', 'N3'],
                 ['E4', A, E, I, 'N4', 'N2'],
                 ['E5', A, E, I, 'N6', 'N2'],
                 ['E6', A, E, I, 'N2', 'N5'],
                 ]

tbl_Nodos = [['N1', 0, 13, 'Libre'],
             ['N2', 3, 6, 'Libre'],
             ['N3', 14, 4, 'Fijo'],
             ['N4', 10, 9, 'Libre'],
             ['N5', 0, 0, 'Fijo'],
             ['N6', 3, 13, 'Libre']]

tbl_Frzas = [[-3, 'N2', 'Dy'],
             [-2, 'N1', 'Dy']]

tbl_Desp = [[0, 'N3', 'Dx'],
            [0, 'N3', 'Dy'],
            [0, 'N5', 'Dx'],
            [0, 'N5', 'Dy']]

print(AnalisisMatricial(tbl_Elementos, tbl_Nodos, tbl_Frzas, tbl_Desp))