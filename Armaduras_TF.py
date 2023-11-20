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
        data.to_excel('KglobEl.xlsx')
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
        [self.k11, self.k12, self.k21, self.k22] = self.MatrizRigidezGlobParcial()

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
            clfx = self.N[Nf][2]
            clfy = self.N[Nf][3]
            element.append(Elemento(name, area, ME, I, xi, yi, xf, yf, [clix, cliy, clfx, clfy]))
        return element

    def MatrizRigidezGlobal(self):
        k = np.zeros((self.nG1, self.nG1))
        for i in range(self.nE):
            ke_global = self.armad[i].Rigidez_global
            for e in range(6):
                for j in range(6):
                    a = self.pi[i, e]-1
                    b = self.pi[i, j]-1
                    k[a, b] = ke_global[e,j] + k[a, b]
        data = pd.DataFrame(k)
        data.to_excel('KGlob.xlsx')
        return k

    def MatrizRigidezGlobParcial(self):
        k11 = self.Kglob[0:self.nkf, 0:self.nkf]
        k12 = self.Kglob[0:self.nkf, self.nkf:self.nG1]
        k21 = self.Kglob[self.nkf:self.nG1, 0:self.nkf]
        k22 = self.Kglob[self.nkf:self.nG1, self.nkf:self.nG1]
        return k11, k12, k21, k22

    def VecFuerzas(self):
        pass

def Graficar(tbl_Nodos:list, tbl_Elementos:list):
    import matplotlib.pyplot as plt
    import matplotlib.patches as ptch
    def apoyo_rigido(x1: float, y1: float, B: float):
        import matplotlib.pyplot as plt
        #plt.figure()
        plt.plot(B * np.array([0, 1]) + x1, B * np.array([0, 0]) + y1, color='k')
        plt.plot(B * np.array([0, 0.5]) + x1, B * np.array([0, 1]) + y1, color='k')
        plt.plot(B * np.array([0.5, 1]) + x1, B * np.array([1, 0]) + y1, color='k')
        plt.plot(B * np.array([-0.1, 1.1]) + x1, B * np.array([0, 0]) + y1, color='k')
        for i in range(11):
            plt.plot(B * np.array([0.0 + 0.1 * i, -0.1 + 0.1 * i])+x1, B * np.array([0, -0.16])+y1, color='k')

    def apoyo_empotrado(x1, y1, B):
        import matplotlib.pyplot as plt
        #plt.figure()
        plt.plot(B * np.array([-1, 1.1]) + x1, B * np.array([-2, -2]) + y1, color='k')
        for i in range(22):
            plt.plot(B * np.array([-1.0 + 0.1 * i, -1.1 + 0.1 * i]) + x1, B * np.array([-2.5, -2.1]) + y1, color='k')

    data1 = pd.DataFrame(tbl_Nodos)
    data2 = pd.DataFrame(tbl_Elementos)
    x = data1[1]
    y = data1[2]
    plt.figure()
    plt.suptitle('Trabajo Final')
    plt.title('Método de la rigidez')
    plt.minorticks_on()
    plt.grid(ls='-')
    plt.xlabel('[m]')
    plt.ylabel('[m]')
    for e in data2.values:
        ni, nf = e[4], e[5]
        x1, y1 = 0, 0
        x2, y2 = 0, 0
        for i in data1.values:
            nod = i[0]
            if nod == ni:
                x1, y1 =i[1], i[2]
            if nod == nf:
                x2, y2 = i[1], i[2]
        plt.plot([x1,x2],[y1,y2], color='k', ls='-')

    #Apoyo empotrado
    apoyo_rigido(13.69,3.3, 0.7)

    #Apoyo Empotrado
    apoyo_empotrado(0,1,0.5)

    style = "Simple, tail_width=0.5, head_width=4, head_length=5"
    kw = dict(arrowstyle=style, color="r")
    a3 = ptch.FancyArrowPatch((7.3,6.9), (6.8, 8.5),
                                 connectionstyle="arc3,rad=-0.9", **kw)
    a4 = ptch.FancyArrowPatch((0, 12.5), (0, 13.7),
                              connectionstyle="arc3,rad=0.9", **kw)
    plt.gca().add_patch(a3)
    plt.gca().add_patch(a4)

    plt.arrow(0,13, 0, -2, color='r', width=0.04)

    plt.arrow(3.2, 14, -0.2, -0.8, color='r', width=0.05)
    plt.arrow(11, 12, -0.9, -2.7, color='r', width=0.05)
    plt.arrow(-0.5, 0.5, 0.3, -0.3, color='r', width=0.05)
    plt.arrow(1.6, 7, 1.2, -0.85, color='r', width=0.05)

    plt.plot([3.2,11],[14,12], color='r', ls='-')
    plt.plot([-0.5,1.6], [0.5, 7], color='r', ls='-')

    plt.text(-0.4, 13, '1', color='b')
    plt.text(3+0.1, 6-0.5, '2', color='b')
    plt.text(14+0.1, 4+0.1, '3', color='b')
    plt.text(10-0.2, 9-0.7, '4', color='b')
    plt.text(0.3, 0.1, '5', color='b')
    plt.text(3-0.4, 13-0.5, '6', color='b')

    plt.text(11.1, 12, '25 kN/m', size=8)
    plt.text(3.2, 14.2, '5 kN/m', size=8)
    plt.text(0,13.8, '4 kN*m', size=8)
    plt.text(0.2, 11, '10 kN', size=8)
    plt.text(1.1, 7.1, '10 kN/m', size=8)
    plt.text(6.5, 8.6, '14 kN*m', size=8)
    plt.text(-0.9, 0.6, '2 kN/m', size=8, rotation=70)
    plt.show()



#E1 = Elemento('Elemento 1', 1, 1, 1, 0, 0, 4,3, [5,6,1,2])
#print(E1)

I = (0.5**3)*0.3/12
A= 0.3*0.5
E = 4700*np.sqrt(28)

print(I)
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
Graficar(tbl_Nodos, tbl_Elementos)

tbl_Frzas = [[-3, 'N2', 'Dy'],
             [-2, 'N1', 'Dy']]

tbl_Desp = [[0, 'N3', 'Dx'],
            [0, 'N3', 'Dy'],
            [0, 'N5', 'Dx'],
            [0, 'N5', 'Dy']]

print(AnalisisMatricial(tbl_Elementos, tbl_Nodos, tbl_Frzas, tbl_Desp))

E=24870.06
A=0.15
I=3.125*10**-3
L=3*np.sqrt(5)
l = np.array([[E*A/L, 0, 0, -E*A/L, 0, 0],
              [0, 12*E*I/(L**3), 6*E*I/(L**2), 0, -12*E*I/(L**3), 6*E*I/(L**2)],
              [0, 6*E*I/(L**2), 4*E*I/L, 0, -6*E*I/(L**2), 2*E*I/L],
              [-E*A/L, 0, 0, E*A/L, 0, 0],
              [0, -12*E*I/(L**3), -6*E*I/(L**2), 0, 12*E*I/(L**3), -6*E*I/(L**2)],
              [0, 6*E*I/(L**2), 2*E*I/L, 0, -6*E*I/(L**2), 4*E*I/L]])

#print(l)