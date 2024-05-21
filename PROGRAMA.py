import csv
import sys
import matplotlib.pyplot as plt

def disenho(X, Y, xy):
    """
    Verificar si el punto xy con coordenadas (x, y) está dentro o fuera de la 
    del polígono definido por los vértices X, Y.

    Args:
        X (list): coordenadas x del polígono.
        y (list): coordenadas y del polígono.
        xy (list): coordenadas del punto x, y.

    Returns:
        bool: el punto xy está o no dentro del polígono.
    """

    num_vertices = len(X)
    intersecciones = 0

    # recorrer los vértices del polígono
    for i in range(num_vertices):
        # print(f"iteración {i}")
        j = (i + 1) % num_vertices

        # coordenadas del punto i
        xi = X[i]
        yi = Y[i]

        # coordenadas del punto j
        xj = X[j]
        yj = Y[j]

        if abs(xj - xi) < 1e-9:
            # print("vert")
            if xy[0] <= xi:
                intersecciones+=1
            
            # print(intersecciones)
            continue

        if abs(yj - yi) < 1e-9:
            # print("hori")
            if abs(xy[1] - yj) < 1e-9:
                intersecciones+=1
            # print(intersecciones)
            continue
        
        # print("ninguna de las dos")
        pend = (yj - yi) / (xj - xi)

        x = (xy[1] - yi + pend * xi) / pend
        
        if  min(xi, xj) <= x <= max(xi, xj) and xy[0] <= x:
            intersecciones+=1

        # print(intersecciones)
        # input()
        
    if intersecciones % 2 != 0:
        return True
    else:
        return False

def beta_1(f_c):
    """Calcular B1 según norma (ACI) teniendo en cuenta f'c    
    Args:
        f_c (float): (f'c) Resistencia característica del concreto
        
    Returns:
        B_1 (float): B_1
    """
    if 17 <= f_c <= 28:
        B1 = 0.85
    elif 28 < f_c < 55:
        B1 = 0.85 - (0.05 *((f_c) - 28 ) / 7 )
    else: 
        B1 = 0.65

    return B1

def a_m( c, B1, h):  
    """ calcular a  (Altura del bloque de compresión) 
    teniendo en cuenta el condicionante a <= h
    Args:
        c (float)  : Distancia de la cara comprimida del hormigón al eje neutro expresada en (m)
        B1 (float) : (Beta 1) Se calcula según norma
        h (float)  : Altura de la sección transversal de la columna expresada en (m)

    Returns:
        a (float): Altura del bloque de compresión en (m).
    """
    if c * B1 < h:
        a = c * B1
    else:
        a = h

    return a

def fs_( Es, c, d_, fy ):
    """ Calcular f's (Esfuerzo de acero a compresión)

    Args:
        Es (float) : modulo elástico del acero, es constante Es = 200'000 Mpa
        c (float)  : Distancia de la cara comprimida del hormigón al eje neutro expresada en (m)
        d_ (float) : (d') distancia más comprimida de la fibra al eje de aceros a compresión (recubrimiento)
        fy (float) : Resistencia a fluencia del acero expresada en (Mpa)

    Returns:
        f's (float): esfuerzo de acero a compresión en (Mpa).
    """
    if Es * 0.003 * ((c - d_) / c) < -fy:
        fs= -fy
    elif Es * 0.003 * ( (c - d_) / c ) > fy:
        fs = fy
    else: 
        fs = Es * 0.003 * ((c - d_) / c)

    return fs

def f_s(Es, d, c,fy):
    """Calcular fs (Esfuerzo de aceros a tracción)

    Args:
        Es (float) : modulo elástico del acero, es constante Es = 200'000 Mpa
        c (float)  : Distancia de la cara comprimida del hormigón al eje neutro expresada en (m)
        d (float)  : d distancia más comprimida de la fibra al eje de aceros a tracción
        fy (float) : Resistencia a fluencia del acero expresada en (Mpa)

    Returns:
        fs (float): Esfuerzo de aceros a tracción en (Mpa)
    """
    
    if Es * 0.003 * (d - c) /c > fy :
        _fs_ = fy
    elif Es * 0.003 * (d - c) /c < -fy:
        _fs_ = -fy        
    else:
        _fs_ = Es * 0.003 * (d - c) /c 
        
    return _fs_

def momento(r, f_c, a, b, h, A_S, f_s, d_, AS, Fs, d,): 
    """ Calcular Mn (Resistencia nominal a momento)

    Args:
        r (float)   : Gama valor constante que vale 0.85
        f_c (float) : (f's) Resistencia característica del concreto 
        a (float)   : Altura del bloque de compresión 
        b (float)   : Ancho de la columna 
        h (float)   : Altura de la sección transversal de la columna
        A_S (float) : (A's) Área de acero a compresión 
        f_s (float) : (f's) Esfuerzo de acero a compresión
        d_ (float)  : (d') distancia más comprimida de la fibra al eje de aceros a compresión (recubrimiento)
        AS (float)  : Área de acero a tracción
        Fs (float)  : Esfuerzo de acero a tracción
        d (float)   : distancia más comprimida de la fibra al eje de aceros a tracción


    Returns:
        Mn (float): Resistencia nominal a momento en KNm 
    """
    M = (r * f_c * a * b * ((h / 2) - (a / 2 )) + 
         (A_S * 0.0004) * f_s * ((h / 2) - d_) + (AS * 0.0004) * Fs * (d - (h / 2))) * 1000
    return M

def presion0Pn(r, f_c, a, b, f_s, A_s, As, fs):
    """calcular Pn (Resistencia axial)

    Args:
        r (float)   : Gama valor constante que vale 0.85
        f_c (float) : (F'c) Resistencia característica del concreto 
        a (float)   : Altura del bloque de compresión 
        b (float)   : Ancho de la columna 
        f_s (float) : (f's) Esfuerzo de acero a compresión
        A_S (float) : (A's) Área de acero a compresión 
        AS (float)  : Área de acero a tracción
        Fs (float)  : Esfuerzo de acero a tracción
                
    Returns:
        Pn (float): resistencia axial Kn
    """
    Pn0 = (r * f_c * a * b + f_s * (A_s * 0.0001) - (As * 0.0001) * fs) * 1000

    return Pn0

def def_unit_AS(d, c):
    """ Deformación unitaria del acero (AS)

    Args:
        d (float)  : distancia más comprimida de la fibra al eje de aceros a tracción
        c (float)  : Distancia de la cara comprimida del hormigón al eje neutro expresada en (m)
        
    Returns:
        def_unit_AS (float): deformación unitaria del acero.
    """
    _def = 0.003 * (d - c) / c
    
    return _def

def φ(_def,fy,Es):
    """factor de seguridad 

    Args:
        _def (float): Deformación unitaria del acero (AS)
        fy (float) : Resistencia a fluencia del acero expresada en (Mpa)
        Es (float) : modulo elástico del acero, es constante Es = 200'000 Mpa

    Returns:
        φ  (float): factor de seguridad.
    """
    if _def >= 0.005:
        φ = 0.90
    elif _def <= fy / Es:
        φ = 0.65
    elif fy / Es <= _def <= 0.005:
        φ = 0.65 + 0.25 * (_def - (fy / Es)) / (0.005 - (fy / Es))
    return φ

def φMn(Mn,φ):
    """ resistencia a momento con factor de seguridad.

    Args:
        Mn (float) : Resistencia nominal a momento en KNm 
        φ (float)  : factor de seguridad según norma 
        
    Returns:
        φMn (float): resistencia a momento con factor de seguridad.
    """
    φMn = Mn * φ
    
    return φMn

def φPn(φ, pn, pn_0_08):
    """Describir la función

    Args:
        φ (float)       : factor de seguridad
        Pn (float)      : resistencia axial 
        pn_0_08 (float) : factor de seguridad a compresión pura.
        
    Returns:
        fs (float): describir la variable.
    """
    if φ * pn > pn_0_08:
        φPn = pn_0_08
    elif φ * pn < pn_0_08:
        φPn = φ * pn
    return φPn

print()
f_c = 20 #(Mpa)
print("f´c = " + str(f_c) + " MPa")

fy = 410 #(Mpa)
print("fy = " + str(fy) + " MPa")

b = 0.4 #(m)
print("b = " + str(b) + " m")

h =	0.4	#(m)
print("h = " + str(h) + " m")

As = 8.04 #(cm^2)
print("As = " + str(As) + " cm^2")

A_S	= 8.04 #(cm^2)
print("A´s = " + str(A_S) + " cm´2")

recub_eje =	0.05 #(m)
print("recub_eje = " + str(recub_eje) + " m")

Es_mod_ellas_Acero = 200000	# Mpa
print("Es (mod ellas acero) = " + str(Es_mod_ellas_Acero) + " MPa")

d_ = recub_eje
print("d´ = " + str(d_))

d = h - recub_eje
print("d = " + str(d))

B1 = beta_1(f_c)
print("B1 = " + str(B1))

p0 = 3379.28 # ultimo valor de Pn (kN)
print("p0 = " + str(p0) + " KN")

P0_0_08 = p0 * 0.65 * 0.8
print("0.8*0*Po = " + str(P0_0_08) + " KNm")

print()

no_espacios = 300
paso = (3 * h) / no_espacios

print(f"c\t\t a\t \tf's\t\tfs\t\tdef_unit\tMn\t\tPn\t\tφ\t\tφMn\t\tφPn")

P = []
M = []
φP = []
φM = []
for i in range(no_espacios):
    C = (i + 1) * paso
    a = a_m(C, B1, h)
    
    _fs = fs_(Es_mod_ellas_Acero, C, d_, fy)
    
    _fs_= f_s(Es_mod_ellas_Acero, d, C,fy)
    
    def_unit = def_unit_AS(d, C)
    
    MN = momento(0.85, f_c, a, b, h, A_S, _fs, d_, As, _fs_, d)
    M.append(MN)
    Pn = presion0Pn(0.85, f_c, a, b, _fs, A_S, As, _fs_)
    P.append(Pn)
    _0_ = φ(def_unit, fy, Es_mod_ellas_Acero)
    
    _0_Mn = φMn(MN, _0_)
    φM.append(_0_Mn)
    _0_Pn = φPn(_0_, Pn, P0_0_08,)
    φP.append(_0_Pn)
    print(f"{C:10.7f}\t{a:10.7f}\t{_fs:+10.7f}\t{_fs_:+10.7f}\t{def_unit:10.7f}\t"
     f"{MN:10.7f}\t{Pn:10.7f}\t{_0_:10.7f}\t{_0_Mn:10.7f}\t{_0_Pn:10.7f}")

# for M, P in zip(φM, φP):
#      print(f"{M:.5f}\t{P:.5f}")

# Listas para almacenar los datos de A e B
A = []
B = []

# Leer el archivo CSV y extraer los datos
with open('file.csv', newline='') as csvfile:
    csvreader = csv.reader(csvfile)
    # next(csvreader)  # Saltar la primera fila (encabezados)
    for line in csvreader:
        A.append(float(line[0]))
        B.append(float(line[1]))
print(f"\tMn\t\tPn")

for i,j in zip(A, B):
    print(i, j)

# ingresar_puntos = A,B
# print(ingresar_puntos)
# sys.exit()

plt.figure()

plt.plot(φM, φP, label='GRAFICA DE DISEÑO')
plt.legend()
plt.xlabel('Mn (KNm)')
plt.ylabel('Pn (KN) ')
plt.title('CURVA DE ITERACION')
plt.grid()
plt.plot(M, P, label='GRAFICA NOMINAL')
plt.legend()

for punto in zip(A, B):
    plt.plot(*punto, 'o')

respuesta = True
for punto in ingresar_puntos:
    if not disenho(φM, φP, punto):
        respuesta = False
        break

print(respuesta)
plt.show()

