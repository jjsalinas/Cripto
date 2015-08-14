"""
    Jose Antonio Jimenez Salinas
    Practicas CRIP - Cifrado Flujo
    Ejemplos ejecuciones
"""

from cifrado_flujo import *

if __name__ == "__main__":
    #Ejer 1
    #Ejemplos secuencias para postulados de Golomb
    #No cumplen los 3 postulados:
    #l=[1,0,1,0,1,1,0,0]
    #l=[0,1,0,1,0,0,1,1,1,0,1,0,1,1,0,0,1,0,1,1,0]
    #Si cumplen los 3 postulados
    l=[1,1,1,1,0,1,0,1,1,0,0,1,0,0,0]
    #Test sobre salidas de LFSR ejer 2:
    #Polinomio primitivos
    # x^4 + x + 1 , Primitivo, semilla=[1, 0, 0, 1] 
    #l=[1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1]
    # x^4 + x + 1 , Primitivo, semilla[0, 0, 0, 1] 
    #l=[0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1]
    print ("Golomb ( ", l, " ) :", golomb(l))
    l=[1,0,1,0,1,1,0,0]
    print ("Golomb ( ", l, " ) :", golomb(l))

    #Ejer 2
    #Ejemplos LFSR con distintos polinomios y semillas
    #coeficientes=[0, 1, 0, 1] #x^4 + x^2 + 1
    #coeficientes=[0, 0, 1, 1] #x^4 + x^3 + 1
    #semilla=[1,0,0,1]
    #Polinomios primitivos:  periodo = 2^L - 1
    coeficientes=[1, 0, 0, 1] #x^4 + x + 1
    semilla=[0,0,0,1]
    #coeficientes=[0, 1, 0, 0, 1] #x^5 + x^2 + 1
    #semilla=[0,1,0,1,0]
    #Polinomio reducible
    #coeficientes=[1,1,1] #x^3 + x^2 + x + 1
    #semilla=[0,1,0]
    #Polinomio irreducible
    #coeficientes=[0,1,0,0,0,1] #x^6 + x^2 +1
    #semilla=[0,0,1,0,1,0]
    len_salida=64
    print('\nLFSR\nCoeficientes: ', coeficientes, '\nSemilla: ', semilla, '\nLongitud salida: ', len_salida)
    print(" --: ", LFSR(coeficientes, semilla, len_salida), '\n')

    #Ejer 3
    #Ejemplos NLFSR
    #F(x,y,z,t)=z + xyz + t + 1
    polinomio=[[0,0,1,0], [1,1,1,0], [0,0,0,1], [0,0,0,0]]
    semilla_nlfsr=[1,0,1,1]
    lonNLFSR=18
    print('\nNLFSR\nPolinomio: ', polinomio, '\nSemilla: ', semilla_nlfsr, '\nLongitud salida: ', lonNLFSR )
    print(" --: ", NLFSR(polinomio, semilla_nlfsr, lonNLFSR), '\n')
    
    #Ejer 4
    #Ejemplos Geffe
    m=[1,0,1,0,1,0,1,1,1,1,1,1,1]
    c=cifrado_geffe(m)
    minverso=descifrado_geffe(c)
    print('Geffe\nCifrando m: ', m, '\nYa cifrado: ', c, '\nDescifro y: ', minverso, '\n')
    
    #Ejer 5
    #Ejemplos obtener complejidad lineal de una secuencia via Berlekamp-Massey
    #secuencia=[1, 0, 0, 1, 1, 1, 1, 0] #L=4
    #secuencia=[0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0] # L=6
    #secuencia=[0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1] # L=5
    secuencia=[0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0] #L=6
    #sumando1, L=4, periodo=15
    sumando1=[0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1]
    #sumando2, L=5, periodo 31
    sumando2=[0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1]
    secuencia_suma=[(i+j)%2 for i,j in zip(sumando1, sumando2)] #L=L1+L2 periodo=mcm(L1,L2)
    secuencia_prod=[(i*j)%2 for i,j in zip(sumando1, sumando2)] #L=L1*L2
    print('\nCalculo complejidad Lineal\nSecuencia: ', secuencia)
    print('Su complejidad Lineal: L =', Berlekamp_Massey(secuencia), '\nSu complejidad Lineal (v2): L =', Berlekamp_Massey_v2(secuencia))
    print('\nComplejidad Lineal\nSecuencia (resultado Suma): \n', secuencia_suma)
    print('Su complejidad Lineal: L =', Berlekamp_Massey(secuencia_suma), '\nSu complejidad Lineal (v2): L =', Berlekamp_Massey_v2(secuencia_suma))
    print('\nComplejidad Lineal\nSecuencia (resultado Producto): ', secuencia_prod)
    print('Su complejidad Lineal: L =', Berlekamp_Massey(secuencia_prod), '\nSu complejidad Lineal (v2): L =', Berlekamp_Massey_v2(secuencia_prod))
    nlfsr=[1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1]
    print("Complejidad lineal del lnfsr ejer 3: ", Berlekamp_Massey_v2(nlfsr))
    
