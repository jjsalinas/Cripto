"""
    Jose Antonio Jimenez Salinas
    Practcias CRIP - Cifrado Flujo
"""

#Ejercicio 1.
#Escribe una función que determine si una secuencia de bits cumple los postulados de Golomb
# En todo periodo la diferencia entre el nº de unos y ceros debe ser a lo sumo 1
# En un periodo el numero de rachas de longitud 1 debe ser el doble al numero de rachas de len 2 y este el doble de rachas de len 3 etc
# La distancia de hamming* entre dos secuencias diferentes, obtenidas mediante desplazamientos circulares de un periodo debe ser const
# *Dist de hamming. nº de bits que tiene que cambiarse para transformar una palabra valida en otra palabra valida. Si dos palabras
# difieren en una distancia d, se necesitan d errores para convertir una en la otra
def golomb(lista):
    n=len(lista)
    #Comprobar primer postulado
    if (lista.count(0)-lista.count(1)) in (-1, 0, 1):
        print ("Postulado #1: True")
        rachas={}
        rachas[1]= rachas[11]= rachas[111]= rachas[1111]= 0
        #Calculo de las rachas    
        rachas[1]=lista.count(1)
        for i in range(n):
            if lista[i:i+2]==[1,1]:
                rachas[11]=rachas[11]+1
            if lista[i:i+3]==[1,1,1]:
                rachas[111]=rachas[111]+1
            if lista[i:i+4]==[1,1,1,1]:
                rachas[1111]=rachas[1111]+1
        #Comprobar segundo postulado
        if rachas[1]>=len(rachas)/2 and (rachas[11]>0 and rachas[11]>=len(rachas)/4):
            if rachas[111]>0 or rachas[1111]>0:
                if rachas[111]>=len(rachas)/8 or rachas[111]>=len(rachas)/4:
                    print("Postulado #2: True")
                if rachas[1111]>=len(rachas)/16 or rachas[1111]>=len(rachas)/8 or rachas[1111]>=len(rachas)/4:
                    print("Postulado #2: True")
            else:
                print("Postulado #2: True")
            #Calculo de las distancias de hamming
            distancias=[]
            for i in range(1, n+1):
                laux=[e for e in lista[-i:]+lista[:-i]]
                d=0
                for i in range(n):
                    if lista[i]!=laux[i]:
                        d+=1
                distancias.append(d)
            #Comprobar tercer postulado
            if all(x==distancias[0] for x in distancias[0:len(distancias)-1]):
                print ("Postulado #3: True")
                return True
    return False

#Ejercicio 2.
#Implementa registros lineales de desplazamiento con retroalimentación (LFSR, [1, Chapter6]).
#La entrada son los coeficientes del polinomio de conexión, la semilla, y la longitud de la secuencia de salida.
#Ilustra con ejemplos la dependencia del periodo de la semilla en el caso de polinomios reducibles, la
#independencia en el caso de polinomios irreducibles, y la maximalidad del periodo en el caso de polinomios primitivos.
#Comprueba que los ejemplos con polinomios primitivos satisfacen los postulados de Golomb (en [1, 4.5.3] hay tablas de polinomios primitivos).
def LFSR(coef, semilla, lenSalida):
    l_coef=len(coef)
    l_semilla=len(semilla)
    if l_coef==1:
        return 1

    s=semilla[:]
    res=semilla[:]
    while(len(res)<lenSalida):
        sj=0
        for i in range(l_semilla):
            sj+=coef[i]*s[(len(res)-1)-i]
        sj=sj%2
        res.append(sj)
        s.append(sj)
    return res

#Ejercicio 3
#Escribe una función que toma como argumentos una función polinómica f , una semilla s y un entero positivo k,
#y devuelve una secuencia de longitud k generada al aplicar a s el registro no lineal de
#desplazamiento con retroalimentación asociado a f .
#Encuentra el periodo de la NLFSR ( (x ∧ y) ∨ z ) ⊕ t con semilla 1011.
#F(x,y,z,t)= z + xyz + t + 1
def NLFSR(polinomio, semilla, k):
    res=semilla[:]
    while len(res)<k:        
        l=len(res)-1
        l_func=len(polinomio)
        nvar=len(polinomio[0])
        aux=0
        for i in range(l_func):
            for j in range(nvar):
                aux+=res[l-(nvar-1)+j]*polinomio[i][j]
        aux=aux%2
        res.append(aux)
    return res
    
#Ejercicio 4
#Implementa el generador de Geffe ([1, 6.50]).
#Encuentra ejemplos donde el periodo de la salida es p 1 p 2 p 3 , con p 1 , p 2 y p 3
#los periodos de los tres LFSRs usados en el generador de Geffe.
#Usa este ejercicio para construir un cifrado en flujo. Con entrada un mensaje m, construye una llave k
#con la misma longitud que m, y devuelve m ⊕ k (donde ⊕ significa suma componente a componente en Z 2 ).
#El descifrado se hace de la misma forma: c ⊕ k (nótese que c ⊕ k = (m ⊕ k) ⊕ k = m ⊕ (k ⊕ k) = m, ya que x ⊕ x = 0 en Z 2 ).
def geffe(L1, L2, L3):
    res=[]
    for i in range(len(L1)):
        aux=L1[i]*L2[i]+L2[i]*L3[i]+L3[i]
        aux=aux%2
        res.append(aux)
    return res

def cifrado_geffe(m):
    lonm=len(m)
    L1=LFSR([0,0,1,1], [1,0,0,0], lonm)
    L2=LFSR([1,0,1,1], [1,1,0,1], lonm)
    L3=LFSR([0,1,0,1], [0,1,0,1], lonm)
    k=geffe(L1, L2, L3)
    res=[]
    for i in range(lonm):
        e=m[i]+k[i]
        e=e%2
        res.append(e)
    return res

def descifrado_geffe(c):
    lonc=len(c)
    L1=LFSR([0,0,1,1], [1,0,0,0], lonc)
    L2=LFSR([1,0,1,1], [1,1,0,1], lonc)
    L3=LFSR([0,1,0,1], [0,1,0,1], lonc)
    k=geffe(L1, L2, L3)
    res=[]
    for i in range(lonc):
        e=c[i]+k[i]
        e=e%2
        res.append(e)
    return res

#Ejercicio 5
#Dada una sucesión de bits periódica, determina la complejidad lineal de dicha sucesión, y el polinomio de conexión que la genera.
#Para esto, usa el algoritmo de Berlekamp-Massey ([1, Algorithm 6.30]).
#Haz ejemplos con sumas y productos de secuencias para ver qué ocurre con la complejidad lineal.
#Implementacion via Handbook of applied cryptography
def Berlekamp_Massey(secuencia):
    #s=secuencia[:]
    n=len(secuencia)
    c=[0 for i in range(n)]
    b=c[:]
    c[0]=1
    b[0]=1
    L=N=0
    m=-1
    while N<n:
        d=secuencia[N]
        for i in range (1, L+1):
            d+=c[i]*secuencia[N-i]
        d=d%2
        if d==1:
            t=c[:]
            j=0
            while(j+N-m)<n:
               c[(j+N-m)]^=b[j]
               j+=1
            
            if L<=N/2:
                L=N+1-L
                m=N
                b=t[:]
        N+=1
    return L

#Version apuntes Algoritmo Berlekamp-Massey en material adicional en swad
def Berlekamp_Massey_v2(secuencia):
    n=len(secuencia)
    for i in range(n):
        if secuencia[i]==0:
            k=i
            break
    
    f=[0 for i in range(n)]
    g=f[:]
    f[0]=1
    g[0]=1
    l=k+1
    a=k
    b=0
    r=k+1
    
    while r<n:
        d=0
        for i in range(l):
            d+=f[i]*secuencia[i+r-l]
        d=d%2
        if d==0:
            b+=1
        if d==1:
            if 2*l>r:
                for i in range(l):
                    f[i]=(f[i]+g[i+b-a])%2
                b+=1
            else:
                aux=f[:]
                for i in range(r+l-1):
                    f[i]=(aux[i+a-b]+g[i])%2
                l=r-l+1
                g=aux[:]
                a=b
                b=r-l+1
        r+=1
    #Devuelve +1 sobre la complejidad L real de la secuencia
    l-=1
    return l, f
                    
