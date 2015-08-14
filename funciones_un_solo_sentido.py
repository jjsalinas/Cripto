"""
    Jose Antonio Jimenez Salinas
    Practicas CRIP - Funciones de un solo sentido
"""
import random
import fractions
import aritmetica_modular_bis as am
import numpy as np
import time
import math
import hashlib as hl

"""
Ejercicio 1. Sea (a 1 , . . . , a k ) una secuencia super-creciente de números positivos (la suma de los términos P
que preceden a a i es menor que a i , para todo i ). Elige n > sum(a i) , y u un entero positivo tal que gcd(n, u) =1.
Define a i ∗ = u * ai mod n. La función mochila (knapsack) asociada a (a 1 ∗ , . . . , a k ∗ ) es
f : Z k 2 → N, f (x 1 , . . . , x k ) = sum (xi  ai* ). desde i=1 hasta k
"""
def greedy(val, lista):
    l=lista[:]
    suma=0
    indices=[]
    while suma!=val and len(l)>0:
        aux=max(l)
        if suma+aux<=val:
            indices.append(lista.index(aux))
            suma+=aux
        l.remove(aux)
    return indices

def knapsack(n, u, lista, m):
    assert(n>sum(lista))
    assert(am.euclides_extendido(n, u)[0]==1)
    nlista=[(a*u)%n for a in lista] #LLave publica
    c=sum([e*a for e,a in zip(m, nlista)])
    return nlista, c

def inv_knapsack(n, u, nlista, c):
    inverso=am.inverso_modular(u, n)
    l_orig=[(a*inverso)%n for a in nlista]
    val=(c*inverso)%n
    indices=greedy(val, l_orig)
    res=[0 for i in l_orig]
    for i in indices:
        res[i]=1
    return res

"""
Ejercicio 2. Sea p un (pseudo-)primo mayor o igual que vuestro número de identidad.
Encuentra un elemento primitivo, α, de Z ∗ p (se puede usar [3, 2.132 (iv)];
para facilitar el criterio, es bueno escoger p de forma que (p − 1)/2 sea también primo, y para ello usamos Miller-Rabin).
Definimos f : Z p → Z p , x → α x .
Calcula el inverso de tu fecha de nacimiento con el formato AAAAMMDD.
"""

#clave privada x
global clave_x
def funcion_ax(p):
    if am.es_primo(p, 30):
        q=(p-1)//2
        alfa=random.randint(1, q)
        temp=am.potencia_modular(alfa, q, p)
        while(temp==1):
            alfa=random.randint(2, p-1)
            temp=am.potencia_modular(alfa, q, p)
        temp=random.randint(1, q)
        while am.euclides_extendido(temp, p-1)[0]!=1:
            temp+=1
        alfa=am.potencia_modular(alfa, temp, p)
        x=random.randint(2, p-2)
        global clave_x
        clave_x=x
        y=am.potencia_modular(alfa, x, p)
        #Clave publica: (p, alfa, y)
        res=(p, alfa, y)
        return res
    else:
        return("P: ", p, " no es primo")

def cifra_ax(m, c_pub):
    p=c_pub[0]
    alfa=c_pub[1]
    y=c_pub[2]
    assert (m<p)
    k=random.randint(2, p-2)
    while am.euclides_extendido(k, p-1)[0]!=1:
        k+=1
    r=am.potencia_modular(alfa, k, p)
    s=am.potencia_modular(y, k, p)
    s=s*m
    return (r,s)

def descifra_ax(rs, c_pub):
    r=rs[0]
    s=rs[1]
    p=c_pub[0]
    global clave_x
    x=clave_x
    m=s*am.potencia_modular(r, (p-1-x), p)
    m=m%p
    return m

def test_gammal(primo, mensaje):
    clave_pub=funcion_ax(primo)
    print("Clave publica generada: ", clave_pub)
    rs=cifra_ax(mensaje, clave_pub)
    desc=descifra_ax(rs, clave_pub)
    print("Valor a cifrar: ", mensaje, '\nCifrado: ', rs)
    print ("Mensaje descifrado: ", desc)

#**************************************************#
"""
global alfa
#Ejer 2 V2
def calculo_alfa(p):
    q=(p-1)//2
    alfa=random.randint(2, p-2)
    while(am.potencia_modular(alfa, q, p)==1):
        alfa=random.randint(2, p-2)
    return alfa

def funcion_ax(alfa, val , p):
    res=am.potencia_modular(alfa, val, p)
    return res

def descifrar_ax(alfa, fax, inv, p):
    a=am.potencia_modular(alfa, inv, p)
    return funcion_ax(alfa, a, p)
    

def inversa_ax(alfa, val, p):
    inv=am.inverso_modular(val, p)
    return inv

def test_fax():
    val=19920608
    #p=gen_safeprime()
    p=3363260039
    alfa=calculo_alfa(p)
    fax=funcion_ax(alfa, val, p)
    invfax=inversa_ax(alfa, val, p)
    desc=descifrar_ax(alfa, fax, invfax, p)
    print("p= ", p, "val=", val, " fax=", fax)
    print("inverso para val: ", invfax)
    print("descifrado: ", desc)
    print(funcion_ax(alfa, desc, p))
"""

#**************************************************#
"""
Ejercicio 3.
n=p*q tq  p y q son primos
Sea f : Z n → Z n la función de Rabin: f (x) = x^2 . Sea n = 48478872564493742276963.
Sabemos que f (12) = 144 = f (37659670402359614687722).
Usando esta información, calcula p y q (mira lademostración de [1, Lemma 2.43])
"""
#Muy Rápido esta descomposición
def descomposicion_ejer3(n, f1, f2):
    r1=f1+f2
    if f1<f2:
        r2=f2-f1
    else:    
        r2=f1-f2
    mcd1=am.euclides_extendido(r1, n)[0]
    if mcd1!=1:
        return (mcd1, n//mcd1)
    else:
        mcd2=am.euclides_extendido(r2, n)[0]
        if mcd2!=1:
            return (mcd2, n//mcd2)
    if am.es_primo(n, 50):
        return (n, 1)
    
#Muy muy lento el encontrar raices iguales
def raices_iguales(n):
    f = lambda x:am.potencia_modular(x, 2, n)
    a=random.randint(2, n//2)
    fa=f(a)
    b=random.randint(n//2, n-1)
    while f(b)!=fa:
        b+=1
        b=b%n
    if b==a:
        return "ERROR: No existen a!=b tal que f(a)==f(b)", -1
    return (a,b)

def descomposicion_ejer3_v2(n):
    if am.es_primo(n, 500):
        return(n ,1)
    if n%2==0:
        return(n//2, 2)
    a,b=raices_iguales(n)
    if type(a)==str:
        return "ERROR: No existen a!=b tal que f(a)==f(b)"
    #print("a=", a, "fa=", fa, "b=", b, "fb=", fb)
    return descomposicion_ejer3(n, a, b)
#**************************************************#
"""
Ejercicio 4.
Elige a 0 y a 1 dos cuadrados arbitrarios módulo n (n como en el Ejercicio 3).
Sea h : Z 2 × (Z n ) ∗ → (Z n ) ∗ , h(b, x) = x 2 a 0 b a 1 1−b .
Usa la construcción de Merkle-Damg ̇ard para implementar una función resumen tomando h como función
de compresión (esta h fue definida por Goldwasser, Micali y Rivest). Los parámetros a 0 , a 1 y n se
hacen públicos (la función debería admitir un parámetro en el que venga especificado el vector inicial).
"""
"""
a0=1
a1=1
x=1


def funcion_hash(v, m):
    a0=v[0]
    a1=v[1]
    n=v[2]
    def h(b,x):
        res=x**2
        res=res*am.potencia_modular(a0, b, n)*am.potencia_modular(a1, 1-b, n)
        return res
    bloques=[]
    res=[]
    mb=bin(m)
    temp=[]
"""

#**************************************************#
"""
Ejercicio 5.
Sea p el menor primo entero mayor o igual que tu número de identidad, y sea q el primer
primo mayor o igual que tu fecha de nacimiento (AAAAMMDD).
Selecciona e tal que gcd(e, (p−1)(q−1)) =1.
Define la función RSA  f : Z n → Z n , x 7→ x e .
Calcula el inverso de 1234567890.
"""
def calculo_e(p, q):
    n=(p-1)*(q-1)
    e=random.randint(2, n-1)
    while am.euclides_extendido(e, n)[0]!=1:
        e=random.randint(2, n-1)
    return e

def calculo_d(p, q, e):
    n=(p-1)*(q-1)
    d=am.inverso_modular(e, n)
    return d

def cifrado_RSA(m, e, n):
    if am.euclides_extendido(m, n)[0]==1:
        c=am.potencia_modular(m, e, n)
        return c
    else:
        return "Mensaje m no es primo relativo con n"
    
def descifrado_RSA(c, d, n):
    inv=am.potencia_modular(c, d, n)
    return inv
#**************************************************#
"""
Ejercicio 6.
Sea n = 50000000385000000551, y que sabemos que una inversa de Z n → Z n , x 7→ x 5 es x 7→
x 10000000074000000101 (esto es, conoces tanto la llave pública como la privada de la función RSA).
Encuentra p y q usando el método explicado en [2, Page 92]. Compara este procedimiento con el algoritmo de Miller-
Rabin y el Ejercicio 3.
"""
def factoriza_n(e, d, n):
    assert(e!=0 and d!=0)
    factores=[]
    a=0
    de=d*e
    de=de-1
    a=am.potencia2(de)
    b=de//am.potencia_modular(2, a, n)
    x=random.randint(1, n-1)
    divisor=am.euclides_extendido(x, n)[0]
    if divisor!=1:
        return (divisor, n//divisor)
    
    y=am.potencia_modular(x, b, n)
    r1=(y-1)%n
    r2=(y+1)%n
    if r1==0 and r2==0:
        return "Fallo"
    while am.euclides_extendido(y, n)[0]==1:
        z=y
        y=am.potencia_modular(y, 2, n)
        r1=(y-1)%n
        r2=(y+1)%n
        if r2==0:
            return "Fallo"
        if r1==0:
            f1=am.euclides_extendido(z-1, n)[0]
            f2=am.euclides_extendido(z+1, n)[0]
            factores.append(f1)
            factores.append(f2)
            return factores
            break
        
    return factores
#**************************************************#
"""
Ejercicio 7.
En este ejercicio se pide implementar un sistema de firma digital y verificación de la firma.
Se puede elegir entre firma RSA o DSS.
Al igual que antes, debe realizar tres tareas: generación de claves (ejercicios anteriores), generación de
firma y verificación de firma.
Para la generación de la firma, se le introducirá un mensaje a cifrar (fichero) y el fichero con la clave
(privada), y deberá generar una firma, que se guardará en un fichero de texto.

Puesto que lo que realmente se firma no es el mensaje, sino un resumen del mensaje, hay que generar
un resumen de dicho mensaje. Para esto emplearemos la función SHA1 (se pueden añadir otras funciones
resumen). Cualquiera de las implementaciones de esta función que hay en la red puede ser usada.

Para la verificación de la firma, se introduce el mensaje (fichero) que se ha firmado, un fichero con la
firma (con el mismo formato que el generado en el apartado anterior) y un fichero con la clave (pública).
Deberá responder si la firma es o no válida.
"""
import hashlib
def genera_claves(fichero_cpublica, fichero_cprivada):
    try:
        publica=open(fichero_cpublica, 'r+')
        privada=open(fichero_cprivada, 'r+')
    except:
        print("Error abriendo los ficheros")
    p=gen_primo(1000, 10000000000)
    q=gen_primo(1000, 10000000000)
    n=p*q
    e=calculo_e(p, q)
    d=calculo_d(p, q, e)
    publica.write(n)
    publica.write(e)
    privada.write(d)
    

def firmado(fichero_m, fichero_cpublica, fichero_cprivada, fichero_salida):
    try:
        mensaje=open(fichero_m, 'r')
        cpublica=open(fichero_cpublica, 'r')
        cprivada=open(fichero_cprivada, 'r')
        f_salida=open(fichero_salida, 'r+')
    except:
        print("Error abriendo los ficheros")
    m=int(mensaje.read())
    resumen=hashlib.sha1(m)
    #RSA
    n=int(cpublica.readline())
    d=int(cprivada.read())
    firma=cifrado_RSA(resumen, d, n)
    f_salida.write(firma)
    mensaje.close()
    cpublica.close()
    cprivada.close()
    f_salida.close()
    
def verificar_firma(fichero_m, fichero_cpublica, fichero_firma):
    try:
        mensaje=open(fichero_m, 'r')
        cpublica=open(fichero_cpublica, 'r')
        f_firma=open(fichero_firma, 'r')
    except:
        print("Error abriendo los ficheros")
    m=int(mensaje.read())
    resumen=hashlib.sha1(m)
    #RSA
    n=int(cpublica.readline())
    e=int(cpublica.readline())
    firma=int(f_firma.read())
    verif=cifrado_RSA(resumen, e, n)
    if verif==firma:
        print("Firma valida")
        return True
    else:
        print("Firma Incorrecta")
        return False

#**************************************************#
def gen_primo(a, b):
    x=random.randint(a, b)
    while am.es_primo(x, 30)==False:
        x=random.randint(a, b)
    return x

def gen_safeprime(a, b):
    p=random.randint(a, b)
    q=(p-1)//2
    while (not am.es_primo(p, 30) or not am.es_primo(q, 30)):
        p=random.randint(a, b)
        q=(p-1)//2
    return p
#**************************************************#
if __name__ == "__main__":
    """
    #Ejemplo ejecuciones ejer 1: Knapsack
    print("Ejercicio 1: knapsack")
    l_numeros=[1, 3, 7, 15, 31, 63, 127, 255]
    n=557
    u=323
    #m=[0,1,1,0,0,1,0,1] #1228
    m=[1,1,1,1,0,0,0,0]
    llave, c= knapsack(n, u, l_numeros, m)
    print("Knapsack para el mensaje: ", m, '\nResultado mensaje cifrado: ', c, " Llave", llave)
    print("Descifrado:", inv_knapsack(n, u, llave, c))
    
    #Ejemplo ejecuciones ejercicio 2
    print("\nEjercicio 2: ElGammal")
    p=gen_safeprime(20078699, 10000000000)
    m=19920608
    test_gammal(p, m)
    
    """
    #Ejercicio 3
    print("\nEjercicio 3 - Descomposicion primos")
    #ejer3()
    #n=48478872564493742276963
    #f1=12
    #f2=37659670402359614687722
    p=gen_primo(10000, 100000)
    q=gen_primo(1000, 10000)
    n=p*q
    print("n=", n)
    #pq=descomposicion_ejer3(n, f1, f2)
    tini=time.time()
    pq=descomposicion_ejer3_v2(n)
    print("Tiempo= %.3f" % (time.time()-tini))
    print("p= ", pq[0], "q=", pq[1])
    #print("p, q = ", pq)

    """
    #Ejercicio 4
    print("\nEjercicio 4: hash")
    
    
    #Ejercicio 5
    print('\nEjercicio 5 - RSA')
    p=20078699
    q=19920608
    m=1234567890

    while(not am.es_primo(p, 30)):
        p+=1
    while (not am.es_primo(q, 30)):
        q+=1
    n=p*q
    print("N=", n)
    e=calculo_e(p, q)
    d=calculo_d(p, q, e)
    print ("e=", e, " d=", d)
    cif=cifrado_RSA(m, e, n)
    descif=descifrado_RSA(cif, d, n)
    print("Cifrado RSA para : ", m, "=", cif)
    print("Descifrado RSA: ", descif)
    
    #Ejercicio 6
    print('\nEjercicio 6')
    n=50000000385000000551
    e=5
    d=10000000074000000101
    factores=factoriza_n(e, d, n)
    print("n=", n, '\nFactores de n: p=', factores[0], " q=", factores[1])
    """
