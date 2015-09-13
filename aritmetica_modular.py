#*************************************#
# Jose Ant. Jimenez Salinas           #
# Practicas CRIP - Aritmetica Modular #
#*************************************#
import numpy
def potencia2(n):
    pot=0
    while n%2==0:
        n=n//2
        pot+=1
    return pot
#Ejercicio 1.
#Implementa el algoritmo extendido de Euclides para el cálculo del máximo común divisor: dados
#dos enteros a y b, encuentra u, v ∈ Z tales que au + bv es el máximo común divisor de a y b
def euclides_extendido(a, b):
    x=0
    y=1
    u=1
    v=0
    while a != 0:
        cociente=b//a
        resto=b%a
        m=x-u*cociente
        n=y-v*cociente

        b=a
        a=resto
        x=u
        y=v
        u=m
        v=n
    mcd = b
    return mcd, x, y

#Ejercicio 2.
#Usando el ejercicio anterior, escribe una función que calcule a^−1 mod b para cualesquiera a, b
#enteros que sean primos relativos.
def inverso_modular(a, b):
    t=euclides_extendido(a,b)
    if t[0]==1:
        return t[1]%b
    else:
        print("No existe el inverso")

#Ejercicio 3.
#Escribe una función que calcule a^b mod n para cualesquiera a, b y n enteros positivos.
#La implementación debería tener en cuenta la representación binaria de b
def potencia_modular(a, b, n):
    temp=bin(b)
    temp=temp[2:]
    temp=temp[::-1]
    k=[int(i) for i in temp]
    #print("b=", k)
    b=1
    if all(i==0 for i in k):
        return b
    A=a
    for i in k:
        if i==1:
            b=(A*b)%n
        A=(A**2)%n
    return b

#Ejercicio 4.
#Dado un entero p, escribe una función para determinar si p es (probablemente) primo usando el
#método de Miller-Rabin
import random
def es_primo(p, iteraciones):
    if p==2:
        return True
    if p%2==0:
        return False
    else:
        esPrimo=True

        n=p-1
        potencia2=0
        while n%2==0:
            n=n//2
            potencia2+=1
        for i in range(iteraciones):
            a=random.randint(2, p-2)
            #x=pow(a, n, p)
            x=potencia_modular(a, n, p)
            if x!=1 and x!=p-1:
                for i in range(1, potencia2):
                    ant=x
                    x=potencia_modular(x, 2, p)
                    if x==1 and ant!=p-1:
                        esPrimo=False
                    elif x==p-1:
                        esPrimo=True
                #return False
                esPrimo=False
            else:
                #return True
                esPrimo=True
        return esPrimo

def test_primos():
    test_primos=[46381, 768479, 9476407, 36780481, 562390847, 1894083629, 65398261921, 364879542899, 8590365927553, 28564333765949, 123456789101119]
    for i in test_primos:
        print('\n', "Es primo ", i, " -: ")
        print(es_primo(i, 50))

#Ejercicio 5.
#Implementa el algoritmo paso enano-paso gigante para el cálculo de algoritmos discretos en Z p
import math
def pasoEnanoGigante(a, b, p):
    m=math.ceil(math.sqrt(p))
    lista=[]
    for i in range(m):
        lista.append([i, potencia_modular(a, i, p)])
    lista.sort(key=lambda x:x[1])
    inv_am=potencia_modular(inverso_modular(a, p), m, p)
    y=b%p
    listah=[]
    for i in range(m-1):
        aux=[i[1] for i in lista]
        if y in aux:
            j=aux.index(y)
            return (i*m+j)%p
        y=(y*inv_am)%p

#Ejercicio 6.
#Sea n = pq, con p y q enteros primos positivos.
#Escribe una función que, dado un entero a y un primo p con (p|a) = 1, devuelve r tal que r^2 ≡ a mod p
# primero te hará falta implementar el símbolo de Jacobi [1, 2.149]).
# Sea a un entero que es residuo cuadrático módulo p y q. Usa el teorema chino de los restos para calcular
# todas las raíces cuadradas de a mod n a partir de las raíces cuadradas de a módulo p y q.

#Simbolo de Jacobi
#Devuelve 0: si n divide a a, 1 si a es residuo cuadratico modulo n, -1 si a no es residuo cuadratico modulo n
def jacobi(a, n):
    #assert(n>2 and n%2==1)
    if a>n:
        a=a%n
    if a==0:
        return 0
    elif a==1:
        return 1
    if a%8==1 or a%8==7:
        return 1
    elif n%8==3 or n%8==5:
        return -1
    if a%2==0:
        return jacobi(a//2, n)
    else:
        return jacobi(n%a, a)

    """
    Mas simple
    Casos base a=0 return 0
               a=1 return 1
    Si a par comprobar modulo 8 y jacobi a//2, n
    Si no par  comprobar modulo 8 y dar la  vuelta: jacobi(n%a, a)

    """

#Función que, dado un entero a y un primo p con (p|a) = 1, devuelve r tal que r^2 ≡ a mod p
def ejer6_1(a, p):
    res=0
    if jacobi(a, p)==1:
        n=random.randint(1, p-1)
        while jacobi(n, p)!=-1:
            n+=1
        potencia=potencia2(p-1)
        aux=(p-1)//potencia_modular(2, potencia, p)
        if potencia==1:
            res=potencia_modular(a, (p+1)//4, p)
            return res
        elif potencia>=2:
            res=potencia_modular(a, (aux+1)//2, p)
            b=potencia_modular(n, aux, p)
            j=0
            c=(inverso_modular(a, p)*res**2)%p
            while j<=(potencia-2):
                if potencia_modular(c, pow(2, potencia-2-j), p)==p-1:
                    res=res*b%p
                b=(b**2)%p
                j+=1
            return (res, p-res)
    else:
        print("jacobi (", a, p, ")!=1")

#Teorema chino del resto
# @input:
#   n - numero ecuaciones
#   a - coeficientes de las ecuaciones
#   m - modulos de las ecuaciones
def teorema_chino_resto(n, a, m):
    if len(a) != len(m) or len(a)!=n or len(m)!=n:
        return "Error en los parametros dados"
    aux=m[:]
    for i in m:
        for j in m:
            if i!=j and euclides_extendido(i,j)[0]!=1:
                return "MCD de los modulos dados != 1"
    x=a[0]
    k=m[0]
    for i in range(1, n):
        ant_x=x
        inv_k=inverso_modular(k, m[i])
        x=((a[i]-ant_x)*inv_k)%m[i]
        r=ant_x+k*x
    return r

# Sea a un entero que es residuo cuadrático módulo p y q. Usa el teorema chino de los restos para calcular
# todas las raíces cuadradas de a mod n a partir de las raíces cuadradas de a módulo p y q.
def ejer6_2(a, p, q):
    c=[]
    n=2
    m=[p, q]

    r=ejer6_1(a, p)
    while type(r)!=int:
        r=ejer6_1(a, p)
    c.append(r)
    r=ejer6_1(a, q)
    while type(r)!=int:
        r=ejer6_1(a, q)
    c.append(r)
    res=teorema_chino_resto(n, c, m)
    return res

#Ejercicio 7.
# Implementa el Método de Fermat para factorización de enteros.
# Implementa el algoritmo de factorización ρ de Pollard
def fermat(p):
    if(p%2==0):
        exp2=potencia2(p)
        return (2, exp2, p//2**exp2)
    a=math.ceil(math.sqrt(p))
    #a=ejer6_1(0, p)+1
    b=int(math.sqrt(p))
    #b=ejer6_1(0, p)
    c=a*a-p
    while b*b != c:
        a=a+1
        c=a*a-p
        #Buscar Mejor Forma logar Raiz Cuadrada (Bottleneck aqui)
        #b=int(math.sqrt(c))
        b=int(numpy.sqrt(c))
        #b=ejer6_1(c, p)[0]
    r1=(a+b)
    r2=(a-b)
    return [r1, r2]

def pollard(p):
    if(p%2==0):
        exp2=potencia2(p)
        return (2, exp2, p//2**exp2)
    c=random.randint(0, 2000)
    m=random.randint(1,15000)
    f=lambda x:potencia_modular(x, 2, p)+c%m
    x=2
    y=2
    d=1
    while(d==1):
        x=f(x)
        y=f(f(y))
        d=euclides_extendido(abs(x-y),p)[0]
    if d==m:
        return None
    else:
        return (d, p//d)

#Ejercicio 8.
#Compara los tiempos de ejecución de tus implementaciones con las de tus compañeros y con las
#primitivas de algunos paquetes de cálculo simbólico como ( GAP , M ATHEMATICA , maxima , . . . ).
def comparativa():
    import time
    ini=time.time()
    for i in range(1000):
        #euclides_extendido(393, 267)
        #inverso_modular(391, 542)
        #potencia_modular(86, 72, 145)
        #for j in range(10): es_primo(123456789101119)
        #pasoEnanoGigante(9, 11, 19)
        #jacobi(5, 299)
        #ejer6_1(19, 31)
        #ejer6_2(19, 31, 5)
        #fermat(132)
        pollard(40259)
    t=time.time()-ini
    return t
#*************************************#
def gen_primo(a, b):
    x=random.randint(a, b)
    while es_primo(x,20)==False:
        x=random.randint(a, b)
    return x
#*************************************#
if __name__ == "__main__":
    """
    #Ejercicio 1 - Algoritmo Euclides Extendido
    #print('\nEjer 1. Euclides Extendido')
    #a=gen_primo(1, 100)
    #b=random.randint(1, 100)
    a=2
    b=0
    c=random.randint(1, 100)
    #print ("Euclides extendido (",a , b, ")=", euclides_extendido(a, b))
    #print ("Euclides extendido (",b , c, ")=", euclides_extendido(b, c))

    #Ejecuciones ejer 2 - Inverso Modular
    #print('\nEjer 2. Inverso Modular')
    p=gen_primo(100, 1000)
    val=random.randint(101, 999)
    #print("Inverso de ", val, " modulo p=", p, " ->", inverso_modular(val, p))

    #Ejecuciones ejer 3 - a^b mod n para cualesquiera a, b, n
    print('\nEjer 3. Potencia Modular')
    #p=gen_primo(100, 1000)
    #a=random.randint(1, 1000)
    #b=random.randint(1, 1000)
    a=134580987679025822862582825852825828
    b=110975366417606590833967199024636437460214709654546954443361638176229883140910302915147996598210506836614251268964446338737195911448213015377273363573074752120243818732263682864461580384738021460717565715722502
    p=110975366417606590833967199024636437460214709654546954443361638176229883140910302915147996598210506836614251268964446338737195911448213015377273363573074752120243818732263682864461580384738021460717565715722503
    print("a=", a, "^ b=", b, " mod p=", p, " resultado=", potencia_modular(a, b, p))
    #print("(a**b)%p =", (a**b)%p)

    #Ejercicio 4
    print('\nEjer 4. Test primos Miller-Rabin')
    #test_primos()
    print (es_primo(110975366417606590833967199024636437460214709654546954443361638176229883140910302915147996598210506836614251268964446338737195911448213015377273363573074752120243818732263682864461580384738021460717565715722505, 50))

    #Ejercicio 5
    print('\nEjer 5. Paso enano - Paso gigante')
    a=3
    b=6
    p=31
    print("Paso enano - paso gigante (", a, ",", b, ",", p, ")=", pasoEnanoGigante(a, b, p))

    #Ejercicio 6
    print('\nEjer 6.')
    #a=646
    #p=809
    #while jacobi(a, p)!=1:
    #    a=random.randint(1, 100)

    a=319
    p=gen_primo(1, 1000)
    r=ejer6_1(a, p)
    print("r tq r^2 ≡ ", a, " mod ", p, "  r=", r)
    if r!=None:
        if type(r)!=int:
            r=r[0]
        print("Comprobacion:", r, "**2", "-", a, "%", p, "=", ((r**2)-a)%p)

    #Prueba Teorema chino restos
    n=2
    a=[5495, 7569]
    m=[7643, 8765]
    print("Teorema chino restos. Coef:", a, " Mod:", m, "Resultado=", teorema_chino_resto(n, a , m))

    a=20
    p=11
    q=19
    r=ejer6_2(a, p , q)
    print("Raices cuadradas de a=", a, "en mod p=", p, "y en mod q=", q, '\nResultado=', r)
    print(r, "**2-", a, "%", p*q, "=", (r**2-a)%(p*q))
    """
    #Ejercicio 7
    print('\nEjer 7. Factorizacion Fermat y Pollard')
    #n=random.randint(1, 100000)
    p=gen_primo(1, 1000000)
    q=gen_primo(1, 100000)
    n=p*q
    if n%2==0:
        n+=1
    print("p:", p, "q:", q, "\np*q=n=", n)
    #print("Es primo:", es_primo(n, 20))
    print('Factorizacion Fermat: ', fermat(n))
    print('Factorizacion Pollard: ', pollard(n))
