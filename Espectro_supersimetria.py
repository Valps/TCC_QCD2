import numpy as np
import matplotlib.pyplot
import time

K_min = 10
K_max = 25

y = 1               # pi*(massa_adjunta**2)/(N*(g**2))

class tupla:      
    def __init__(self, n_index_list, coeff):
        self.n_index_list = n_index_list
        self.coeff = coeff

    def print_tupla(self):
        print("coef = ", self.coeff, ", upla = (", end = "")
        for j in len(self.n_index_list):
            print(self.n_index_list[j], ", ", end = "")
        print(")")

def kr(a,b):   # Kronecker
    if(a != b):
        return 0
    return 1

def inverse_coeff(n):
    return 1/n

def sum_coeff(n):
    m = 1
    sum = 0
    while(m < n):
        sum += 1/((n-m)**2)
        m += 2
    return sum

def A_ni(n1,n2,n3,n4):
    if(n1-n3 != 0):
        return (1/((n1-n3)**2)) - (1/((n1+n2)**2))
    else:
        return -(1/((n1+n2)**2))        # renormalização?

def B_ni(n1,n2,n3,n4):
    return (1/((n2+n3)**2)) - (1/((n1+n2)**2))

def sum_ni_1(n3,n4):
    sum = 0
    n1 = 1
    while(n1 < n3 + n4):
        n2 = n3 + n4 - n1
        sum += A_ni(n1,n2,n3,n4)
        n1 += 2
    return sum



def operating_two(coefficient, state_operator_list, length):
    i = 0
    sum = 0
    while(i < length):
        sum += coefficient(state_operator_list[i])
        i += 1

    return sum

def sum_core_1(state, q1, q2, tuple_list):      # B^dagger B^dagger B B
    p = len(state)
    
    n1 = 1
    mq = state[q1]
    mq1 = state[q2]

    while(n1 < mq + mq1):     #   n1 < n3 + n4
        n2 = mq + mq1 - n1
        coef = -2*A_ni(n1,n2,mq,mq1)
        if(coef != 0):
            n_index_list = []
            i = 0
            while(i < p):
                n_index_list.append( state[i] + kr(i,q1)*( n1 - state[i] ) + kr(i,q2)*( n2 - state[i] ) )
                i += 1
            tuple_list.append(tupla(n_index_list, coef))
        n1 += 2
    
    return 

def sum_core_2(state, q, tuple_list):       # B^dagger B^dagger B^dagger B
    p = len(state)

    n1 = 1
    n2 = 1
    mq = state[q]

    while(mq - n1 - n2 >= 1):
        while(mq - n1 - n2 >= 1):
            n3 = mq - n1 - n2
            coef = -2*B_ni(n1,n2,n3,mq)
            if(coef != 0):
                n_index_list = []
                i = 0

                while(i < p):
                    if(i != q):
                        n_index_list.append( state[i] )
                    else:
                        n_index_list.append( n1 )
                        n_index_list.append( n2 )
                        n_index_list.append( n3 )
                    
                    i += 1

                tuple_list.append(tupla(n_index_list, coef))

            n2 += 2
        n2 = 1          # resetar n2
        n1 += 2         # incrementar n1
    
    return 


def sum_core_3(state, q1, q2, q3, tuple_list, parity):
    p = len(state)
    
    mq = state[q1]
    mq1 = state[q2]
    mq2 = state[q3]

    n4 = mq + mq1 + mq2

    if(p % 2 != 0):     # número ímpar de operadores
        coef = -2*B_ni(mq,mq1,mq2,n4)
    else:               # número par de operadores
        coef = -2*parity*B_ni(mq,mq1,mq2,n4)

    if(coef != 0):
        n_index_list = []
        n_index_list.append( n4 )

        i = q3+1
        counter = 0

        while(counter < p-3):
            if(i >= p):
                i = 0
            n_index_list.append( state[i] )
            i += 1
            counter += 1
        
        tuple_list.append(tupla(n_index_list, coef))
    return


def verify(list_1,list_2):      # verificar se duas sequências são idênticas
    for l in range(len(list_1)):
        if(list_1[l] != list_2[l]):
            return False
    return True


def cycle(list1):       # fazer permutação cíclica
    k = list1[0]
    list1.pop(0)
    list1.append(k)
    return

# verificar se duas sequências são permutações cíclicas
# e caso positivo, indicar a paridade pós permutação de traço
def verify_cyclic(array_1, array_2):   
    is_cyclic = False
    result = -1    #  0: is_cyclic = False,   -1: is_cyclic = True & parity = -1,   1: is_cyclic = True & parity = 1
    cycles = 0
    concatene = array_1 + array_1
    p = len(array_1)
    for i in range(p):
        if(verify(concatene[i:i+p], array_2)):
            is_cyclic = True
            break
        cycles += 1

    parity = (-1)**( (len(array_1)+1)*cycles )

    result = is_cyclic*parity
    return result


def matriz_hamiltoniana_T(set_states_list):
    i = 0
    lenght = len(set_states_list)
    matrix_A = np.zeros((lenght,lenght))

    

    while(i < len(set_states_list)):
        state = set_states_list[i]

        tuple_list = []

        # Termo B^dagger B^dagger B B

        p = len(state)
        q = 0
        while(q < p-1):           
            sum_core_1(state, q, q+1, tuple_list)
            q += 1
        sum_core_1(state, p-1, 0, tuple_list)

        # termo B^dagger B^dagger B^dagger B

        q = 0
        while(q < p):           
            sum_core_2(state, q, tuple_list)
            q += 1

        # termo B^dagger B B B
        
        q = 0
        parity = 1
        while(q < p-2):           
            sum_core_3(state, q, q+1, q+2, tuple_list, parity)
            parity *= 1
            q += 1
        sum_core_3(state, p-2, p-1, 0, tuple_list, parity)
        sum_core_3(state, p-1, 0 , 1, tuple_list, (-1)*parity)

        # ordenar e identificar as tuplas
        
        for k in range(len(tuple_list)):
            number_op = len(tuple_list[k].n_index_list)
            if(number_op == 1):
                continue
            for j in range(len(set_states_list)):
                if(number_op == len(set_states_list[j])):
                    order_test_1 = (tuple_list[k].n_index_list).copy()
                    order_test_2 = (set_states_list[j]).copy()
                    order_test_1.sort()
                    order_test_2.sort()
                    if( not (verify(order_test_1, order_test_2)) ):
                        continue
                    
                    result = verify_cyclic(tuple_list[k].n_index_list, set_states_list[j])

                    if(abs(result)):    # se é uma permutação cíclica
                        matrix_A[j,i] += (tuple_list[k].coeff)*result
                        break


        # termo B^dagger B
        matrix_A[i,i] += 4*operating_two(sum_coeff, set_states_list[i], p)

        i += 1
    return matrix_A


def matriz_hamiltoniana_V(set_states_list):
    i = 0
    lenght = len(set_states_list)
    matrix_A = np.zeros((lenght,lenght))

    while(i < len(set_states_list)):
        state = set_states_list[i]
        p = len(state)
        matrix_A[i,i] += operating_two(inverse_coeff, state, p)
        i += 1

    return matrix_A

def partitions(n, I=1):
    yield (n,)
    for i in range(I, n//2 + 1,2):
        for p in partitions(n-i, i):
            yield (i,) + p

def substitute(necklace,value_array):
    final_array = []
    for i in range(len(necklace)):
        final_array.append(value_array[necklace[i]])
    return final_array

def SimpleFix(t, p,    n, necklace_list, aa, bb, lenght):
    if t > n:
        if n % p == 0:
            necklace_list.append(aa.copy())
    else:
        for j in range(aa[t - p - 1], lenght):
            if bb[j] > 0:
                aa[t - 1] = j
                bb[j] -= 1
                if j == aa[t-p-1]:
                    SimpleFix(t+1, p, n, necklace_list, aa, bb, lenght)
                else:
                    SimpleFix(t+1, t, n, necklace_list, aa, bb, lenght)
                bb[j] += 1

def array_counter(array):
    value_list = []
    count_list = []
    
    p = len(array)
    count_list.append(1)
    value_list.append(array[0])
    i = 1
    while(i < p):
        j = 0
        pj = len(value_list)
        match = False
        while(j < pj):
            if(array[i] == value_list[j]):
                count_list[j] += 1
                match = True
                break
            j += 1
        if (not (match)):
            count_list.append(1)
            value_list.append(array[i])
            pj += 1
        i += 1
    return value_list, count_list

def necklace_it(array):
    value, bb = array_counter(array)
    n = len(array)
    aa = [0] * n
    lenght = len(bb)
    necklace_list = []
    SimpleFix(1,1,n,necklace_list,aa,bb,lenght)
    return necklace_list, value

def create_basis(k):
    lista = list(partitions(k))
    lista2 = []
    m = len(lista)

    if (k % 2):     # K ímpar
        for j in range(1,m):
            if(len(lista[j]) % 2):
                lista2.append(lista[j])
    else:
        for j in range(1,m):
            if(len(lista[j]) % 2 == 0):
                lista2.append(lista[j])

    # lista de partições pré-pronta

    basis_list = []
    for i in range(len(lista2)):
        necklaces, values = necklace_it(lista2[i])
        p = len(necklaces)
        j = 0

        if (k % 2):     # K ímpar
            while(j < p):
                basis_list.append( substitute(necklaces[j], values) )
                j += 1
        else:           # K par
            while(j < p):
                temp = substitute(necklaces[j], values)
                if(len(set(temp)) == 1):
                    break
                
                # verificar estados com traço nulo
                match = False
                lenght = len(temp)
                half_lenght = lenght//2
                cycles_counter = 0
                temp_cycled = temp.copy()
                while(cycles_counter < lenght):
                    cycle(temp_cycled)
                    cycles_counter += 1
                    if( verify(temp_cycled, temp ) and (cycles_counter % 2 != 0) ):
                        match = True
                        break
                    if(cycles_counter > half_lenght):
                        break
                if(not (match)):
                    basis_list.append( temp )

                j += 1
    
    return basis_list





k = K_min
counter = 0

lightest_hadron = []
light_2_hadron = []
inverse_k_axis = []

inicio = time.time()
while(k <= K_max):
    if(k >= 24):
        print("Gerando base K =", k)
    basis = create_basis(k)
    if(k >= 17):
        print("Calculando matrizes para K =", k)
    
    Matriz_V = matriz_hamiltoniana_V(basis)
    Matriz_T = matriz_hamiltoniana_T(basis)
    Hamiltoniana = y*Matriz_V + Matriz_T
    
    eigenvalues = np.linalg.eigvals(Hamiltoniana)

    eigenvalues = k*eigenvalues

    eigenvalues.sort()
    lightest_hadron.append(eigenvalues[0])
    light_2_hadron.append(eigenvalues[1])
    inverse_k_axis.append(1/k)

    p = len(eigenvalues)
    x_axis = []
    for i in range(p):
        x_axis.append(1/k)

    if(counter < 2):
        if(k % 2 != 0):     # fermions
            matplotlib.pyplot.plot(x_axis, eigenvalues, marker='o',color='red', linestyle='', label = 'Férmions', markersize=3)
        else:               # bosons
            matplotlib.pyplot.plot(x_axis, eigenvalues, marker='o',color='blue', linestyle='', label = 'Bósons', markersize=3)
    else:
        if(k % 2 != 0):     # fermions
            matplotlib.pyplot.plot(x_axis, eigenvalues, marker='o',color='red', linestyle='', markersize=3)
        else:               # bosons
            matplotlib.pyplot.plot(x_axis, eigenvalues, marker='o',color='blue', linestyle='', markersize=3)

    counter += 1
    k += 1




# Extrapolação


if(K_min != K_max):
    n = 3            # número de parâmetros
    n1 = n - 1       # jmax

    def g(j,x):      # funções de parâmetro
        if(j == 0):
            return 1
        if(j == 1):
            return x
        if(j == 2):
            return x**2
        else:
            raise Exception("Não há funções suficientes.")

    m = len(lightest_hadron)      # número de pontos / amostras
    m1 = m - 1       # imax, kmax

    gjl = np.zeros((n,n))
    gammal = np.zeros(n)

    def gjl_matrix(gjl, m1, n1, xi):
        j = 0
        while(j <= n1):
            l = 0
            while(l <= n1):
                i = 0
                sum = 0
                while(i <= m1):
                    sum += g(j,xi[i])*g(l,xi[i])
                    i += 1
                gjl[j,l] = sum
                l += 1
            j += 1

    def gamma_column(gammal, m1, n1, xi, yi):
        #gammal = np.zeros(n)
        l = 0
        while(l <= n1):
            sum = 0
            i = 0
            while(i <= m1):
                sum += yi[i]*g(l,xi[i])
                i += 1
            gammal[l] = sum
            l += 1


    # funcao de aproximação
    # f(x) = c0*g(0,x) + c1*g(1,x) + c2*g(2,x) + ... + c[n-1]*g(n-1,x)
    def f(x,n1,cj):
        j = 0
        sum = 0
        while(j <= n1):
            sum += cj[j]*g(j,x)
            j += 1
        return sum

    # lightest hadron states

    gjl_matrix(gjl, m1, n1, inverse_k_axis)
    gamma_column(gammal, m1, n1, inverse_k_axis, lightest_hadron)
    cj = np.linalg.solve(gjl, gammal) 

    # 2nd lightest hadron states
    gammal = np.zeros(n)
    gamma_column(gammal, m1, n1, inverse_k_axis, light_2_hadron)
    dj = np.linalg.solve(gjl, gammal)






    # intervalo de plottagem
    a = 0
    b = 0.14
    pontos = 200
    dx = (b-a)/pontos

    # eixos
    ex = []
    ey_lightest = []
    ey_lightest_2 = []

    x = a
    p = 0
    while(p <= pontos):
        ex.append(x)
        ey_lightest.append(f(x,n1,cj))
        ey_lightest_2.append(f(x,n1,dj))
        p += 1
        x += dx

    matplotlib.pyplot.plot(ex, ey_lightest, linestyle='--')
    matplotlib.pyplot.plot(ex, ey_lightest_2, linestyle='--')


    print("Extrapolação: ")
    print("Hadron mais leve: ", f(0,n1,cj))
    print("Segundo hadron mais leve: ", f(0,n1,dj))
    print(" ")
fim = time.time()
print("Tempo gasto:", fim - inicio, " segundos")

matplotlib.pyplot.title(label='Espectro para y = ' + str(y))
matplotlib.pyplot.xlabel("1/K")
#matplotlib.pyplot.ylabel("Z")
matplotlib.pyplot.xlim([0,0.11])
matplotlib.pyplot.ylim([0,113])
matplotlib.pyplot.legend()
matplotlib.pyplot.show()

