#from sys import getsizeof
K = 37

def partitions(n, I=1):         # produzir partições
    yield (n,)
    for i in range(I, n//2 + 1,2):
        for p in partitions(n-i, i):
            yield (i,) + p

def verify(list_1,list_2):      # verificar se duas sequências são idênticas
    for l in range(len(list_1)):
        if(list_1[l] != list_2[l]):
            return False
    return True

# contador
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


####################

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

def substitute(necklace,value_array):
    final_array = []
    for i in range(len(necklace)):
        final_array.append(value_array[necklace[i]])
    return final_array


def necklace_it(array):
    value, bb = array_counter(array)
    n = len(array)
    aa = [0] * n
    lenght = len(bb)
    necklace_list = []
    SimpleFix(1,1,n,necklace_list,aa,bb,lenght)
    return necklace_list, value

def cycle(list1):       # fazer permutação cíclica
    k = list1[0]
    list1.pop(0)
    list1.append(k)
    return

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

basis = create_basis(K)
print(len(basis))


#memory_used_bytes = getsizeof(basis)
#memory_used_kb = memory_used_bytes / (1024)
#print(f'Memory used : {memory_used_kb:.2f} KB')

#for u in range(10, K+1):
#    basis = create_basis(u)
#    print(len(basis))


