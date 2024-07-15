import numpy as np
import pandas as pd
import itertools
import random
import igraph as ig
import math
from mpmath import nsum, exp, inf, fac


def make_flat(l): # A complementary function
    f=[x for s in l for x in s]
    return f

def clean_odd_hole(g): # A function that detects clean odd holes
    for v1, v2, v3 in itertools.combinations(g.vs,3):
        p12=g.get_shortest_paths(v1,to=v2)[0]
        p13=g.get_shortest_paths(v1,to=v3)[0]
        p23=g.get_shortest_paths(v2,to=v3)[0]
        if(len(set(p12+p13+p23))>=5 and len(set(p12+p13+p23))%2 == 1):
            f = g.induced_subgraph(list(set(p12+p13+p23)),implementation='auto')
            if(all(np.array(f.vs.degree())==2)):
                #print(set(p12+p13+p23))
                return True
    return False


def clean_hole_free_generation(n,p): # Clean odd hole free generation function that takes order and density as input and returns a igraph graph object
    
    def clean_cycle_check(g,v1,v2): # Function checks for clean odd holes
        vlist=list(g.vs)
        vlist.remove(v1)
        vlist.remove(v2)
        for v3 in vlist:
            p13=g.get_shortest_paths(v1,to=v3)[0]
            p23=g.get_shortest_paths(v2,to=v3)[0]
            if(len(set(p13+p23))>=5 and len(set(p13+p23))%2 == 1):
                indg=g.induced_subgraph(list(set(p13+p23)),implementation='auto')
                if(all(np.array(indg.vs.degree())==2)):
                    return True
        return False
    
    def advanced_procedure(g): # When the process stuck this function is called and make some changes to continue to the generation process
        if(random.uniform(0,1)>0.2):
            if(len(max(g.maximal_cliques(),key=len))>2):
                lar_cliques = []
                for i in g.maximal_cliques():
                    if(len(i)>2): lar_cliques.append(i)
                a_lar_clique = random.choice(lar_cliques)
                new_vertex=g.add_vertex()
                for i in a_lar_clique:
                    g.add_edge(new_vertex.index,i)
                g.delete_vertices(int(random.uniform(0,g.vcount()-1)))

            else:
                a_lar_clique = random.choice(g.maximal_cliques())
                new_vertex=g.add_vertex()
                for i in a_lar_clique:
                    g.add_edge(new_vertex.index,i)
                g.delete_vertices(int(random.uniform(0,g.vcount()-1)))
        
        else:
            a_lar_clique = random.choice(g.maximal_cliques())
            new_vertex=g.add_vertex()
            for i in a_lar_clique:
                g.add_edge(new_vertex.index,i)
            g.delete_vertices(int(random.uniform(0,g.vcount()-1)))
    
    g = ig.Graph(n)
    e=round(((n*(n-1))/2)*p)
    k=0
    while(g.ecount()<e):

        if(k>50):
            advanced_procedure(g)
            k=0
            continue

        number_of_edge = np.random.poisson(lam=((0.40*(g.ecount()/e))), size=None)+1
        s=[]
        for i in range(number_of_edge):
            s.append(g.vs[random.sample(range(g.vcount()),2)])
        
        for i in range(number_of_edge):
            g.add_edge(s[i][0],s[i][1])
                    
        edge_result=[]
        for i in range(number_of_edge):
            edge_result.append(clean_cycle_check(g,s[i][0],s[i][1]))
             
        if(any(edge_result)):  
            for i in range(number_of_edge):
                g.delete_edges([(s[i][0],s[i][1])])
            k=k+1
            continue
        
        k=0

        g.simplify()

    return g


def c5_free_generation(n,p): # C5 free generation function that takes order and density as input and returns a igraph graph object

    def eliminate_c5(g,c5_list):
        two_vertex_combination=list(itertools.combinations(c5_list,2))
        random.shuffle(two_vertex_combination)
        for j in two_vertex_combination:
            if(g.are_connected(j[0],j[1])):
                g.delete_edges([(j[0],j[1])])
                if(c5_check_elim(g,j[0],j[1])):
                    g.add_edge(j[0],j[1])
                    continue
                else:
                    break
            else:
                g.add_edge(j[0],j[1])
                if(c5_check_elim(g,j[0],j[1])):
                    g.delete_edges([(j[0],j[1])])
                    continue
                else:
                    break
    
    def c5_check_elim(g,v1,v2):
        nhood=g.neighborhood([v1, v2], order=2)
        nhood[0].remove(v1)
        nhood[1].remove(v2)
        if(len(nhood[0])<len(nhood[1])):
            v1list=nhood[0]
            v2list=nhood[1]
        else:
            v1list=nhood[1]
            v2list=nhood[0]
        
        for c1,c2 in itertools.combinations(v1list,2):
            for c3 in v2list:
                indg=g.induced_subgraph([v1,v2,c1,c2,c3],implementation='auto')
                if(all(np.array(indg.vs.degree())==2)):
                    return True
        return False
    
    def c5_check(g,v1,v2):
        nhood=make_flat(g.neighborhood([v1, v2], order=2))
        nhood.remove(v1)
        nhood.remove(v2)
        
        for c1,c2,c3 in itertools.combinations(nhood,3):
            indg=g.induced_subgraph([v1,v2,c1,c2,c3],implementation='auto')
            if(all(np.array(indg.vs.degree())==2)):
                return [v1,v2,c1,c2,c3]
        return False
    
    g = ig.Graph(n)
    e=round(((n*(n-1))/2)*p)
    k=list(itertools.combinations(g.vs,2))

    while(g.ecount()<e):
        random.shuffle(k)

        for i in k:
            if(not g.are_connected(i[0],i[1])):
                g.add_edge(i[0],i[1])
            else:
                continue

            c5_list=c5_check(g,i[0].index,i[1].index)
            if(c5_list!=False):
                eliminate_c5(g,c5_list)
                break
            else:
                break

        print(g.density())
        
    return g



def bell_number(n): # Function that calculates Bell number for a set of size n
    return int((float(nsum(lambda k: ((k**n)/fac(k)), [0, inf])))/math.e)

def random_partition(vertex_list): # Function that uniformly samples from all possible partitions of a set, used during random GSG generation
    n=len(vertex_list)
    cum_list=[]
    csum=0
    for i in range(1,n+1): 
        k=i
        csum=csum+(k**n)/(math.factorial(k)*math.e*bell_number(n))
        cum_list.append(csum)
    
    rn=random.uniform(0, 1)
    K=0
    for i in range(len(cum_list)):
        if(cum_list[i]>rn):
            K=i+1
            break
    
    color=[]
    for i in range(n):
        color.append(int(random.uniform(1, K+1)))
        
    partition=[]
    for i in list(set(color)):
        same_color=[]
        for j in range(len(color)):
            if(color[j]==i):
                same_color.append(vertex_list[j])
        partition.append(same_color)
    return partition

def central_clique_size(n): # Function that randomly calculates the size of the central clique during random GSG generation
    L_nk=0
    for i in range(n+1):
        L_nk = L_nk + math.comb(n,i)*(2**(i*(n-i)))*bell_number(n-i)
    
    c_sum=0
    prob_list=[]
    for i in range(n):
        c_sum=c_sum + ((math.comb(n,i)*(2**(i*(n-i)))*bell_number(n-i))/L_nk)
        prob_list.append(c_sum)
        
    rn=random.uniform(0, 1)
    k=0
    for i in range(len(prob_list)):
        if(prob_list[i]>rn):
            k=i
            break

    return k


def gsg_generator(n,is_counipolar): # Function generates a random GSG that takes order and a binary number as input and returns igraph graph object
    g = ig.Graph(n)
    k=central_clique_size(n)
    partition=random_partition(list(np.arange(k,n)))
    for i in range(k):
        for j in range(k):
            if(i!=j):
                if(not g.are_connected(i,j)):
                    g.add_edge(i,j)
        
    for i in partition:
        for j in i:
            for l in i:
                if(j!=l):
                    if(not g.are_connected(j,l)):
                        g.add_edge(j,l)
    
    flat_partition=make_flat(partition)
    
    for i in np.arange(k,n):
        for j in range(k):
            if(0.5<random.uniform(0, 1)):
                g.add_edge(i,j)
    
    g.simplify()
    if (is_counipolar==True):
        h=g.complementer(loops=False)
        return h
    else:
        return g
    

def gsg_modifier(g): # Function that make perturbations on GSG, takes igraph graph object as input and make perturbations on location
    def c5_check(g,v1,v2):
        nhood=g.neighborhood([v1, v2], order=2)
        nhood[0].remove(v1.index)
        nhood[1].remove(v2.index)
        if(len(nhood[0])<len(nhood[1])):
            v1list=nhood[0]
            v2list=nhood[1]
        else:
            v1list=nhood[1]
            v2list=nhood[0]
        
        for c1,c2 in itertools.combinations(v1list,2):
            for c3 in v2list:
                indg=g.induced_subgraph([v1,v2,c1,c2,c3],implementation='auto')
                if(all(np.array(indg.vs.degree())==2)):
                    return True
        return False
        
    def clean_cycle_check(g,v1,v2):
        vlist=list(g.vs)
        vlist.remove(v1)
        vlist.remove(v2)
        for v3 in vlist:
            p12=g.get_shortest_paths(v1,to=v2)[0]
            p13=g.get_shortest_paths(v1,to=v3)[0]
            p23=g.get_shortest_paths(v2,to=v3)[0]
            if(len(set(p12+p13+p23))>=5 and len(set(p12+p13+p23))%2 == 1):
                indg=g.induced_subgraph(list(set(p12+p13+p23)),implementation='auto')
                if(all(np.array(indg.vs.degree())==2)):
                    return True
        return False
    
    vertex_list=list(g.vs)
    number_of_iterations=0
    number_of_modifications=0
    
    while(number_of_modifications<250):
        number_of_iterations=number_of_iterations+1
        v1,v2=random.sample(vertex_list,2)
        
        if(g.are_connected(v1,v2)):
            g.delete_edges([(v1,v2)])
            #if(clean_cycle_check(g,v1,v2)):
                #g.add_edge(v1,v2)
                #break
            if(c5_check(g,v1,v2)):
                g.add_edge(v1,v2)
            else:
                number_of_modifications=number_of_modifications+1
        else:
            g.add_edge(v1,v2)
            #if(clean_cycle_check(g,v1,v2)):
                #g.delete_edges([(v1,v2)])
                #break
            if(c5_check(g,v1,v2)):
                g.delete_edges([(v1,v2)])
            else:
                number_of_modifications=number_of_modifications+1
        
    print("number_of_modifications: ",number_of_modifications)
    print("number_of_iterations: ",number_of_iterations)
    return True    


def clean_hole_repair(graph_list): # Function repairs non-perfect graphs, takes list of non-perfect igraph graph objects and returns list of repaired igraph graph objects
    def c5_check(g,v1,v2):
        nhood=g.neighborhood([v1, v2], order=2)
        nhood[0].remove(v1)
        nhood[1].remove(v2)
        if(len(nhood[0])<len(nhood[1])):
            v1list=nhood[0]
            v2list=nhood[1]
        else:
            v1list=nhood[1]
            v2list=nhood[0]
        
        for c1,c2 in itertools.combinations(v1list,2):
            for c3 in v2list:
                indg=g.induced_subgraph([v1,v2,c1,c2,c3],implementation='auto')
                if(all(np.array(indg.vs.degree())==2)):
                    return True
        return False
    
    def clean_odd_hole_counter(g):
        s=list()
        for v1, v2, v3 in itertools.combinations(g.vs,3):
            p12=g.get_shortest_paths(v1,to=v2)[0]
            p13=g.get_shortest_paths(v1,to=v3)[0]
            p23=g.get_shortest_paths(v2,to=v3)[0]
            if(len(set(p12+p13+p23))>=5 and len(set(p12+p13+p23))%2 == 1):
                f = g.induced_subgraph(list(set(p12+p13+p23)),implementation='auto')
                if(all(np.array(f.vs.degree())==2)):
                    s.append(list(set(p12+p13+p23)))
        return s
    
    s_list=[]
    for i in (range(len(graph_list))):
        unique_list=[]
        a=clean_odd_hole_counter(graph_list[i])
        for j in a:
            if(j in unique_list):
                continue
            else:
                unique_list.append(j)
        s_list.append(unique_list)

    len_list=[]
    for i in s_list:
        len_list.append(len(i))

    thold=100
    
    small_cycles=np.array(s_list)[np.array(len_list)<thold]
    small_cycle_graphs=np.array(graph_list)[np.array(len_list)<thold]

    for i in (range(len(small_cycle_graphs))):
        for l in small_cycles[i]:
            two_vertex_combination=list(itertools.combinations(l,2))
            random.shuffle(two_vertex_combination)
            for j in two_vertex_combination:
                if(small_cycle_graphs[i].are_connected(j[0],j[1])):
                    small_cycle_graphs[i].delete_edges([(j[0],j[1])])
                    if(c5_check(small_cycle_graphs[i],j[0],j[1])):
                        small_cycle_graphs[i].add_edge(j[0],j[1])
                        continue
                    else:
                        break
                else:
                    small_cycle_graphs[i].add_edge(j[0],j[1])
                    if(c5_check(small_cycle_graphs[i],j[0],j[1])):
                        small_cycle_graphs[i].delete_edges([(j[0],j[1])])
                        continue
                    else:
                        break
    return small_cycle_graphs