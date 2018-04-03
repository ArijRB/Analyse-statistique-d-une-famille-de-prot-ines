# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 20:04:24 2018

@author: amrou
"""
import numpy as np
import matplotlib.pyplot as plt
import math


"""######### A LA FIN DE CHAQUE PARTIE, UN ENSEMBLE DE TEST DES FONCTIONS COMMENTES, IL FAUT DECOMMENTER POUR TESTER ET TRACER LES COURBES###########"""""

#Initialisation des variables non changées
M=5643
L=48
q=21
alpha_acide=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-']  

#Fonction de lecture du fichier

def readFile(filepath):
    file=open(filepath,'r')
    Mat=[]
    for l in file:
        if l[0]!= '>' :
            list=[]
            for char in l:
                if char!= '\n':
                    list.append(char)
            Mat.append(list)
    return np.array(Mat)
    


  
#fonction qui retourne un dictionnaire qui pour un tuple (colonne,acide) associe le nombre d'occurence de l'acide dans cette colonne  
def nombre_occurences(matrice):
    dic={}
    #on initialise le dictionnaire
    for i in range(0,L):
        for char in alpha_acide:
            dic[(i,char)]=0
    
    for i in range(0,L):
        colonne=matrice[:,i]
        for char in colonne:
            dic[(i,char)]=dic[(i,char)]+1
    
    return dic
    

#fonction poids qui etant donné le dictionnaire des nb d'occurence retourne la matrice de poids
def poids(dicNbOcc):
    w=[] #aura pour taille q*L
    for i in range(0,L):
        temp=[]
        for char in alpha_acide:
            temp.append((dicNbOcc[(i,char)]+1.0)/(M+q))
        w.append(temp)
    return w



#print(len(w[0]))


#fonction qui renvoie une liste de toute les entropies de chaque colonne
#prend en parametre la liste des poids des colonne  de la matrice
def entropie(w):
    S=[] #contient toutes les entropies
    for i in range(0,L):
        somme=0
        wCol=[] #le poid d'une colonne on l'extrait de w 
        for j in range (0,q):
            wCol.append(w[i][j])#
        for k in wCol:
            somme=somme+(k*math.log(k,2))
        
        S.append(math.log(q,2)+somme)
    
    return S
    

##s_triee= sorted(s, key=lambda t:t[1],reverse=True)
##print(len(s_triee))
##s_triee.getValue()
##print(s_triee)
##s_triee=s.sort(reverse=True) #on trie la liste des entropie de manière décroissante pour afficher les 3 plus grandes
##il faut trouver les colonne associées





#retourne les acide amines des 3 position les plus conservés
def les3plusGrandeEntropie(w,s):
    s_sorted=sorted(s, reverse=True) #on trie d'abord la liste des entropie
    res=[] #liste que lon va renvoyer et qui va contenir le tuple (position,acide_amine)
    list=[] #servira pour recuperer les position avec les 3 plus grande entropie
    
    for i in range(0,3):
        list.append(s.index(s_sorted[i]))
    
        
    for j in list:
        wj=[]
        for k in range (0,q):
            wj.append(w[j][k])
    
        res.append((j,alpha_acide[wj.index(max(wj))]))
    
    return res
    





def fonctionFmodeleNul(wb,alphab):
    param=[]
    for i in range (0,q):
        somme=0
        for j in range(len(wb)):
            somme= somme+(wb[j][alpha_acide.index(alphab[i])])
        somme=somme*(1.0/L)
        param.append(somme)
    
    return param


# retourne la liste des vraisemblace 
def LogVraisemblance(nom,Wb):
     test=readFile(nom)
     liste=[]
     f0=fonctionFmodeleNul(Wb,alpha_acide)
     for i in range(0,len(test[0])-len(Wb)):
         
         summ=0
         
         for j in range(i,i+L):
             ind=alpha_acide.index(test[0][j])
             
             val=(Wb[j-i][ind]*1.0)/f0[ind] # application de la formule
             summ+=math.log(val,2)
         liste.append(summ)
     return liste
    
def moyenne_liste(L):
    s=0 # pour calcule la somme des éléments de la liste L
    for elm in L:
        s+=elm
    return s/len(L)


"""########################################################################"""
"""#################TEST Premiere partie###################################"""

##Lecture du fichier
#
#Mat=readFile("Dtrain.txt")  
#print(Mat)
#print(Mat[1])#Ca marche bien




#dicti=nombre_occurences(Mat)
#w= poids(dicti)
#print(w[0][20]) #test ω0(‘ − ‘) ≃ 0.31  marche bien
#s=entropie(w)
#print(len(s))
#print(s[0]) # le test S0 ≃ 1.85 marche

###Graphe pour l'entropie
##plt.plot(s)
##plt.xlabel('Position i')
##plt.ylabel('Entropie Si')
##plt.show()   




#Les acide amines des trois plus grande entropie avec leur position dans notre cas
#pl=les3plusGrandeEntropie(w,s)
#print (pl) # affiche que W P et G sont les plus conservé!!! affiche [(31, 'W'), (46, 'P'), (43, 'G')]



#Decommenter pour utiliser! 
#TEST de fonction vraisemblance
#Matr=readFile("Dtrain.txt")
#Dico=nombre_occurences(Matr)
#W=poids(Dico)
#liste=LogVraisemblance("test_seq.txt",W)
#
#print (liste[0]) # affiche la meme valeur que le test -115.7964212127327 a peu près -116

###Graphe pour la vraisemblance
#plt.plot(liste)
#plt.xlabel('')
#plt.ylabel('LogVraissemblance')
#plt.show()



###########FIN PARTIE 1#########################





"""##########################################################################"""



"""##############################2EME PARTIE######################################"""


### a utiliser avec le poids
## fonction qui calcule le nombre d’occurence nij (a, b)
def nombre_occurence2(mat):
    
    dic={}
    
    for i in range (1,L):
        for j in range(i):
            c1=mat[:,j]
            c2=mat[:,i]
            
            for k in range(len(c1)):
                if(j,i,c1[k],c2[k]) in dic: # si combinaison deja presente on incrémente son occurence
                    dic[(j,i,c1[k],c2[k])] +=1
                else: # sinon on l'ajoute au dictionnaire
                    dic[(j,i,c1[k],c2[k])] =1
    return dic
    
    
    


## fonction qui calcule le poids    
def poid2(Nij):
    Pij= {} ##dictionnair de poid pour chaque quadruplet
    for i in range (1,L):  #combi
        for j in range(i):
            for char1 in alpha_acide: #toutes les lettres
                for char2 in alpha_acide:
                    if (j,i,char1,char2) in Nij: # si on a deja la combinaison dans le dictionnaire
                        Pij[(j,i,char1,char2)]=(Nij[(j,i,char1,char2)]+1.0/q)/(M+q)
                    else:  #sinon
                        Pij[(j,i,char1,char2)]=(1.0/q)/(M+q)
    return Pij
    
    



#fonction qui calcule l'information mutuelle selon la formule donnée
# renvoie une mat de tuple (triplet)
#Wij poids 2eme partie, W poid 1ere partie
def info_mutuelle(Wij,W):
    mutu=[]
    liste=[]
    for i in range(1,L):
        l=[]
        for j in range(i):
            
            somme=0
            for char1 in alpha_acide:
                for char2 in alpha_acide:
                    quoti=Wij[(j,i,char1,char2)]/(W[j][alpha_acide.index(char1)]*W[i][alpha_acide.index(char2)])
                    somme=somme + Wij[(j,i,char1,char2)]*math.log(quoti,2)
            if (j,i) not in liste:
                liste.append((j,i))
                l.append((j,i,somme))
        
        mutu.append(l)
    return mutu

    
    

#fonction qui va calculer les 50 paire de dist et les fraction qui vont avec
def Distances(Mij,file):
    TuplesMij=[]    
    for lignes in range(1,L):
        for ind in range(0,len(Mij[lignes-1])):
            TuplesMij.append(Mij[lignes-1][ind])
    TuplesMij=sorted(TuplesMij, key=lambda L: L[2],reverse=True)
    TuplesMij=TuplesMij[0:50]
    Listes50=[]
    for i in TuplesMij:
        Listes50.append((i[0],i[1]))
    fic=open(file,'r')
    ListesD=[]    
    for ligne in fic:
        
        temp=ligne[0:len(ligne)].split()
        if float(temp[2])<8:
            ListesD.append((int(temp[0]),int(temp[1])))
     
    fic.close()
    vect=[]
    cpt=0
    for i in range(0,50):
        
        if Listes50[i] in ListesD:
            cpt+=1
        vect.append((cpt*1.0)/(i+1))
    return vect    


"""########################################################################"""
## TEST POUR LA DEUXIEME PARTIE
# Decommenter pour utiliser
    
#Matr=readFile("Dtrain.txt") 
#Dico=nombre_occurences(Matr)
#W=poids(Dico)   
#Dico2=nombre_occurence2(Matr)
#Wij=poid2(Dico2)
#Mij=info_mutuelle(Wij,W)
#vect=Distances(Mij,"distances.txt")
#
#print(Mij[0]) ##affiche [(0, 1, 0.4040474999905438)] qui est bien 
#
#
###TRACER COURBE FRACTION
#plt.plot(vect)
#plt.xlabel('nombre de paires considerées')
#plt.ylabel('fract')
#plt.show()
































