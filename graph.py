from collections.abc import Iterable
import math


class DeBrujinGraph:

    # initialise la structure de données
    def __init__(self, str, k=21):
        self.noeuds = []
        self.nbNoeuds = 0

        for kmer in str:
            #On cree un objet Noeud pour chaque kmer
            self.unPack(kmer)

        self.Map = HashMap(self.nbNoeuds,k)
        for i in self.noeuds:
            #Et on ajoute tous les noeuds dans le HashMap
            self.Map.__setitem__(i.display(),i.display())


    #Methode qui cree un Noeud et l'ajoute dans la liste de noeuds
    def unPack(self, kmer):
        nouvNoeud = Noeud(kmer)
        self.add(nouvNoeud)
        self.nbNoeuds+=1

    #Méthode qui verifie si le graphe contient le noeud en question
    def __contains__(self, str) -> bool:
        if self.Map.__getitem__(str) is not None:
            return True
        else:
            return False

    # retourne un itérable sur les noeuds du graphe
    def __iter__(self):
        return self.nodes()

    # calcule le facteur de charge de la table de hachage sous-jacente
    def load_factor(self) -> float:
        return (self.Map._n / len(self.Map))

    #Ajoute le noeud au graphe
    def add(self, noeud):
        self.noeuds.append(noeud)

    #Enleve un noeud du graphe
    def remove(self, str):
        self.Map.__delitem__(str)

    #retourne un itérable sur les noeuds du graphe
    def nodes(self):
        return self.Map.display()

    #retourne tous les predecesseurs du noeud
    def predecessors(self, str):

        tempStr = str[:-1]
        predecesseurs = []
        for i in "ACGT":
            tempStr2 = i + tempStr
            if self.__contains__(tempStr2):
                predecesseurs.append(tempStr2)

        return predecesseurs

    #retourne tous les successeurs du noeud
    def successors(self, str):
        tempStr = str[1:]
        successeurs = []
        for i in "ACGT":
            tempStr2 = tempStr + i
            if self.__contains__(tempStr2):
                successeurs.append(tempStr2)
        return successeurs

    #parcours le graphe pour obtenir des segments contigus
    def kmer_walk(self):

        #On se fait une liste de tous les debuts de sequence
        noPredecessors = [i for i in self.noPredecessors()]

        #Au cas ou il y en ait en double
        startSet = set()
        for i in noPredecessors:
            startSet.add(i)
        start = [i for i in startSet]

        #Et on passe a travers chaque starter
        for i in start:
            successeursAVisiter = []
            successeursAVisiter.append(i)

            #Pendant qu'on a encore des successeurs non visites du start actuel
            while successeursAVisiter:
                alreadyVisited = False
                prochainKmer = successeursAVisiter.pop()
                contig = None
                closed = set()

                # Tant qu'on arrive pas deux fois au meme kmers
                while alreadyVisited is False:
                    kmerActuel = prochainKmer


                    if contig is None:
                        contig = kmerActuel
                    else:
                        contig += kmerActuel[-1]
                    #Si le kmer a deja ete visite
                    if kmerActuel in closed:
                        yield contig
                        alreadyVisited = True

                    else:
                        closed.add(kmerActuel)


                    #On cherche les successeurs du noeud actuel
                    successeurs = self.successors(kmerActuel)
                    #S'il y a plus qu'un successeur
                    if len(successeurs) > 1:
                        for i in range(len(successeurs)):
                            if i==0:
                                #On selectionne un pour le prochain
                                prochainKmer = successeurs[0]

                            else:
                                #Et on s'occupera des autres apres cette boucle while
                                successeursAVisiter.append(successeurs[i])

                    elif len(successeurs) == 1:
                        #Sinon on va directement au prochain successeur
                        prochainKmer = successeurs.pop()

                    #S'il n'y a plus rien a visiter
                    elif len(successeurs) == 0:
                        yield contig
                        #On passe au prochain
                        alreadyVisited = True



    #Cherche pour tous les noeuds qui ont aucun predecesseur
    def noPredecessors(self):
        for i in self.noeuds:
            if self.predecessors(i.display()) == []:
                yield i.display()



class Noeud:
    def __init__(self, string):
        self.kmer = string


    def display(self):
        return self.kmer





class HashMap:

    def __init__(self, l, k):
        self.k = k
        #on definit la taille intiale
        self.size = int(math.ceil((l - k + 1)))
        self.T = self.size * [None]
        self._n = 0     #le nombre d'elements dnas le HashMap



    # fonction de hachage, environ 20% de collision
    #reference a ce code: https://gist.github.com/mengzhuo/180cd6be8ba9e2743753
    def _hash_function(self, k):
        hash = 5381
        for x in k:
           hash = ((hash << 5) + hash) + ord(x)
        return (hash & 0xFFFFFFFF)%self.size

    # taille (nombre d'elements) dans la table
    def __len__(self):
        return self.size

    #methode pour verifier si un item est dans le hashmap
    def __getitem__(self, key):
        # on calcule l'index du bucket
        key_hash = self._hash_function(key)

        if self.T[key_hash] is not None:
            # on traverse le bucket
            for item in self.T[key_hash]:
                # s'il s'agit bien de l'item
                if item[0] == key:
                    # on le retourne
                    return item[1]
        return None

    # insertion d'un item
    def __setitem__(self, k, v):
        # on calcule l'index du bucket
        j = self._hash_function(k)
        item = [k, v]

        # si le bucket existe pas encore
        if self.T[j] is None:
            # on le cree et on ajoute l'item
            self.T[j] = list([item])
            self._n += 1
            self.resizeCheck()
            return
        else:
            # si l'index existe deja, on overright la valeur
            for i in self.T[j]:
                if i[0] == k:
                    i[1] = v
                    return

            # sinon on l'ajoute dans le bucket
            self.T[j].append(item)
            self._n += 1
            # on verifie a chaque fois si le tableau se rempli pas trop
            self.resizeCheck()
            return




    def resizeCheck(self):
        # si le nombre d'elements dpasse 75% de la taille de la table
        # on la double
        if self._n > len(self.T) * 0.75:
            self._resize(2 * len(self.T))

    # methode pour supprimer un item
    def __delitem__(self, key):
        # on calcule l'index
        key_hash = self._hash_function(key)
        # s'il y a rien a cet index
        if self.T[key_hash] is None:
            return False
        # sinon on cherche l'item dans le bucket
        for i in range(0, len(self.T[key_hash])):
            if self.T[key_hash][i][0] == key:
                self.T[key_hash].pop(i)
                return True
        return False

    # on change ta taille du hashmap
    def _resize(self, c):
        # c est la nouvelle taille
        newT = c * [None]
        self.size = c
        # on rehash tous les items
        for i in self.T:
            if i is not None:
                for u in i:
                    newT = self.setItemResize(newT,u[0],u[1])
        self.T = newT

    # fonction de setItem specifique pour le resizing
    # ( puisqu'on travaille pas sur le hashmap de vieille taille
    # mais sur un temporaire de nouvelle taille)
    def setItemResize(self,tempT, k, v):
        # on calcule l'index du bucket
        tempHash = self._hash_function(k)
        item = [k, v]

        # on insre l'ment de clk dans le bucket
        if tempT[tempHash] is None:
            tempT[tempHash] = list([item])
            return tempT
        else:
            for i in tempT[tempHash]:
                if i[0] == k:
                    i[1] = v
                    return tempT

            tempT[tempHash].append(item)
            return tempT

    # afficher tout les items non nulls du hashmap
    def display(self):
        tCopy = []
        for i in self.T:
            if i is not None:
                tCopy.append(i)
        return tCopy






