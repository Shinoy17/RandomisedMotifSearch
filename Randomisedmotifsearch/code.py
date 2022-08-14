'''
RANDOMIZEDMOTIFSEARCH(Dna, k, t)
        randomly select k-mers Motifs = (Motif1, …, Motift) in each string
            from Dna
        BestMotifs ← Motifs
        while forever
            Profile ← Profile(Motifs)
            Motifs ← Motifs(Profile, Dna)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
            else
                return BestMotifs
'''
import random

#dna --- array of strings of DNA
#k --- length of k-merMotifs to be searched
def iterate_randomized_motif_search(dna, k):
    i = 0
    best_motifs = []
    for i in range(len(dna)):
        best_motifs.append(dna[i][:k])
    while i < 1000:
        motifs = randomized_motif_search(dna, k)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
        i = i+1
        #print(i)
    return best_motifs


def randomized_motif_search(dna, k):
    best_motifs = random_motifs(dna, k)
    while True:
        profile = create_profile(best_motifs)
        motifs = create_motifs(profile, dna)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs


# dna --- array of DNA strings
# profile -- profile matrix
# returns an array profile most probable kmers of each dna string
def create_motifs(profile, dna):
    motifs = []
    for i in dna:
        motifs.append(profile_most_probable(i, profile, len(profile)))
    return motifs


# calculate probability of a motif candidate (string) given the profile matrix(profile --- 2d array) 
def calculate_probablity(profile, string):
    probablity = 1
    for i in range(0, len(string)):
        if string[i] == 'A':
            probablity = probablity * profile[i][0]
        elif string[i] == 'C':
            probablity = probablity * profile[i][1]
        elif string[i] == 'G':
            probablity = probablity * profile[i][2]
        elif string[i] == 'T':
            probablity = probablity * profile[i][3]
    return probablity

# dna_string --- a single string of DNA
# profile --- profile matrix (2d array)
# k --- lenght of Kmer
# returns profile most probable kmer of the given DNA string according to given profile matrix
def profile_most_probable(dna_string, profile, k):
    best_pattern = dna_string[0:0 + k]
    best_probability = 0
    for i in range(len(dna_string) - k + 1):
        string = dna_string[i:i + k]
        new_probablity = calculate_probablity(profile, string)
        if new_probablity > best_probability:
            best_pattern = string
            best_probability = new_probablity
    return best_pattern


# motifs is the array of selected motifs
# creates profile matrix as an array of arrays
def create_profile(motifs):
    profile =[]
    for i in range(len(motifs[0])):
        # laplace smoothing for escaping zero probabilities
        count_A, count_C, count_G, count_T = 1, 1, 1, 1
        for motif in motifs:
            if motif[i] == 'A':
                count_A += 1
            elif motif[i] == 'C':
                count_C += 1
            elif motif[i] == 'G':
                count_G += 1
            elif motif[i] == 'T':
                count_T += 1
        # append an array of probabilities of A,C,G,T to profile array
        profile.append([count_A / (len(motifs) + 4), count_C/ (len(motifs) + 4),
                        count_G / (len(motifs) + 4), count_T / (len(motifs) + 4)])
    return profile


# calculates and returns the hamming distance between two strings of same length
def hamming_distance(str1, str2):
    counter = 0
    for s1, s2 in zip(str1, str2):
        if s1 != s2:
            counter += 1
    return counter


# motifs is the array of selected motifs
# retuns the consensus string of the given set of motifs 
def find_consensus(motifs):
    consensus = ''
    for i in range(len(motifs[0])):
        count_A, count_C, count_G, count_T = 0, 0, 0, 0
        for motif in motifs:
            if motif[i] == 'A':
                count_A += 1
            elif motif[i] == 'C':
                count_C += 1
            elif motif[i] == 'G':
                count_G += 1
            elif motif[i] == 'T':
                count_T += 1
        if count_A >= max(count_C, count_G, count_T):
            consensus += "A"
        elif count_C >= max(count_A, count_G, count_T):
            consensus += "C"
        elif count_G >= max(count_C, count_A, count_T):
            consensus += "G"
        elif count_T >= max(count_C, count_G, count_A):
            consensus += "T"
    return consensus


# motifs is the array of selected motifs
# returns the score of set of motifs
def score(motifs):
    consensus = find_consensus(motifs)
    score = 0
    for motif in motifs:
        score += hamming_distance(consensus, motif)
    return score


# dna is array of DNA strings
# selects a random motif from each dna string and appends it to an array then return the array 
def random_motifs(dna,k):
    motifs = []
    for i in dna:
        point = random.randint(0, len(i)-k)
        motifs.append(i[point:point + k])
    return motifs

if __name__ == "__main__":
    data = "".join(open('C:\\SEM-2 All Assignments\\Bio-2\\rosalind_ba2f (1).txt')).split()
    best_motifs = iterate_randomized_motif_search(data[2:], int(data[0]))
    for i in best_motifs:
        print(i)
        
        
        