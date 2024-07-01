import textdistance as td
import nltk as nk
nk.download('punkt')

def hamming(s1, s2):
    return td.hamming.normalized_similarity(s1, s2)

def levenshtein(s1, s2):
    return td.levenshtein.normalized_similarity(s1, s2)

def jaro_winkler(s1, s2):
    return td.jaro_winkler(s1, s2)

def jaccard(s1, s2):
    set1 = set(s1.split(']'))
    set2 = set(s2.split(']'))
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union
    
def n_gram(s1, s2, n):
    ngrams = []
    for x in [s1, s2]:
        tokens = nk.word_tokenize(x)
        grams = list(nk.ngrams(tokens, n))
        ngrams.append(grams)
    unique_grams = set(ngrams[0]).union(set(ngrams[1]))
    shared_grams = len(ngrams[0]) + len(ngrams[1]) - len(list(unique_grams))
    return (shared_grams / unique_grams)