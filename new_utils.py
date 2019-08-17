import random
from typing import List


def weighted_die(probabilities):
    kmer = ''  # output variable
    # your code here
    n = random.uniform(0, 1)
    for p in probabilities:
        n -= probabilities[p]
        if n <= 0:
            return p
    return kmer


def random_motifs(dna, k, t) -> List[str]:
    # random.randint(1, len(Dna))
    rand_motifs = []
    for i in range(t):
        r = random.randint(1, len(dna[0]) - k)  # 1 is not added as it is inclusive of last element also
        rand_motifs.append(dna[i][r: r + k])
    return rand_motifs


def count_with_pseudocounts(motifs):
    count = {}
    for nucl in ['A', 'T', 'G', 'C']:
        count[nucl] = [1] * len(motifs[0])
        for seq in motifs:
            for n, char in enumerate(seq):
                count[nucl][n] = count[nucl][n] + int(char == nucl)
    return count


def profile_with_pseudocounts(motifs):
    count = count_with_pseudocounts(motifs)
    profile = {nucl: [v / (len(motifs) + 4) for v in val] for nucl, val in count.items()}
    return profile


def pr(text, profile):
    # print(Profile)
    ans = 1
    for n, nucl in enumerate(text):
        ans = ans * profile[nucl][n]
    return ans


def profile_most_probable_kmer(text, k, profile):
    n_max = 0
    pos = 0
    for start in range(len(text) - k + 1):
        local = pr(text[start: start + k], profile)
        if local > n_max:
            n_max = local
            pos = start
    return ''.join(text[pos: pos + k])


def get_motifs(profile, dna):
    k = len(profile['A'])
    rv = []
    for sting in dna:
        rv.append(profile_most_probable_kmer(text=sting, k=k, profile=profile))
    return rv


def get_consensus(motifs: List[str]) -> str:
    consensus: List[str] = []
    profile = profile_with_pseudocounts(motifs)
    for n in range(len(motifs[0])):
        max_on_this_pos = max([profile[char][n] for char in ['A', 'T', 'G', 'C']])
        consensus.append(''.join([key for key, value in profile.items() if value[n] == max_on_this_pos][0]))
    return ''.join(consensus)


def get_score(motifs):
    score = 0
    consensus_motif = get_consensus(motifs)
    count = count_with_pseudocounts(motifs)
    for n in range(len(motifs[0])):
        score += sum([count[nucl][n] for nucl in ['A', 'T', 'G', 'C'] if nucl != consensus_motif[n]])
    return score


def randomized_motif_search(dna, k, t):
    rand_motif = random_motifs(dna, k, t)
    best_motifs = rand_motif
    while True:
        profile = profile_with_pseudocounts(rand_motif)
        rand_motif = get_motifs(profile, dna)
        if get_score(rand_motif) < get_score(best_motifs):
            best_motifs = rand_motif
        else:
            return best_motifs


def normalize(probabilities):
    coef = sum(list(probabilities.values()))
    for ch, val in probabilities.items():
        probabilities[ch] = val / coef
    return probabilities


def profile_generated_string(text, profile, k):
    n = len(text)
    probabilities = {}
    for i in range(0, n - k + 1):
        probabilities[text[i:i + k]] = pr(text[i:i + k], profile)
    probabilities = normalize(probabilities)
    return weighted_die(probabilities)


def GibbsSampler(dna, k, t, N):
    # Motifs = (Motif1, …, Motift) in each string from Dna
    motifs = random_motifs(dna, k, t)
    # BestMotifs ← Motifs
    best_motifs = motifs
    for j in range(1, N):
        # i ← randomly generated integer between 1 and t
        i = random.randint(0, t - 1)
        # Profile ← profile matrix formed from all strings in Motifs except for Motifi
        profile = profile_with_pseudocounts(list(set(motifs) - set(motifs[i])))
        motifs[i] = profile_generated_string(text=motifs[i], profile=profile, k=k)
        if get_score(motifs) < get_score(best_motifs):
            best_motifs = motifs
    return best_motifs
