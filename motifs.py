# -*- coding: utf-8 -*-
from typing import List, Dict


def calc_nucl_count(motifs: List[List[str]]) -> Dict[str, List[int]]:
    r"""Function Count(Motifs) that takes a list of strings Motifs as input
    and returns the count matrix of  Motifs (as a dictionary of lists)
    :param motifs:
    :return: dict
    """
    count = {}
    for nucl in ['A', 'T', 'G', 'C']:
        count[nucl] = [0] * len(motifs[0])
        for seq in motifs:
            for n, char in enumerate(seq):
                count[nucl][n] = count[nucl][n] + int(char == nucl)
    return count


def profile_motifs(motifs):
    r"""
    that takes Motifs as input and returns their profile matrix as a dictionary of lists
    :param motifs:
    :return: dict
    """
    count = calc_nucl_count(motifs)
    profile = {}
    for key, val in count.items():
        profile[key] = val / len(motifs)
    return profile
