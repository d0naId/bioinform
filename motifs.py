# -*- coding: utf-8 -*-
from typing import List, Dict
import random


def calc_nucl_count(motifs: List[List[str]]) -> Dict[str, List[int]]:
    r"""function count(motifs) that takes a list of strings motifs as input
    and returns the count matrix of  motifs (as a dictionary of lists)
    :param motifs:
    :return: dict
    """
    count = {}
    for nucl in ['a', 't', 'g', 'c']:
        count[nucl] = [0] * len(motifs[0])
        for seq in motifs:
            for n, char in enumerate(seq):
                count[nucl][n] = count[nucl][n] + int(char == nucl)
    return count


def profile_motifs(motifs):
    r"""
    that takes motifs as input and returns their profile matrix as a dictionary of lists
    :param motifs:
    :return: dict
    """
    count = calc_nucl_count(motifs)
    profile = {}
    for nucl, counts in count.items():
        profile[nucl] = [sum / len(motifs) for sum in counts]
    return profile





