# =====================================================================================
# drug_synergy_1.py  --  ANNOTATED COPY FOR LEARNING
# The code lines are IDENTICAL to drug_synergy_1.py; only comment lines were added.
# [PAPER] = how the line relates to Cheng, Kovacs & Barabasi, Nat Commun 2019;10:1197
# [FLAW n] = a problem confirmed by testing on your real network (details in the review doc)
#
# BIG PICTURE (what / why / how)
#   WHAT: for a disease and a pair of drugs, measure
#         (1) how close each drug's targets are to the disease genes on the protein network
#         (2) how close the two drugs' targets are to EACH OTHER
#         then label the combination (Overlapping, Complementary, Single, ...).
#   WHY : a combination whose drugs hit the disease module from different (non-overlapping)
#         network neighbourhoods is a candidate for complementary / synergistic action.
#   HOW : distances come from a pre-computed table (shortest_paths_cache.pkl). Each real
#         distance is compared with 1000 random target sets of similar connectivity.
#         z-score = how many standard deviations closer than random. Negative = closer.
# =====================================================================================
#imports
# [NOTE] SECTION 1 - IMPORTS: bring in tools written by others.
# `import X` loads library X. `import X as Y` gives it a short nickname Y.
# pickle = built-in tool to save/load Python objects (here: the huge distance table)
import pickle
# random = built-in random number tools (random.choice, random.seed)
import random
# logging = prints timestamped messages to screen AND file (better than print)
import logging
# pandas (pd) = tables called DataFrames; reads Excel / tab-separated files
import pandas as pd
# numpy (np) = fast maths on arrays (mean, std)
import numpy as np
# networkx (nx) = graph library; used to build the protein network and find neighbours
import networkx as nx
# typing = 'type hints' such as Set[str] (a set of strings). They only DOCUMENT what goes in/out;
# Python does not enforce them. Optional[int] means 'an int or None'.
from typing import Set, List, Tuple, Dict, Optional

# [NOTE] SECTION 2 - LOGGING SETUP (runs once when the file loads).
# basicConfig sets global rules:
#   level=INFO        -> show INFO, WARNING, ERROR (hide DEBUG)
#   format=...        -> each line looks like: time - LEVEL - message
#   handlers=[...]    -> two destinations: a file (drug_synergy_analysis.log, appended
#                        to each run) and the screen (StreamHandler)
#Logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('drug_synergy_analysis.log'),
        logging.StreamHandler()])
# logger = the object we call logger.info / .warning / .error on.
# __name__ is the module name ('__main__' when this file is run directly).
logger = logging.getLogger(__name__)

# [NOTE] SECTION 3 - THE CLASS.
# A class bundles DATA (self.shortest_paths) with FUNCTIONS that use it (called methods).
# You create one object from it later:  analyzer = NetworkProximityAnalyzer()
class NetworkProximityAnalyzer:
    # __init__ runs automatically when the object is created. `self` = the object itself.
    def __init__(self):
        # Start with no distance table; load_shortest_paths fills it. None = 'nothing yet'.
        self.shortest_paths = None
        
    #Network
    #Network
    # load_shortest_paths(self, cache_file: str) -> Dict
    # `cache_file: str` = hint that the file name is text; `-> Dict` = returns a dictionary.
    def load_shortest_paths(self, cache_file: str) -> Dict:
        """
        Load precomputed shortest paths from file
        """
        # f"...{x}..." is an f-string: {x} is replaced by the value of x.
        logger.info(f"Loading precomputed shortest paths from {cache_file}...")
        # try/except = 'attempt this; if it crashes, jump to except instead of dying immediately'.
        try:
            # open(file, 'rb') = open for Reading in Binary mode (pickle files are binary).
            # `with ... as f` closes the file automatically afterwards, even on error.
            with open(cache_file, 'rb') as f:
                # pickle.load rebuilds the saved dictionary: {gene_A: {gene_B: distance, gene_C: distance, ...}, ...}
                self.shortest_paths = pickle.load(f)
            # len(dict) = number of genes (keys) in the table
            logger.info(f"Loaded shortest paths for {len(self.shortest_paths)} nodes")
            # return = hand the result back to whoever called the method
            return self.shortest_paths
        # `except Exception as e` catches any error and names it e
        except Exception as e:
            logger.error(f"Failed to load shortest paths: {e}")
            # raise = re-throw the same error so the script stops (we log it first, but do not hide it)
            raise
    
    # get_network_nodes: returns every gene that has a row in the table.
    def get_network_nodes(self) -> Set[str]:
        """
        Get all nodes in the network from shortest paths data
        """
        # .keys() = the gene names; set(...) = unordered collection without duplicates (fast `in` tests)
        return set(self.shortest_paths.keys())
    
    # get_degree_binning: groups nodes by degree (number of connections) into bins of >= min_bin_size (100) nodes.
    # [PAPER] the reference distribution must use random nodes of 'matching size and degree distribution'.
    # This function is how the script tries to do that - but see the flaw on the degree-counting line below.
    def get_degree_binning(self, min_bin_size: int = 100) -> List[Tuple[int, int, List[str]]]:
        """
        Bin nodes by degree for degree-preserving randomization
        Uses the shortest paths data to calculate degrees
        """
        #Calculate degrees from shortest paths
        # degrees = {} -> empty dictionary {gene: degree}
        degrees = {}
        # `for node, connections in dict.items()` loops over key AND value together
        for node, connections in self.shortest_paths.items():
            # [FLAW 1] DEGREE IS COUNTED WRONGLY. The cache is built from unweighted hops (precompute uses
            # single_source_shortest_path_length), so the ONLY entry equal to 0 in a node's row is the node itself.
            # Result on your real network (14,301 nodes): every node gets degree 1 -> ONE bin holding all nodes ->
            # every node is a 'look-alike' of every other node. Measured: BRAF cluster mean true degree = 89, but the
            # random sets drawn have mean degree = 17. Hubs sit close to everything, so the comparison is unfair.
            # Minimal fix (valid because the cache stores hops): count entries with distance == 1 instead of 0.0
            # (the same change is needed in get_degree_equivalents). Not applied here - the code is unchanged.
            direct_neighbors = [n for n, dist in connections.items() if dist == 0.0]
            # degree = how many items the list above found
            degrees[node] = len(direct_neighbors)
        # degree_to_nodes = {degree: [genes with that degree]}
        degree_to_nodes = {}
        for node, degree in degrees.items():
            # setdefault(key, []) = 'if key is missing create an empty list there'; .append adds the gene to it
            degree_to_nodes.setdefault(degree, []).append(node)
        
        # sorted unique degree values, e.g. [1, 2, 3, ...]
        values = sorted(degree_to_nodes.keys())
        # bins = final result list, filled below
        bins = []
        # i = position pointer in `values`
        i = 0
        
        # OUTER loop: walks through the degree values from low to high
        while i < len(values):
            # low = first degree of the new bin
            low = values[i]
            # val = genes in this bin so far. .copy() so we don't change the original list.
            val = degree_to_nodes[values[i]].copy()
            # INNER loop: while bin is smaller than min_bin_size AND more degrees remain, absorb the next degree.
            # (`and` = both conditions must be true)
            while len(val) < min_bin_size and i < len(values) - 1:
                # move pointer to the next degree
                i += 1
                # extend adds all genes of that degree into this bin
                val.extend(degree_to_nodes[values[i]])
            # high = last degree included in this bin
            high = values[i]
            # move on to start the next bin
            i += 1
            # Leftover problem: the last bin may still be too small. If a previous bin exists...
            if len(val) < min_bin_size and bins:
                # unpack the previous bin (a tuple of 3) into three variables
                low_prev, high_prev, val_prev = bins[-1]
                # ...merge the leftover genes into the previous bin (keep its low, extend its high)
                bins[-1] = (low_prev, high, val_prev + val)
            # otherwise (bin is big enough, or it is the very first bin)...
            else:
                # ...store this bin as a tuple (low, high, genes)
                bins.append((low, high, val))
        # return the list of bins
        return bins
    
    # get_degree_equivalents: for each of YOUR genes, list the other genes that share its bin
    # (its 'look-alikes' in connectivity). These are the candidates used to build random sets.
    def get_degree_equivalents(self, nodes: Set[str], bins: List[Tuple[int, int, List[str]]]) -> Dict[str, List[str]]:
        """
        Find degree-equivalent nodes for randomization
        """
        # [FLAW 1] same wrong degree rule as in get_degree_binning (counts only the node itself).
        #Calculate degrees from shortest paths
        degrees = {}
        # for each gene in the input set
        for node in nodes:
            # only if the gene exists in the table
            if node in self.shortest_paths:
                direct_neighbors = [n for n, dist in self.shortest_paths[node].items() if dist == 0.0]
                degrees[node] = len(direct_neighbors)
            else:
                # unknown gene -> treat degree as 0
                degrees[node] = 0
        # node_to_equivalents = {gene: [look-alike genes]}
        node_to_equivalents = {}
        for node in nodes:
            d = degrees[node]
            # loop over the bins as 3 values: low, high, genes
            for low, high, bin_nodes in bins:
                # chained comparison: True when low <= d AND d <= high
                if low <= d <= high:
                    # equivalents = bin members except the gene itself
                    equivalents = [n for n in bin_nodes if n != node]
                    node_to_equivalents[node] = equivalents
                    # break = leave the `for` loop early (bin found, no need to check the rest)
                    break
        # return the dictionary
        return node_to_equivalents
    
    # get_random_nodes_matching_degree: builds n_random (1000) random replacements of a gene set.
    # Each real gene is swapped for a random look-alike. seed=None means different results every run;
    # a number would make it reproducible.
    def get_random_nodes_matching_degree(self, original_nodes: Set[str], n_random: int = 1000, min_bin_size: int = 100, seed: Optional[int] = None) -> List[Set[str]]:
        """
        Generate random node sets with degree distribution matching original nodes
        """
        # if a seed was given, fix the random generator so runs can be repeated exactly
        if seed is not None:
            random.seed(seed)
        # step 1: make the bins
        bins = self.get_degree_binning(min_bin_size)
        # step 2: find look-alikes for each original gene
        node_to_equivalents = self.get_degree_equivalents(original_nodes, bins)
        # random_sets = list that will hold the 1000 random sets
        random_sets = []
        # all genes in the network (used as fallback)
        network_nodes = self.get_network_nodes()
        # repeat n_random times; `_` = loop variable we do not need
        for _ in range(n_random):
            # start a new empty random set
            random_set = set()
            # go through each real gene
            for node in original_nodes:
                # if it has at least one look-alike...
                if node in node_to_equivalents and node_to_equivalents[node]:
                    # ...pick one at random (random.choice)
                    equivalent = random.choice(node_to_equivalents[node])
                    # set.add: sets ignore duplicates, so a random set can end up SMALLER than the original
                    random_set.add(equivalent)
                # else (no look-alike)...
                else:
                    #Fallback: random node from network
                    # ...fallback: any random gene from the whole network (list(...) is rebuilt every time = slow)
                    random_set.add(random.choice(list(network_nodes)))
            # store this random set
            random_sets.append(random_set)
        # return the list of 1000 sets
        return random_sets
    
    #Proximity
    #Proximity
    # [PAPER] Eq.(1): d(X,Y) = 1/||Y|| * SUM over DISEASE nodes y of ( min over DRUG nodes x of d(x,y) ).
    # [FLAW 2] This function averages over set_a. calculate_proximity calls it as (drug, disease), so it averages over
    # the DRUG nodes - the reverse of Eq.(1). To follow the paper, the call should be (disease, drug).
    def calculate_closest_distance(self, set_a: Set[str], set_b: Set[str]) -> float:
        """
        Calculate mean closest distance between two sets of nodes using precomputed paths
        """
        # distances = list of one minimum per gene in A
        distances = []
        for node_a in set_a:
            # skip genes missing from the table
            if node_a not in self.shortest_paths:
                # `continue` = jump straight to the next gene in the loop
                continue
                
            # float('inf') = infinity, a starting value larger than any real distance
            min_dist = float('inf')
            # scan every gene in B
            for node_b in set_b:
                # check that a distance A->B is recorded
                if node_b in self.shortest_paths[node_a]:
                    # look it up: shortest_paths[A][B] (dictionary inside a dictionary)
                    dist = self.shortest_paths[node_a][node_b]
                    # if smaller than the best so far...
                    if dist < min_dist:
                        # ...remember it (this is the 'running minimum' pattern)
                        min_dist = dist
            
            # if at least one B gene was reachable...
            if min_dist != float('inf'):
                # ...keep the minimum for this A gene
                distances.append(min_dist)
        # Conditional expression `X if condition else Y`: return the mean, or infinity if nothing was found
        return np.mean(distances) if distances else float('inf')
    
    # [PAPER] z = (d - mu) / sigma, with mu and sigma from random sets: matches the paper (random drug AND disease sets).
    # Note: the z-score is calculated here from one set of 1000 random draws; the p-value later uses a different draw.
    def calculate_proximity(self, drug_targets: Set[str], disease_genes: Set[str], n_random: int = 1000, min_bin_size: int = 100,seed: Optional[int] = None) -> Tuple[float, float, float, float]:
        """
        Calculate drug-disease proximity (z-score)
        Returns: (z_score, actual_distance, mean_random, std_random)
        """
        #Filter nodes that exist in network
        #Filter nodes that exist in network
        network_nodes = self.get_network_nodes()
        # `&` = set intersection: keep only genes that exist in the network.
        # !! This is where disease genes missing from STRING are silently dropped (main() only warns).
        drug_targets = drug_targets & network_nodes
        # same filter for disease genes
        disease_genes = disease_genes & network_nodes
        
        # `not X` is true for an empty set
        if not drug_targets or not disease_genes:
            # raise ValueError(...) = stop with an error message (caught in main() and written as an ERROR row)
            raise ValueError("No overlapping nodes between input sets and network")
        
        #Calculate actual distance
        #Calculate actual distance
        # the real drug->disease closest distance
        actual_distance = self.calculate_closest_distance(drug_targets, disease_genes)
        
        #Generate random sets
        #Generate random sets
        # 1000 random versions of the drug target set
        random_drug_sets = self.get_random_nodes_matching_degree(drug_targets, n_random, min_bin_size, seed)
        # 1000 random versions of the disease gene set
        random_disease_sets = self.get_random_nodes_matching_degree(disease_genes, n_random, min_bin_size, seed)
        
        #Calculate random distances
        #Calculate random distances
        random_distances = []
        # zip pairs the i-th random drug set with the i-th random disease set
        for rand_drugs, rand_diseases in zip(random_drug_sets, random_disease_sets):
            # distance for this random pair
            dist = self.calculate_closest_distance(rand_drugs, rand_diseases)
            # ignore infinite results
            if dist != float('inf'):
                random_distances.append(dist)
        
        #Calculate z-score
        #Calculate z-score
        # average of the 1000 random distances
        mean_random = np.mean(random_distances)
        # standard deviation of them (np.std = population SD)
        std_random = np.std(random_distances)
        
        # protect against division by zero (happens if all random distances are equal)
        if std_random == 0:
            # z = 0 means 'no information'
            z_score = 0.0
        # otherwise...
        else:
            # z = (real - random mean) / random SD. NEGATIVE = real distance shorter than chance = 'close'.
            z_score = (actual_distance - mean_random) / std_random
        # return 4 values at once (a tuple)
        return z_score, actual_distance, mean_random, std_random
    
    # [PAPER] Eq.(2): s_AB = <d_AB> - (<d_AA> + <d_BB>)/2.  Aligned. s_AB < 0 -> overlapping, >= 0 -> separated.
    def calculate_separation(self, set_a: Set[str], set_b: Set[str]) -> float:
        """
        Calculate separation score between two sets of nodes (s_AB)
        s_AB = <d_AB> - (<d_AA> + <d_BB>) / 2
        """
        #Calculate d_AA (within set A)
        #Calculate d_AA (within set A)
        d_aa = self.calculate_within_set_distance(set_a)
        
        #Calculate d_BB (within set B)
        #Calculate d_BB (within set B)
        d_bb = self.calculate_within_set_distance(set_b)
        
        #Calculate d_AB (between sets A and B)
        #Calculate d_AB (between sets A and B)
        d_ab = self.calculate_between_set_distance(set_a, set_b)
        # the formula; `2.0` forces decimal division
        separation = d_ab - (d_aa + d_bb) / 2.0
        return separation
    
    # [PAPER] d_AA = average over nodes of the nearest OTHER node in the same set (aligned).
    # A set with one node returns 0.0. The paper only used drugs with >= 2 known targets, so this case is outside it.
    def calculate_within_set_distance(self, node_set: Set[str]) -> float:
        """
        Calculate mean closest distance within a set of nodes
        """
        # a set with 0 or 1 gene has no internal distance -> 0.0
        if len(node_set) <= 1:
            return 0.0
        
        # distances = one minimum per gene
        distances = []
        # convert set to list so we can use positions i and j
        nodes_list = list(node_set)
        
        # enumerate gives (position, value) pairs: i, node_a
        for i, node_a in enumerate(nodes_list):
            if node_a not in self.shortest_paths:
                continue
                
            min_dist = float('inf')
            # second loop over the same list
            for j, node_b in enumerate(nodes_list):
                # `i != j` = don't compare a gene with itself (otherwise the minimum would always be 0)
                if i != j and node_b in self.shortest_paths[node_a]:
                    dist = self.shortest_paths[node_a][node_b]
                    if dist < min_dist:
                        min_dist = dist
            
            if min_dist != float('inf'):
                distances.append(min_dist)
        return np.mean(distances) if distances else float('inf')
    
    # [PAPER] Eq.(13) 'closest': pooled A->B and B->A nearest distances divided by |A|+|B| - the pooled mean here is the
    # same thing (aligned). Shared nodes have distance 0 by definition (also as in the paper).
    def calculate_between_set_distance(self, set_a: Set[str], set_b: Set[str]) -> float:
        """
        Calculate mean closest distance between two sets of nodes
        """
        distances_ab = []
        
        #From A to B
        #From A to B
        for node_a in set_a:
            if node_a not in self.shortest_paths:
                continue
                
            min_dist = float('inf')
            for node_b in set_b:
                if node_b in self.shortest_paths[node_a]:
                    dist = self.shortest_paths[node_a][node_b]
                    if dist < min_dist:
                        min_dist = dist
            
            if min_dist != float('inf'):
                distances_ab.append(min_dist)
        
        #From B to A
        #From B to A
        for node_b in set_b:
            if node_b not in self.shortest_paths:
                continue
                
            min_dist = float('inf')
            for node_a in set_a:
                if node_a in self.shortest_paths[node_b]:
                    dist = self.shortest_paths[node_b][node_a]
                    if dist < min_dist:
                        min_dist = dist
            
            if min_dist != float('inf'):
                distances_ab.append(min_dist)
        
        # mean of all pooled minima (infinity if none)
        return np.mean(distances_ab) if distances_ab else float('inf')
    
    #Significance
    #Significance
    # [PAPER] The paper does not use a p-value for drug-disease proximity; it uses the sign of z (z < 0 = close).
    # These p-values are computed but never change the label (see the classification block).
    def calculate_empirical_p_value(self, actual_distance: float, drug_targets: Set[str],  disease_genes: Set[str], n_random: int = 1000) -> float:
        """
        Calculate empirical p-value for drug-disease proximity
        """
        # random_distances = [] list to fill
        random_distances = []
        
        #Generate random distances
        # (note: new random sets are drawn here, separate from those used for the z-score)
        #Generate random distances
        random_drug_sets = self.get_random_nodes_matching_degree(drug_targets, n_random)
        random_disease_sets = self.get_random_nodes_matching_degree(disease_genes, n_random)
        
        for rand_drugs, rand_diseases in zip(random_drug_sets, random_disease_sets):
            dist = self.calculate_closest_distance(rand_drugs, rand_diseases)
            if dist != float('inf'):
                random_distances.append(dist)
        
        # p-value: (random_distances <= actual) makes a True/False array; np.mean of True/False
        # = the FRACTION of True. Example: 12 of 1000 random distances <= real -> p = 0.012.
        # (No '+1' correction, so p can be exactly 0.)
        p_value = np.mean(np.array(random_distances) <= actual_distance)

        # return the p-value
        return p_value
    
    # [PAPER] The paper explicitly says randomisation is NOT appropriate for s_AB because drugs have few targets
    # (non-Gaussian). This p-value, and the 'Significant ...' prefix built on it, is an addition by the script's author,
    # not part of the paper's method. The six classes themselves use only signs.
    def calculate_separation_p_value(self, set_a: Set[str], set_b: Set[str], actual_separation: float, n_random: int = 1000) -> float:
        """
        Calculate empirical p-value for drug-drug separation
        """
        # list of random s_AB values
        random_separations = []
        
        #Generate random separation scores
        # random replacements of drug A's targets
        random_sets_a = self.get_random_nodes_matching_degree(set_a, n_random)
        # random replacements of drug B's targets
        random_sets_b = self.get_random_nodes_matching_degree(set_b, n_random)
        
        # pair them up
        for rand_a, rand_b in zip(random_sets_a, random_sets_b):
            # s_AB for this random pair
            separation = self.calculate_separation(rand_a, rand_b)
            if separation != float('inf'):
                random_separations.append(separation)
        
        #For separation: more negative values are significant for overlap (s_AB < 0)
        #More positive values are significant for separation (s_AB > 0)
        #For separation: more negative values are significant for overlap (s_AB < 0)
        #More positive values are significant for separation (s_AB > 0)
        # one-sided test chosen by the sign of the real value:
        if actual_separation < 0:  #Overlap case
            # overlap case -> fraction of random values at least as negative
            p_value = np.mean(np.array(random_separations) <= actual_separation)
        # else (separation case)...
        else:  #Separation case
            # ...fraction of random values at least as positive
            p_value = np.mean(np.array(random_separations) >= actual_separation)
        
        # return the p-value
        return p_value
    
    
    #Main
    #Main
    # [PAPER] Six classes (Fig. 2): P1 Overlapping (z_A<0, z_B<0, sAB<0); P2 Complementary (z_A<0, z_B<0, sAB>=0);
    # P3 Indirect (one z<0, sAB<0); P4 Single (one z<0, sAB>=0); P5 Non-exposure (no z<0, sAB<0); P6 Independent (no z<0, sAB>=0).
    # The thresholds below (z < 0, s_AB < 0) match the paper.
    def classify_drug_combination_with_confidence(self, drug_a_targets: Set[str], drug_b_targets: Set[str], disease_genes: Set[str], n_random: int = 1000, z_threshold: float = 0.0, alpha: float = 0.05) -> Dict:
        """
        Classify drug combination with statistical confidence metrics
        Returns dictionary with classification, p-values, and confidence indicators
        """
        #Calculate drug-disease proximities with p-values
        # z-score for drug A and drug B (the `_` placeholders ignore the 3rd and 4th returned values)
        #Calculate drug-disease proximities with p-values
        z_da, actual_da, _, _ = self.calculate_proximity(drug_a_targets, disease_genes, n_random)
        z_db, actual_db, _, _ = self.calculate_proximity(drug_b_targets, disease_genes, n_random)
        
        #Calculate empirical p-values (this redraws a fresh random set: each drug-disease test is run twice)
        #Calculate empirical p-values for drug-disease proximity
        p_da = self.calculate_empirical_p_value(actual_da, drug_a_targets, disease_genes, n_random)
        p_db = self.calculate_empirical_p_value(actual_db, drug_b_targets, disease_genes, n_random)
        
        #Calculate drug-drug separation with p-value
        #Calculate drug-drug separation with p-value
        s_ab = self.calculate_separation(drug_a_targets, drug_b_targets)
        p_separation = self.calculate_separation_p_value(drug_a_targets, drug_b_targets, s_ab, n_random)
        
        #Determine significance (True/False)
        #Determine significance
        drug_a_significant = p_da < alpha
        drug_b_significant = p_db < alpha
        separation_significant = p_separation < alpha
        
        # etermine the three binary conditions (typo in the original: 'Determine')
        # etermine the three binary conditions
        # is drug A 'close' to the disease? (z below threshold)
        drug_a_hits_disease = z_da < z_threshold
        # is drug B 'close'?
        drug_b_hits_disease = z_db < z_threshold
        # do the drugs overlap? (s_AB negative)
        drugs_overlap = s_ab < 0

        # {z_da:.4f} = show 4 decimals;  {p_separation:.4e} = scientific notation
        logger.info(f"Proximity of drug A - Z score: {z_da:.4f}, pvalue: {p_da}")
        logger.info(f"Proximity of drug B - Z score: {z_db:.4f}, pvalue: {p_db}")
        logger.info(f"Separation between drug A & drug B - Sab score: {s_ab:.4f}, P value: {p_separation:.4e}")
        
        #Classify based on the combination of conditions
        # DECISION TABLE (A close?, B close?, overlap?):
        #   both close  + overlap -> Overlapping Exposure     both close  + separate -> Complementary Exposure
        #   one close   + overlap -> Indirect Exposure        one close   + separate -> Single Exposure
        #   none close  + overlap -> Non-exposure             none close  + separate -> Independent Action
        # 'Significant ...' is added only if p_separation < alpha, else '... but not significant'.
        # !! Observation: in every branch below, the `if drug_a_significant and drug_b_significant` part and
        # its `else` part assign IDENTICAL labels (only the log message differs), so the drug-disease
        # p-values never change the final label.
        #Classify based on the combination of conditions
        # CASE 1: both drugs close to the disease
        if drug_a_hits_disease and drug_b_hits_disease:
            # (both drug-disease p-values significant -> only changes the log text)
            if drug_a_significant and drug_b_significant:
                logger.info("Significant- Both drug targets in proximity with the disease module")
                # overlap or separate? then significant or not? -> picks one of 4 labels
                if drugs_overlap:
                    if separation_significant:
                        classification = "Significant Overlapping Exposure"
                    else:
                        classification = "Overlapping Exposure but not significant"
                else:
                    if separation_significant:
                        classification = "Significant Complementary Exposure"
                    else:
                        classification = "Complementary Exposure but not significant"
            # (else: not both significant) -> same four labels again, as noted above
            else:
                logger.info("Not Significant- Both drug targets in proximity with the disease module")
                if drugs_overlap:
                    if separation_significant:
                        classification = "Significant Overlapping Exposure"
                    else:
                        classification = "Overlapping Exposure but not significant"
                else:
                    if separation_significant:
                        classification = "Significant Complementary Exposure"
                    else:
                        classification = "Complementary Exposure but not significant"
        
        # CASE 2: only drug A close
        elif drug_a_hits_disease and not drug_b_hits_disease:
            if drug_a_significant and drug_b_significant:
                logger.info("Significant- Only drug target A in proximity with the disease module")
                if drugs_overlap:
                    if separation_significant:
                        classification = "Significant Indirect Exposure"
                    else:
                        classification = "Indirect Exposure but not significant"
                else:
                    if separation_significant:
                        classification = "Significant Single Exposure"
                    else:
                        classification = "Single Exposure but not significant"
            else:
                logger.info("Not Significant- Only drug target A in proximity with the disease module")
                if drugs_overlap:
                    if separation_significant:
                        classification = "Significant Indirect Exposure"
                    else:
                        classification = "Indirect Exposure but not significant"
                else:
                    if separation_significant:
                        classification = "Significant Single Exposure"
                    else:
                        classification = "Single Exposure but not significant"
  
        # CASE 3: only drug B close
        elif not drug_a_hits_disease and drug_b_hits_disease:
            if drug_a_significant and drug_b_significant:
                logger.info("Significant- Only drug target B in proximity with the disease module")
                if drugs_overlap:
                    if separation_significant:
                        classification = "Significant Indirect Exposure"
                    else:
                        classification = "Indirect Exposure but not significant"
                else:
                    if separation_significant:
                        classification = "Significant Single Exposure"
                    else:
                        classification = "Single Exposure but not significant"
            else:
                logger.info("Not Significant- Only drug target B in proximity with the disease module")
                if drugs_overlap:
                    if separation_significant:
                        classification = "Significant Indirect Exposure"
                    else:
                        classification = "Indirect Exposure but not significant"
                else:
                    if separation_significant:
                        classification = "Significant Single Exposure"
                    else:
                        classification = "Single Exposure but not significant"
        
        # CASE 4 (else): neither drug close
        else:
            if drug_a_significant and drug_b_significant:
                logger.info("Significant- None of the drug target in proximity with the disease module")
                if drugs_overlap:
                    if separation_significant:
                        classification = "Significant Non-exposure"
                    else:
                        classification = "Non-exposure but not significant"
                else:
                    if separation_significant:
                        classification = "Significant Independent Action"
                    else:
                        classification = "Independent Action but not significant"
            else:
                logger.info("Not Significant- None of the drug target in proximity with the disease module")
                if drugs_overlap:
                    if separation_significant:
                        classification = "Significant Non-exposure"
                    else:
                        classification = "Non-exposure but not significant"
                else:
                    if separation_significant:
                        classification = "Significant Independent Action"
                    else:
                        classification = "Independent Action but not significant"
        
        #Return all metrics in a dictionary: {name: value}. main() reads them by name, e.g. result['z_da']
        #Return all metrics in a dictionary
        return {
            "classification": classification,
            "z_da": z_da,
            "pvalue_da":p_da,
            "z_db": z_db,
            "pvalue_db":p_db,
            "s_ab": s_ab,
            "p_separation": p_separation,
            "drug_a_significant": drug_a_significant,
            "drug_b_significant": drug_b_significant,
            "separation_significant": separation_significant
        }


# SECTION 4 - VALIDATION HELPERS (they only WARN; they do not remove anything).
# validate_moduleprotein: counts disease genes that appear in neither the source nor the target column of the STRING table.
def validate_moduleprotein(module, string):
    not_string = []
    for i in module:
        # `not in` on a numpy array scans linearly = slow, but works. Names of missing genes are NOT logged, only the count.
        if i not in string.source.values and i not in string.target.values:
            not_string.append(i)
    if len(not_string) > 0:
        logger.warning(f"{len(not_string)} out of {len(module)} proteins in the disease module are not in string")

# validate_targetprotein: same check for drug targets
def validate_targetprotein(drug_targets, string):
    not_string = []
    for i in drug_targets:
        if i not in string.source.values and i not in string.target.values:
            not_string.append(i)
    if len(not_string) > 0:
        logger.warning(f"{len(not_string)} out of {len(drug_targets)} drug targets are not in string!")


# Main function
# SECTION 5 - main(): the script's workflow
# Main function
def main():
    # first log line
    logger.info("Starting drug synergy analysis")
    
    # try ... except at the very end catches any fatal error
    try:
        # create the analyzer object from the class
        analyzer = NetworkProximityAnalyzer()
        # load the distance table (the .pkl must be in the same folder)
        analyzer.load_shortest_paths("shortest_paths_cache.pkl")
        
        logger.info("Loading STRING PPI network")
        # read STRING_PPI_filtered.txt; sep='\t' = tab-separated
        string = pd.read_csv("STRING_PPI_filtered.txt", sep="\t")
        # [PAPER] The paper's interactome = experimentally supported physical PPIs (Y2H, AP-MS, literature, structures, kinase-
        # substrate, signalling); computationally inferred links were EXCLUDED. STRING min900 also includes text-mining,
        # co-expression etc. channels (about 4% of your edges have no experimental or database evidence). Ask your lead.
        G = nx.from_pandas_edgelist(string, source='source', target='target', edge_attr=True)
        
        logger.info("Loading drug target data")
        # read sheet 'Drug_combination' (Disease, Drug_target1, Drug_target2)
        combination = pd.read_excel('Drug_target.xlsx', sheet_name='Drug_combination')
        # read sheet 'Module_genes' (one column per disease)
        modules = pd.read_excel('Drug_target.xlsx', sheet_name='Module_genes')
        
        # Create a list to store results
        # results = list that collects one dictionary (= one output row) per combination
        results = []
        
        # `for module in modules` over a DataFrame loops over COLUMN NAMES (Melanoma, CRC...)
        for module in modules:
            logger.info("=" * 50)
            logger.info(f"Running proximity analysis for {module}")
            # that disease's genes: .dropna() removes blanks (columns have different lengths), .to_list() makes a list
            protein = modules[module].dropna().to_list()
            # warn if some disease genes are not in STRING
            validate_moduleprotein(protein, string)
            
            # boolean filter: keep rows of the combination table whose Disease equals this column name (must match exactly)
            module_combinations = combination[combination["Disease"] == module]
            # .empty is True when no rows matched
            if module_combinations.empty:
                logger.warning(f"No drug combinations found for disease: {module}")
                # skip to the next disease
                continue
                
            # iterrows() gives (row number, row) for each combination
            for index, row in module_combinations.iterrows():
                # 'MAP2K1,MAP2K2' -> .strip() trims outer spaces -> .split(',') -> ['MAP2K1', 'MAP2K2'].
                # !! inner spaces are NOT removed, so 'MAP2K1, MAP2K2' would produce ' MAP2K2' and fail.
                drugA = row["Drug_target1"].strip().split(",")
                # same for drug B
                drugB = row["Drug_target2"].strip().split(",")
                # log which combination is running
                logger.info(f"Running combination: {drugA} + {drugB}")
                
                # seed_genes = all target genes of both drugs
                seed_genes = []
                # extend adds the items of a list (vs append, which adds one item)
                seed_genes.extend(drugA)
                seed_genes.extend(drugB)
                # warn if any target is missing from STRING
                validate_targetprotein(seed_genes, string)
                
                # clusters = {seed gene: its expanded gene list}
                clusters = {}
                # enumerate gives i (0,1,2...) and the gene
                for i, seed_gene in enumerate(seed_genes):
                    # `seed_gene in G` = is this gene a node in the network graph?
                    if seed_gene in G:
                        clusters[seed_gene] = {
                            # [FLAW 3] NEIGHBOUR EXPANSION. [PAPER] a drug's target module = its experimentally validated targets (average 3).
                            # Here every target is replaced by itself + ALL its direct neighbours. Consequences on your real network:
                            # - sets get large (BRAF 47 nodes, EGFR 239 nodes);
                            # - d_AA and d_BB become 1.0 by construction (every neighbour is 1 step from its hub), so s_AB depends mainly on d_AB;
                            # - Melanoma pair: s_AB = -0.458 (overlap) with expansion, +0.500 (separate) with targets only -> the label changes.
                            'nodes': [seed_gene] + [n for n in G.neighbors(seed_gene) if n in G],
                            # 'name' is just a label; it is never used later
                            'name': f'Cluster {i+1} ({seed_gene})'
                        }
                    # else: gene not in network...
                    else:
                        # warn
                        logger.warning(f"Seed gene {seed_gene} not found in network")
                        # ...give it an empty list so the code below doesn't crash
                        clusters[seed_gene] = {'nodes': [], 'name': f'Cluster {i+1} ({seed_gene})'}
                
                # lists that will hold each drug's full expanded target genes
                drug_a_targets = []
                drug_b_targets = []
                # for each gene of drug A (e.g. MAP2K1, MAP2K2) add its cluster
                for i in drugA:
                    drug_a_targets.extend(clusters[i]["nodes"])
                # same for drug B
                for i in drugB:
                    drug_b_targets.extend(clusters[i]["nodes"])

                # set() removes duplicates (genes shared between clusters)
                drug_a_targets = set(drug_a_targets)
                drug_b_targets = set(drug_b_targets)
                # the disease genes as a set (genes missing from STRING are dropped later, in calculate_proximity)
                disease_genes = set(protein)
                
                # try: run the full analysis for this combination (n_random=1000 randomisations, alpha=0.05)
                try:
                    result = analyzer.classify_drug_combination_with_confidence(drug_a_targets, drug_b_targets, disease_genes, n_random=1000, alpha=0.05)
                    
                    # store one result row (a dictionary) in the list. Column names here become the CSV headers.
                    results.append({
                        "Disease": module,
                        "Drug_target1": row["Drug_target1"],
                        "Drug_target2": row["Drug_target2"],
                        "Z_score_drugA": result["z_da"],
                        # (commented-out lines starting with # are ignored: the author chose not to output the drug-disease p-values)
                        #"Pvalue_drugA":result["pvalue_da"],
                        "Z_score_drugB": result["z_db"],
                        #"Pvalue_drugB":result["pvalue_db"],
                        "Sab_score": result["s_ab"],
                        "P_value_separation": result["p_separation"],
                        "Classification": result["classification"]})
                    # log the final label
                    logger.info(f"Classification: {result['classification']}")
                # if anything failed (e.g. no overlapping nodes) record an ERROR row instead of stopping everything
                except Exception as e:
                    logger.error(f"Error processing combination {drugA} + {drugB} for {module}: {e}")
                    results.append({
                        "Disease": module,
                        "Drug_target1": row["Drug_target1"],
                        "Drug_target2": row["Drug_target2"],
                        "Z_score_drugA": "ERROR",
                        #"Pvalue_drugA": "NA",
                        "Z_score_drugB": "ERROR",
                        #"Pvalue_drugB":"NA",
                        "Sab_score": "ERROR",
                        "P_value_separation": "ERROR",
                        "Classification": f"Error: {str(e)}"})
        
        # turn the list of dictionaries into a table
        results_df = pd.DataFrame(results)
        # write Drug_combination_results.csv; index=False = don't write row numbers
        results_df.to_csv("Drug_combination_results.csv", index=False)
        logger.info(f"Results saved to Drug_combination_results.csv with {len(results_df)} records")

        # count rows without 'Error' in the label. The part `not in results_df.columns` is always true
        # (a label is never a column name) - harmless but pointless.
        successful_runs = len([r for r in results if r["Classification"] not in results_df.columns and "Error" not in r["Classification"]])
        # errors = total rows - successful rows
        errors = len(results) - successful_runs
        logger.info(f"Analysis completed: {successful_runs} successful runs, {errors} errors")
    # outer except: log the fatal error and re-raise it
    except Exception as e:
        logger.error(f"Fatal error in main execution: {e}")
        raise

# `if __name__ == \"__main__\":` = run main() only when this file is executed directly
# (python drug_synergy.py), NOT when another script imports it as a library.
if __name__ == "__main__":
    main()