#imports
# =====================================================================================
# drug_synergy_2.py  --  corrected version of drug_synergy_1.py
#
# Every place that differs from drug_synergy_1.py is marked with a comment block that
# starts with  "# [MODIFIED n]"  (n = number in the review document):
#   [MODIFIED 1]  degree counting fixed (distance == 1, not 0.0)
#   [MODIFIED 2]  Eq.(1) direction: average over DISEASE nodes (as in the paper)
#   [MODIFIED 3]  drug module = validated targets only (no neighbour expansion)
#   [MODIFIED 4]  randomisation p-value for s_AB and the "Significant" prefix removed
#   [MODIFIED 5]  names of missing / dropped / unreachable genes are now printed
#   [MODIFIED 6]  borderline drugs (z < 0 but p >= alpha) are now WARNED about; the paper's rule
#                 (close = z < 0) is the default, a stricter z<0 AND p<alpha rule is switchable
#   [MODIFIED 7]  z-score and p-value come from the SAME random sets; reproducible seed
#   [MODIFIED 8]  spaces inside target lists ("MAP2K1, MAP2K2") are handled
#   [MODIFIED 9]  pointless condition in the "successful runs" counter removed
#   [MODIFIED 10] results also saved as Drug_combination_results.xlsx
#   [MODIFIED 11] optional drug names (Drug1_name, Drug2_name) carried into the output
# Everything else is unchanged.
# =====================================================================================
import pickle
import random
import logging
import pandas as pd
import numpy as np
import networkx as nx
from typing import Set, List, Tuple, Dict, Optional

#Logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('drug_synergy_analysis.log'),
        logging.StreamHandler()])
logger = logging.getLogger(__name__)

# [MODIFIED 7] one fixed seed so that repeated runs give identical numbers.
# Change the number to get a different (but again reproducible) random draw.
RANDOM_SEED = 42

# [MODIFIED 6] FINAL DECISION: follow the paper -> a drug is "close" to the disease when z < 0.
#   False (default, = paper): close if z < 0. Borderline cases (z < 0 but p >= 0.05) are still
#                             reported with a WARNING and their p-values are written to the output.
#   True  (stricter)        : close only if z < 0 AND p < 0.05 (not used in the paper).
REQUIRE_P_FOR_CLOSE = False


class NetworkProximityAnalyzer:
    def __init__(self):
        self.shortest_paths = None

    #Network
    def load_shortest_paths(self, cache_file: str) -> Dict:
        """
        Load precomputed shortest paths from file
        """
        logger.info(f"Loading precomputed shortest paths from {cache_file}...")
        try:
            with open(cache_file, 'rb') as f:
                self.shortest_paths = pickle.load(f)
            logger.info(f"Loaded shortest paths for {len(self.shortest_paths)} nodes")
            return self.shortest_paths
        except Exception as e:
            logger.error(f"Failed to load shortest paths: {e}")
            raise

    def get_network_nodes(self) -> Set[str]:
        """
        Get all nodes in the network from shortest paths data
        """
        return set(self.shortest_paths.keys())

    def get_degree_binning(self, min_bin_size: int = 100) -> List[Tuple[int, int, List[str]]]:
        """
        Bin nodes by degree for degree-preserving randomization
        Uses the shortest paths data to calculate degrees
        """
        #Calculate degrees from shortest paths
        degrees = {}
        for node, connections in self.shortest_paths.items():
            # [MODIFIED 1] was:  if dist == 0.0
            # The cache stores HOP counts, so the only 0 in a row is the node itself
            # (every node got degree 1 -> one single bin). A direct neighbour is at
            # distance 1, so we count entries equal to 1 = the true degree.
            direct_neighbors = [n for n, dist in connections.items() if dist == 1]
            degrees[node] = len(direct_neighbors)
        degree_to_nodes = {}
        for node, degree in degrees.items():
            degree_to_nodes.setdefault(degree, []).append(node)

        values = sorted(degree_to_nodes.keys())
        bins = []
        i = 0

        while i < len(values):
            low = values[i]
            val = degree_to_nodes[values[i]].copy()
            while len(val) < min_bin_size and i < len(values) - 1:
                i += 1
                val.extend(degree_to_nodes[values[i]])
            high = values[i]
            i += 1
            if len(val) < min_bin_size and bins:
                low_prev, high_prev, val_prev = bins[-1]
                bins[-1] = (low_prev, high, val_prev + val)
            else:
                bins.append((low, high, val))
        return bins

    def get_degree_equivalents(self, nodes: Set[str], bins: List[Tuple[int, int, List[str]]]) -> Dict[str, List[str]]:
        """
        Find degree-equivalent nodes for randomization
        """
        #Calculate degrees from shortest paths
        degrees = {}
        for node in nodes:
            if node in self.shortest_paths:
                # [MODIFIED 1] was:  if dist == 0.0   (same fix as in get_degree_binning)
                direct_neighbors = [n for n, dist in self.shortest_paths[node].items() if dist == 1]
                degrees[node] = len(direct_neighbors)
            else:
                degrees[node] = 0
        node_to_equivalents = {}
        for node in nodes:
            d = degrees[node]
            for low, high, bin_nodes in bins:
                if low <= d <= high:
                    equivalents = [n for n in bin_nodes if n != node]
                    node_to_equivalents[node] = equivalents
                    break
        return node_to_equivalents

    def get_random_nodes_matching_degree(self, original_nodes: Set[str], n_random: int = 1000, min_bin_size: int = 100, seed: Optional[int] = None) -> List[Set[str]]:
        """
        Generate random node sets with degree distribution matching original nodes
        """
        if seed is not None:
            random.seed(seed)
        bins = self.get_degree_binning(min_bin_size)
        node_to_equivalents = self.get_degree_equivalents(original_nodes, bins)
        random_sets = []
        network_nodes = self.get_network_nodes()
        # [MODIFIED 7] sorted list built once (was list(network_nodes) inside the loop).
        # A Python set of strings has a different order in every run, so without sorting
        # the same seed would still give different random sets.
        network_nodes_list = sorted(network_nodes)
        for _ in range(n_random):
            random_set = set()
            # [MODIFIED 7] was:  for node in original_nodes   (set order differs per run)
            for node in sorted(original_nodes):
                if node in node_to_equivalents and node_to_equivalents[node]:
                    equivalent = random.choice(node_to_equivalents[node])
                    random_set.add(equivalent)
                else:
                    #Fallback: random node from network
                    random_set.add(random.choice(network_nodes_list))
            random_sets.append(random_set)
        return random_sets

    #Proximity
    def calculate_closest_distance(self, set_a: Set[str], set_b: Set[str]) -> float:
        """
        Calculate mean closest distance between two sets of nodes using precomputed paths
        (averages over the nodes of set_a; for each one the nearest node of set_b)
        """
        distances = []
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
                distances.append(min_dist)
        return np.mean(distances) if distances else float('inf')

    # [MODIFIED 7] now returns FIVE values (z, actual, mean_random, std_random, p_value).
    # The empirical p-value is computed from the SAME random sets as the z-score
    # (before: a second, separate set of random draws was used for the p-value).
    def calculate_proximity(self, drug_targets: Set[str], disease_genes: Set[str], n_random: int = 1000, min_bin_size: int = 100,seed: Optional[int] = None) -> Tuple[float, float, float, float, float]:
        """
        Calculate drug-disease proximity (z-score)
        Returns: (z_score, actual_distance, mean_random, std_random, p_value)
        """
        #Filter nodes that exist in network
        network_nodes = self.get_network_nodes()

        # [MODIFIED 5] genes that are dropped here were dropped silently before.
        # Now the names are printed (counts alone were printed in main()).
        dropped_drug = sorted(drug_targets - network_nodes)
        dropped_disease = sorted(disease_genes - network_nodes)
        if dropped_drug:
            logger.warning(f"{len(dropped_drug)} drug target(s) not in the distance table, dropped from the calculation: {dropped_drug}")
        if dropped_disease:
            logger.warning(f"{len(dropped_disease)} disease gene(s) not in the distance table, dropped from the calculation: {dropped_disease}")

        drug_targets = drug_targets & network_nodes
        disease_genes = disease_genes & network_nodes

        if not drug_targets or not disease_genes:
            raise ValueError("No overlapping nodes between input sets and network")

        # [MODIFIED 5] disease genes that ARE in the network but have no path to any drug
        # target (other connected component). calculate_closest_distance skips them silently.
        unreachable = sorted(g for g in disease_genes
                             if not any(g in self.shortest_paths[t] for t in drug_targets))
        if unreachable:
            logger.warning(f"{len(unreachable)} disease gene(s) unreachable from these drug targets (different network component), ignored in the distance: {unreachable}")

        #Calculate actual distance
        # [MODIFIED 2] was:  calculate_closest_distance(drug_targets, disease_genes)
        # Paper Eq.(1): average over the DISEASE nodes of the distance to the nearest drug target,
        # so the disease set must be the first argument (the one that is averaged over).
        actual_distance = self.calculate_closest_distance(disease_genes, drug_targets)

        #Generate random sets
        random_drug_sets = self.get_random_nodes_matching_degree(drug_targets, n_random, min_bin_size, seed)
        random_disease_sets = self.get_random_nodes_matching_degree(disease_genes, n_random, min_bin_size, seed)

        #Calculate random distances
        random_distances = []
        for rand_drugs, rand_diseases in zip(random_drug_sets, random_disease_sets):
            # [MODIFIED 2] was:  calculate_closest_distance(rand_drugs, rand_diseases)
            # The random distances must use the same direction as the real distance.
            dist = self.calculate_closest_distance(rand_diseases, rand_drugs)
            if dist != float('inf'):
                random_distances.append(dist)

        # [MODIFIED 7] clear error instead of a numpy warning + NaN if nothing was usable
        if not random_distances:
            raise ValueError("No finite random distances could be calculated")

        #Calculate z-score
        mean_random = np.mean(random_distances)
        std_random = np.std(random_distances)

        if std_random == 0:
            # [MODIFIED 6] was:  z_score = 0.0   (silently meant "not close")
            # Now NaN + warning: a NaN is never "< 0", so the drug is still NOT called close,
            # but the output shows that the z-score could not be calculated.
            logger.warning("Random distances have zero spread: z-score cannot be calculated (set to NaN)")
            z_score = float('nan')
        else:
            z_score = (actual_distance - mean_random) / std_random

        # [MODIFIED 7] empirical p-value from the SAME random distances as the z-score.
        # "+1" correction: p can never be exactly 0 (smallest value = 1/(n_random+1)).
        p_value = (1 + np.sum(np.array(random_distances) <= actual_distance)) / (1 + len(random_distances))
        return z_score, actual_distance, mean_random, std_random, p_value

    def calculate_separation(self, set_a: Set[str], set_b: Set[str]) -> float:
        """
        Calculate separation score between two sets of nodes (s_AB)
        s_AB = <d_AB> - (<d_AA> + <d_BB>) / 2
        """
        #Calculate d_AA (within set A)
        d_aa = self.calculate_within_set_distance(set_a)

        #Calculate d_BB (within set B)
        d_bb = self.calculate_within_set_distance(set_b)

        #Calculate d_AB (between sets A and B)
        d_ab = self.calculate_between_set_distance(set_a, set_b)
        separation = d_ab - (d_aa + d_bb) / 2.0
        return separation

    def calculate_within_set_distance(self, node_set: Set[str]) -> float:
        """
        Calculate mean closest distance within a set of nodes
        """
        # NOTE: a single-target drug has no internal distance -> 0.0 (the paper only used drugs
        # with >= 2 targets; here single-target drugs are accepted on purpose, see review doc)
        if len(node_set) <= 1:
            return 0.0

        distances = []
        nodes_list = list(node_set)

        for i, node_a in enumerate(nodes_list):
            if node_a not in self.shortest_paths:
                continue

            min_dist = float('inf')
            for j, node_b in enumerate(nodes_list):
                if i != j and node_b in self.shortest_paths[node_a]:
                    dist = self.shortest_paths[node_a][node_b]
                    if dist < min_dist:
                        min_dist = dist

            if min_dist != float('inf'):
                distances.append(min_dist)
        return np.mean(distances) if distances else float('inf')

    def calculate_between_set_distance(self, set_a: Set[str], set_b: Set[str]) -> float:
        """
        Calculate mean closest distance between two sets of nodes
        """
        distances_ab = []

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

        return np.mean(distances_ab) if distances_ab else float('inf')

    #Significance
    # [MODIFIED 7] NO LONGER CALLED. The p-value is now calculated inside calculate_proximity from
    # the same random sets as the z-score. Kept only so nothing else breaks; Eq.(1) direction
    # fixed here too ([MODIFIED 2]) in case someone calls it.
    def calculate_empirical_p_value(self, actual_distance: float, drug_targets: Set[str],  disease_genes: Set[str], n_random: int = 1000) -> float:
        """
        Calculate empirical p-value for drug-disease proximity (UNUSED, see note above)
        """
        random_distances = []

        #Generate random distances
        random_drug_sets = self.get_random_nodes_matching_degree(drug_targets, n_random)
        random_disease_sets = self.get_random_nodes_matching_degree(disease_genes, n_random)

        for rand_drugs, rand_diseases in zip(random_drug_sets, random_disease_sets):
            dist = self.calculate_closest_distance(rand_diseases, rand_drugs)
            if dist != float('inf'):
                random_distances.append(dist)

        p_value = np.mean(np.array(random_distances) <= actual_distance)

        return p_value

    # [MODIFIED 4] NO LONGER CALLED. The paper states that randomisation is not suitable for s_AB
    # (drugs have only a few targets, so the random distribution is not Gaussian). The paper
    # classifies with the SIGN of s_AB only. Function kept for reference, never run.
    def calculate_separation_p_value(self, set_a: Set[str], set_b: Set[str], actual_separation: float, n_random: int = 1000) -> float:
        """
        Calculate empirical p-value for drug-drug separation (UNUSED, see note above)
        """
        random_separations = []

        #Generate random separation scores
        random_sets_a = self.get_random_nodes_matching_degree(set_a, n_random)
        random_sets_b = self.get_random_nodes_matching_degree(set_b, n_random)

        for rand_a, rand_b in zip(random_sets_a, random_sets_b):
            separation = self.calculate_separation(rand_a, rand_b)
            if separation != float('inf'):
                random_separations.append(separation)

        #For separation: more negative values are significant for overlap (s_AB < 0)
        #More positive values are significant for separation (s_AB > 0)
        if actual_separation < 0:  #Overlap case
            p_value = np.mean(np.array(random_separations) <= actual_separation)
        else:  #Separation case
            p_value = np.mean(np.array(random_separations) >= actual_separation)

        return p_value

    #Main
    # [MODIFIED 6] new parameter require_p_for_close (default False = paper).
    #   False: a drug is "close" if z < z_threshold only                 (exactly as in the paper)
    #   True : a drug is "close" only if z < z_threshold AND p < alpha   (stricter than the paper)
    def classify_drug_combination_with_confidence(self, drug_a_targets: Set[str], drug_b_targets: Set[str], disease_genes: Set[str], n_random: int = 1000, z_threshold: float = 0.0, alpha: float = 0.05, require_p_for_close: bool = False) -> Dict:
        """
        Classify drug combination into the six classes of the paper (P1-P6)
        Returns dictionary with classification, z-scores, p-values and s_AB
        """
        #Calculate drug-disease proximities with p-values
        # [MODIFIED 7] p-value now returned by calculate_proximity (same random sets as the z-score);
        # the two separate calls of calculate_empirical_p_value were removed.
        z_da, actual_da, _, _, p_da = self.calculate_proximity(drug_a_targets, disease_genes, n_random)
        z_db, actual_db, _, _, p_db = self.calculate_proximity(drug_b_targets, disease_genes, n_random)

        #Calculate drug-drug separation (no randomisation, see [MODIFIED 4])
        s_ab = self.calculate_separation(drug_a_targets, drug_b_targets)
        # [MODIFIED 4] removed:  p_separation = self.calculate_separation_p_value(...)
        #                        separation_significant = p_separation < alpha

        #Determine the three binary conditions
        drug_a_hits_disease = z_da < z_threshold
        drug_b_hits_disease = z_db < z_threshold
        # [MODIFIED 6] borderline check: z < 0 (called "close" under the paper's rule) but p >= alpha
        # means there is little evidence of real closeness (e.g. z = -0.01, p = 0.45). It is now
        # reported instead of passing unnoticed.
        for label, z, p in (("A", z_da, p_da), ("B", z_db, p_db)):
            if z < z_threshold and not p < alpha:
                logger.warning(f"Drug {label}: z = {z:.4f} < {z_threshold} but p = {p:.4f} >= {alpha} -> borderline, weak evidence of closeness")
        # Optional stricter rule (NaN z-scores are never "< 0", so they stay not close):
        if require_p_for_close:
            drug_a_hits_disease = drug_a_hits_disease and (p_da < alpha)
            drug_b_hits_disease = drug_b_hits_disease and (p_db < alpha)
        drugs_overlap = s_ab < 0

        logger.info(f"Proximity of drug A - Z score: {z_da:.4f}, pvalue: {p_da}")
        logger.info(f"Proximity of drug B - Z score: {z_db:.4f}, pvalue: {p_db}")
        logger.info(f"Separation between drug A & drug B - Sab score: {s_ab:.4f}")

        #Classify based on the combination of conditions (six classes, Fig. 2 of the paper)
        # [MODIFIED 4 + 6] The old block had four near-identical if/else branches that added
        # "Significant ..." / "... but not significant" and gave the same label whether the
        # drug-disease p-values were significant or not. Now: one label per class, nothing else.
        #   both close  + overlap -> Overlapping Exposure     both close  + separate -> Complementary Exposure
        #   one close   + overlap -> Indirect Exposure        one close   + separate -> Single Exposure
        #   none close  + overlap -> Non-exposure             none close  + separate -> Independent Action
        if drug_a_hits_disease and drug_b_hits_disease:
            logger.info("Both drug targets in proximity with the disease module")
            classification = "Overlapping Exposure" if drugs_overlap else "Complementary Exposure"
        elif drug_a_hits_disease and not drug_b_hits_disease:
            logger.info("Only drug target A in proximity with the disease module")
            classification = "Indirect Exposure" if drugs_overlap else "Single Exposure"
        elif not drug_a_hits_disease and drug_b_hits_disease:
            logger.info("Only drug target B in proximity with the disease module")
            classification = "Indirect Exposure" if drugs_overlap else "Single Exposure"
        else:
            logger.info("None of the drug targets in proximity with the disease module")
            classification = "Non-exposure" if drugs_overlap else "Independent Action"

        #Return all metrics in a dictionary
        return {
            "classification": classification,
            "z_da": z_da,
            "pvalue_da":p_da,
            "z_db": z_db,
            "pvalue_db":p_db,
            "s_ab": s_ab,
            "drug_a_close": drug_a_hits_disease,
            "drug_b_close": drug_b_hits_disease,
        }


# [MODIFIED 5] both helpers now print the NAMES of the missing genes (before: only the count).
# The lookup also uses a set (fast) instead of scanning the numpy arrays for every gene.
def validate_moduleprotein(module, string):
    string_nodes = set(string.source.values) | set(string.target.values)
    not_string = [i for i in module if i not in string_nodes]
    if len(not_string) > 0:
        logger.warning(f"{len(not_string)} out of {len(module)} proteins in the disease module are not in string")
        logger.warning(f"Disease genes missing from string: {sorted(not_string)}")

def validate_targetprotein(drug_targets, string):
    string_nodes = set(string.source.values) | set(string.target.values)
    not_string = [i for i in drug_targets if i not in string_nodes]
    if len(not_string) > 0:
        logger.warning(f"{len(not_string)} out of {len(drug_targets)} drug targets are not in string!")
        logger.warning(f"Drug targets missing from string: {sorted(not_string)}")


# Main function
def main():
    logger.info("Starting drug synergy analysis")

    try:
        analyzer = NetworkProximityAnalyzer()
        analyzer.load_shortest_paths("shortest_paths_cache.pkl")

        logger.info("Loading STRING PPI network")
        string = pd.read_csv("STRING_PPI_filtered.txt", sep="\t")
        G = nx.from_pandas_edgelist(string, source='source', target='target', edge_attr=True)

        logger.info("Loading drug target data")
        combination = pd.read_excel('Drug_target.xlsx', sheet_name='Drug_combination')
        modules = pd.read_excel('Drug_target.xlsx', sheet_name='Module_genes')

        # Create a list to store results
        results = []

        for module in modules:
            logger.info("=" * 50)
            logger.info(f"Running proximity analysis for {module}")
            # [MODIFIED 8] also strip spaces around each disease gene name (cell values like "BRAF ")
            protein = [str(g).strip() for g in modules[module].dropna().to_list()]
            validate_moduleprotein(protein, string)

            module_combinations = combination[combination["Disease"] == module]
            if module_combinations.empty:
                logger.warning(f"No drug combinations found for disease: {module}")
                continue

            for index, row in module_combinations.iterrows():
                # [MODIFIED 8] was:  row["Drug_target1"].strip().split(",")
                # that left a space on the second gene of "MAP2K1, MAP2K2" (' MAP2K2' is not found).
                # Now every gene is stripped and empty entries (e.g. trailing comma) are dropped.
                drugA = [g.strip() for g in str(row["Drug_target1"]).split(",") if g.strip()]
                drugB = [g.strip() for g in str(row["Drug_target2"]).split(",") if g.strip()]
                logger.info(f"Running combination: {drugA} + {drugB}")

                # [MODIFIED 7] same seed for every combination: result of one pair does not depend on
                # which pairs ran before it, and a repeat run gives identical numbers.
                random.seed(RANDOM_SEED)

                seed_genes = []
                seed_genes.extend(drugA)
                seed_genes.extend(drugB)
                validate_targetprotein(seed_genes, string)

                clusters = {}
                for i, seed_gene in enumerate(seed_genes):
                    if seed_gene in G:
                        clusters[seed_gene] = {
                            # [MODIFIED 3] was:  'nodes': [seed_gene] + [n for n in G.neighbors(seed_gene) if n in G],
                            # The paper's drug module = the drug's validated targets only. Adding all
                            # neighbours made d_AA = d_BB = 1 by construction and could flip overlap/separate.
                            # (Single-target drugs are accepted on purpose for this user story.)
                            'nodes': [seed_gene],
                            'name': f'Cluster {i+1} ({seed_gene})'
                        }
                    else:
                        logger.warning(f"Seed gene {seed_gene} not found in network")
                        clusters[seed_gene] = {'nodes': [], 'name': f'Cluster {i+1} ({seed_gene})'}

                drug_a_targets = []
                drug_b_targets = []
                for i in drugA:
                    drug_a_targets.extend(clusters[i]["nodes"])
                for i in drugB:
                    drug_b_targets.extend(clusters[i]["nodes"])

                drug_a_targets = set(drug_a_targets)
                drug_b_targets = set(drug_b_targets)
                disease_genes = set(protein)

                # [MODIFIED 3] reminder in the log: with one target, d_AA (or d_BB) is 0 by convention
                for label, tg in (("A", drug_a_targets), ("B", drug_b_targets)):
                    if len(tg) == 1:
                        logger.warning(f"Drug {label} has a single target {sorted(tg)}: d_{label}{label} = 0 by convention (the paper used drugs with >= 2 targets)")

                try:
                    result = analyzer.classify_drug_combination_with_confidence(drug_a_targets, drug_b_targets, disease_genes, n_random=1000, alpha=0.05, require_p_for_close=REQUIRE_P_FOR_CLOSE)

                    # [MODIFIED 4 + 7] P_value_separation column removed (no randomisation of s_AB);
                    # the drug-disease p-value columns (commented out before) are switched on.
                    # [MODIFIED 11] optional drug-name columns (Drug1_name / Drug2_name in the Drug_combination
                    # sheet) are copied to the output so the result rows read "Dabrafenib + Trametinib".
                    # If the columns are absent the output is empty there; nothing else depends on them.
                    results.append({
                        "Disease": module,
                        "Drug_A": row.get("Drug1_name", ""),
                        "Drug_B": row.get("Drug2_name", ""),
                        "Drug_target1": row["Drug_target1"],
                        "Drug_target2": row["Drug_target2"],
                        "Z_score_drugA": result["z_da"],
                        "Pvalue_drugA":result["pvalue_da"],
                        "Z_score_drugB": result["z_db"],
                        "Pvalue_drugB":result["pvalue_db"],
                        "Sab_score": result["s_ab"],
                        "Classification": result["classification"]})
                    logger.info(f"Classification: {result['classification']}")
                except Exception as e:
                    logger.error(f"Error processing combination {drugA} + {drugB} for {module}: {e}")
                    results.append({
                        "Disease": module,
                        "Drug_A": row.get("Drug1_name", ""),
                        "Drug_B": row.get("Drug2_name", ""),
                        "Drug_target1": row["Drug_target1"],
                        "Drug_target2": row["Drug_target2"],
                        "Z_score_drugA": "ERROR",
                        "Pvalue_drugA": "NA",
                        "Z_score_drugB": "ERROR",
                        "Pvalue_drugB":"NA",
                        "Sab_score": "ERROR",
                        "Classification": f"Error: {str(e)}"})

        results_df = pd.DataFrame(results)
        results_df.to_csv("Drug_combination_results.csv", index=False)
        # [MODIFIED 10] results are now ALSO written as an Excel file (user story asks for Excel output).
        # The CSV line above is unchanged, so the old output still exists.
        results_df.to_excel("Drug_combination_results.xlsx", index=False)
        logger.info(f"Results saved to Drug_combination_results.csv and Drug_combination_results.xlsx with {len(results_df)} records")

        # [MODIFIED 9] was:  ... if r["Classification"] not in results_df.columns and "Error" not in r["Classification"]
        # The first condition was always true (a label is never a column name), so it is removed.
        successful_runs = len([r for r in results if "Error" not in r["Classification"]])
        errors = len(results) - successful_runs
        logger.info(f"Analysis completed: {successful_runs} successful runs, {errors} errors")
    except Exception as e:
        logger.error(f"Fatal error in main execution: {e}")
        raise

if __name__ == "__main__":
    main()
