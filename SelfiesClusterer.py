import json
import random
import numpy as np
import selfies as sf

class MySelfiesClusterer:
    """
    A class for clustering SELFIES strings using a K-Means algorithm with a custom edit distance metric.
    
    Attributes:
        k (int): The number of clusters.
        max_iterations (int): The maximum number of iterations for the K-Means clustering algorithm.
    """

    def __init__(self, k=3, max_iterations=100):
        """
        Initializes the SelfiesClusterer with the specified number of clusters and maximum iterations.
        
        Parameters:
            k (int): The number of clusters.
            max_iterations (int): The maximum number of iterations for the clustering algorithm.
        """
        self.k = k
        self.max_iterations = max_iterations

    def load_selfies_from_json(self, file_path):
        """
        Loads SELFIES strings from a JSON file.
        
        Parameters:
            file_path (str): The path to the JSON file containing SELFIES strings.
        
        Returns:
            list: A list of SELFIES strings loaded from the file.
        """
        with open(file_path, 'r') as file:
            data = json.load(file)
            selfies_strings = [entry['selfies'] for entry in data]
        print(f"Loaded {len(selfies_strings)} SELFIE strings from the JSON file.")
        return selfies_strings

    @staticmethod
    def selfies_edit_operations(s1, s2):
        """
        Calculates the number of edit operations (insertions, deletions, substitutions) required to transform
        one SELFIES string into another.

        Parameters:
        - s1, s2 (str): SELFIES strings representing two molecular structures.

        Returns:
        - tuple: A tuple containing counts of deletions, insertions, and substitutions (in that order).
        """
        # Split both SELFIES strings into lists of their symbols.
        tokens1 = list(sf.split_selfies(s1))
        tokens2 = list(sf.split_selfies(s2))
        
        # Get the lengths of the tokenized strings.
        m, n = len(tokens1), len(tokens2)
        
        # Initialize a matrix to store the counts of operations required to match substrings.
        # Each cell will store a tuple (deletions, insertions, substitutions).
        operations = [[(0, 0, 0) for _ in range(n + 1)] for _ in range(m + 1)]
        
        # Initialize the first row and column of the matrix with the number of operations
        # needed to match an empty string: all deletions for the first column, all insertions for the first row.
        for i in range(1, m + 1):
            operations[i][0] = (i, 0, 0)  # i deletions
        for j in range(1, n + 1):
            operations[0][j] = (0, j, 0)  # j insertions
        
        # Iterate through the matrix, filling in the counts of operations.
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if tokens1[i - 1] == tokens2[j - 1]:
                    operations[i][j] = operations[i - 1][j - 1]  # No operation needed for a match.
                else:
                    # Calculate costs for each operation from the perspective of operations needed.
                    delete = (operations[i - 1][j][0] + 1, operations[i - 1][j][1], operations[i - 1][j][2])  # Deletion
                    insert = (operations[i][j - 1][0], operations[i][j - 1][1] + 1, operations[i][j - 1][2])  # Insertion
                    substitute = (operations[i - 1][j - 1][0], operations[i - 1][j - 1][1], operations[i - 1][j - 1][2] + 1)  # Substitution
                    
                    # Choose the operation with the minimum total cost, prioritizing substitutions for ties.
                    operations[i][j] = min([delete, insert, substitute], key=lambda x: (sum(x), -x[2]))
                    
        # The final tuple in the matrix contains the counts of each operation required to transform s1 into s2.
        return operations[m][n]

    @classmethod
    def edit_distance_wrapper(cls, s1, s2):
        """
        Calculates the total edit distance between two SELFIES strings.
        
        Parameters:
            s1 (str): The first SELFIES string.
            s2 (str): The second SELFIES string.
        
        Returns:
            int: The total edit distance between the two SELFIES strings.
        """
        deletions, insertions, substitutions = cls.selfies_edit_operations(s1, s2)
        return deletions + insertions + substitutions

    def find_closest_centroid(self, string, centroids):
        """
        Finds the index of the closest centroid to a given SELFIES string.
        
        Parameters:
            string (str): The SELFIES string.
            centroids (list): A list of centroid SELFIES strings.
        
        Returns:
            int: The index of the closest centroid.
        """
        min_distance = float('inf')
        closest_centroid_index = 0
        for i, centroid in enumerate(centroids):
            distance = self.edit_distance_wrapper(string, centroid)
            if distance < min_distance:
                min_distance = distance
                closest_centroid_index = i
        return closest_centroid_index

    def calculate_new_centroid(self, cluster):
        """
        Calculates the new centroid for a cluster based on the minimum sum of edit distances.
        
        Parameters:
            cluster (list): A list of SELFIES strings in the cluster.
        
        Returns:
            str: The new centroid SELFIES string for the cluster.
        """
        min_sum_distance = float('inf')
        new_centroid = cluster[0]
        for string in cluster:
            sum_distance = sum(self.edit_distance_wrapper(string, other) for other in cluster)
            if sum_distance < min_sum_distance:
                min_sum_distance = sum_distance
                new_centroid = string
        return new_centroid

    def calculate_ssd(self, clusters, centroids):
        """
        Calculates the sum of squared distances (SSD) for all clusters.
        
        Parameters:
            clusters (list): A list of clusters, each containing SELFIES strings.
            centroids (list): A list of centroid SELFIES strings.
        
        Returns:
            float: The sum of squared distances for all clusters.
        """
        ssd = 0
        for cluster, centroid in zip(clusters, centroids):
            for string in cluster:
                distance = self.edit_distance_wrapper(string, centroid)
                ssd += distance ** 2
        return ssd

    def k_means_plus_plus_initialization(self, strings):
        """
        Initializes centroids using the K-Means++ algorithm.
        
        Parameters:
            strings (list): A list of SELFIES strings to cluster.
        
        Returns:
            list: A list of initial centroid SELFIES strings.
        """
        centroids = [random.choice(strings)]
        while len(centroids) < self.k:
            distances = [min(self.edit_distance_wrapper(string, centroid) for centroid in centroids) for string in strings]
            distances_sum = sum(distances)
            probabilities = [distance / distances_sum for distance in distances]
            cumulative_probabilities = np.cumsum(probabilities)
            r = random.random()
            for i, cp in enumerate(cumulative_probabilities):
                if r < cp:
                    centroids.append(strings[i])
                    break
        return centroids

    def k_means_clustering(self, strings):
        """
        Performs K-Means clustering on the given SELFIES strings.
        
        Parameters:
            strings (list): A list of SELFIES strings to cluster.
        
        Returns:
            tuple: A tuple containing the clusters and their centroids.
        """
        centroids = self.k_means_plus_plus_initialization(strings)
        for iteration in range(self.max_iterations):
            print(f"Iteration: {iteration+1}/{self.max_iterations}")
            clusters = [[] for _ in range(self.k)]
            for string in strings:
                closest_centroid_index = self.find_closest_centroid(string, centroids)
                clusters[closest_centroid_index].append(string)
            
            new_centroids = []
            for cluster in clusters:
                if cluster:  # Ensure the cluster is not empty
                    new_centroid = self.calculate_new_centroid(cluster)
                    new_centroids.append(new_centroid)
                else:
                    new_centroids.append(random.choice(strings))  # Handle empty cluster

            ssd = self.calculate_ssd(clusters, centroids)
            print(f"Iteration {iteration+1}, SSD: {ssd}")

            if new_centroids == centroids:
                break
            else:
                centroids = new_centroids

        return clusters, centroids

    def run(self, file_path):
        """
        Loads SELFIES strings from a file and performs K-Means clustering.
        
        Parameters:
            file_path (str): The path to the JSON file containing SELFIES strings.
        """
        selfies_strings = self.load_selfies_from_json(file_path)
        clusters, centroids = self.k_means_clustering(selfies_strings)

        # Optionally, print the results
        for i, cluster in enumerate(clusters):
            print(f"Cluster {i+1} (size {len(cluster)}):", cluster)
        print("Centroids:", centroids)


# Example usage
file_path = '/home/james/Workspace/Data Dump/SMILES/SELFIES_table_sampled.json'
clusterer = MySelfiesClusterer(k=3, max_iterations=100)
clusterer.run(file_path)
