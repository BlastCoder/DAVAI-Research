# DAVAI Research
The field of drug discovery is constantly evolving, with researchers seeking new and innovative methods to accelerate the development of life-saving medications. Computational drug discovery offers a promising approach by harnessing the power of artificial intelligence to streamline the process. DAVAI is a research project exploring one such avenue: the use of large language models (LLMs) specifically trained for generating novel drug candidates. The ultimate objective is to generate random graphs of protein interactions and analyze them for drug repurposing and other machine learning applications. By leveraging natural language models to extract essential drug properties from scientific literature, machine learning algorithms can subsequently generate potential chemicals and proteins that exhibit these properties. This approach aims to harness the synergy between natural language processing and machine learning to innovate and expedite the drug discovery process.

## My Research
My role within this research project is to cluster small molecules within a Cartesian space using various similarity metrics. By grouping molecules with similar properties, we can train an LLM to recognize these patterns and predict the creation of novel molecules based on a vector input of desired characteristics. This methodology relies heavily on molecular fingerprinting techniques, which are computerized representations of molecules. These usually contain a unique sequence of bits or tokens that encodes certain properties, structural and chemical, of a molecule. I am currently exploring two particularly promising fingerprinting methods: SELFIES (Self- Referencing Embedded Strings) and network graphs. I am currently exploring the development of these techniques, analyze the string edit distance and maximal subgraph isomorphism similarity metrics used in conjunction with them, and discuss the distinct benefits and challenges associated with each approach. By leveraging molecular fingerprinting techniques, researchers can assign similarity metrics to cluster molecules with specific attributes, streamlining the drug discovery process. These results are discussed within the similarity measures paper I have written.

## Data Acquisition
String and graph representations of the molecules were obtained through a two-step process. First, bulk downloads of molecule data in .mol file format were retrieved from the PubChem database. These .mol files represent the molecular structure in a standardized format.

Next, using Python’s RDKit library, the .mol files were converted into SMILES strings. SMILES strings offer a compact, text-based representation of the molecule’s structure. However, for this study, we aimed for a string representation with guaranteed bijectivity (one-to-one correspondence between string and molecule). Therefore, we further con- verted the SMILES strings into SELFIES strings using a custom conversion model. This ensures that every SELFIES string uniquely corresponds to a valid molecule.

For the graph representation, the same .mol files were used to generate adjacency matrices. These matrices capture the connectivity information between atoms in the molecule. Finally, the NetworkX library was employed to convert the adjacency matrices into network planar graphs. In these graphs, nodes represent individual atoms, and edges represent the bonds connecting them. This network representation allows us to leverage graph-based algorithms for further analysis and exploration of the molecular properties.

## Molecular Fingerprints
While this is discussed more comprehensively within the paper, the two promising fingerprinting methods I have researched are SELFIES (Self- Referencing Embedded Strings) and network graphs. SELFIES, an extension of SMILES, offers a robust compact string representation which is both bijective and chemically valid. As such, it allows the opportunity to introduce random mutations to build a latent space of molecules for machine learning models. Network graphs, while being more memory-intensive, captures deeper chemical properties by providing a realistic representation of compounds, offering a wider range of graph theory algorithms.

## Similarity Metrics
While this is discussed more comprehensively within the paper, here is a brief overview of the various promising similarity metrics:
SELFIES:
  1. Hamming Distance
  2. Levenshtein Distance
  3. Jaro-Winker Metric
  4. Jaccard Similarity
  5. n-gram Analysis

Network Graphs:
  1. Tanimoto Metric
  2. QSBR Coefficient
  3. RASCAL Value

## Python Programs
Within this research, I built various Python programs to compare the different similarity metrics and molecular fingerprinting techniques. Here is a brief description of all the code I have developed:
Format Converters:
  1. molecular name to .mol file
  2. .mol file to SMILES string
  3. .mol file to SELFIES string
  4. .mol file to InChl key
  5. .mol file to adjacency matrix
  6. adjacency matrix to NetworkX graph

Similarity Measures:
  1. Hamming Distance (SELFIES)
  2. Levenshtein Distance (SELFIES)
  3. Jaro-Winker Metric (SELFIES)
  4. Jaccard Similarity (SELFIES)
  5. n-gram Analysis (SELFIES)
  6. Maximum common subgraph (graph)
  7. Tanimoto metric (graph)
  8. QSBR coefficient (graph)
  9. RASCAL value (graph)


