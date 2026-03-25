**Genome evolution in plant pathogenic bacteria**

**Signficance**
Bacterial plant pathogens can devastate the production of crops used for food, fiber, or fuel.  New emerging pathogens capable of causing epidemics due to host range expansion, acquired pathogenicity, overcoming host resistance, or antibiotic resistance, is an ever-present threat in today’s world of agriculture and global trade. Understanding the drivers and mechanisms of evolution of these traits, especially pathogenicity, is important for informing predictive disease models, creating new disease resistant plant cultivars, developing biological control methods, and even informing public policy for surveillance and quarantine guidelines.  Additionally, much of the information gained about bacterial plant pathogen evolution may also be translatable and informative of pathogens causing infectious disease in humans.  



**Materials and Methods**

1. Reference Genome Selection
The reference genome species were selected from Agrios 2005 and Young et al. 1996 (Review of Plant Pathology) and then searched in the NCBI database. In order to be selected, multiple assembled genomes for that species or pathovar needed to be present, and the official reference genomes were selected when available. Metadata was primarily gathered from CABI Compendium Datasheets and other literature as needed.

2. Phylogenetic Tree
A phylogenetic tree for comparative genomic analysis was constructed using OrthoFinder (v2.4.0) (Emms and Kelly 2019), which infers orthogroups from protein sequences using default settings. Orthogroups were identified via DIAMOND-based method. The unrooted gene and species trees are inferred using the DendroBLAST  and  STAG (Species Tree from All Genes) algorithm, respectively. The STAG species tree is rooted with the STRIDE (Species Tree Root Inference from Duplication Events) algorithm by identifying high-confidence gene duplication events. The resulting rooted species tree from OrthoFinder output (SpeciesTree_rooted_node_labels) was employed for subsequent comparative analyses. The tree annotation was added using ggtree and gheatmap functions in R. Genome sizes were determined from genome files using samtools. The presence of Type III secretion system (T3SS) was identified if present in the respective genome used through NCBI genome annotations, otherwise asterisk (*) is added if presence is not consistent across species.
3. SecReT6
To find presence of T6SS, the local version of SecReT6 v3 was download to the Alabama Supercomputer Authority HPC and genome assemblies downloaded from ascension number (Supplementary Table 1) from NCBI. The SecReT6 database was used and the pipeline was modified to run on loop and be compatible with multipartite assemblies. Eight genomes returned no result (neither present or absence of T6SS) and are listed as NA.

4. Correlation
Correlation coefficients were calculated in Microsoft Excel using the CORREL function with the metadata Supplementary Table 2. Present was converted to 1 and absent to -1. Vascular and biotroph was converted to 1, non-vascular and necrotroph was converted to  -1, and both or hemibiotroph was converted to 0. 

5. T2SS
We screened the genomes of plant pathogenic bacteria strains for the presence of various genes involved in breakdowns Carbohydrate Esterases (CEs), Polysaccharide Lyases (PLs), Glycoside Hydrolases (GHs),  GlycosylTransferases (GTs) of carbohydrates, Auxiliary Activities (AAs), and Carbohydrate-Binding Modules (CBMs). Protein sequences (.faa) previously annotated via Prokka (v1.14.5) were subjected to CAZyme domain assignment using the run_dbcan3 tool (https://github.com/linnabrown/run_dbcan3). Searches were performed against the HMMER, DIAMOND, and eCAMI databases using default settings. To ensure high-confidence assignments, final CAZyme domains were retained only if supported by the best hits from at least two of the three databases.


**References**



