Cancer Systems Biology, Technical University of Denmark, 2800, Lyngby, Denmark
Cancer Structural Biology, Danish Cancer Institute, 2100, Copenhagen, Denmark
Repository associated to the publication:

```
MAVISp: Multi-layered Assessment of VarIants by Structure for proteins
Matteo Arnaudi, Ludovica Beltrame, Kristine Degn, Mattia Utichi, Pablo Sánchez-Izquierdo, Simone Scrima, Francesca Maselli, Karolina Krzesińska, Terézia Dorčaková, Jordan Safer, Katrine Meldgård, Julie Bruun Brockhoff, Amalie Drud Nielsen, Alberto Pettenella, Jérémy Vinhas, Peter Wad Sackett, Claudia Cava, Sumaiya Iqbal, View ORCID ProfileMatteo Lambrughi, Matteo Tiberti, Elena Papaleo
biorxiv, doi: https://doi.org/10.1101/2022.10.22.513328
```

# CABS-flex pipeline

## Overview

This snakemake[^Mölder2021]-based pipeline is designed to perform in-silico modeling of protein flexibility with CABSflex [^Kurcinski2019]. 
In particular, it is designed to run CABSflex on a set of structural variants of a given protein, each one stored in a separate PDB file.

### CABSflex

CABSflex is a tool that creates fast simulations of protein structure flexibility based on simulations of protein dynamics using coarse-grained protein models with an optional reconstruction to an all-atom representations. CABSflex is computationally lighter than a classic all-atom molecular dynamics which is the merit of its use.

The input to the CABSflex tool is, as a minimum, a PDB file starting structure. Optionally a set of restraints and simulation parameters can be chosen. Once started, the CABS simulation generates a conformational trajectory of 1000 models. These 1000 models are reduced to k clusters and each cluster is assigned a representative conformation (medoid). The k representatives are then reconstructed to all-atom representation and are the final output.

Optionally, our pipeline allows to add a ligand in the final output structure by rigid superimposition on a reference structure in which the ligand is present.

Once all the models have been generated, the THESEUS software[^Theobald2006] is used to generate a single trajectory and best-superimpose all the models for each run.

Finally, DSSP [^Kabsch1983] is used to assign secondary structure definition for all the models and the starting structure, and the SOV_refine score [^Liu2018], which is a score designed to quantify the secondary structure overlap between conformations, is calculated between the original structure and the ones generated using CABSflex.

The input to this CABSflex pipeline is a comma-separated file, variants.csv, containing information on one or more directories containing PDB files. The content of this file is described below in the design of the input files section. Additionally, a configuration file and the snake file itself are required. 

## Requirements

The user must have `Snakemake` [^Mölder2021], `python` v3.7 or higher, a distribution of Perl, `CABS-flex` [^Kurcinski2019], `MODELLER` [^Webb2016], `DSSP` [^Kabsch1983] and `THESEUS` [^Theobald2006] installed. This pipeline also uses the BioPython, MDAnalysis and Pandas Python packages. Additionally, the Perl script [`SOV_refine.pl`](http://dna.cs.miami.edu/SOV/) [^Liu2018] should also be present in the directory the pipeline will be run from.  

The user will also need to provide:

- A directory containing one or more PDB files. PDB files should be named as follows:

`{identifier}_{chain}_{WT residue type}{residue number}{variant residue type}.pdb`

so for instance:

`2XWR_A_ALA129ASP.pdb`

the corresponding wild-type structure is also identified and used as:

`{identifier}_{chain}_WT.pdb`

- An appropriately formatted variants.csv and config.yaml file (see below)

- if the user wishes to include a ligand in the final structure, they need to provide a PDB file containing
the complex to be used as reference to place the ligand in the model structures. This can have any name
and must be specified in the configuration file (see below).

## Usage

1. Clone this repository where you would like to run to the pipeline:

`git clone https://github.com/ELELAB/cabs_production.git`
`cd cabs_production/cabsflex`

2. Edit the config.yaml and variants.csv files appropriately for your system.

3. run snakemake:

`snakemake --cores NCORES`

### Options

Some widely-used snakemake options that might come useful:

|Option  |Meaning  |
|---|---|
|`--cores`  |Number of cores to be used|
|`-p`    |Print out all actions that are being performed|
|`--dry-run` |Generate DAG but don't run rules - for testing purposes|



### Design of input files

#### snakefile

It is the file that contains all the instructions to be executed. Usually, it doesn't need to be changed.

#### variants.csv
It is a CSV file that uses a comma as a column delimiter which specifies the structures that should undergo flexibility modeling. Each line specifies a set of mutants (i.e. a set of PDB files). For each line, all the run types defined in the config file will be performed (see below for more details).

The following information must be added:

|Column name|Meaning|Example|
|---|---|---|
|`source_structure`|The identifier of the original structure, in this case PDB ID and original chain identifier|2XWR|
|`chain_in_source`|It is the chain of interest in the source structure|A|
|`aa_num_start`|It is the numbering of the first amino acid in the structure|91|
|`aa_num_end`|It is the numbering of the last amino acid in the structure|289|
|`source_type`|The type of structure in question, e.g., X-ray, NMR or model|model|
|`pdb_dir`|the directory in which the structure can be found. This path starts from the `pdb_input_dir` directory specified in the config file|2XWR_91-289/test/|

Example:
source_structure,chain_in_source,aa_num_start,aa_num_end;source_type,pdb_dir
2XWRa,A,91,289,model,2XWRa_91-289/test/

 N.B. The header line within the table must be kept.

In this example, we are considering chain A on the X-ray structure with PDB id 2XWR, which has undergone computational mutations. Each model of a mutation is available in `2XWRa_91-289/test/`

#### config.yaml
A YAML file containing the script configuration. 

The following options and parameters must be set within the configuration file:

|Generic options|Meaning|
|---|---|
|`variant_csv`|It is a file “,” separated containing all the input names, as described above|
|`pdb_input_dir`|Path of the directory containing the folders with the PDB structures (as specified in the `pdb_dir` column of the csv file)|
|`out_dir`|Path of the folder where the outputs will be written, per default is “./results”|

|CABS-flex options|Meaning|
|---|---|
|`env`|This is the command used to set up the environment right before running CABSflex, so that the CABSflex executable is available. Most of the times, this will be activating a virtual environment (see example)|
|`k-medoids`|Number of medoids in k-medoids clustering algorithm|
|`dssp_location`|Path for the DSSP program|
|`verbose`|Controls how explicit the program output is. It ranges from 0 (only critical messages) to 4 (maximum verbosity)|

The `runtypes` section allows to define different runtypes, i.e. set of parameters per run. All of
them will be considered and run in combination with the PDB sets defined in the input csv file. Options
for each runtype are:

|runtype options|Meaning|
|---|---|
|`mc-cycles`| This option,(-y, --mc-cycles NUM) sets the number of Monte Carlo cycles to NUM (NUM>0, default value = 50).|
|`mc-steps`| This option, (-s, --mc-steps NUM) sets number of Monte Carlo cycles between trajectory frames to NUM (NUM > 0, default value = 50). |
|`mc-annealing`| This option, (-a, --mc-annealing NUM) sets number of Monte Carlo temperature annealing cycles to NUM (NUM > 0, default value = 20, changing default value is recommended only for advanced users). |
|`rebuild-models`| Rebuild final models to all-atom representation (requires MODELLER installed), the default value is true |
|`random_seed`| Setting a random seed for the run. As a default all random seeds are set to 10, which ensures reproducibility of runs. |
|`modelling_restraint`| This option allow the user to restrain how flexible residues in secondary structures should be, the default value is "ss2 3 3.8 8.0". The options are ss2; where residues in contact are restrained if they are both in structured SS. ss1; where a residue in a secondary structure is restrained to stay in that secondary structure. all; where all residues are restrained to remain in their original secondary structure assignment. The numerical values describe the gap between restrained pairs (3), the minimum distance of restraining (3.8 Å) the maximim distance of restraint (8.0 Å)|

|File Availability|Meaning|
|`theseus_location`| The placement and name of the theseus program: e.g. "/usr/local/theseus-3.3.0/theseus"
|`dssp_location`| The placement and name of the dssp program: e.g. "/usr/local/dssp-3.0.10/bin/mkdssp"|
|`sov_location`| The placement and name of the SOV script, e.g. "./SOV_refine.pl"|

#### Specifying run types

A run type is a specific combination of CABSflex restraints and input secondary structure definition. One or more runtypes
must be specified in the config.yaml configuration file.

Run types can be specified a specific structure, for instance:

  run_types:
    default:
      ss_def: null
      restraints: null
      ligand:
        place_ligand: true
        ligand_structure: '2XWR.pdb'
        model_fit_selection: "protein and name CA and segid A"
        ligand_fit_selection: "protein and name CA and segid A"
        ligand_selection: "segid A and resname ZN"
    metal_bound:
      ss_def: null
      restraints: 
        - '--ca-rest-add 179:A 176:A 5.5 1.0'
        - '--ca-rest-add 179:A 238:A 6.6 1.0'
        - '--ca-rest-add 179:A 242:A 8.1 1.0'
        - '--ca-rest-add 176:A 238:A 7.1 1.0'
        - '--ca-rest-add 176:A 242:A 6.5 1.0'
        - '--ca-rest-add 238:A 242:A 6.7 1.0'
      ligand:
        place_ligand: true
        ligand_structure: '2XWR.pdb'
        model_fit_selection: "protein and name CA and segid A"
        ligand_fit_selection: "protein and name CA and segid A"
        ligand_selection: "segid A and resname ZN"

For each run type (in this example, “default and “metal_bound"), both a secondary structure definition and restraints can be defined.

Secondary structure definition (ss_def) can be either a string of letters, as specified by CABSflex, or null to indicate no secondary structure defined.

Restraints can either be null if none need to be applied, or a list of command-line options to define restraints. These are passed directly to the CABSflex command line.

additionaly, each runtype has a `ligand` section which controls whether to add a ligand by rigid superimposition with a reference structure.
The section has the following options:

|ligand options|Meaning|
|---|---|
|place_ligand|can be either `true` or `false`. Controls whether to place ligand by rigid superimposition; if this is false, the following options will be ignored|
|ligand_structure|structure from which the atoms and position of the ligand will be taken. Can have any name, but the path refers to the directory the `snakefile` is in|
|model_fit_selection|this is an MDAnalysis selection string to select the part of the model that will be used to superimpose model and reference ligand structure|
|ligand_fit_selection|this is an MDAnalysis selection string to select the part of the reference ligand structure that will be used for the superimposition|
|ligand_selection|this is an MDAnalysis selection string to select the residue in the reference ligand structure that corresponds to the ligand|

### Outputs

Directories of typical CABSflex outputs, one per run type, per model are expected for each line within the table, variants.csv.
The folders containing the output will be built and named with the information provided in the table.

Example:

[source_structure]_[chain_in_source]_[aa_num_start]-[aa_num_end]_[WT_amino_acid][AA_num][Mutated_amino_acid]
└── model
    └── metal_bound
        ├── 2XWRa_A_ALA129ASP.pdb
        ├── CABS.log
        ├── config.ini
        ├── input.pdb -> 2XWRa_A_ALA129ASP.pdb
        ├── output_data
        │   ├── all_rmsds_A.txt
        │   ├── DSSP_output_input.txt
        │   └── ... 
        ├── output_pdbs
        │   ├── cluster_0.pdb
        │   ├── cluster_10.pdb
        │   ├── cluster_11.pdb
        │   ├── ...
        │   ├── model_0.pdb
        │   ├── model_10.pdb
        │   ├── model_11.pdb
        │   ├── model_12.pdb        
        │   ├── ...
        │   ├── model_ligand_0.pdb
        │   ├── model_ligand_10.pdb
        │   ├── model_ligand_11.pdb
        │   └── ...
        ├── model_quality
        │   ├── model_0.sov
        │   ├── model_0.ss
        │   ├── ...
        │   ├── model_19.sov
        │   ├── model_19.ss
        │   ├── reference_structure.ss
        │   ├── model_SS_summary.csv
        │   └── ...
        └── theseus
            ├── theseus_ave.pdb
            ├── theseus_LS_tree.nxs
            ├── theseus_ML_tree.nxs
            ├── theseus_residuals.txt
            ├── theseus_sup.pdb
            ├── theseus_transf.txt
            └── theseus_variances.txt

Notice that there is only "model" in this example, however, experimental inputs are also applicable defined by the method, e.g., xray.

### Example

See the directory example_run to see an example of input files, one type of configuration and the resulting output. This directory also 
contains a small readme to explain the example run. 

### References

[^Mölder2021]: Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.

[^Kurcinski2019]: Mateusz Kurcinski, Tymoteusz Oleniecki, Maciej Pawel Ciemny, Aleksander Kuriata, Andrzej Kolinski, Sebastian Kmiecik, CABS-flex standalone: a simulation environment for fast modeling of protein flexibility, Bioinformatics, Volume 35, Issue 4, 15 February 2019, Pages 694–695, https://doi.org/10.1093/bioinformatics/bty685

[^Webb2016]: Webb, B.; Sali, A. Comparative protein structure modeling using MODELLER. Curr. Protoc. Bioinforma.2016, 54, 5.6.1-5.6.37.

[^Kabsch1983]: Kabsch W, Sander C. Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers. 1983 Dec;22(12):2577-637.

[^Theobald2006]: Theobald, Douglas L, and Deborah S Wuttke. “THESEUS: maximum likelihood superpositioning and analysis of macromolecular structures.” Bioinformatics (Oxford, England) vol. 22,17 (2006): 2171-2. doi:10.1093/bioinformatics/btl332

[^Liu2018]: Liu, Tong, and Zheng Wang. "SOV_refine: a further refined definition of segment overlap score and its significance for protein structure similarity." Source code for biology and medicine 13.1 (2018): 1-10.