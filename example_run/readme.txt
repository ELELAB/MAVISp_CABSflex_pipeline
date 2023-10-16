#Readme for the example run. 

module load python/3.7/modulefile
module load cabs/0.9.18/modulefile

snakemake --cores 1

#This example is run on three variants and the wildtype of 2XWR chain A. A structure covering the 
#DNA binding domain of the protein p53. 

#Only no restrains are used. 

#Here the variant and the wildtype is present in the directory input_structures. 
#This placement is reflected in the variant.csv file, where it is evident in 
#the last column "input_structures" and in config.yaml, where the run directory is defined. 


