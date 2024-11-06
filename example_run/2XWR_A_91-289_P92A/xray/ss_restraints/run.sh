        
        ln -s 2XWR_A_P92A.pdb input.pdb
        set +u; module unload python; module load cabs; set -u
        
        CABSflex            -i input.pdb            -k 20            -v 3            -s 50            -y 50            -a 20            -A            -z 10            -g all 3 3.8 8.0            -C            --dssp-command /usr/local/dssp-3.0.10/bin/mkdssp                        --log
        
        set +u; module unload cabs; module load python; set -u
