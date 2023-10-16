        
        ln -s 2XWR_A_WT.pdb input.pdb
        set +u; source /usr/local/envs/cabsflex/bin/activate; set -u

        CABSflex            -i input.pdb            -k 20            -v 3            -s 50            -y 50            -a 20            -A            -z 10            -C             --dssp-command /usr/local/dssp-3.0.10/bin/mkdssp                        --log
