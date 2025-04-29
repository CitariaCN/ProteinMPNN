#!/bin/bash
#SBATCH -n 1                              # request one thread for each job in the array
#SBATCH -t 0-03:00                        # runtime in D-HH:MM format
#SBATCH -p short                   # partition to run in
#SBATCH --mem-per-cpu=4G
#SBATCH --array=0-1
#SBATCH -o  /runlogs/%A_%a.out                     # where to write stdout
#SBATCH -e  /runlogs/%A_%a.err              # where to write stderr

source /net/software/miniforge3/etc/profile.d/conda.sh
conda activate ligandmpnn

python /net/software/protein_design/LigandMPNN/run.py \
        --model_type "protein_mpnn" \
        --temperature 0.05 \
        --pdb_path "/pdbsref/MHET/MHET.pdb" \
        --batch_size 12 \
        --file_ending "proteinmpnn" \
        --out_folder "/protein_mpnn/outputs/MHET/heyfrank" \
        --save_stats 1 \
        --omit_AA "C" \
        --fixed_residues "A92 A93 A94 A185 A186 A187 A215 A218 A219 A358 A372 A376 A377 A378 A380 A381 A385 A453 A455 A456 A489 A490 A1 A2 A12 A15 A34 A39 A53 A55 A78 A79 A82 A84 A85 A90 A92 A93 A117 A119 A123 A126 A127 A166 A173 A174 A180 A181 A184 A185 A186 A188 A189 A190 A191 A192 A196 A197 A200 A204 A205 A208 A209 A210 A211 A212 A264 A265 A268 A269 A272 A273 A278 A281 A292 A330 A340 A434 A435 A438 A442 A443 A444 A449 A450 A453 A464 A465 A478 A480 A482 A485 A486 A489 A490 A493 A495 A505 A508 A509 A510 A512 A514 A515 A535 A538 A541 A548 A550 A560 A12 A53 A185 A264 A281 A301 A309 A490 A538 A560 A186 A453 A489" \
        --checkpoint_protein_mpnn "/net/software/protein_design/LigandMPNN/model_params/proteinmpnn_v_48_020.pt" \
