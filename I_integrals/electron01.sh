#!/bin/bash
#SBATCH --job-name=omega_job
#SBATCH --array=0-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=50:00:00
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

# Load the Python module and install necessary packages in your user environment
module load python/3.12
python -m pip install --user numpy fitsio astropy

# Set parameters
x=1840        # Example value for x


# Set y to 20.0 (overriding any calculation)
END_Y=1999
STEP=1

# Generate the array of omega scales (from 0.01 to END_Y in steps of STEP)
mapfile -t y < <(seq 0 $STEP $END_Y)
num_scales=${#y[@]}

# Set block size and compute total number of computation blocks
block_size=400
TOTAL_BLOCKS=$(( (num_scales + block_size - 1) / block_size ))
echo "Total omega scales: $num_scales"
echo "Block size: $block_size"
echo "Total computation blocks: $TOTAL_BLOCKS"

BASE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
RESULTDIR="${BASE_DIR}/results"
mkdir -p "$RESULTDIR"

if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    # MERGE TASK
    echo "Merge task: Checking for completed blocks in $RESULTDIR ..."
    while true; do
        missing_blocks=()
        for i in $(seq 1 $TOTAL_BLOCKS); do
            start_idx=$(( (i - 1) * block_size ))
            end_idx=$(( i * block_size - 1 ))
            if [ $end_idx -ge $num_scales ]; then
                end_idx=$(( num_scales - 1 ))
            fi
            block_complete=true
            for (( j = start_idx; j <= end_idx; j++ )); do
                omega_int=${y[$j]}                       # 10, 20, … 2000  (already int)
                file_name="${RESULTDIR}/output_x${x}_omega${omega_int}.fits"

                if [ ! -f "$file_name" ]; then           # file missing → block incomplete
                    block_complete=false
                    break
                fi
            done

            if [ "$block_complete" = false ]; then
                missing_blocks+=($i)
            fi
        done

        if [ ${#missing_blocks[@]} -eq 0 ]; then
            echo "All blocks complete. Proceeding with merge..."
            break
        else
            echo "Still missing blocks: ${missing_blocks[*]}"
            sleep 60
        fi
    done

    # Run the merge script, which should combine all the output files from RESULTDIR.
    python -u merge_fits.py "$x" "$RESULTDIR"
    cp merged_output_x${x}.fits "$SLURM_SUBMIT_DIR"
    cp $TMPDIR/*.npy "$SLURM_SUBMIT_DIR"
else
    # COMPUTATION TASKS: SLURM_ARRAY_TASK_ID is 1 to TOTAL_BLOCKS
    task_id=$SLURM_ARRAY_TASK_ID
    start_idx=$(( (task_id - 1) * block_size ))
    end_idx=$(( task_id * block_size - 1 ))
    if [ $end_idx -ge $num_scales ]; then
         end_idx=$(( num_scales - 1 ))
    fi

    echo "Computation task $task_id: initially assigned indices $start_idx to $end_idx (total scales: $num_scales)."
    
    # Determine the minimal missing index in this block.
    new_start_idx=$start_idx
    missing_found=false
    for (( i = start_idx; i <= end_idx; i++ )); do
        omega_int=${y[$i]}                               # 10, 20, … 2000
        file_name="${RESULTDIR}/output_x${x}_omega${omega_int}.fits"

        if [ ! -f "$file_name" ]; then                   # file absent → we must run it
            new_start_idx=$i
            missing_found=true
            echo "First missing ω found at index $i  (omega=${omega_int})."
            break
        fi
    done

    # If no missing file was found, set new_start_idx to end_idx + 1.
    if [ "$missing_found" = false ]; then
        new_start_idx=$(( end_idx + 1 ))
    fi

    if [ $new_start_idx -gt $end_idx ]; then
        echo "Block $task_id is already complete. Skipping."
        exit 0
    fi
    
    echo "Computation task $task_id: processing missing omega scales from index $new_start_idx to $end_idx."

    # Change directory to RESULTDIR so that output files are saved there.
    cd "$RESULTDIR" || { echo "Cannot change to $RESULTDIR"; exit 1; }
    
    for (( i = new_start_idx; i <= end_idx; i++ )); do
         current_omega="${y[$i]}" 
         echo "Running makeintegrand.py for x=$x and omega_scale=$current_omega"
         python -u "$SLURM_SUBMIT_DIR/makeintegrand.py" "$x" "$current_omega"
    done

    # Wait briefly to ensure all file operations are complete
    sleep 5

    # Optionally, if this block is now complete, create a block-level flag.
    touch "${RESULTDIR}/output_block_${task_id}.fits"
fi

