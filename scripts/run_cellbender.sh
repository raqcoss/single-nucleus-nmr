tmux
raquelcr@aleph:~$ conda env create -f /home/raquelcr/seurat/cellbender/celbender-env.yml
conda activate cellbender-env
mkdir -p /home/raquelcr/scanpy/denoised

pendingfiles=ls /home/raquelcr/scanpy/raw_data

for file in "${pendingfiles[@]}"; do
    infile="/home/raquelcr/scanpy/raw_data/$file"
    outfile="/home/raquelcr/scanpy/denoised/${file/_raw.h5ad/_denoised.h5ad}"
    logfile="${outfile%.h5ad}.log"

    echo "Processing $infile → $outfile"
    cellbender remove-background \
        --input "$infile" \
        --output "$outfile" \
        --total-droplets-included 20000 \
        > "$logfile" 2>&1
done