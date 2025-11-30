# Parallel ChEMBL Retrieval on USC CARC

This guide explains how to run the parallel ChEMBL retrieval on USC CARC computing nodes.

## Overview

The parallel version (`scripts/retrieve_chembl_parallel.py`) is optimized for HPC environments with:

- **Multiprocessing**: Processes multiple UniProt IDs simultaneously
- **Resume capability**: Automatically skips already-processed UniProt IDs
- **Memory efficient**: Can process in chunks to manage memory
- **Logging**: Comprehensive logs for debugging and monitoring
- **SLURM integration**: Ready-to-use job script for CARC

## Performance

**Sequential (original)**:
- 1 UniProt ID per API call
- ~0.3s rate limit between calls
- For 1000 UniProt IDs: ~300 seconds (~5 minutes)

**Parallel (16 workers)**:
- 16 UniProt IDs processed simultaneously
- For 1000 UniProt IDs: ~20 seconds (15x speedup)

**Parallel (32 workers)**:
- 32 UniProt IDs processed simultaneously  
- For 1000 UniProt IDs: ~10 seconds (30x speedup)

## Setup on CARC

### 1. Transfer Files

```bash
# From your local machine
scp scripts/retrieve_chembl_parallel.py <username>@discovery.usc.edu:~/PLATE-VS/scripts/
scp slurm/run_chembl_parallel.slurm <username>@discovery.usc.edu:~/PLATE-VS/slurm/
scp data/annotation_table_with_uniprot.parquet <username>@discovery.usc.edu:~/PLATE-VS/data/
```

### 2. Setup Python Environment

```bash
# SSH into CARC
ssh <username>@discovery.usc.edu

# Navigate to project directory
cd PLATE-VS

# Load Python module
module load python/3.11.3

# Create/activate conda environment
conda create -n chembl_env python=3.11
conda activate chembl_env

# Install required packages
pip install pandas numpy pyarrow chembl-webresource-client requests
```

### 3. Edit SLURM Script

Edit `run_chembl_parallel.slurm`:

```bash
#SBATCH --account=<YOUR_ACCOUNT>          # Replace with your CARC account
#SBATCH --mail-user=<YOUR_EMAIL>          # Replace with your email
```

Adjust resources if needed:
```bash
#SBATCH --cpus-per-task=32                # More CPUs = faster (16, 32, 64)
#SBATCH --mem=64GB                        # More memory if needed
#SBATCH --time=48:00:00                   # Adjust time limit
```

## Running the Job

### Basic Usage

```bash
# Submit the job
sbatch run_chembl_parallel.slurm

# Check job status
squeue -u $USER

# Check job output (while running)
tail -f logs/chembl_retrieval_<JOB_ID>.out

# Cancel job if needed
scancel <JOB_ID>
```

### Custom Parameters

```bash
# Use more workers for faster processing
sbatch --cpus-per-task=32 --mem=64GB run_chembl_parallel.slurm

# Process in chunks (for very large datasets)
python scripts/retrieve_chembl_parallel.py \
    --input data/annotation_table_with_uniprot.parquet \
    --output data/chembl_parallel \
    --workers 32 \
    --chunk-size 100

# Skip final merge (useful for testing)
python scripts/retrieve_chembl_parallel.py \
    --input data/annotation_table_with_uniprot.parquet \
    --output data/chembl_parallel \
    --workers 16 \
    --skip-merge
```

## Output Structure

```
data/chembl_parallel/
├── logs/
│   └── chembl_retrieval_<timestamp>.log
├── uniprot_P12345/
│   ├── P12345_chembl_activities.parquet
│   └── P12345_metadata.json
├── uniprot_Q67890/
│   ├── Q67890_chembl_activities.parquet
│   └── Q67890_metadata.json
├── ...
└── retrieval_summary.json

data/
└── chembl_activities_merged.parquet       # Final merged output
```

## Resume Capability

The script automatically skips already-processed UniProt IDs, so you can:

1. **Stop and resume**: If job times out, just resubmit - it will continue
2. **Process incrementally**: Add new systems and rerun - only new UniProt IDs are processed
3. **Recover from failures**: Failed UniProt IDs can be retried without reprocessing successful ones

## Monitoring Progress

```bash
# Watch log file in real-time
tail -f logs/chembl_retrieval_<JOB_ID>.out

# Check how many UniProt IDs processed so far
ls data/chembl_parallel/uniprot_* | wc -l

# View summary (after completion)
cat data/chembl_parallel/retrieval_summary.json | python -m json.tool
```

## Troubleshooting

### Job runs out of memory
- Increase `--mem` in SLURM script
- Use `--chunk-size` to process in smaller batches
- Reduce `--workers`

### Job times out
- Increase `--time` in SLURM script
- Just resubmit - resume capability will skip completed work

### API rate limiting issues
- ChEMBL API is usually generous, but if you hit limits:
- Reduce `--workers` (fewer simultaneous requests)
- Add delays in the code if needed

### Connection errors
- Usually transient - the script will log errors and continue
- Check logs to see which UniProt IDs failed
- Rerun to retry failed ones

## Performance Tips

1. **Start with fewer workers**: Test with `--workers 8` first
2. **Use chunks for large datasets**: `--chunk-size 200` for 1000+ UniProt IDs
3. **Monitor memory usage**: `sstat -j <JOB_ID> --format=MaxRSS`
4. **Request appropriate time**: Estimate ~1 second per UniProt ID per worker

## Example Workflow

```bash
# 1. Submit job
sbatch run_chembl_parallel.slurm

# 2. Monitor
watch -n 30 'ls data/chembl_parallel/uniprot_* | wc -l'

# 3. After completion, check results
python -c "
import pandas as pd
df = pd.read_parquet('data/chembl_activities_merged.parquet')
print(f'Total activities: {len(df)}')
print(f'Unique compounds: {df[\"molecule_chembl_id\"].nunique()}')
print(f'Unique UniProt IDs: {df[\"mapped_uniprot_id\"].nunique()}')
"
```

## Questions?

Contact CARC support: hpc@usc.edu
