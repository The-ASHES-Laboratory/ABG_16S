# Scripts

## demultiplex_summarize.sh

Summarize demultiplexed reads in QIIME 2 and optionally import sequencing data from a manifest file.

### Usage

```bash
# Summarize an existing demultiplexed artifact
./demultiplex_summarize.sh -i demux.qza -o output_dir -l runA

# Import from a manifest then summarize
./demultiplex_summarize.sh -M manifest.csv -o output_dir -l runB
```

The script requires an active QIIME 2 environment and produces a summary visualization (`<label>-summary.qzv`). If a manifest is supplied, the imported demultiplexed artifact (`<label>.qza`) is also generated.
