# sample-input

rnaseq-quant takes FASTQ files (`reads`) and an index from rnaseq-index (`index`). The demo
data of the bulk RNA-seq workflows is Discover dataset 427 (RE-JOIN; CC BY 4.0): 12
single-end libraries of mouse trigeminal ganglion neurons (female/male x naive/tongue tumour,
3 pools each), with the authors' gene counts and SPARC `samples.xlsx` / `subjects.xlsx`.
`scripts/fetch-sample-data.sh` in `workflow-bulk-rnaseq` downloads it with GENCODE M39:

```sh
cd ../workflow-bulk-rnaseq && AWS_PROFILE=<any> scripts/fetch-sample-data.sh
ls data/tg-tumour/primary/        # pool-FN1 ... pool-FT3, pool-N1..N3, pool-T1, pool-T2, pool-MT
```

Each `pool-<id>/` holds one library `<id>_S<n>_L007_R1.fastq.gz` (~20-30 M reads) and the
authors' `<id>_Genes_ReadCount.txt`. Each library pools the trigeminal ganglia of two mice;
rnaseq-quant joins both mice's SPARC rows by the `pool-` folder (values that differ, such as the
sample and subject ids, are joined with `;`).

## 1 M-read subsample (what the workflow builders test on)

```sh
D=../workflow-bulk-rnaseq/data
mkdir -p $D/work/tg-tumour/reads-1m && cp $D/tg-tumour/{samples,subjects}.xlsx $D/work/tg-tumour/reads-1m/
for f in $D/tg-tumour/primary/pool-*/*_R1.fastq.gz; do
  o=$D/work/tg-tumour/reads-1m/${f#$D/tg-tumour/primary/}; mkdir -p $(dirname $o)
  seqtk sample -s42 $f 1000000 | gzip -c > $o
done
```

Then, with the index from rnaseq-index in `$D/work/gencode-mouse-M39/index/`:

```sh
docker run --rm --platform linux/amd64 --memory 6g --cpus 4 --network none --read-only --tmpfs /tmp \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output \
  -v $D/work/tg-tumour/reads-1m:/data/input/reads:ro \
  -v $D/work/gencode-mouse-M39/index:/data/input/index:ro \
  -v $D/work/tg-tumour/quant:/data/output rnaseq-quant:dev
```

Outputs land in `data/work/tg-tumour/quant/{counts,transcripts}/`; `counts/` feeds bulk-de.

## Your own data

Put the FASTQ files anywhere under the reads input. Illumina names group themselves
(`<sample>_S1_L001_R1_001.fastq.gz`, `<sample>_R1.fastq.gz`, `<sample>_1.fq.gz`; lanes of a
sample are combined). Otherwise, or to add metadata, add `samples.csv`:

```csv
sample_id,fastq_1,fastq_2,condition,batch
ctrl_1,run1/ctrl1_L001_R1.fq.gz,run1/ctrl1_L001_R2.fq.gz,control,b1
ctrl_1,run1/ctrl1_L002_R1.fq.gz,run1/ctrl1_L002_R2.fq.gz,control,b1
treated_1,run1/t1_R1.fq.gz,run1/t1_R2.fq.gz,treated,b1
```

(one row per sample, or per lane; `fastq_2` empty for single-end; paths equal to, or the end
of, the file's path under the input).

## Synthetic fixture (no download; needs salmon on PATH to build the fixture's index)

```sh
make fixture PYTHON=.venv/bin/python    # data/fixture/{reads,index}/
make run PYTHON=.venv/bin/python        # Docker over data/fixture/ -> data/output/
```
