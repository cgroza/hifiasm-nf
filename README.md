# hifiasm Nextflow pipeline on Rorqual

This pipeline assembles PacBio HiFi reads with hifiasm. It can run either an
unphased HiFi assembly or a trio assembly using parental reads and yak. Existing
haplotype assemblies can also be supplied for the optional variant-calling
steps.

The instructions below target the Alliance of BiCompute Canada **Rorqual**
cluster and its Slurm scheduler.

## Requirements

- An Alliance account with [access to Rorqual](https://docs.alliancecan.ca/wiki/Rorqual)
- A Slurm account such as `def-pi`
- Nextflow and Apptainer modules on Rorqual
- A prebuilt Apptainer image made from [`container.def`](container.def)
- Input data in `/project` or `/scratch`

Rorqual compute nodes do not have internet access. Clone the pipeline and
transfer the container and input files before submitting the workflow. Do not
place large inputs, containers, or Nextflow work directories in `$HOME`.

## 1. Build and transfer the container

Build the image on a Linux system where you have root or working Apptainer
fakeroot support:

```bash
git clone https://github.com/cgroza/hifiasm_nextflow.git
cd hifiasm_nextflow
sudo apptainer build hifiasm-nextflow.sif container.def
```

Transfer the repository and image to Rorqual with Globus, `rsync`, or `scp`.
For example, after logging in to Rorqual:

```bash
mkdir -p /project/def-pi/$USER/hifiasm-nextflow
cd /project/def-pi/$USER/hifiasm-nextflow
git clone https://github.com/cgroza/hifiasm_nextflow.git pipeline
mkdir -p containers
# Transfer hifiasm-nextflow.sif into containers/ outside this shell.
```

Building off-cluster is recommended because the definition downloads and
compiles software, while Rorqual compute nodes cannot access the internet.

## 2. Prepare a design CSV

Use absolute paths so that Nextflow can stage every input reliably.

### HiFi assembly

```csv
sample,hifi
sample1,/project/def-pi/data/sample1.hifi.fastq.gz
sample2,/project/def-pi/data/sample2.hifi.fastq.gz
```

### Trio assembly with paired-end parental reads

```csv
sample,hifi,maternal_1,maternal_2,paternal_1,paternal_2
sample1,/project/def-pi/data/sample1.hifi.fastq.gz,/project/def-pi/data/mother_R1.fastq.gz,/project/def-pi/data/mother_R2.fastq.gz,/project/def-pi/data/father_R1.fastq.gz,/project/def-pi/data/father_R2.fastq.gz
```

### Trio assembly with HiFi parental reads

```csv
sample,hifi,maternal,paternal
sample1,/project/def-pi/data/sample1.hifi.fastq.gz,/project/def-pi/data/mother.hifi.fastq.gz,/project/def-pi/data/father.hifi.fastq.gz
```

Sample names become directory and file-name components. Keep them unique and
use only letters, numbers, periods, underscores, and hyphens.

## 3. Add Rorqual scheduler settings

Create `rorqual.config` beside the pipeline. Replace `def-pi` with the Slurm
account shown under **My Projects > My Resources and Allocations** in the
Alliance account portal.

```groovy
process {
    clusterOptions = '--account=def-pi'
}

executor {
    queueSize = 100
    submitRateLimit = '60/1min'
}
```

The repository's `cluster` profile selects the Slurm executor and runs each
task in `$SLURM_TMPDIR`. The extra configuration supplies the account and
respects the Alliance-recommended submission limits.

Processes request 40 CPUs and 160 GB by default. Override both for all tasks
with `--cpus` and `--memory`. Rorqual's standard CPU nodes provide 192 cores
and up to 750 GB of allocatable memory.

## 4. Create the Nextflow launcher

Run the lightweight Nextflow controller on a Rorqual login node inside
`tmux`. The controller remains on the login node and submits the computational
processes to Slurm through the `cluster` profile. Save the following as
`run-hifiasm.sh` and update the paths:

```bash
#!/bin/bash

set -euo pipefail

module purge
module load nextflow
module load apptainer

export NXF_OPTS='-XX:ActiveProcessorCount=2 -Dleveldb.mmap=false'
export NXF_HOME="${HOME}/.nextflow"

PIPELINE_DIR="/project/def-pi/${USER}/hifiasm-nextflow/pipeline"
RUN_DIR="/project/def-pi/${USER}/hifiasm-nextflow/run"
CONTAINER="/project/def-pi/${USER}/hifiasm-nextflow/containers/hifiasm-nextflow.sif"
DESIGN="${RUN_DIR}/design.csv"
WORK_DIR="${HOME}/scratch/hifiasm-nextflow-work"
export NXF_TEMP="${HOME}/scratch/hifiasm-nextflow-tmp"

mkdir -p "${RUN_DIR}" "${WORK_DIR}" "${NXF_TEMP}"
cd "${RUN_DIR}"

nextflow run "${PIPELINE_DIR}/main.nf" \
    -profile cluster \
    -c "${PIPELINE_DIR}/rorqual.config" \
    -with-apptainer "${CONTAINER}" \
    -work-dir "${WORK_DIR}" \
    --design "${DESIGN}" \
    --cpus 40 \
    --memory 160G
```

Start a named `tmux` session on the login node and run the controller:

```bash
tmux new -s hifiasm
bash run-hifiasm.sh 2>&1 | tee hifiasm-nextflow.log
```

Detach from the session with <kbd>Ctrl-b</kbd>, then <kbd>d</kbd>. The
controller continues running after the SSH connection closes. Reattach and
monitor the submitted Slurm jobs with:

```bash
tmux attach -t hifiasm
sq
```

Do not submit `run-hifiasm.sh` with `sbatch`: that would place the Nextflow
controller on a compute node. The `tmux` session must remain running until all
pipeline jobs finish.

## Run modes

Change the final `nextflow run` arguments in the launcher for the desired mode.

**HiFi assembly:**

```text
--design /absolute/path/design.csv
```

**Trio assembly with paired-end parental reads:**

```text
--design /absolute/path/design.csv --trio true
```

**Trio assembly with HiFi parental reads:**

```text
--design /absolute/path/design.csv --trio true --hifi_parents true
```

**Use existing haplotype assemblies instead of assembling reads:**

```text
--assemblies /absolute/path/assemblies.csv
```

The assemblies CSV must have this form:

```csv
sample,hap1,hap2
sample1,/project/def-pi/data/sample1_hap1.fa.gz,/project/def-pi/data/sample1_hap2.fa.gz
```

**Run dipcall after assembly:**

```text
--design /absolute/path/design.csv --dipcall true --ref /absolute/path/reference.fa
```

**Run SVIM-asm after assembly:**

```text
--design /absolute/path/design.csv --svim_asm true --ref /absolute/path/reference.fa
```

`--dipcall` and `--svim_asm` are mutually exclusive in the workflow. Variant
calling requires a reference FASTA whose contig names match the expected
genome build. The SNV portion of the SVIM-asm mode currently expects human
contigs `chr1` through `chr22` and `chrX`.

## Outputs and working data

Outputs are published relative to the launch directory:

- `assemblies/`: bgzip-compressed haplotype FASTA files
- `variants/`: compressed variant calls when a variant mode is enabled
- `.nextflow.log`: controller log
- `hifiasm-nextflow.log`: controller console log

The Nextflow cache and intermediate task data are stored under the path passed
to `-work-dir`. Rorqual scratch storage is large but is not backed up and is
purged, so copy or publish final results under `/project`.

## Resume or restart

If the controller or a task fails, fix the cause and run the same command with
`-resume`:

```text
nextflow run ... -resume
```

Keep the same launch directory and `-work-dir`; both the `.nextflow` cache and
work files are needed to reuse completed tasks. Do not delete the work
directory until the results have been checked and backed up.

## Rorqual notes

- Run only the lightweight Nextflow controller on the login node. The
  `cluster` profile submits all computational processes to Slurm; never run
  those processes directly on the login node.
- Rorqual jobs can run for at most 7 days, and each user may have at most 1,000
  running plus pending jobs.
- `$SLURM_TMPDIR` provides node-local task scratch. If an assembly needs more
  than the default local capacity, append an appropriate request such as
  `--tmp=1000G` to `process.clusterOptions` in `rorqual.config`. Rorqual accepts
  values from 370 GB through 3360 GB.
- If jobs remain pending, use `sq` to inspect the scheduler reason. Large
  memory or long time requests can increase queue time.
- Review a failed task's `.command.err`, `.command.log`, and `.command.sh`
  inside its Nextflow work directory before resuming.

See the Alliance documentation for
[Rorqual](https://docs.alliancecan.ca/wiki/Rorqual),
[Nextflow](https://docs.alliancecan.ca/wiki/Nextflow), and
[running Slurm jobs](https://docs.alliancecan.ca/wiki/Running_jobs).
