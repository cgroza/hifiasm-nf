# hifiasm Nextflow pipeline on Rorqual

The `hifi-trios` branch assembles PacBio HiFi BAMs with hifiasm. It converts each
BAM to FASTQ, assembles maternal and paternal samples independently, and
assembles probands in trio-binning mode using parental yak k-mers from the same
family. Existing haplotype assemblies can also be supplied for the optional
variant-calling steps.

The instructions below target the Alliance of BiCompute Canada **Rorqual**
cluster and its Slurm scheduler. For Azure CycleCloud-managed Slurm clusters,
see [CycleCloud Slurm](#cyclecloud-slurm).

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

### Docker alternative

The [`Dockerfile`](Dockerfile) builds a task image containing hifiasm, yak,
samtools, minimap2, bcftools, bgzip, GNU Parallel, SVIM-asm, the SNV parsing
script and its Python dependencies, and the complete dipcall toolkit. Build
from the repository root:

```bash
docker build --platform linux/amd64 -t hifiasm-nextflow:local .
```

The amd64 platform is required by the official dipcall release binaries;
Docker on Apple Silicon must use amd64 emulation. The image does not include
the Nextflow controller. On a machine with Docker and Nextflow installed, run
the pipeline with Apptainer explicitly disabled:

```bash
printf 'apptainer.enabled = false\n' > docker.config
nextflow run main.nf -profile standard \
    -c docker.config \
    -with-docker hifiasm-nextflow:local \
    --design /absolute/path/design.csv
```

Docker is not used on Rorqual compute nodes; use the Apptainer instructions
above for cluster execution.

## 2. Prepare a design CSV

Use absolute paths so that Nextflow can stage every input reliably.

### HiFi families

```csv
BAM,sample_id,sample_family,role
/project/def-pi/data/mother.hifi.bam,mother,family1,mat
/project/def-pi/data/father.hifi.bam,father,family1,pat
/project/def-pi/data/child.hifi.bam,child,family1,proband
/project/def-pi/data/sibling.hifi.bam,sibling,family1,proband
```

Provide one HiFi BAM per sample. `sample_id` must be globally unique.
`sample_family` connects each proband to its parents, and `role` must be `mat`,
`pat`, or `proband`. Every family containing a proband must have exactly one
`mat` and one `pat`. Multiple probands can share the same parents; their parental
k-mer databases are computed once and reused. Parents are always assembled
independently, even in families without a proband. Duplicate parent roles,
missing parents, invalid roles, missing BAMs, and duplicate sample IDs fail
validation before tasks are launched.

Sample and family identifiers must start with a letter or number and contain
only letters, numbers, periods, underscores, and hyphens. Output directory and
file names use `sample_id`, not `sample_family`. In trio assemblies, haplotype
1 is paternal (`hifiasm -1`) and haplotype 2 is maternal (`hifiasm -2`).

This branch replaces the previous FASTQ design formats and the `--trio` and
`--hifi_parents` switches. Roles in the CSV now select the assembly mode.

## 3. Add Rorqual scheduler settings

Create `rorqual.config` beside the pipeline. Replace `def-pi` with the Slurm
account shown under **My Projects > My Resources and Allocations** in the
Alliance account portal.

```groovy
docker.enabled = false
apptainer.enabled = true
apptainer.autoMounts = true
apptainer.runOptions = '--contain --bind $(pwd):/tmp'

process {
    clusterOptions = '--account=def-pi'
}

executor {
    queueSize = 10
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

## CycleCloud Slurm

Use `-profile cyclecloud` on an Azure CycleCloud-managed Slurm cluster. This
profile submits tasks to Slurm, runs them with Docker using
`cgroza/nf-hifiasm:latest`, and uses `/scratch` for node-local scratch.
Ensure that
this filesystem has enough space for BAM conversion and assembly; configure
`process.scratch` to a larger node-local path in the site configuration if
necessary.

Install Nextflow and Java on the submit/login node, and Docker on all compute
nodes. The job user must have permission to access the running Docker daemon.
Slurm commands such as `sbatch`, `squeue`, and `scancel`
must be available to the controller. Use **x86_64 compute nodes** because the
published container includes amd64 dipcall binaries. Keep the repository,
design CSV, BAMs, reference, launch directory, and Nextflow
work directory on storage mounted at the same paths on all nodes.

Create `cyclecloud.config` for your cluster's partition and, if required,
account. Replace the example values with those configured in Slurm:

```groovy
process {
    queue = 'your-partition'
    // Add only if your cluster requires a Slurm account:
    // clusterOptions = '--account=your-account'
}
```

The profile otherwise uses Slurm's default partition and account. Select a
partition whose CycleCloud node configuration provides the requested CPUs,
memory, and local disk. Each process requests 40 CPUs and 160 GB by default.

Compute nodes must be able to pull the image from Docker Hub, or have it
preloaded with `docker load` before launching. The following uses `/shared`
as an example mount; replace it with your cluster's actual shared filesystem:

```bash
mkdir -p /shared/hifiasm/run /shared/hifiasm/work

cd /shared/hifiasm/run
tmux new -s hifiasm
nextflow run /shared/hifiasm/pipeline/main.nf \
    -profile cyclecloud \
    -c /shared/hifiasm/pipeline/cyclecloud.config \
    -work-dir /shared/hifiasm/work \
    --design /shared/hifiasm/run/design.csv \
    --cpus 40 \
    --memory 160G
```

Run the controller on the submit/login node inside `tmux`, not in an `sbatch`
job. Remove any previous `-with-apptainer` overrides from your launcher.
Use `squeue` to monitor jobs and
`-resume` from the same launch and work directories to resume a run.

## Run modes

Change the final `nextflow run` arguments in the launcher for the desired mode.

**Assemble parents and trio-bin probands from HiFi BAMs:**

```text
--design /absolute/path/design.csv
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

- `assemblies/<sample_id>_diploid/<sample_id>_hap{1,2}.fa.gz`:
  bgzip-compressed haplotype FASTA files for every parent and proband
- `variants/`: compressed variant calls when a variant mode is enabled
- `.nextflow.log`: controller log
- `hifiasm-nextflow.log`: controller console log

Set `--assemblies_outdir /path/to/assemblies` to change the assembly output
directory (default: `assemblies`). Relative paths are resolved from the launch
directory. The `<sample_id>_diploid/` subdirectories and FASTA filenames are
unchanged. This also applies to renamed `--assemblies` inputs and does not
change the variant output or Nextflow work directories.

The downstream `rename_contigs_pansn` process renames contigs in both haplotype
FASTAs to the PanSN format `sample_id#haplotype#contig`. For example,
`>h1tg000001l` becomes `>child#1#h1tg000001l` in haplotype 1 and
`>h2tg000001l` becomes `>child#2#h2tg000001l` in haplotype 2. Sequences, contig
order, and header descriptions are preserved; output files remain bgzip
compressed with the same sample-based filenames. Only renamed assemblies are
published, and optional variant calling uses these renamed FASTAs.

Renaming also applies to `--assemblies` inputs, using the CSV's `sample`
column as the sample identifier. Both plain and gzip-compressed input FASTAs
are supported. Sample identifiers follow the same rules as `sample_id` in
the design CSV, and input contig identifiers must not contain `#` (including
already PanSN-named inputs). Original input files are never modified.

The Nextflow cache and intermediate task data are stored under the path passed
to `-work-dir`. Rorqual scratch storage is large but is not backed up and is
purged, so copy or publish final results under `/project`.

## Resume or restart

If the controller or a task fails, fix the cause and run the same command with
`-resume`:

```text
nextflow run ... -resume
```

## Testing

With Nextflow, Python 3, samtools, and bgzip on `PATH`, run:

```bash
NXF_SYNTAX_PARSER=v1 python3 -m unittest discover -s tests -v
NXF_SYNTAX_PARSER=v2 python3 -m unittest discover -s tests -v
```

Run both commands to check compatibility with the legacy and strict parsers
(the strict parser requires Nextflow 25.04 or later). Tests default to the
strict parser when `NXF_SYNTAX_PARSER` is unset.

The integration tests use real BAM-to-FASTQ conversion and compression, mock
hifiasm and yak to check commands and family routing without a large assembly,
and use Nextflow stubs to check optional variant-calling outputs. They cover
multiple families, siblings, parent-only inputs, invalid designs, corrupt BAMs,
existing assemblies, and PanSN renaming with sequence preservation. Set
`NEXTFLOW` to an executable path if it is not on
`PATH`.

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
