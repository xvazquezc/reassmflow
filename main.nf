nextflow.enable.dsl = 2

process CONCAT_REFERENCE {
  tag { id }
  label 'mapping'

  input:
  tuple val(id), path(refdir)

  output:
  tuple val(id), path('reference.fasta'), emit: ref

  script:
  """find -L ${refdir} -maxdepth 1 -type f \\( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \\) -print -quit | grep -q . || { echo "No FASTA files found in ${refdir}" >&2; exit 1; }
  find -L ${refdir} -maxdepth 1 -type f \\( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \\) -print0 | sort -z | xargs -0 cat | seqkit rename -n > reference.fasta"""

  stub:
  """printf '>stub_reference\\nACGT\\n' > reference.fasta"""
}
process MINIMAP2_MAP {
  tag { id }
  label 'mapping'

  input:
  tuple val(id), path(reads), path(ref)

  output:
  tuple val(id), path('*.bam'), emit: bam

  script:
  """minimap2 -t ${task.cpus} -ax lr:hq ${ref} ${reads} | samtools view -@ ${task.cpus} -bF 4 - | samtools sort -@ ${task.cpus} -o ${id}.${task.index}.bam -"""

  stub:
  """touch ${id}.${task.index}.bam"""
}
process MERGE_LONG_BAMS {
  tag { id }
  label 'mapping'

  input:
  tuple val(id), path(bams)

  output:
  tuple val(id), path('long_recovered.fastq.gz'), emit: reads

  script:
  """samtools merge -@ ${task.cpus} - ${bams.join(' ')} | samtools fastq -@ ${task.cpus} - | pigz -p ${task.cpus} > long_recovered.fastq.gz"""

  stub:
  """touch long_recovered.fastq.gz"""
}
process FLYE_META_ASSEMBLY {
  tag { id }
  label 'lr_assm'
  publishDir { "${params.outdir}/${id}/flye" }, mode: 'copy'

  input:
  tuple val(id), path(reads)

  output:
  path 'flye_out'

  script:
  """state_dir='${params.outdir}/.resume/${id}/flye'
  mkdir -p "\$state_dir"
  if [[ -f "\$state_dir/params.json" ]]; then
    flye -t ${task.cpus} --meta --nano-corr ${reads} --resume -o "\$state_dir"
  else
    flye -t ${task.cpus} --meta --nano-corr ${reads} -o "\$state_dir"
  fi
  find "\$state_dir" -type f -name '*dump' -delete
  cp -a "\$state_dir" flye_out"""

  stub:
  """mkdir -p flye_out; touch flye_out/assembly.fasta"""
}
process MYLOASM_ASSEMBLY {
  tag { id }
  label 'lr_assm'
  publishDir { "${params.outdir}/${id}/myloasm" }, mode: 'copy'

  input:
  tuple val(id), path(reads)

  output:
  path 'mylo_out'

  script:
  """state_dir='${params.outdir}/.resume/${id}/myloasm'
  mkdir -p "\$(dirname "\$state_dir")"
  if [[ -d "\$state_dir" ]]; then
    myloasm exist -t ${task.cpus} -o "\$state_dir"
  else
    myloasm ${reads} -t ${task.cpus} -o "\$state_dir"
  fi
  rm -rf "\$state_dir/binary_temp"
  mylotools annotate-gfa --gfa "\$state_dir/final_contig_graph.gfa" --fasta "\$state_dir/assembly_primary.fa" --output "\$state_dir/final_contig_graph.annot.gfa"
  cp -a "\$state_dir" mylo_out"""

  stub:
  """mkdir -p mylo_out; touch mylo_out/assembly_primary.fa mylo_out/final_contig_graph.annot.gfa"""
}
process BOWTIE2_BUILD {
  tag { id }
  label 'mapping'

  input:
  tuple val(id), path(ref)

  output:
  tuple val(id), path('reference.fasta.*'), emit: index

  script:
  """bowtie2-build ${ref} reference.fasta"""

  stub:
  """touch reference.fasta.1.bt2"""
}
process BOWTIE2_MAP {
  tag { "${id}:${lib}" }
  label 'mapping'

  input:
  tuple val(id), val(lib), path(r1), path(r2), path(single), path(index)

  output:
  tuple val(id), path('*.bam'), emit: bam

  script:
  """bowtie2 -p ${task.cpus} -x reference.fasta -1 ${r1} -2 ${r2} -U ${single} | samtools view -@ ${task.cpus} -bF 12 - | samtools sort -@ ${task.cpus} -o ${lib}.bam -"""

  stub:
  """touch ${lib}.bam"""
}
process MERGE_SHORT_BAMS {
  tag { id }
  label 'mapping'

  input:
  tuple val(id), path(bams)

  output:
  tuple val(id), path('sr_R1.fastq.gz'), path('sr_R2.fastq.gz'), path('sr_singles.fastq.gz'), emit: reads

  script:
  """samtools merge -@ ${task.cpus} - ${bams.join(' ')} | samtools sort -@ ${task.cpus} -n | samtools fastq -@ ${task.cpus} -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -s sr_singles.fastq.gz -0 /dev/null -"""

  stub:
  """touch sr_R1.fastq.gz sr_R2.fastq.gz sr_singles.fastq.gz"""
}
process MEGAHIT_DEFAULT_ASSEMBLY {
  tag { id }
  label 'sr_assm'
  publishDir { "${params.outdir}/${id}/megahit-default" }, mode: 'copy'

  input:
  tuple val(id), path(r1), path(r2), path(single)

  output:
  path 'megahit-default'

  script:
  """state_dir='${params.outdir}/.resume/${id}/megahit-default'
  mkdir -p "\$(dirname "\$state_dir")"
  single_arg=()
  if [[ -n "\$(zcat -f ${single} 2>/dev/null | head -c1)" ]]; then
    single_arg=(-r ${single})
  fi
  if [[ -d "\$state_dir" ]]; then
    megahit -t ${task.cpus} -1 ${r1} -2 ${r2} "\${single_arg[@]}" --continue -o "\$state_dir"
  else
    megahit -t ${task.cpus} -1 ${r1} -2 ${r2} "\${single_arg[@]}" -o "\$state_dir"
  fi
  kmer=\$(find "\$state_dir/intermediate_contigs" -name 'k*.contigs.fa' -printf '%f\\n' | sed -E 's/^k([0-9]+)\\.contigs\\.fa\$/\\1/' | sort -n | tail -1)
  megahit_toolkit contig2fastg "\$kmer" "\$state_dir/intermediate_contigs/k\${kmer}.contigs.fa" > "\$state_dir/k\${kmer}.contigs.fastg"
  cp -a "\$state_dir" megahit-default"""

  stub:
  """mkdir -p megahit-default; touch megahit-default/final.contigs.fa"""
}
process SPADES_META_HYBRID {
  tag { id }
  label 'sr_assm'
  publishDir { "${params.outdir}/${id}/spades-meta-hybrid" }, mode: 'copy'

  input:
  tuple val(id), path(r1), path(r2), path(s), path(lr)

  output:
  path 'spades-meta-hybrid'

  script:
  """state_dir='${params.outdir}/.resume/${id}/spades-meta-hybrid'
  mkdir -p "\$(dirname "\$state_dir")"
  if [[ -d "\$state_dir" ]]; then
    spades.py --continue -o "\$state_dir"
  else
    spades.py -t ${task.cpus} -m ${(task.memory.toGiga() * 0.9) as int} --meta -1 ${r1} -2 ${r2} --nanopore ${lr} -o "\$state_dir"
  fi
  cp -a "\$state_dir" spades-meta-hybrid"""

  stub:
  """mkdir -p spades-meta-hybrid; touch spades-meta-hybrid/contigs.fasta"""
}
process SPADES_META_HYBRID_KEXT {
  tag { id }
  label 'sr_assm'
  publishDir { "${params.outdir}/${id}/spades-meta-hybrid-kext" }, mode: 'copy'

  input:
  tuple val(id), path(r1), path(r2), path(s), path(lr)

  output:
  path 'spades-meta-hybrid-kext'

  script:
  """state_dir='${params.outdir}/.resume/${id}/spades-meta-hybrid-kext'
  mkdir -p "\$(dirname "\$state_dir")"
  if [[ -d "\$state_dir" ]]; then
    spades.py --continue -o "\$state_dir"
  else
    spades.py -t ${task.cpus} -m ${(task.memory.toGiga() * 0.9) as int} --meta -1 ${r1} -2 ${r2} -k 21,33,55,77,101,127 --nanopore ${lr} -o "\$state_dir"
  fi
  cp -a "\$state_dir" spades-meta-hybrid-kext"""

  stub:
  """mkdir -p spades-meta-hybrid-kext; touch spades-meta-hybrid-kext/contigs.fasta"""
}
process SPADES_META_KEXT {
  tag { id }
  label 'sr_assm'
  publishDir { "${params.outdir}/${id}/spades-meta-kext" }, mode: 'copy'

  input:
  tuple val(id), path(r1), path(r2), path(s), path(lr)

  output:
  path 'spades-meta-kext'

  script:
  """state_dir='${params.outdir}/.resume/${id}/spades-meta-kext'
  mkdir -p "\$(dirname "\$state_dir")"
  if [[ -d "\$state_dir" ]]; then
    spades.py --continue -o "\$state_dir"
  else
    spades.py -t ${task.cpus} -m ${(task.memory.toGiga() * 0.9) as int} --meta -1 ${r1} -2 ${r2} -k 21,33,55,77,101,127 --nanopore ${lr} -o "\$state_dir"
  fi
  cp -a "\$state_dir" spades-meta-kext"""

  stub:
  """mkdir -p spades-meta-kext; touch spades-meta-kext/contigs.fasta"""
}
process SPADES_META_KEXT_TRUSTED {
  tag { id }
  label 'sr_assm'
  publishDir { "${params.outdir}/${id}/spades-meta-kext-trusted" }, mode: 'copy'

  input:
  tuple val(id), path(r1), path(r2), path(s), path(lr), path(trusted)

  output:
  path 'spades-meta-kext-trusted'

  script:
  """state_dir='${params.outdir}/.resume/${id}/spades-meta-kext-trusted'
  mkdir -p "\$(dirname "\$state_dir")"
  if [[ -d "\$state_dir" ]]; then
    spades.py --continue -o "\$state_dir"
  else
    spades.py -t ${task.cpus} --meta -m ${(task.memory.toGiga() * 0.9) as int} -1 ${r1} -2 ${r2} --trusted-contigs ${trusted} -k 21,33,55,77,101,127 -o "\$state_dir"
  fi
  cp -a "\$state_dir" spades-meta-kext-trusted"""

  stub:
  """mkdir -p spades-meta-kext-trusted; touch spades-meta-kext-trusted/contigs.fasta"""
}

workflow {
  if (!params.input_csv) {
    error('Specify input_csv in a parameter file or with --input_csv.')
  }

  def raw = params.raw_libraries ?: [:]
  def globalLibraries = params.global_libraries ?: [:]
  if (!raw && !globalLibraries) {
    error('Specify at least one library under raw_libraries or global_libraries.')
  }

  def samples = channel.fromPath(params.input_csv, checkIfExists: true).splitCsv(header: true).map { row -> tuple(row.id as String, file(row.path, checkIfExists: true)) }
  def ref = CONCAT_REFERENCE(samples)
  def globalLongPaths = globalLibraries.collectMany { name, library -> library.long_reads ?: [] }
  def globalShortLibraries = globalLibraries.collectMany { name, library -> library.short_reads ?: [] }

  def longOut = null
  if (raw.values().any { it.long_reads } || globalLongPaths) {
    def longReads = samples.flatMap { id, refdir ->
      def paths = (raw[id]?.long_reads ?: []) + globalLongPaths
      def resolved = paths.collectMany { path -> files(path) }
      def key = groupKey(id, resolved.size())
      resolved.collect { read -> tuple(id, key, read) }
    }
    def mappedLongReads = longReads.join(ref.ref).map { id, key, read, reference -> tuple(key, read, reference) }
    longOut = MERGE_LONG_BAMS(MINIMAP2_MAP(mappedLongReads).bam.groupTuple())
    FLYE_META_ASSEMBLY(longOut.reads)
    MYLOASM_ASSEMBLY(longOut.reads)
  }

  def shortOut = null
  if (raw.values().any { it.short_reads } || globalShortLibraries) {
    def shortReads = samples.flatMap { id, refdir ->
      def libraries = (raw[id]?.short_reads ?: []) + globalShortLibraries
      def key = groupKey(id, libraries.size())
      libraries
        .withIndex()
        .collect { library, index ->
          tuple(id, key, "${id}_short_${index + 1}", file(library.r1, checkIfExists: true), file(library.r2, checkIfExists: true), library.singletons ? file(library.singletons, checkIfExists: true) : file('/dev/null'))
        }
    }
    def mappedShortReads = shortReads.join(BOWTIE2_BUILD(ref.ref).index).map { id, key, library, r1, r2, singletons, index -> tuple(key, library, r1, r2, singletons, index) }
    shortOut = MERGE_SHORT_BAMS(BOWTIE2_MAP(mappedShortReads).bam.groupTuple())
    MEGAHIT_DEFAULT_ASSEMBLY(shortOut.reads)
  }

  if (longOut && shortOut) {
    def hybridReads = shortOut.reads.join(longOut.reads)
    SPADES_META_HYBRID(hybridReads)
    SPADES_META_HYBRID_KEXT(hybridReads)
    SPADES_META_KEXT(hybridReads)
    SPADES_META_KEXT_TRUSTED(hybridReads.join(ref.ref))
  }
}
