nextflow.enable.dsl=2

process CONCAT_REFERENCE {
  tag { id }
  label 'mapping'
  input:
  tuple val(id), path(refdir)
  output:
  tuple val(id), path('reference.fasta'), emit: ref
  script:
  """find -L ${refdir} -maxdepth 1 -type f \\( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \\) -print0 | xargs -r -0 cat | seqkit rename -n > reference.fasta""" }
process MINIMAP2_MAP {
  tag { id }
  label 'mapping'
  input:
  tuple val(id), path(reads), path(ref)
  output:
  tuple val(id), path('*.bam'), emit: bam
  script:
  """minimap2 -t ${task.cpus} -ax lr:hq ${ref} ${reads} | samtools view -@ ${task.cpus} -bF 4 - | samtools sort -@ ${task.cpus} -o ${id}.bam -""" }
process MERGE_LONG_BAMS {
  tag { id }
  label 'mapping'
  input:
  tuple val(id), path(bams)
  output:
  tuple val(id), path('long_recovered.fastq.gz'), emit: reads
  script:
  """samtools merge -@ ${task.cpus} -f merged.bam ${bams.join(' ')}; samtools fastq -@ ${task.cpus} merged.bam | pigz -p ${task.cpus} > long_recovered.fastq.gz""" }
process FLYE_META_ASSEMBLY {
  tag { id }
  label 'lr_assm'
  input:
  tuple val(id), path(reads)
  output:
  path('flye_out')
  script:
  """flye -t ${task.cpus} --meta --nano-corr ${reads} -o flye_out; find flye_out -type f -name '*dump' -delete""" }
process MYLOASM_ASSEMBLY {
  tag { id }
  label 'lr_assm'
  input:
  tuple val(id), path(reads)
  output:
  path('mylo_out')
  script:
  """myloasm ${reads} -t ${task.cpus} -o mylo_out; rm -rf mylo_out/binary_temp; mylotools annotate-gfa --gfa mylo_out/final_contig_graph.gfa --fasta mylo_out/assembly_primary.fa --output mylo_out/final_contig_graph.annot.gfa""" }
process BOWTIE2_BUILD {
  tag { id }
  label 'mapping'
  input:
  tuple val(id), path(ref)
  output:
  tuple val(id), path('reference.fasta.*'), emit:index
  script:
  """bowtie2-build ${ref} reference.fasta""" }
process BOWTIE2_MAP {
  tag { "${id}:${lib}" }
  label 'mapping'
  input:
  tuple val(id), val(lib), path(r1), path(r2), path(single), path(index)
  output:
  tuple val(id), path('*.bam'), emit:bam
  script:
  """bowtie2 -p ${task.cpus} -x reference.fasta -1 ${r1} -2 ${r2} -U ${single} | samtools view -@ ${task.cpus} -bF 12 - | samtools sort -@ ${task.cpus} -o ${lib}.bam -""" }
process MERGE_SHORT_BAMS {
  tag { id }
  label 'mapping'
  input:
  tuple val(id), path(bams)
  output:
  tuple val(id), path('sr_R1.fastq.gz'), path('sr_R2.fastq.gz'), path('sr_singles.fastq.gz'), emit:reads
  script:
  """samtools merge -@ ${task.cpus} -f merged.bam ${bams.join(' ')}; samtools sort -@ ${task.cpus} -n merged.bam -o n.bam; samtools fastq -@ ${task.cpus} n.bam -1 sr_R1.fastq.gz -2 sr_R2.fastq.gz -s sr_singles.fastq.gz -0 /dev/null -""" }
process MEGAHIT_DEFAULT_ASSEMBLY {
  tag { id }
  label 'sr_assm'
  input:
  tuple val(id), path(r1), path(r2), path(single)
  output:
  path('megahit-default')
  script:
  """megahit -t ${task.cpus} -1 ${r1} -2 ${r2} -r ${single} -o megahit-default; kmer=\$(find megahit-default/intermediate_contigs -name 'k*.contigs.fa' -printf '%f\\n' | sed -E 's/^k([0-9]+)\\.contigs\\.fa\$/\\1/' | sort -n | tail -1); megahit_toolkit contig2fastg "\$kmer" "megahit-default/intermediate_contigs/k\${kmer}.contigs.fa" > "megahit-default/k\${kmer}.contigs.fastg""" }
process SPADES_META_HYBRID {
  tag { id }
  label 'sr_assm'
  input:
  tuple val(id), path(r1), path(r2), path(s), path(lr)
  output:
  path('spades-meta-hybrid')
  script:
  """spades.py -t ${task.cpus} -m 100 --meta -1 ${r1} -2 ${r2} -s ${s} --nanopore ${lr} -o spades-meta-hybrid""" }
process SPADES_META_HYBRID_KEXT {
  tag { id }
  label 'sr_assm'
  input:
  tuple val(id), path(r1), path(r2), path(s), path(lr)
  output:
  path('spades-meta-hybrid-kext')
  script:
  """spades.py -t ${task.cpus} -m 100 --meta -1 ${r1} -2 ${r2} -s ${s} -k 21,33,55,77,101,127 --nanopore ${lr} -o spades-meta-hybrid-kext""" }
process SPADES_META_KEXT {
  tag { id }
  label 'sr_assm'
  input:
  tuple val(id), path(r1), path(r2), path(s), path(lr)
  output:
  path('spades-meta-kext')
  script:
  """spades.py -t ${task.cpus} -m 100 --meta -1 ${r1} -2 ${r2} -s ${s} -k 21,33,55,77,101,127 --nanopore ${lr} -o spades-meta-kext""" }
process SPADES_META_HYBRID_KEXT_TRUSTED {
  tag { id }
  label 'sr_assm'
  input:
  tuple val(id), path(r1), path(r2), path(s), path(lr), path(trusted)
  output:
  path('spades-meta-hybrid-kext-trusted')
  script:
  """spades.py -t ${task.cpus} -m 100 -1 ${r1} -2 ${r2} -s ${s} --nanopore ${lr} --trusted-contigs ${trusted} -k 21,33,55,77,101,127 -o spades-meta-hybrid-kext-trusted""" }

workflow {
  def samples=channel.fromPath(params.input_csv,checkIfExists:true).splitCsv(header:true).map { row -> tuple(row.id as String,file(row.path,checkIfExists:true)) }
  def ref=CONCAT_REFERENCE(samples)
  def raw=params.raw_libraries ?: [:]
  def longSpecs=raw.collectMany { id,x -> (x.long_reads ?: []).collect { p -> tuple(id as String,p as String) } }
  def longReads=channel.fromList(longSpecs).flatMap { id,p -> files(p).collect { f -> tuple(id,f) } }
  def longOut=MERGE_LONG_BAMS(MINIMAP2_MAP(longReads.join(ref.ref)).bam.groupTuple())
  FLYE_META_ASSEMBLY(longOut.reads); MYLOASM_ASSEMBLY(longOut.reads)
  def shortSpecs=raw.collectMany { id,x -> (x.short_reads ?: []).withIndex().collect { z,n -> tuple(id as String,"${id}_short_${n+1}",z) } }
  def shortReads=channel.fromList(shortSpecs).map { id,lib,x -> tuple(id,lib,file(x.r1,checkIfExists:true),file(x.r2,checkIfExists:true),x.singletons ? file(x.singletons,checkIfExists:true) : file('/dev/null')) }
  def shortOut=MERGE_SHORT_BAMS(BOWTIE2_MAP(shortReads.join(BOWTIE2_BUILD(ref.ref).index)).bam.groupTuple())
  MEGAHIT_DEFAULT_ASSEMBLY(shortOut.reads)
  def h=shortOut.reads.join(longOut.reads)
  SPADES_META_HYBRID(h); SPADES_META_HYBRID_KEXT(h); SPADES_META_KEXT(h); SPADES_META_HYBRID_KEXT_TRUSTED(h.join(ref.ref))
}
