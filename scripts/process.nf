params.baseDir = "/home/nicolai/tnSeqQC/scripts/"
params.repeatSeq = "TACGAAGACCGGGGACTTATCATCCAACCTGT"
params.barcodeReplacementSeq = "TGTATGAGACGTCAGAATTGGTTAATTCGTCTCT"
params.spacerSeq = "CAGAATTGGGAGTCTACGAAGACCGGGGACTTATCATCCAACCTGT"
params.identifier = "ACTG"
params.genomePath = "/home/nicolai/for_carlos_3/seq/GCF_000154205.1_ASM15420v1_genomic_ATCC_8492.fna"
params.reads = "/home/nicolai/for_carlos_3/seq/seq/*/*_1_sequence.fastq"
// For testing
//params.reads = "/home/nicolai/for_carlos_3/seq/seq2/*/*_1_sequence.fastq"
params.readTrim5 = 110
params.pathToIndex = "/home/nicolai/for_carlos_3/seq/atcc_8492_concatenated"
params.publishDir = "/home/nicolai/testNF"

Channel.fromPath(params.reads).map{ str ->
  def sampleID = str.name.replaceAll("_1_sequence.fastq", "")
  return tuple(sampleID, str)
}.into{read_pairs_ch; read_pairs_ch2; read_pairs_ch3; read_pairs_ch4; read_pairs_ch5; read_pairs_ch6}
genomePath = Channel.value(params.genomePath)
trim5 = Channel.value(params.readTrim5)
rs = Channel.value(params.repeatSeq)
brs = Channel.value(params.barcodeReplacementSeq)
i = Channel.value(params.pathToIndex)
spacerSeq = Channel.value(params.spacerSeq)
identifier = Channel.value(params.identifier)
baseDir = Channel.value(params.baseDir)

process map_reads {

   input:
   tuple val(s), path (sampleID) from read_pairs_ch
   val index from i
   val trim from trim5

   output:
   tuple val(s), file("sampleIDbam") into bamFile

   """
   bowtie2 -x ${index} --trim5 ${trim} --very-sensitive-local -U ${sampleID} -S sampleIDbam
   """
}

process sort_bam {

   input:
   tuple val(s), file(b) from bamFile

   output:
   tuple val(s), file(bamFileSorted) into bamFileSorted
   tuple val(s), file(bamFileSorted) into bamFileSorted2
   tuple val(s), file(bamFileSorted) into bamFileSorted3

   """
   samtools sort ${b} > bamFileSorted
   """

}

bamFileSorted = bamFileSorted.join(read_pairs_ch2)
bamFileSorted3 = bamFileSorted3.join(read_pairs_ch6)

process get_depth_at_pos {


  input:
  tuple val(s), file(bamFileSorted), file(sampleID) from bamFileSorted3

  output:
  file "*_depthAtPos" into depthAtPos

  """
  samtools depth -d 0 ${bamFileSorted} | sed "s/\t//" | sed "s/^/${s}\t/" > ${sampleID}_depthAtPos
  """
}

process get_read_info {

  input:
  tuple val(s), file(bamFileSorted), file(sampleID) from bamFileSorted

  publishDir params.publishDir

  output:
  tuple val(s), file("*_readInfo") into ri

  """
  samtools view ${bamFileSorted} | cut -f1,2,4 | sed "s/^/${s}\t/"  > ${sampleID}_readInfo
  """
}

bamFileSorted2 = bamFileSorted2.join(read_pairs_ch3)

process get_depth_info {

  input:
  tuple val(s), file(bamFileSorted), path(sampleID) from bamFileSorted2

  publishDir params.publishDir

  output:
  tuple val(s), file("*_depth") into de

  """
  samtools depth ${bamFileSorted} | sed "s/^/${s}\t/" > ${s}_depth
  """

}

process reads_to_fasta {

  input:
  tuple val(s), path(sampleID) from read_pairs_ch4

  publishDir '/home/nicolai/testNF'

  output:
  tuple val(s), file("*_fastaRead.fa") into fastaReads

  """
  paste <(cat $sampleID | paste - - - - | cut -f1 | sed "s/^@/>/") <(cat $sampleID | paste - - - - | cut -f2) | tr "\t" "\n" > ${sampleID}_fastaRead.fa
  """

}

/*
ri = ri.map{ str ->
  def sampleID = str.name.replaceAll("_1_sequence.fastq_readInfo", "")
  return tuple(sampleID, str)
}


fastaReads = fastaReads.map{ str ->
  def sampleID = str.name.replaceAll("_1_sequence.fastq_fastaRead.fa", "")
  return tuple(sampleID, str)
}

read_pairs_ch5 = read_pairs_ch5.map{ str ->
  def sampleID = str.name.replaceAll("_1_sequence.fastq", "")
  return tuple(sampleID, str)
}
*/

all = read_pairs_ch5.join(fastaReads).join(ri)

/*
Channel
.fromPath("${params.input_dir}/{fna,faa,hmm}/*.{fna,faa,hmm}")
.map { file ->
def sampleId = file.name.replaceAll(/\.[hf][mna][ma]$/, "")
return tuple(sampleId, file)
}
.groupTuple()
.set { markers_ch }

Channel
.fromPath(params.input_dir + "/" + params.file_pattern)
.map { file ->
def sample = file.name.replaceAll(suffix_pattern, "")
sample = sample.replaceAll(/\.$/, "")
return tuple(sample, file)
}
.groupTuple()
.set { samples_ch }
input:
set sample, file(bamfile) from samples_ch
*/

process parse_read_stats {
  input:
  val rs from rs
  val brs from brs
  val spacerSeq from spacerSeq
  val identifier from identifier
  //path sampleID from read_pairs_ch5
  // Get rid of this channel; get fasta from fastq.
  //path fastaFile from fastaReads
  //path readInfo from ri
  tuple val(s), path(sampleID), path(fastaFile), path(readInfo) from all
  val baseDir from baseDir

  output:
  file "*_readInfoParsed" into readInfoParsed
  file "*_readInfoParsed2" into readInfoParsed2

  """
  python ${baseDir}analyse_barcode_distribution_2021_03_13.py ${fastaFile} ${readInfo} ${sampleID}_readInfoParsed ${rs} ${brs}
  python ${baseDir}analyse_barcode_distribution.py ${spacerSeq} ${identifier} ${fastaFile} ${sampleID}_readInfoParsed2 ${readInfo}
  """
}

/*
Channel
.fromPath("${params.input_dir}/{fna,faa,hmm}/*.{fna,faa,hmm}")
.map { file ->
def sampleId = file.name.replaceAll(/\.[hf][mna][ma]$/, "")
return tuple(sampleId, file)
}
.groupTuple()
.set { markers_ch }
Channel
.fromPath(params.input_dir + "/" + params.file_pattern)
.map { file ->
def sample = file.name.replaceAll(suffix_pattern, "")
sample = sample.replaceAll(/\.$/, "")
return tuple(sample, file)
}
.groupTuple()
.set { samples_ch }
input:
set sample, file(bamfile) from samples_ch
*/

readInfoParsedConcat = readInfoParsed.collectFile(name: 'o')
readInfoParsedConcat2 = readInfoParsed2.collectFile(name: 'o2')
depthAtPosConcat = depthAtPos.collectFile(name: "depthAtPos")

readInfoParsedConcat = readInfoParsedConcat.last()
readInfoParsedConcat2 = readInfoParsedConcat2.last()
depthAtPosConcat = depthAtPosConcat.last()

process get_barcode_distrib_plots {

 input:
 file r from readInfoParsedConcat
 file r2 from readInfoParsedConcat2
 file depthAtPos from depthAtPosConcat
 val baseDir from baseDir

 publishDir params.publishDir

 output:
 file "barcode_distrib_plots.png" into devNull
 file "testPlot*" into devNull2

 """
 Rscript ${baseDir}analyse_read_situation_final.r ${r} ${r2} depthAtPos barcode_distrib_plots.png testPlot
 """

}
