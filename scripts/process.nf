/*
 * pipeline input parameters
 */
params.repeatSeq = "TACGAAGACCGGGGACTTATCATCCAACCTGT"
params.barcodeReplacementSeq = "TGTATGAGACGTCAGAATTGGTTAATTCGTCTCT"
params.genomePath = "/home/nicolai/for_carlos_3/seq/GCF_000154205.1_ASM15420v1_genomic_ATCC_8492.fna"
params.speciesAndStrainName = "BacteroidesUnformisXXX"
params.reads = "/home/nicolai/for_carlos_3/seq/seq/*/*_1_sequence.fastq"
params.readTrim5 = 110
params.pathToIndex = "/home/nicolai/for_carlos_3/seq/atcc_8492_concatenated"
params.publishDir = "/home/nicolai"

Channel.fromPath(params.reads).map{ str ->
  def sampleID = str.name.replaceAll("_1_sequence.fastq", "")
  return tuple(sampleID, str)
}.into{read_pairs_ch; read_pairs_ch2; read_pairs_ch3; read_pairs_ch4; read_pairs_ch5}
genomePath = Channel.value(params.genomePath)
genomeName = Channel.value(params.speciesAndStrainName)
trim5 = Channel.value(params.readTrim5)
rs = Channel.value(params.repeatSeq)
brs = Channel.value(params.barcodeReplacementSeq)
i = Channel.value(params.pathToIndex)

process concatenate_genome {
    input:
    file genomePath from genomePath
    file genomeName from genomeName

    output:
    file "genome_concatenated" into genomePathConcatenated

    """
    cat $genomePath | grep -v "^>" > tmp; echo ">$genomeName" >> genome_concatenated; cat tmp >> genome_concatenated
    """
}

/*
process index_genome {
   input:
   file genomeFile from genomePathConcatenated

   output:
   path "testIndexName*" into genomeIndexName

   """
   bowtie2-build $genomeFile testIndexName && touch testIndexName
   """
}
*/

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

   """
   samtools sort ${b} > bamFileSorted
   """

}

bamFileSorted = bamFileSorted.join(read_pairs_ch2)

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

  publishDir '/home/nicolai/testNF'

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

// TODO: Figure out how to modify channels. Specifically, understand how to turn a list of strings
// into a list of two element tuples, where the first element is an ID that can be used to join channels.

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
  //path sampleID from read_pairs_ch5
  // Get rid of this channel; get fasta from fastq.
  //path fastaFile from fastaReads
  //path readInfo from ri
  tuple val(s), path(sampleID), path(fastaFile), path(readInfo) from all

  publishDir '/home/nicolai/testNF'

  output:
  file "*_readInfoParsed" into readInfoParsed

  """
  python /home/nicolai/for_carlos_3/seq/analyse_barcode_distribution_2021_03_13.py ${fastaFile} ${readInfo} ${sampleID}_readInfoParsed ${rs} ${brs}
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

//process collect_parsed_read_stats {
// input:
//  file r from readInfoParsed
//  publishDir '/home/nicolai/testNF'
//  output:
//  file o into readInfoParsedConcat
//  """
//  cat ${r} >> o
//  """
//}

readInfoParsedConcat = readInfoParsed.collectFile(name: 'o')

// We've combined everything into o above, so make sure channel has only one 'member'.
readInfoParsedConcat = readInfoParsedConcat.last()

process get_barcode_distrib_plots {

 input:
 file r from readInfoParsedConcat

 publishDir '/home/nicolai/testNF'

 output:
 file "barcode_distrib_plots.png" into devNull

 """
 Rscript /home/nicolai/for_carlos_3/seq/analyse_read_situation_final.r ${r} barcode_distrib_plots.png
 """

}
