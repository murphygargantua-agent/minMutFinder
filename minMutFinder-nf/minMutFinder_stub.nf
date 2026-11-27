// minMutFinder_stub.nf — stub mimicking the minMutFinder parameter interface.
// Produces the expected output directory structure + merged-result files so that
// run_multi.sh orchestration (parallel per-sample runs + merge) can be validated
// without the real (tool-dependent) bioinformatics steps.
params {
  help       = false
  h          = false
  version    = false
  v          = false
  out_path   = false
  r1         = false
  r2         = false
  ref_seq    = false
  annotate   = false
  syn_muts   = "no"
  manifest   = ""
}

process writeSampleFiles {
  tag "${sample_name}"
  input:
  path out_dir
  path ref_seq
  path qc_csv
  path muts_csv
  val sample_name

  output:
  path "${out_dir}/qc/qc_metrics.csv"
  path "${out_dir}/mutations/${sample_name}_mutations.csv"
  path "${out_dir}/assembly/ref_seq.fasta"

  script:
  """
  echo 'sample,test,score' > ${out_dir}/qc/qc_metrics.csv
  echo '${sample_name},reads_total,1000' >> ${out_dir}/qc/qc_metrics.csv
  echo '${sample_name},reads_filt_qc,800' >> ${out_dir}/qc/qc_metrics.csv

  echo 'SampleID,Gene,Mutation_type,Aa_change,Type_of_aa_change,Nt_mutation,Mutation_frequency' > ${out_dir}/mutations/${sample_name}_mutations.csv
  echo '${sample_name},UL1,missense,A,gain,A,1.0' >> ${out_dir}/mutations/${sample_name}_mutations.csv
  echo '${sample_name},UL2,synonymous,S,loss,S,0.5' >> ${out_dir}/mutations/${sample_name}_mutations.csv

  cp ${ref_seq} ${out_dir}/assembly/ref_seq.fasta
  """
}

def versionMessage() {
  log.info "minMutFinder_stub v1.0.0 (stub)"
}

def helpMessage() {
  log.info "Stub pipeline for run_multi.sh validation."
}

workflow {
  if (params.version || params.v) { versionMessage(); exit 0 }
  if (params.help || params.h || params.isEmpty()) { helpMessage(); exit 0 }

  def out_path = params.out_path
  def out_dir = file(out_path); out_dir.mkdirs()

  // Manifest mode: one row per sample -> create the per-sample output structure.
  if (params.manifest) {
    Channel.fromPath(params.manifest)
        .splitCsv(header: true, quote: '"')
        .map { row -> [row.sample_name, new File(row.r1).getCanonicalPath(), new File(row.r2).getCanonicalPath()] }
        .set { ch_samples }

    // Pre-create the per-sample output dirs and submit the file-writing process.
    ch_samples.map { name, r1, r2 ->
      def sample_dir = file("$out_path/$name"); sample_dir.mkdirs()
      ["fastq","assembly","variant_calling","mutations","plots","qc"].each { subdir ->
        file("$sample_dir/$subdir").mkdirs()
      }
      sample_dir
    }.set { ch_out_dirs }

    def ch_ref_seqs = ch_samples.map { name, r1, r2 -> params.ref_seq }
    def ch_qc_csvs = ch_samples.map { name, r1, r2 -> "$out_path/$name/qc/qc_metrics.csv" }
    def ch_muts_csvs = ch_samples.map { name, r1, r2 -> "$out_path/$name/mutations/${name}_mutations.csv" }

    writeSampleFiles(ch_out_dirs, ch_ref_seqs, ch_qc_csvs, ch_muts_csvs, ch_samples.map { name, r1, r2 -> name })

    println "stub: manifest mode with ${ch_samples.count()} sample(s)"
  } else {
    // Legacy single-sample mode: out_dir is the user-supplied out_path so the
    // per-sample naming behaviour (out_path.baseName) is preserved.
    def r1 = params.r1 && !params.r1.equals(false) ? new File(params.r1).getCanonicalPath() : null
    def r2 = params.r2 && !params.r2.equals(false) ? new File(params.r2).getCanonicalPath() : null
    Channel.of([null, r1, r2]).set { ch_samples }

    def SAMPLE = new File(out_path).getName()
    ["fastq","assembly","variant_calling","mutations","plots","qc"].each { subdir ->
      file("$out_path/$subdir").mkdirs()
    }
    // single sample: out_dir is the user-supplied out_path
    writeSampleFiles(
      file(out_path), params.ref_seq, "$out_path/qc/qc_metrics.csv",
      "$out_path/mutations/${SAMPLE}_mutations.csv", SAMPLE)
  }

  if (params.syn_muts == "yes") {
    println "stub: syn_muts=yes (would append UL3 mutation)"
  }
  if (params.annotate) {
    def ann = file(params.annotate)
    if (ann.exists()) {
      println "stub: annotate=$ann"
    }
  }
}
