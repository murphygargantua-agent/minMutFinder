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
  ["fastq","assembly","variant_calling","mutations","plots","qc"].each { subdir ->
    file("$out_path/$subdir").mkdirs()
  }

  // copy reference + index
  def ref = file(params.ref_seq)
  "cp ${ref} ${out_path}/assembly/ref_seq.fasta".execute()
  "samtools faidx ${out_path}/assembly/ref_seq.fasta".execute()

  def SAMPLE = new File(out_path).getName()
  def muts_csv = "$out_path/mutations/${SAMPLE}_mutations.csv"
  def qc_csv = "$out_path/qc/qc_metrics.csv"

  def bash = ["/bin/bash", "-c"]
  def qcmd = "echo 'sample,test,score' > $qc_csv && echo '${SAMPLE},qc_test,0.99' >> $qc_csv"
  (bash + [qcmd]).execute()
  def mcmd = "echo 'SampleID,Gene,Mutation_type,Aa_change,Type_of_aa_change,Nt_mutation,Mutation_frequency' > $muts_csv && echo '${SAMPLE},UL1,missense,A,gain,A,1.0' >> $muts_csv && echo '${SAMPLE},UL2,synonymous,S,loss,S,0.5' >> $muts_csv"
  (bash + [mcmd]).execute()

  if (params.syn_muts == "yes") {
    (bash + ["echo '${SAMPLE},UL3,missense,G,gain,G,0.8' >> $muts_csv"]).execute()
  }

  if (params.annotate) {
    def ann = file(params.annotate)
    if (ann.exists()) {
      (bash + ["echo annotated >> $muts_csv"]).execute()
    }
  }
}
