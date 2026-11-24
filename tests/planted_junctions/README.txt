The junction that Modules 2 and 3 report, tested on planted junctions
=====================================================================

A test against known answers of what Tandem reports for a junction found in reads: its
position, its microhomology, and whether a direct repeat flanks it. Duplications with
known junctions are planted in a reference genome, reads are simulated from a population
that carries each one, and the scan and Module 3 are run on the reads. Nothing in the
reference folders is changed.

Run, from the folder that holds genome/pop_ncbi_ref (one <sp>/<sp>.fna per genome):

  THREADS=8 bash tests/planted_junctions/run_junction_check.sh

Other genomes are added with EXTRA_REFS="a.fasta b.fasta". Results go to
junction_check_v<version>/ (OUT=... chooses another folder): summary.txt, and for each
genome mutants/truth.tsv, reads_planted.tsv and random_junctions.tsv; the .tar.gz beside
the folder holds the same without genomes or reads. SEED=2 plants a different set;
RANDOM_N=0 skips step 3. Needs tandem on the PATH.

What the run does for every genome
  1. make_mutants.py plants ten duplications of 1.5-48 kb, one per mutant genome, in the
     longest sequence; both ends of every unit lie in sequence found once in the genome.
       m01      no microhomology (two bases are also inserted at the junction in the
                mutant genome; the simulated reads do not carry them)
       m02      no microhomology
       m03-m05  1, 2 and 4 bp of microhomology already present in the genome
       m06-m08  7, 12 and 18 bp of microhomology, planted before the unit
       m09-m10  a junction at the middle of a pair of direct repeats of 60 and 80 bp at
                95% identity (the two copies differ at one base in every 20, so they
                read the same for 19 bp at the junction)
     In every design the base beyond the microhomology differs, so the planted length is
     exact.
  2. reads_planted.py cuts the genome 120 kb either side of the unit, simulates reads from
     a population with the duplication in 30% of cells (tests/simulate_reads.py: 150 bp
     pairs from 350 bp fragments, 100x, 0.1% errors), and runs the scan (no coordinates)
     and Module 3 (coordinates given) with this excerpt as the reference. Every junction
     has its own read seed (read_seed in reads_planted.tsv). With one seed for all, the
     simulator places its fragments at the same positions in every excerpt of the same
     length, and the read counts of different genomes repeat each other.
  3. random_junctions.py draws 100 duplications per size class (1.5-10 kb, 10-50 kb,
     50 kb-1 Mb) whose two ends are unrelated, and calls the junction description at
     their exact coordinates on the whole genome (no reads needed).

Reading summary.txt
  scan: junction found       a junction at exactly the position Module 3 reports
  scan: other calls          junctions reported anywhere else in the excerpt
  Module 3: placed leftmost  the junction moved left by the microhomology, as designed
  microhom. exact            the reported length equals the planted one
  repeat where none planted  a flanking repeat reported at m01-m08
  planted repeat called      a flanking repeat reported at m09-m10
The reported repeat is shorter than the planted one: the junction lies inside the repeat,
and only the part of the repeat on one side of it is compared.
  Planted repeats            at m09-m10 the junction can lie anywhere along the bases
                             that the two copies of the repeat share around the planted
                             crossover; it is expected at the leftmost of these
                             positions, with their number as its microhomology
                             (expected_start and expected_microhomology in the table)
  Dosage                     every junction is planted in 30% of cells, so the dosage
                             should be 0.30: mean, spread and how many of the 95%
                             intervals hold 0.30, for Module 3 and for the scan
