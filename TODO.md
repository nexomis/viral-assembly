# TO DO list

List of things to implement in the pipeline.

## Consensus from Alignment
- Obtain a consensus from the alignment using a local re-assembly software like Platypus
- remapping vs the consensus 
- pairwise align consensus vs reference

make_consensus_from bam with make_consensus | realign_with_indel | pairwise_align_consensus_v_reference

## Evaluate ElasticBLAST: in-depth characterization of unassigned reads (conta, large viral derivation, ...) by blasting of a subset of unmapped reads without managing db nt?

## add other QCs
- BUSCO with viral databases
- Other quality metrics REAPR, FRCBAM (check redundancy with QUAST)
- Add diagnostic plot or metrics such as dot-plots

## Iterative reassembly

(see ICORN2)

## Genome Finishing with Pilon or alternatives pypolca, ICORN2 (iterative !).
  - Gap closer
  - Genome missassembly breaker
  - identify missasemblie: palindromic sequence search (something less manual than dotplot) - or bandage in relation with icarus plot ??

