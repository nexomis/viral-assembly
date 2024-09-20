# TO DO list

List of things to implement in the pipeline
## Reads corretion
TODO: check that read correction in spades is enabled for all interest modes (meta?), otherwise add an explicit step for read correction

## ICARUS/QUAST
- possibilities to add annotation
- do report at batch level: useful for Contig alignment viewer of ICARUS
  
## Indel Realign
- Manage BAM realignments with ABRA, ABRA2 or GATK3.
- Add diagnostic plot around indels.
~~If not InDels-realignement: cut reads extreminity on bam, but very archaic approches !~~

## Subworkflow map_with_indels
Build a subworkflow to map | sort | indel realign | resort (| remove duplicate ?)

## Consensus from Alignment
- Obtain a consensus from the alignment using a local re-assembly software like Platypus
- remapping vs the consensus 
- pairwise align consensus vs reference

make_consensus_from bam with make_consensus | realign_with_indel | pairwise_align_consensus_v_reference

## Selection of contigs with reference with reordering but without scaffolding.

Use Abacas or alternative to reorder contigs/scaffolds based on a reference genome without scaffolding.
*Note: If not used abacas contigs: realign all candidate contigs in the same orientation as the reference*
*TODO: understand that Abacas aggregates contigs according to coverage? Add N? If not manageable, compare with scaffolding tools*

## Replace krake2n with Centrifuge, at least for the classification of reads unmapped to the assembly.

## Evaluate ElasticBLAST: in-depth characterization of unassigned reads (conta, large viral derivation, ...) by blasting of a subset of unmapped reads without managing db nt?

## Selection of contigs based on proteome. 

- Build an incidence matrix between contigs and input proteins
- Input based on tblastn, diamond or Miniprot.
 - use protein coverage as incidence values
 - use sign to indicate strand 
- Exploit this incidence matrix to select and re-orient contigs/scaffolds

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

## Integrate with hannot annotation pipeline

Annotate the final assemblies based on hannot pipeline and an input proteom
