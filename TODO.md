# TO DO list

List of things to implement in the pipeline

## Indel Realign
- Manage BAM realignments with ABRA, ABRA2 or GATK3.
- Add diagnostic plot around indels.

## Subworkflow map_with_indels
Build a subworkflow to align | sort | indel realign | resort

## Consensus from Alignment
- Obtain a consensus from the alignment using a local re-assembly software like Platypus
- realign vs the consensus
- pairwise align consensus vs reference

make_consensus_from bam with make_consensus | realign_with_indel | pairwise_align_consensus_v_reference

## Selection of contigs with reference with reordering but without scaffolding.

Use Abacas or alternative to reorder contigs/scaffolds based on a reference genome without scaffolding.

## Replace krake2n with Centrifuge, at least for the classification of reads unmapped to the assembly.

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

## Integrate with hannot annotation pipeline

Annotate the final assemblies based on hannot pipeline and an input proteom
