# Simple Coverage

A very simple, easy to digest way of determining coverage of primers/probes/random sequences.

## Input

All input is in FASTA format and can only contain the four standard nucleotides (A, C, G, T).  

A list of primer/probe/whatever sequences  
Target sequence(s)

## Output

The percentage the target sequence(s) which is covered by the list of input sequences

### Possible future outputs

- mapping of individual input seqs to their target seq alignment spot
- mismatch metrics (counts, positions, composition, etc.)
- coverage of target sequence by individual input seqs

## Usage

```bash
python scrun.py -i input.fasta -t target.fasta
```
