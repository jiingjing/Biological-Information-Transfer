# Biological-Information-Transfer

## Overview

This work demonstrates the processes of information transfer that makes up the central dogma of molecular biology:

- DNA replicates
- DNA is transcribed into RNA
- RNA is translated into protein

Additional utilities are included for reading-frame analysis and sequence-level calculations.

## Contents

- [Repository structure](#repository-structure)
- [Example Gene Expression](#example-gene-expression)
- [Usage](#usage)
- [Biological Background Information](#biological-background-information)
  - [DNA Replication](#dna-replication)
  - [Transcription](#transcription)
  - [Translation](#translation)
  - [Reading frames](#reading-frames)
  - [Nucleotide rules](#nucleotide-rules)
- [Notes](#notes)
- [References](#references)

## Repository structure

- `information_transfer.py` - core processes for DNA replication, RNA transcription, and protein translation
- `reading_frame.py` - helper class for three-frame and six-frame translation analysis
- `sequence_utils.py` - sequence utilities such as GC content calculation

## Example Gene Expression

Given a DNA template strand:

- Template: `3' TACTTACACACG 5'`
- Coding strand: `5' ATGAATGTGTGC 3'`

The process produces:

- DNA replication: `5' ATGAATGTGTGC 3'`
- RNA transcription: `5' AUGAAUGUGUGC 3'`
- Protein translation: `N-Met-Asn-Val-Cys-C`

## Usage

```
    git clone https://github.com/jiingjing/Biological-Information-Transfer.git
```

In the `Biological-Information-Transfer` directory, create and activate the virtual env:

```
    python3 -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt
```

Example usage in Python:

```python
import pandas as pd
from information_transfer import BioInfoTrans
from reading_frame import ReadingFrameAnalysis
from sequence_utils import gc_content

# codon table DataFrame
# with columns incl. Codon, Short Name, Letter
codon_df = pd.read_csv(...)

template_strand = "3' TACTTACACACG 5'"
bio = BioInfoTrans(template_strand, codon_df)

dna_copy, rna, protein = bio.info_transfer()

analysis = ReadingFrameAnalysis(bio, rna)
analysis.six_frames_translation()

print(gc_content(analysis.rna_bases))
```

```
output:

DNA Template Strand: 3'TACTTACACACG5'
DNA Replication: 5'ATGAATGTGTGC3'
RNA Transcription: 5'AUGAAUGUGUGC3'
Protein Translation: N-Met-Asn-Val-Cys-C
+3    E  C  V
+2   *  M  C
+1  M  N  V  C
→   AUGAAUGUGUGC
←   UACUUACACACG
-1    H  I  H  A
-2      F  T  H
-3     S  H  T
41.66666666666667
```

Recommended codon table

- Preferred: download the RNA codon table CSV from https://codontable.org/codon-table
- Alternative: use the table at `codon_tables/codon-table-standard-rna.csv`
- Format: CSV (UTF-8) with a header row. Expected columns: `Codon`, `Short Name`, `Letter`

## Biological Background Information

### DNA Replication

1. The two strands of the original dsDNA separate
2. RNA primers attach to the separated DNA strands
3. The DNA is replicated off of the 3’ end of the RNA primers
4. Gives two identical DNA molecules

---

- Template strand = strand being 'filled in'
- New strand is synthesised in the 5' → 3' direction
- New strand is anti-parallel (has opposite orientation) and has complementary base pairing to the template strand

  ![The steps of DNA Replication](https://openlearninglibrary.mit.edu/assets/courseware/v1/02dbcc2f4ff227938abdf49e01f463bd/asset-v1:OCW+Pre-7.01+1T2020+type@asset+block/m0_3.png)

### Transcription

1. A promoter on the DNA tells RNA polymerase where to start transcription
2. RNA polymerase binds to one DNA strand, called the template strand
3. RNA polymerase reads the template strand 3' → 5' and builds RNA 5' → 3'
4. RNA nucleotides are joined together to make an RNA strand
5. When RNA polymerase reaches a stop signal, transcription ends
6. The RNA transcript detaches

---

- RNA uses uracil (U) instead of thymine (T)

  ![Transcription](https://openlearninglibrary.mit.edu/assets/courseware/v1/70f34ce5ec809f12c0a99d03a0251460/asset-v1:OCW+Pre-7.01+1T2020+type@asset+block/n0_newBW.png)

### Translation

1. mRNA is read in groups of three bases called codons
2. Translation starts at the start codon `AUG`, which codes for methionine
3. Each codon specifies a particular amino acid using the genetic code
4. Amino acids are joined together to form a growing protein chain
5. Translation ends when a stop codon such as `UAA`, `UAG`, or `UGA` is reached

---

- A codon is a set of 3 consecutive nucleotides that together code for an amino acid
- Protein chains grow from the amino (N) terminus to the carboxyl (C) terminus

  ![A Codon Chart](https://openlearninglibrary.mit.edu/assets/courseware/v1/522a130818e51e7d580bc9834ad314b6/asset-v1:OCW+Pre-7.01+1T2020+type@asset+block/p0.png)

  _A Codon Chart_

### Reading frames

- DNA and RNA can be read in codons of three bases each
- There are three forward reading frames for a single strand
- There are three reverse reading frames for the complementary strand
- Frame analysis is useful when the exact start position is unknown

  ![Six different reading frames give six different ways for the same fragment of DNA to be transcribed and translated, three from each strand. The top three amino acid strings are read from left to right, whereas the bottom three strings are read from right to left. Stop codons are represented by XXX.](https://rosalind.info/media/problems/ba4b/reading_frames.png)

### Nucleotide rules

- DNA bases: `A`, `T`, `G`, `C`
- RNA bases: `A`, `U`, `G`, `C`
- Base pairing:
  - `A` pairs with `T` in DNA, `A` pairs with `U` in RNA
  - `G` pairs with `C`

## Notes

- `sequence_utils.py` isolates general sequence calculations from the central-dogma logic
- `reading_frame.py` keeps reading-frame analysis separate from the core `BioInfoTrans` class
- Transcription is applied to the entire input DNA template sequence

## References

For Biological Background Information:

- MIT OpenCourseWare - Fundamentals of Biology: https://ocw.mit.edu/courses/7-01sc-fundamentals-of-biology-fall-2011/pages/molecular-biology/ (accessed 20 May 2026)
- MIT Open Learning Library - Pre-7.01 course: https://openlearninglibrary.mit.edu/courses/course-v1:OCW+Pre-7.01+1T2020/course/#block-v1:OCW+Pre-7.01+1T2020+type@chapter+block@b6ee7fedeabc46de8167c1b954a0b8fb (accessed 20 May 2026)
