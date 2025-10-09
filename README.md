# gfa2gff

`gfa2gff` is a command-line tool that converts a **GFA (Graphical Fragment
Assembly)** file representing a **compacted de Bruijn graph** into a **GFF3**
annotation file.  
Each segment in the GFA graph becomes one feature in the GFF file.
This makes it possible for downstream tools to view and analyze the graph directly in genomic coordinates.

It is compatible with graphs produced by tools such as **Bifrost**, **GGCAT**, **TwoPaCo**, **BCALM2**, **Cuttlefish2** and any other that generate GFA representations of compacted de Bruijn graphs.

------------------------------------------------------------

## What this means:

Say you have a genome with three chromosomes.  
Two of them (**chr1** and **chr2**) differ by a single base in the middle (`A`  ↔ `T`),
while **chr3** is the same as **chr1** but with an `N` inserted.

`example.fa`
```
>chr1
GATACAT
>chr2
GATTCAT
>chr3
GATNACAT
```

---

Its compacted de Bruijn graph (without considering canonical unitigs) is the following:

```
         ┌───(2) ATACA───┐
(1) GAT──┤               ├──(4) CAT
         └───(3) ATTCA───┘
```

---

The GFF output by `gfa2gff` is shown below.
Here:

* `ID=` is the ID of the node in the GFA.
* `genome=` is the name of the FASTA file it comes from (`example.fa`).
* `substr=` is present when the string is a substring of the node.  
  In this case, the string `ACA` is a substring of node `2`, from position `3` to `5`.
  (These positions are 1-based index and both inclusive, like in the GFF).

In **chr3** the 4th base (`N`) is not represented, since it is not an
A/C/G/T character and most compacted de Bruijn graph tools ignore such regions.

```
##gff-version 3.1.26
##sequence-region chr1 1 7
##sequence-region chr2 1 7
##sequence-region chr3 1 8
chr1	gfa2gff	SO:0000856	1	3	.	+	.	ID=1;genome=example
chr1	gfa2gff	SO:0000856	2	6	.	+	.	ID=2;genome=example
chr1	gfa2gff	SO:0000856	5	7	.	+	.	ID=4;genome=example
chr2	gfa2gff	SO:0000856	1	3	.	+	.	ID=1;genome=example
chr2	gfa2gff	SO:0000856	2	6	.	+	.	ID=3;genome=example
chr2	gfa2gff	SO:0000856	5	7	.	+	.	ID=4;genome=example
chr3	gfa2gff	SO:0000856	1	3	.	+	.	ID=1;genome=example
chr3	gfa2gff	SO:0000856	5	7	.	+	.	ID=2;genome=example;substr=(3,5)
chr3	gfa2gff	SO:0000856	6	8	.	+	.	ID=4;genome=example
```

---

## Install

### Dependencies

- **C++17**   
- **OpenMP**   
- **zlib** 

Clone the repository and build:

```bash
git clone https://github.com/lucaparmigiani/gfa2gff
cd gfa2gff
make
```

## Run

Basic usage:

```bash
gfa2gff <k> <graph.gfa> <genome.fasta> [more.fasta ...] [options]
```

Options:

* `-t, --threads <N>`: number of threads (default: all available cores)

---

## Example

Run with three genomes:

```bash
gfa2gff 11 example/graph.gfa example/a.fa example/b.fa example/c.fa > graph.gff
```

Output:
```
##gff-version 3.1.26
##sequence-region A 1 77
##sequence-region B 1 137
##sequence-region C 1 206
A	gfa2gff	SO:0000856	1	11	.	+	.	ID=10;genome=a
A	gfa2gff	SO:0000856	2	12	.	+	.	ID=11;genome=a
A	gfa2gff	SO:0000856	3	13	.	+	.	ID=12;genome=a
A	gfa2gff	SO:0000856	4	14	.	-	.	ID=12;genome=a
A	gfa2gff	SO:0000856	5	15	.	+	.	ID=13;genome=a
A	gfa2gff	SO:0000856	6	16	.	-	.	ID=13;genome=a
...
B	gfa2gff	SO:0000856	1	11	.	-	.	ID=15;genome=b
B	gfa2gff	SO:0000856	2	14	.	+	.	ID=3;genome=b
B	gfa2gff	SO:0000856	5	16	.	+	.	ID=2;genome=b
B	gfa2gff	SO:0000856	7	17	.	-	.	ID=14;genome=b
B	gfa2gff	SO:0000856	8	18	.	-	.	ID=15;genome=b
...
C	gfa2gff	SO:0000856	126	137	.	+	.	ID=2;genome=c
C	gfa2gff	SO:0000856	128	138	.	-	.	ID=14;genome=c
C	gfa2gff	SO:0000856	129	139	.	-	.	ID=15;genome=c
C	gfa2gff	SO:0000856	130	142	.	+	.	ID=3;genome=c
C	gfa2gff	SO:0000856	133	144	.	+	.	ID=2;genome=c
C	gfa2gff	SO:0000856	135	145	.	-	.	ID=14;genome=c
C	gfa2gff	SO:0000856	136	146	.	-	.	ID=15;genome=c
C	gfa2gff	SO:0000856	137	156	.	+	.	ID=4;genome=c
C	gfa2gff	SO:0000856	165	182	.	+	.	ID=9;genome=c;substr=(3,20)
C	gfa2gff	SO:0000856	189	206	.	-	.	ID=9;genome=c;substr=(3,20)
```

