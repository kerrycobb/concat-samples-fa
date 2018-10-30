# concat-samples-fa

Generate concatenated consensus sequence for each individual from Stacks populations fasta_samples output. Effectively generating the phylip_var_all output.

#### Installation
Nim must be installed in order to compile. For instructions on installing Nim, see the documentation [here](https://nim-lang.org/install.html)

Download script `wget https://raw.githubusercontent.com/kerrycobb/concat-samples-fa/master/concat_samples_fa.nim`

Compile script with `nim --run c -o:concat-samples-fa concat_samples_fa.nim`

Place compiled binary in your path with `mv concat-samples-fa /usr/local/bin`


#### Usage
```
Usage:
  concat-fasta-samples INPUT [OUTPUT | -f FORMAT]

Arguments:
  INPUT     Input file path
  OUTPUT    Output file path

Options:
  -h, --help   Show this screen
  -f FORMAT    Output file format, not needed if using OUTPUT argument,
                expects: 'nexus' or 'phylip'
```

#### Examples
`concat-samples-fa input.fa` outputs input.phy and input.part files.
`concat-samples-fa input.fa -nexus` outputs input.nex containing a partition block
`concat-samples-fa input.fa output.phy` outputs output.phy and output.part
