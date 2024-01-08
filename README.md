# BEguider

**Welcome to the BEguider!**

**BEguider is a deep-learning model, could design sgRNA for multiple gene loci in batches to predict the efficiency and products of each sgRNA.**

- Base editors with Strict Limits of NGG PAM: `ABE7.10-NGG`, `BE4-NGG`
- PAM-less base editors: `ABEmax-SpRY`, `ABE8e-SL-SpRY`, `ABE8e-NL-SpRY`,`BE4max-SpRY`, `FNLS-YE1-SpRY`, `YE1-SpRY`
- The pairs of Base Editors and SNVs:
  - For **ABE**: _pos_  A > G ,  _neg_  T > C
  - For **CBE**: _pos_  C > T ,  _neg_  G > A

We offer a variety of ways to use it:

- Our [online website](http://beguider.bmicc.org/)
  - Easy to use
- Docker of local Web tool
  - Fast, no environment configuration is required
- A local Web tool / command line tool
  - Customizable

Data from ClinVar were predicted by PAM-less base editor for reference, please download `datas/*.csv`.

## Installation

### Docker

The docker image is available at [Docker Hub](https://hub.docker.com/repository/docker/moloch0/beguider/general).

Please use `docker pull xxx/xxx`to download the image.

```shell
docker pull moloch0/beguider
```

Then, you can run the image with the following command:

```shell
docker run -p 8501:8501 beguider
```

The website will be available at `http://localhost:8501`.

### Local Web tool / command line tool

First thing first, clone this repository.

```shell
git clone https://github.com/Wangxiaoyue-lab/BEguider.git
```

[Miniconda](<https://docs.conda.io/en/latest/miniconda.html>) is strongly recommended to run the local Web tool / command line tool.

We provide a `environment.yml` file for you to create the environment.

```shell
conda env create -f environment.yml
```

Then, you can activate the environment with the following command:

```shell
conda activate beguider
```

`pip` is still needed to install some dependencies.

```shell
pip install requirements.txt
```

**_! Please download the release of [model weight](<https://github.com/Wangxiaoyue-lab/BEguider/archive/refs/tags/model.tar.gz>) and put it in the `BEguider_v3\saved_model` folder._**

Finally, you can run the local Web tool with the following command:

```shell
cd webs_main
streamlit run Intro.py
```

The website will be available at `http://localhost:8501`.

Also, you can run the command line tool with the following command:

```shell
cd BEguider_V3
python BEguider.py [options]
```

## Usage of  the command line tool

The options are as follows:

- `-g` or `--genes` : A txt file including genes and sequences separated by comma(,).
- `-c` or `--chromosome` : A txt file including chromosomes, coordinates and genetic type separated by comma(,). Genetic type: `r` means editing wild genes, `a` means editing mutant genes.
- `-s` or `--rsID` : A txt file including rsID and genetic type separated by comma(,).
- `-b` or `--BaseEditor` : Base Editors: `ALL` / A Specific BE name. Design sgRNA for a specific BE in optional BEs or for all BEs.
  
- _If you have installed [casoffinder](<http://www.rgenome.net/cas-offinder/>), the following two options are available._
  - `-f` or `--offtarget` : True / False. True: predict off-target sites in hg38 genome. (Default = False)
  - `-m` or `--mismatch` : Allowed maximum mismatch site between sgRNA and genome while searching for off-target sites. (Default = 3)
- `-o` or `--output` : Output directory. (Default = current directory)

## Examples

1. use genes as input file:

   ```csv
   Genes,Seqs
   TP53,CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCT
   BRCA1,TGGCTGAAGAATTTGCTAAGCAATCAGGAAAGCTGGTGG
   ARID1A,CTCCACCGAAGGGAGGACCCACTGCCCCCAGCCGGGGTCTCG
   ```

2. use chromosomes and coordinates as input file:

   ```csv
   Chrom,Coordinate,Type
   chr1,145634,r
   chrX,87632,a
   ```

3. use rsIDs as input file:

   ```csv
   SNP,Type
   rs5297,r
   rs12603332,a
   rs75456785,r
   ```
