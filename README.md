# Co-Tracer
An Assembly-Based Pipeline for Colocalization Analysis and Host Tracing

## 🏢 Institution

This project was developed at the:  
Hunan Provincial University Key Laboratory for Environmental and Ecological Health, Hunan Provincial University Key Laboratory for Environmental Behavior and Control Principle of New Pollutants, College of Environment and Resources, Xiangtan University, Xiangtan 411105, China

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language](https://img.shields.io/badge/Language-Bash%20%26%20Python-blue.svg)]()

---

### Overview

**Co-Tracer** is a bioinformatics pipeline designed for advanced **colocalization analysis** within metagenomic data. Developed at the Laboratory of Environmental Ecology and Health, Xiangtan University, this tool moves beyond simple co-occurrence statistics by employing an assembly-based approach to identify physical linkages between genes of interest—specifically Antibiotic Resistance Genes (ARGs) and Mobile Genetic Elements (MGEs)—and trace their taxonomic origins.

The core methodology of Co-Tracer is centered on analyzing genomic context from assembled contigs. It reconstructs genomes from raw reads, annotates genetic features, and then precisely measures the physical distance between ARGs and MGEs to confirm their colocalization. This robust approach not only provides strong evidence for potential horizontal gene transfer events but also enables the direct identification of the microbial hosts carrying these linked genetic elements.

---

### Core Features

*   **Assembly-Centric Colocalization**: Utilizes metagenomic assembly (`MEGAHIT`) as its foundation, enabling the analysis of gene order and genomic context, which is crucial for accurate colocalization assessment.
*   **Physical Proximity Mapping**: Calculates the precise base-pair distance between annotated ARGs and MGEs on the same contig, using a defined distance threshold to reliably identify colocalization events.
*   **Host Tracing**: Implements taxonomic classification (`CAT/BAT`) on contigs that harbor ARG-MGE pairs, directly linking mobile resistance determinants to their potential microbial hosts.
*   **Comprehensive Annotation**: Integrates `Prodigal` for gene prediction, `Diamond` (vs. CARD) for ARG annotation, and `HMMER` (vs. Pfam) for identifying MGE-associated protein domains.
*   **Automated End-to-End Workflow**: A series of modular scripts automates the entire process, from raw sequencing data to publication-quality visualizations and statistical reports.
*   **In-depth Statistical Analysis**: Employs Fisher's Exact Test to determine the statistical significance of observed ARG-MGE pairings, with correction for multiple testing.

---

### Analysis Workflow

---

### Dependencies

**1. Bioinformatics Tools:**
*   `fastp`
*   `MEGAHIT`
*   `QUAST`
*   `Prodigal`
*   `DIAMOND`
*   `HMMER`
*   `CAT/BAT` (Optional, for host tracing)

**2. Python 3 and Libraries:**
*   `pandas`, `numpy`, `scipy`, `statsmodels`
*   `matplotlib`, `seaborn`, `matplotlib-venn`, `adjustText`
*   `openpyxl`
*   `networkx` (Optional)

**3. Databases:**
*   **CARD**: [Comprehensive Antibiotic Resistance Database](https://card.mcmaster.ca/)
*   **Pfam-A**: [Protein families database](http://pfam.xfam.org/)
*   **CAT/BAT Database**: Required for taxonomic classification.

---

### Usage

1.  **Clone the repository**
    ```bash
    git clone https://github.com/RaySinclair01/Co-Tracer.git
    cd Co-Tracer
    ```

2.  **Prepare the environment**
    We recommend using `conda` to manage dependencies.
    ```bash
    # Example: create a conda environment
    conda create -n cotracer -c bioconda fastp megahit quast prodigal diamond hmmer cat
    conda activate cotracer
    # Install Python libraries
    pip install pandas numpy matplotlib seaborn openpyxl scipy statsmodels matplotlib-venn adjustText networkx
    ```

3.  **Prepare your data**
    *   Place your raw sequencing files (`*.R1.raw.fastq.gz` and `*.R2.raw.fastq.gz`) into the `./raw_data` directory.
    *   Run the database preparation script (`scripts/prepare_databases.sh`) or ensure pre-built databases are available at the paths specified in the scripts.

4.  **Run the pipeline**
    Execute the scripts sequentially as laid out in `Analysis_MGEs.txt` to run the full workflow.

---

### How to Cite

If you use Co-Tracer in your research, please cite this repository. A manuscript is in preparation, and this section will be updated with a formal citation upon publication.

*   **Project:** Co-Tracer
*   **Institution:** Laboratory of Environmental Ecology and Health, Xiangtan University
*   **GitHub:** https://github.com/RaySinclair01/Co-Tracer

---

### Contact

For questions, bug reports, or suggestions, please contact **lifeng6220@xtu.edu.cn** or open an issue on this GitHub repository.

---

### License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
