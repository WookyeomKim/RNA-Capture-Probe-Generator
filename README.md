# Hybridization Capture Probe Designer for RNA Extraction

This repository provides a Python tool for designing sequence-specific hybridization capture probes. It is particularly effective for designing a minimal number of DNA capture probes, such as one or two, for direct RNA extraction from plasma and whole blood, as described in our research.

The algorithm performs a single-pass screening of genomic sequences to rank candidate probes based on mismatch tolerance, GC content, and thermodynamic binding stability calculated using NUPACK. Accessibility, defined as the probability that contiguous 10 nt within the capture site are fully unpaired, is also calculated using the RNAplfold module in ViennaRNA.

A detailed description of the algorithm, its underlying rationale, and experimental validation will be provided in our manuscript, currently under submission.


## 1. Key Features
- **Single-Pass Ranking:** Evaluates all potential sites once and ranks them by performance.
- **Memory Optimization:** Efficiently handles large FASTA datasets using batch processing and incremental history saving.
- **Thermodynamic Validation:** Integrates NUPACK for precise binding energy (Delta G) calculations.
- **Accessibility Scan:** Uses `RNAplfold` to ensure probes target accessible regions of the target RNA.

## 2. Prerequisites
Before running the script, ensure the following are installed on your system:

### **External Software**
- **ViennaRNA Package (RNAplfold):** This is required for accessibility scanning. 
  - Install via Conda: `conda install -c bioconda viennarna`
  - Ensure `RNAplfold` is in your system's PATH.

### **Python Packages**
Install the necessary Python libraries using the provided `requirements.txt`:
```bash
pip install -r requirements.txt


## 2. Installation & Prerequisites

### Step 1: Download the Code
Download this repository as a ZIP file, or clone it using Git:
```bash
git clone https://github.com/WookyeomKim/RNA-Capture-Probe-Generator.git
cd [RNA-Capture-Probe-Generator]

###Step 2: Install External Software (ViennaRNA)
This tool requires ViennaRNA for accessibility scanning. It must be installed on your system (not via pip). The easiest way to install it is using Conda:

Bash

conda install -c bioconda viennarna


###Step 3: Install Python Packages
Install the required Python dependencies using the provided requirements.txt:

Bash

pip install -r requirements.txt

##3. Usage
To design custom capture probes for your specific targets:

Prepare Data: Place your target genomic sequences in a FASTA format file (e.g., HIV_2021_M_with_CRFs.fasta) in the same folder as the script.

Configure Parameters: Open Capture_Probe_Finder.py in a text editor. Modify the settings inside the ProbeParameters class at the top of the file to match your experimental conditions (e.g., sequence_file, mismatch_tolerance, temp, deltaG_standard).

Run the Script: Execute the Python file from your terminal:

Bash

python Capture_Probe_Finder.py

##4. Output Data
The tool generates an Excel file containing the ranked probes with the following key columns:

Rank: Overall performance ranking.

Initial_Matches_(Mismatch): Number of targets matched under the mismatch criteria.

Captured_Targets_(DeltaG): Number of targets captured under thermodynamic criteria.

Coverage_Percent: Percentage of input sequences captured.

Accessibility_Mean & Percentile(%): Relative accessibility score (higher percentile means more accessible).

Thermodynamic Values: Mean, Min, Max, and Std of Delta G.

##5. Citation
A detailed description of the algorithm and experimental validation will be provided in our manuscript (currently under submission). Please cite our work if you use this tool in your research.

##6. License
This project is licensed under the MIT License.