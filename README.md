# Intrusion Detection System (IDS) Research Project

## Project Overview

This repository contains the R scripts, datasets, intermediate data, processed results, figures, and reports used in our Intrusion Detection System (IDS) research project.

The project includes data preparation, feature processing, dataset management, analysis, and evaluation-related workflows. It also maintains separate safe and leaky dataset splits for research and experimental analysis.

## Repository Structure

```text
2.Final_codes/
│
├── data/
│   ├── interim/
│   │   ├── cor_matrix.rds
│   │   ├── inventory.rds
│   │   ├── keep_features.rds
│   │   └── working_set.rds
│   │
│   └── processed/
│       ├── arms.rds
│       ├── balanced_idx.rds
│       ├── macro_f1.rds
│       ├── overlap_curve.rds
│       ├── random_baseline.rds
│       ├── selection.rds
│       ├── split_leaky.rds
│       ├── split_safe.rds
│       ├── substitute_quality.rds
│       └── uncertainty.rds
│
├── figures/
├── reports/
│
├── 01_inventory.R
├── 02_*.R
├── ...
├── 12_substitute_quality.R
│
├── Final_codes.Rproj
├── .gitignore
├── .gitattributes
└── README.md
```

## Dataset Description

The repository contains intermediate and processed datasets in RDS format.

| File | Description |
|---|---|
| `working_set.rds` | Intermediate working dataset |
| `split_safe.rds` | Safe dataset split |
| `split_leaky.rds` | Leaky dataset split |
| `arms.rds` | Processed data object |
| `balanced_idx.rds` | Balanced dataset indices |
| `keep_features.rds` | Retained feature information |
| `cor_matrix.rds` | Correlation matrix |

Other processed files contain experiment-related results and analysis outputs.

## Technologies Used

- R Programming Language
- RStudio
- RDS data storage
- Git and GitHub
- Git Large File Storage (Git LFS)

## How to Clone the Repository

First, install Git and Git LFS.

Clone the repository:

```bash
git clone YOUR_GITHUB_REPOSITORY_URL
```

Move into the project directory:

```bash
cd 2.Final_codes
```

Download the actual large dataset files:

```bash
git lfs pull
```

Open the project using the `Final_codes.Rproj` file in RStudio.

## Running the Project

1. Clone the repository and download the LFS files.
2. Open `Final_codes.Rproj` in RStudio.
3. Check the required R packages and dependencies.
4. Run the R scripts in their intended execution order.
5. Review the generated figures, reports, and processed outputs.

## Large File Management

Large `.rds` files are managed using Git LFS.

The `.gitattributes` file specifies which files are tracked by Git LFS, allowing large datasets to be stored separately from the regular Git repository objects.

## Research Team

This repository contains the implementation and data-related resources for our IDS research project.

## Note

Dataset source, research references, package versions, and detailed experimental configurations should be documented here as the research documentation is finalized.
