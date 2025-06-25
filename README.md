# Project Programming Course

This file contains informations about the Project for the exam Programming Approaches for Bioinformatics. AY 2025, Matteo Poncina

## Getting Started

Follow these steps to set up the project and run the analysis.

### 1. Build the Docker Image

```bash
docker build -t final_project_image .
This will create a Docker image with all the required software and packages for the analysis.
```

### 2. Run the Docker Container
To execute the analysis, run the following command. This will mount your local directory to the container:

```bash
docker run -it -v "${PWD}:/usr/src/app" final_project_image /bin/bash
Rscript Final_project.R
```


### Files included
fullMatrix: an R package to perform the first step of the analysis, which is to build the full Matrix starting from barcodes, features and matrix files. Vignette is included.
Dockerfile: The configuration to build the Docker image with the necessary software and libraries.
Final_project.R: the entire R script used to perform analysis and create plots.
Project_markdown.html: The main R Markdown script used to perform the analysis and generate the report.

