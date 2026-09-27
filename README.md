# Stanford Multi-omics Hackathon 2026 Track 2

## Omic Discordance Explained

*Omic layers within the MoTrPAC data are often discordant. Can this discordance be modeled or explained?*

### Challenge

Compare two or more compatible layers — such as RNA, total protein, and post-translational modifications — within a defined tissue, and identify a model that explains when they agree and when they disagree.

### Data

Compatible MoTrPAC transcriptomic, proteomic, and PTM measurements with protein, site, pathway, and interaction annotations.

### Potential Outputs

- A catalog of concordant and discordant events
- A predictive model of discordance

> [!IMPORTANT]
> PTM signal is not automatically a measure of modification occupancy, enzyme activity, or functional consequence.

### Project Snapshot 
#### Title: MoTrPAC Explorer 
#### One-line Purpose: Comprehensive query of MoTrPAC stack based on a disease signatures
#### Team:
##### Sheng-Ye Wang: AI Engineer  
##### Christopher Thai: Team Lead
##### Bradley Haraguchi: Chief Developer 
##### Daniyal Rahman: Project Manager
##### Nur-Taz Rahman: Quality Assurance/Resident Mom

### Intended Users
Clinicians & researchers looking for comprehensive OMICS data based on their disease signature of interest 

### Why It Matters
Tools are required to reliably query publicly available OMICS datasets to accelerate hypothesis generation & target validation. 

### Research Question:
Is there a relationship between exercise and disease outcomes?

#### Problem:
, hypothesis or objective, scope, and 

##### success criteria
Visualizations for the OMICS layers 
Side-by-side comparison of transcriptomics & proteomics 

### Workflow 
#### See WORKFLOWS.md

### Setup 
#### Prerequisites + versions
#### Install commands
#### Data downloads + pthas 
#### Containers or notebooks 


### Inputs & Outputs
Input is a list of disease signatures 
Output is visualization (heatmaps, line charts, dot plots) from MoTrPAC data. 

#### Quick start
#### See BOOTSTRAP.md

#### FastAPI  
#### Dashboards: See DESIGN_PROVENANCE.md 

### Methods
#### FastAPI to connect to MoTrPAC; live query, Codex, Kiro, Claude, & ChatGPT. 
#### License: MIT License 

### Validation 
#### Small test dataset/data subset, Expected Output, screenshots or plots, 

##### Failure Modes 
No output generated for certain gene lists 


