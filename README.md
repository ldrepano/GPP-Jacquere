# GPP Jacquere Manuscript Repository

This repository features the source data and code used in the development of the CRISPick Aggregate CFD metric and Jacquere Library (including validation screen data and analysis scripts) as presented in "Balancing off-target and on-target considerations for optimized Cas9 CRISPRko library design" (Drepanos et al., *Cell Genomics*, 2025). 

External source data that exceeds 100MB is omitted, and information on the origin of such data is provided in the corresponding jupyter notebooks.

For calculation of Aggregate CFD scores and other guide selection criteria implemented in Jacquere, please refer to the GPP library design webtool, [CRISPick](broad.io/crispick)

To reproduce this code, it is recommended to clone this repository, and download all required packages within a new virtual environment. 

```
git clone https://github.com/ldrepano/GPP-Jacquere
cd GPP-Jacquere
python3 -m venv run_jacquere_code
source run_jacquere_code/bin/activate
pip install -r requirements.txt 
```
