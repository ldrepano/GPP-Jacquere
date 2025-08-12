# GPP Jacquere Manuscript Repository

Code to accompany the manuscript "Balancing off-target and on-target considerations for optimized Cas9 CRISPRko library design" 

This repository contains the code and data required for the development of Aggregate CFD and design of the Jacquere Library. Data that exceeds 100MB is omitted, and information on the origin of such data is provided in the corresponding jupyter notebooks.

For calculation of Aggregate CFD scores and other guide selection criteria implemented in Jacquere, please refer to the GPP library design webtool, [CRISPick](broad.io/crispick)

To reproduce this code, it is recommended to clone this repository, and download all required packages within a new virtual environment. 

```
git clone https://github.com/ldrepano/GPP-Jacquere
cd GPP-Jacquere
python3 -m venv run_jacquere_code
source run_jacquere_code/bin/activate
pip install -r requirements.txt 
```
