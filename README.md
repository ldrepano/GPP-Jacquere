# GPP Jacquere Manuscript Repository

This repository contains the code and data required for the development of Aggregate CFD and design of the Jacquere Library. Data that exceeds 100MB is omitted, and information on the origin of such data is provided in the corresponding jupyter notebooks.

For calculation of Aggregate CFD scores and other guide selection criteria implemented in Jacquere, please refer to the GPP library design webtool, [CRISPick](broad.io/crispick)

To reproduce this code, it is recommended to clone this repository, start a virtual environment within the directory, and download all required packages. 

```
git clone https://github.com/ldrepano/GPP-Jacquere
cd GPP-Jacquere
python3 -m venv ~/run_jacquere_code
source ~/run_jacquere_code/bin/activate
pip install -r requirements.txt 
```
