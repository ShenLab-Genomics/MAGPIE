## MAGPIE framework
### Basic requirements 
1. Python 3.9.12.
2. Python packages listed in requirements.

### Additional requirements for training from denovo
3. SpliceAI (conda environment created by spliceai.yml recommended).
4. MATLAB CLI.
5. AnnoVar (register required).
6. OMIM database(application required).
### Usage
#### Use pretrained MAGPIE model to predict annotated files
`python magpie.py --mode pred --test_file [filepath] --visualization` e.g. `bash magpie.sh --mode pred --test_file data/datasets/test.csv --visualization`
Results would be saved in `data/result`.

Reminder: running MAGPIE on single CPU may take some time because single process of autoFE.
#### Train MAGPIE model from denovo
1. Install and decompress all required packages and pieces of software. Download and decompress required database for annotating using `bash download.sh`.
2. Apply for OMIM database access, and place 'genemap2.txt' in `data/annotation_database`
3. Run `source magpie.sh --mode train --input_file [filepath]` e.g. `source magpie.sh --mode train --input_file data/datasets/denovo.csv`.

Model would be saved in `output/result/`.

### Technical support
Please feel free to contact me(yichengliu at zju.edu.cn) for technical support.



