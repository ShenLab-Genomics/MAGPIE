## MAGPIE framework
### Basic requirements 
1. Python 3.9.12.

2. Python packages listed in requirements.

Use `conda env create -f magpie.yml` to create the conda environment is recommanded. 

### Additional requirements for training from denovo
3. SpliceAI (Use `conda env create -f spliceai.yml` to create the conda environment is recommanded).

4. MATLAB CLI.
5. AnnoVar (register required).
6. OMIM database(application required).

### Usage
#### Use pretrained MAGPIE model to predict annotated files
1. Install packages listed in requirements or use magpie conda environment.

2. Run `bash magpie.sh --mode pred --test_file [filepath] --visualization` e.g. `bash magpie.sh --mode pred --test_file data/datasets/test.csv --visualization`
Results would be saved in `data/result`.

Reminder: running MAGPIE on single CPU may take some time because single process of autoFE.
#### Train MAGPIE model from denovo
1. Download and decompress required database for annotating using `bash download.sh`.
2. Apply for AnnoVar access, and place all execuatable annotation tools in `./annovar`.
3. Apply for OMIM database access, and place 'genemap2.txt' in `data/annotation_database`.
4. Install packages and pieces of software manually or use `Dockerfile` to create a docker image. 
    ```
    docker build -t magpie .
    docker run -it magpie
    ```

5. Run `source magpie.sh --mode train --input_file [filepath]` e.g. `source magpie.sh --mode train --input_file data/datasets/denovo.csv`.
Model would be saved in `output/result/`.

### Technical support
Please feel free to contact me(yichengliu at zju.edu.cn) for technical support.



