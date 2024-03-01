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
Reminder: running MAGPIE on single CPU may take some time because single process of autoFE.

#### Input format
MAGPIE supports variants in CSV format as input. The input file should contain at least 5 columns in the header as follows. [Sample file](data/datasets/test.csv)

|  Chr  | Start |  End  |  Ref  |  Alt  |  ...  |
| ----- | ----- | ----- | ----- | ----- | ----- |

#### Use pretrained MAGPIE model to predict variants

1. Install packages listed in requirements or use magpie conda environment.

##### Annotated variants

2. Run `source magpie.sh --mode pred --test_file [filepath] --file_state annotated --visualization` e.g. `source magpie.sh --mode pred --test_file data/datasets/test.csv --file_state annotated --visualization`


##### Unannotated variants

2. Run `source magpie.sh --mode pred --test_file [filepath] --file_state unannotated --visualization` e.g. `source magpie.sh --mode pred --test_file data/datasets/test.csv --file_state unannotated --visualization`

Results would be saved in `data/result`.


#### Train MAGPIE model from denovo
1. Download and decompress required database for annotating using `bash download.sh`.
2. Apply for AnnoVar access, and place all execuatable annotation tools in `./annovar`.
3. Apply for OMIM database access, and place 'genemap2.txt' in `data/annotation_database`.
4. Install packages and pieces of software manually or use `Dockerfile` to create and run a docker image. 
    ```
    docker build -t magpie .
    docker run -it magpie
    ```

5. Run `source magpie.sh --mode train --input_file [filepath]` e.g. `source magpie.sh --mode train --input_file data/datasets/denovo.csv`.
Model would be saved in `output/result/`.

### Technical support
Please feel free to contact me(yichengliu at zju.edu.cn) for technical support.