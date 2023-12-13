import argparse
import os.path
import time
import pandas as pd
from predict import predict
from train import train
from data_process import prepare_input, prepare_training_data
from annotation import annotate
from visualization import visualize

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
parser = argparse.ArgumentParser()
parser.add_argument('--input_file', type=str, dest='input_file', default=os.path.join(root, 'data/datasets/input.csv'),
                    help='Filename of input file for denovo training.')
parser.add_argument('--test_file', type=str, dest='test_file', default=os.path.join(root, 'data/datasets/test.csv'),
                    help='Filename of test file.')
parser.add_argument('--model_file', type=str, dest='model_file', default=os.path.join(root, 'data/result/MAGPIE.model'),
                    help='Filename of trained model.')
parser.add_argument('--feature', type=str, dest='feature', default=os.path.join(root, 'data/result/openFE.features'),
                    help='Filename of autoFE generated features.')
parser.add_argument('--selection', type=str, dest='selection', default=os.path.join(root, 'data/result/selection.csv'),
                    help='Filename of autoFE selected features.')
parser.add_argument('--spliceai_output', type=str, dest='spliceai_output', required=False,
                    help='Filename of output file of SpliceAI annotation.')
parser.add_argument('-m', '--mode', type=str, dest='mode', required=True,
                    help='Mode of running magpie, train/pred supported.')
parser.add_argument('--visualization', action='store_true', default = False,
                    help = 'Enable MAGPIE results visualization mode. Files are stored in /output/visualization/')
args = parser.parse_args()
input_file = args.input_file
test_file = args.test_file
model_file = args.model_file
feature = args.feature
selection = args.selection
spliceai_output = args.spliceai_output
mode = args.mode
visualization = args.visualization

if mode == 'pred':
    test = pd.read_csv(test_file, low_memory = False)
    test_pred = predict(test, feature, selection, model_file, os.path.splitext(os.path.basename(test_file))[0])
    if visualization:
        visualize(test_pred, os.path.splitext(os.path.basename(test_file))[0])

elif mode == 'prepare':
    prepare_input(input_file, annovar_dir = os.path.join(root, 'data/output/annovar'), spliceai_dir = os.path.join(root, 'data/output/spliceai'))

elif mode == 'merge':
    data = annotate(input_file, spliceai_output)
    prepare_training_data(data).to_csv(os.path.join(root, f"data/temp/{os.path.splitext(os.path.basename(input_file))[0].replace('.hg38_multianno', '')}.csv"), index = False)


elif mode == 'train':
    train(input_file)

