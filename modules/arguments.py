import argparse
from argparse import RawTextHelpFormatter

def get_args():
    parser = argparse.ArgumentParser(description="A simple argument parser", formatter_class=RawTextHelpFormatter)
    
    parser.add_argument("--xyz", metavar=("file1"), help = "Input a xyz file in standard format")
    
    
    args = parser.parse_args()
    return args