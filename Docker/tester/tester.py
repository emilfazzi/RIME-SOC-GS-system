import os
import argparse
import pathlib

def main():

    my_parser = argparse.ArgumentParser(description='Telemetry to raw', epilog='Example: python script.py -i inputFolder -o outputFolder -l outputFolderLogFiles -vid productVersionID')

    my_parser.add_argument('-d',
                       '--directory',
                       action='store',
                       required=True,
                       help='directory to scan')
    args = (my_parser.parse_args())
    custom_directory = pathlib.Path(args.directory)


    current_directory = os.getcwd()
    print(f"Contents of the current directory '{current_directory}':")
    list_directory_contents(current_directory)

    print(f"Contents of the current directory '{current_directory}':")
    list_directory_contents(custom_directory)


def list_directory_contents(directory):
    try:
        contents = os.listdir(directory)
        for item in contents:
            print(item)
    except OSError as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
